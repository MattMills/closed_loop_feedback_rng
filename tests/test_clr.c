/* Unit tests.  Includes the implementation directly to reach its internal steps. */
#include "../src/clr.c"

#include <stdio.h>
#include <stdlib.h>

static int failures;

#define CHECK(cond, ...)                                     \
    do {                                                     \
        if (!(cond)) {                                       \
            failures++;                                      \
            printf("FAIL %s:%d: ", __FILE__, __LINE__);      \
            printf(__VA_ARGS__);                             \
            printf("\n");                                    \
        }                                                    \
    } while (0)

static uint64_t sm_state = 0x9E3779B97F4A7C15u;
static uint64_t splitmix(void)
{
    uint64_t z = (sm_state += 0x9E3779B97F4A7C15u);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
    return z ^ (z >> 31);
}

static void random_state(clr_state *st)
{
    memset(st, 0, sizeof *st);
    for (unsigned k = 0; k < CLR_LAYERS; k++)
        for (unsigned c = 0; c < CLR_CELLS; c++)
            st->cell[k][c] = (uint8_t)splitmix();
    st->gold_u = splitmix() & M61;
    st->gold_v = splitmix() & M61;
    st->rounds = CLR_DEFAULT_ROUNDS;
    st->buf_pos = CLR_BLOCK;
}

static int same_state(const clr_state *a, const clr_state *b)
{
    return !memcmp(a->cell, b->cell, sizeof a->cell) && a->gold_u == b->gold_u && a->gold_v == b->gold_v;
}

static uint8_t gmul_slow(uint8_t a, uint8_t b)
{
    uint8_t r = 0;
    while (b) {
        if (b & 1) r ^= a;
        a = x2(a);
        b >>= 1;
    }
    return r;
}

static void test_tables(void)
{
    for (unsigned a = 1; a < 256; a++) {
        CHECK(CLR_GF_EXP[CLR_GF_LOG[a]] == a, "exp/log mismatch at %u", a);
        CHECK(gmul_slow((uint8_t)a, 3) == CLR_GF_EXP[CLR_GF_LOG[a] + 1], "alpha != 0x03 at %u", a);
    }
    for (unsigned a = 0; a < 256; a++) {
        uint8_t inv = 0;
        for (unsigned b = 1; b < 256 && a; b++)
            if (gmul_slow((uint8_t)a, (uint8_t)b) == 1) inv = (uint8_t)b;
        uint8_t s = inv;
        for (unsigned r = 1; r < 5; r++)
            s ^= (uint8_t)((inv << r) | (inv >> (8 - r)));
        s ^= 0x63;
        CHECK(CLR_SBOX[a] == s, "S-box entry %u is not affine(inverse)", a);
    }
}

static uint8_t galois_byte_slow(uint64_t *s, uint64_t mask)
{
    uint8_t out = 0;
    for (unsigned i = 0; i < 8; i++) {
        uint64_t o = *s & 1;
        *s >>= 1;
        if (o) *s ^= mask;
        out |= (uint8_t)(o << i);
    }
    return out;
}

static void test_gold(void)
{
    for (int t = 0; t < 10000; t++) {
        clr_state st;
        random_state(&st);
        uint64_t u = st.gold_u, v = st.gold_v;
        uint8_t ref = galois_byte_slow(&u, CLR_GOLD_MASK_U) ^ galois_byte_slow(&v, CLR_GOLD_MASK_V);
        CHECK(gold_byte(&st) == ref, "table-driven Gold byte differs from bit-serial");
        CHECK(st.gold_u == u && st.gold_v == v, "table-driven LFSR step differs from bit-serial");
        uint64_t s = splitmix() & M61, s1 = s;
        galois_byte_slow(&s1, CLR_GOLD_MASK_U);
        for (int i = 0; i < 8; i++) s1 = galois_unstep(s1, CLR_GOLD_MASK_U);
        CHECK(s1 == s, "galois_unstep is not the inverse of a step");
    }
    /* m-sequence sanity: from state 1 the u register must not return to 1 within 2^20 steps */
    uint64_t s = 1;
    for (unsigned i = 1; i <= (1u << 20); i++) {
        uint64_t o = s & 1;
        s = (s >> 1) ^ (o ? CLR_GOLD_MASK_U : 0);
        if (s == 1) {
            CHECK(0, "u register has period %u", i);
            break;
        }
    }
}

static void half_step(uint8_t *c, unsigned color, const uint8_t rc[CLR_DIM])
{
    uint64_t row[CLR_DIM];
    load_rows(row, c);
    half_step_rows(row, color, rc);
    store_rows(c, row);
}

/* Cell-by-cell definition of the half-step; the word-parallel version must match it. */
static void half_step_ref(uint8_t *c, unsigned color, const uint8_t rc[CLR_DIM])
{
    for (unsigned i = 0; i < CLR_DIM; i++) {
        for (unsigned j = (i + color) & 1; j < CLR_DIM; j += 2) {
            uint8_t n = c[AT(i - 1, j)], e = c[AT(i, j + 1)];
            uint8_t s = c[AT(i + 1, j)], w = c[AT(i, j - 1)];
            uint8_t x = (uint8_t)(n ^ x2(e) ^ x2(x2(s)) ^ x2(x2(x2(w))) ^ CLR_GF_EXP[AT(i, j)] ^ rc[i]);
            c[AT(i, j)] ^= CLR_SBOX[x];
        }
    }
}

static void test_half_step_matches_reference(void)
{
    for (int t = 0; t < 10000; t++) {
        clr_state a, b;
        random_state(&a);
        b = a;
        uint8_t rc[CLR_DIM];
        for (unsigned i = 0; i < CLR_DIM; i++) rc[i] = (uint8_t)splitmix();
        half_step(a.cell[0], t & 1, rc);
        half_step_ref(b.cell[0], t & 1, rc);
        CHECK(same_state(&a, &b), "word-parallel half-step differs from the cell-by-cell definition");
    }
}

/* Direct definition of the column RS decode-and-apply step. */
static void correct_ref(const uint8_t *src, uint8_t *dst)
{
    for (unsigned j = 0; j < CLR_DIM; j++) {
        uint8_t s1 = 0, s2 = 0;
        for (unsigned i = 0; i < CLR_DIM; i++) {
            uint8_t a_i = CLR_GF_EXP[i], a_2i = CLR_GF_EXP[2 * i];
            s1 ^= gmul_slow(src[AT(i, j)], a_i);
            s2 ^= gmul_slow(src[AT(i, j)], a_2i);
        }
        if (!s1 || !s2) continue;
        uint8_t x = 0, e = 0;                   /* x = S2/S1 = a^loc, e = S1^2/S2 */
        for (unsigned b = 1; b < 256; b++) {
            if (gmul_slow(s1, (uint8_t)b) == s2) x = (uint8_t)b;
            if (gmul_slow(s2, (uint8_t)b) == gmul_slow(s1, s1)) e = (uint8_t)b;
        }
        dst[AT(CLR_GF_LOG[x], j)] ^= e;
    }
}

static void test_correct_matches_reference(void)
{
    for (int t = 0; t < 2000; t++) {
        clr_state a, b;
        random_state(&a);
        if (t < 8)                              /* include columns with zero syndromes */
            memset(a.cell[0], 0, CLR_CELLS);
        b = a;
        correct(a.cell[0], a.cell[1]);
        correct_ref(b.cell[0], b.cell[1]);
        CHECK(same_state(&a, &b), "Horner-rule RS decode differs from the direct definition");
    }
}

static void test_components_invertible(void)
{
    for (int t = 0; t < 2000; t++) {
        clr_state a, b;
        random_state(&a);
        b = a;
        uint8_t rc[CLR_DIM];
        for (unsigned i = 0; i < CLR_DIM; i++) rc[i] = (uint8_t)splitmix();
        half_step(a.cell[1], t & 1, rc);
        CHECK(!same_state(&a, &b) , "half-step changed nothing");
        half_step(a.cell[1], t & 1, rc);
        CHECK(same_state(&a, &b), "half-step is not an involution");
        correct(a.cell[2], a.cell[3]);
        correct(a.cell[2], a.cell[3]);
        CHECK(same_state(&a, &b), "RS correction coupling is not an involution");
        gold_feedback(&a);
        gold_feedback(&a);
        CHECK(same_state(&a, &b), "Gold feedback is not an involution");
    }
}

static void test_round_invertible(void)
{
    for (int t = 0; t < 10000; t++) {
        clr_state a, b;
        random_state(&a);
        b = a;
        clr_round(&a);
        CHECK(!same_state(&a, &b), "round changed nothing");
        clr_unround(&a);
        CHECK(same_state(&a, &b), "unround(round(s)) != s");
        clr_unround(&a);
        clr_round(&a);
        CHECK(same_state(&a, &b), "round(unround(s)) != s");
    }
    /* many passes forward then the same number back */
    clr_state a, b;
    random_state(&a);
    b = a;
    for (int r = 0; r < 1000; r++) clr_round(&a);
    for (int r = 0; r < 1000; r++) clr_unround(&a);
    CHECK(same_state(&a, &b), "1000 rounds are not undone by 1000 unrounds");
}

static void test_fill(void)
{
    clr_state a, b;
    clr_init(&a, "abc", 3, 0);
    clr_init(&b, "abc", 3, 0);
    uint8_t x[1000], y[1000];
    clr_fill(&a, x, sizeof x);
    for (size_t off = 0, n = 1; off < sizeof y; off += n, n = n % 37 + 1) {
        if (off + n > sizeof y) n = sizeof y - off;
        clr_fill(&b, y + off, n);
    }
    CHECK(!memcmp(x, y, sizeof x), "chunked fill differs from one-shot fill");
}

static void test_seed_avalanche(void)
{
    uint8_t seed[32] = {0}, base[CLR_BLOCK], out[CLR_BLOCK];
    clr_state st;
    clr_init(&st, seed, sizeof seed, 0);
    clr_block(&st, base);
    unsigned long total = 0;
    for (unsigned bit = 0; bit < 256; bit++) {
        seed[bit / 8] ^= (uint8_t)(1u << (bit % 8));
        clr_init(&st, seed, sizeof seed, 0);
        clr_block(&st, out);
        seed[bit / 8] ^= (uint8_t)(1u << (bit % 8));
        unsigned d = 0;
        for (unsigned b = 0; b < CLR_BLOCK; b++) d += (unsigned)__builtin_popcount(base[b] ^ out[b]);
        CHECK(d > 80 && d < 176, "seed bit %u flips only %u of 256 output bits", bit, d);
        total += d;
    }
    double mean = total / 256.0;
    CHECK(mean > 124 && mean < 132, "mean seed avalanche %.2f bits (expected ~128)", mean);
}

static void test_known_answer(void)
{
    static const uint8_t expect[CLR_BLOCK] = {
#include "kat_abc.inc"
    };
    clr_state st;
    uint8_t out[CLR_BLOCK];
    clr_init(&st, "abc", 3, 0);
    clr_block(&st, out);
    CHECK(!memcmp(out, expect, CLR_BLOCK), "known-answer vector for seed \"abc\" changed");
}

int main(void)
{
    CHECK(clr_state_size() == sizeof(clr_state), "state size");
    test_tables();
    test_gold();
    test_half_step_matches_reference();
    test_correct_matches_reference();
    test_components_invertible();
    test_round_invertible();
    test_fill();
    test_seed_avalanche();
    test_known_answer();
    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }
    printf("all tests passed\n");
    return 0;
}
