#include "clr.h"
#include "clr_tables.h"

#include <string.h>

#define M61 ((UINT64_C(1) << 61) - 1)
#define AT(i, j) (((unsigned)((i) & (CLR_DIM - 1)) * CLR_DIM) + (unsigned)((j) & (CLR_DIM - 1)))

#if CLR_LAYERS < 2
#error "the layer coupling is a Feistel step between distinct layers; CLR_LAYERS must be >= 2"
#endif
#if CLR_DIM != 8
#error "the column RS code (length 8), row words and output fold assume an 8x8 lattice"
#endif

/* ------------------------------------------------------------ GF(2^8) */

#define BYTES(x) (UINT64_C(0x0101010101010101) * (x))

static inline uint8_t x2(uint8_t a) { return (uint8_t)((a << 1) ^ ((a >> 7) * 0x1B)); }

/* x2 on eight packed cells at once. */
static inline uint64_t x2w(uint64_t v)
{
    return ((v & BYTES(0x7F)) << 1) ^ (((v >> 7) & BYTES(0x01)) * 0x1B);
}

static inline uint64_t rotr8(uint64_t v) { return (v >> 8) | (v << 56); }
static inline uint64_t rotl8(uint64_t v) { return (v << 8) | (v >> 56); }

/*
 * CLR_LINEAR_CONTROL builds a control for the analysis scripts: the identical
 * architecture with each nonlinear step swapped for a GF(2)-linear one (S-box
 * -> multiply by 0x03, RS decode -> raw syndrome injection).  Never use it.
 */
static inline uint8_t sub(uint8_t x)
{
#ifdef CLR_LINEAR_CONTROL
    return (uint8_t)(x2(x) ^ x);
#else
    return CLR_SBOX[x];
#endif
}

/* ------------------------------------------------------------ Gold pair */

/* Eight steps of both Galois LFSRs; returns 8 bits of the Gold sequence u XOR v. */
static inline uint8_t gold_byte(clr_state *st)
{
    unsigned iu = (unsigned)(st->gold_u & 0xFF), iv = (unsigned)(st->gold_v & 0xFF);
    uint8_t g = (uint8_t)(CLR_GOLD_OU[iu] ^ CLR_GOLD_OV[iv]);
    st->gold_u = (st->gold_u >> 8) ^ CLR_GOLD_TU[iu];
    st->gold_v = (st->gold_v >> 8) ^ CLR_GOLD_TV[iv];
    return g;
}

/* Inverse of one Galois step s' = (s >> 1) ^ (s & 1 ? mask : 0); needs mask bit 60 set. */
static inline uint64_t galois_unstep(uint64_t s, uint64_t mask)
{
    uint64_t o = (s >> 60) & 1;
    return (((s ^ (o ? mask : 0)) << 1) | o) & M61;
}

static void gold_rewind(clr_state *st, unsigned nbytes)
{
    for (unsigned i = 0; i < 8 * nbytes; i++) {
        st->gold_u = galois_unstep(st->gold_u, CLR_GOLD_MASK_U);
        st->gold_v = galois_unstep(st->gold_v, CLR_GOLD_MASK_V);
    }
}

static inline uint64_t load64le(const uint8_t *p)
{
    return (uint64_t)p[0] | (uint64_t)p[1] << 8 | (uint64_t)p[2] << 16 | (uint64_t)p[3] << 24 |
           (uint64_t)p[4] << 32 | (uint64_t)p[5] << 40 | (uint64_t)p[6] << 48 | (uint64_t)p[7] << 56;
}

static inline void store64le(uint8_t *p, uint64_t v)
{
    for (unsigned b = 0; b < 8; b++)
        p[b] = (uint8_t)(v >> (8 * b));
}

/* Outer layer drives the Gold pair.  Gold state is updated from layer 0 only: invertible. */
static void gold_feedback(clr_state *st)
{
    const uint8_t *c = st->cell[0];
    uint64_t a = load64le(c) ^ load64le(c + 16) ^ load64le(c + 32) ^ load64le(c + 48);
    uint64_t b = load64le(c + 8) ^ load64le(c + 24) ^ load64le(c + 40) ^ load64le(c + 56);
    st->gold_u ^= a & M61;
    st->gold_v ^= b & M61;
}

/* ------------------------------------------------------------ spatial self-interaction */

/*
 * Checkerboard half-step on one layer: every cell of one colour is XORed with
 * S(N + 2E + 4S + 8W + position + Gold), where N/E/S/W are its torus
 * neighbours (all of the other colour, so unchanged during the half-step).
 * Each half-step is therefore an involution: applying it twice is identity.
 *
 * A row is one little-endian word (byte j = column j), so the S-box inputs of
 * a whole row are computed at once: E and W are the row rotated by one cell.
 */
static void half_step_rows(uint64_t row[CLR_DIM], unsigned color, const uint8_t rc[CLR_DIM])
{
    uint64_t upd[CLR_DIM];
    for (unsigned i = 0; i < CLR_DIM; i++) {
        uint64_t n = row[(i - 1) & (CLR_DIM - 1)], s = row[(i + 1) & (CLR_DIM - 1)];
        uint64_t e = rotr8(row[i]), w = rotl8(row[i]);
        uint64_t x = n ^ x2w(e) ^ x2w(x2w(s)) ^ x2w(x2w(x2w(w)))
                   ^ load64le(CLR_GF_EXP + CLR_DIM * i) ^ BYTES(rc[i]);
        uint64_t u = 0;
        for (unsigned j = (i + color) & 1; j < CLR_DIM; j += 2)
            u |= (uint64_t)sub((uint8_t)(x >> (8 * j))) << (8 * j);
        upd[i] = u;
    }
    for (unsigned i = 0; i < CLR_DIM; i++)
        row[i] ^= upd[i];
}

static void load_rows(uint64_t row[CLR_DIM], const uint8_t *c)
{
    for (unsigned i = 0; i < CLR_DIM; i++)
        row[i] = load64le(c + CLR_DIM * i);
}

static void store_rows(uint8_t *c, const uint64_t row[CLR_DIM])
{
    for (unsigned i = 0; i < CLR_DIM; i++)
        store64le(c + CLR_DIM * i, row[i]);
}

static void layer_round(clr_state *st, unsigned k)
{
    uint8_t rc[2 * CLR_DIM];
    uint64_t row[CLR_DIM];
    for (unsigned i = 0; i < 2 * CLR_DIM; i++)
        rc[i] = gold_byte(st);
    load_rows(row, st->cell[k]);
    half_step_rows(row, 0, rc);
    half_step_rows(row, 1, rc + CLR_DIM);
    store_rows(st->cell[k], row);
}

static void layer_unround(clr_state *st, unsigned k)
{
    uint8_t rc[2 * CLR_DIM];
    uint64_t row[CLR_DIM];
    gold_rewind(st, 2 * CLR_DIM);
    for (unsigned i = 0; i < 2 * CLR_DIM; i++)
        rc[i] = gold_byte(st);
    gold_rewind(st, 2 * CLR_DIM);
    load_rows(row, st->cell[k]);
    half_step_rows(row, 1, rc + CLR_DIM);
    half_step_rows(row, 0, rc);
    store_rows(st->cell[k], row);
}

/* ------------------------------------------------------------ Galois self-correction */

/*
 * Treat each column of `src` as a word of a length-8 Reed-Solomon code over
 * GF(2^8) and run the single-error decoder: syndromes S1 = sum c_i a^i,
 * S2 = sum c_i a^(2i), error location a^loc = S2/S1, magnitude e = S1^2/S2.
 * The correction is applied to the same column of `dst` (a different layer),
 * so the step is a Feistel update and loses nothing: running it again undoes it.
 *
 * Syndromes of all 8 columns are computed at once by Horner's rule on row
 * words (a = 0x03, a^2 = 0x05).
 */
static void correct(const uint8_t *src, uint8_t *dst)
{
    uint64_t s1w = load64le(src + CLR_DIM * (CLR_DIM - 1)), s2w = s1w;
    for (int i = CLR_DIM - 2; i >= 0; i--) {
        uint64_t r = load64le(src + CLR_DIM * (unsigned)i);
        s1w = x2w(s1w) ^ s1w ^ r;
        s2w = x2w(x2w(s2w)) ^ s2w ^ r;
    }
    for (unsigned j = 0; j < CLR_DIM; j++) {
        unsigned s1 = (uint8_t)(s1w >> (8 * j)), s2 = (uint8_t)(s2w >> (8 * j));
#ifdef CLR_LINEAR_CONTROL
        dst[AT(0, j)] ^= (uint8_t)s1;
        dst[AT(1, j)] ^= (uint8_t)s2;
        continue;
#endif
        if (s1 && s2) {
            unsigned l1 = CLR_GF_LOG[s1], l2 = CLR_GF_LOG[s2];
            unsigned loc = (l2 + 255 - l1) % 255;
            uint8_t e = CLR_GF_EXP[(2 * l1 + 255 - l2) % 255];
            dst[AT(loc, j)] ^= e;
        }
    }
}

/* ------------------------------------------------------------ the loop */

void clr_round(clr_state *st)
{
    for (unsigned k = 0; k < CLR_LAYERS; k++) {
        layer_round(st, k);
        correct(st->cell[k], st->cell[(k + 1) % CLR_LAYERS]);
    }
    gold_feedback(st);
}

void clr_unround(clr_state *st)
{
    gold_feedback(st);
    for (int k = CLR_LAYERS - 1; k >= 0; k--) {
        correct(st->cell[k], st->cell[(k + 1) % CLR_LAYERS]);
        layer_unround(st, (unsigned)k);
    }
}

/* ------------------------------------------------------------ seeding and output */

void clr_absorb(clr_state *st, const void *data, size_t len)
{
    const uint8_t *p = data;
    for (;;) {
        size_t n = len < CLR_CELLS ? len : CLR_CELLS;
        for (size_t b = 0; b < n; b++)
            st->cell[0][b] ^= p[b];
        p += n;
        len -= n;
        if (n < CLR_CELLS) {                    /* final chunk: pad10*1 */
            st->cell[0][n] ^= 0x01;
            st->cell[0][CLR_CELLS - 1] ^= 0x80;
            break;
        }
        for (uint32_t r = 0; r < st->rounds; r++)
            clr_round(st);
    }
    for (uint32_t r = 0; r < CLR_SEED_ROUNDS; r++)
        clr_round(st);
    memset(st->buf, 0, sizeof st->buf);
    st->buf_pos = CLR_BLOCK;
}

void clr_init(clr_state *st, const void *seed, size_t len, uint32_t rounds)
{
    memset(st, 0, sizeof *st);
    st->rounds = rounds ? rounds : CLR_DEFAULT_ROUNDS;
    st->gold_u = 1;
    st->gold_v = 1;
    for (unsigned b = 0; b < 8; b++)            /* domain-separate seeds by length */
        st->cell[CLR_LAYERS - 1][b] = (uint8_t)((uint64_t)len >> (8 * b));
    clr_absorb(st, seed, len);
}

void clr_block(clr_state *st, uint8_t out[CLR_BLOCK])
{
    for (uint32_t r = 0; r < st->rounds; r++)
        clr_round(st);
    for (unsigned b = 0; b < CLR_BLOCK; b++)
        out[b] = st->cell[0][b] ^ st->cell[0][b + CLR_BLOCK];
}

void clr_fill(clr_state *st, void *out, size_t len)
{
    uint8_t *p = out;
    while (len) {
        if (st->buf_pos == CLR_BLOCK) {
            clr_block(st, st->buf);
            st->buf_pos = 0;
        }
        size_t n = CLR_BLOCK - st->buf_pos;
        if (n > len)
            n = len;
        memcpy(p, st->buf + st->buf_pos, n);
        memset(st->buf + st->buf_pos, 0, n);
        st->buf_pos += (uint32_t)n;
        p += n;
        len -= n;
    }
}

uint64_t clr_u64(clr_state *st)
{
    uint8_t b[8];
    clr_fill(st, b, sizeof b);
    return load64le(b);
}

size_t clr_state_size(void) { return sizeof(clr_state); }
