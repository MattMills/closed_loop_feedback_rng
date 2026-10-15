#include "clr.h"

#include <string.h>

#define AT(i, j) (((unsigned)((i) & (CLR_DIM - 1)) * CLR_DIM) + (unsigned)((j) & (CLR_DIM - 1)))

/* Erase a buffer in a way the compiler may not drop as a dead store. */
static void wipe(void *p, size_t n)
{
    volatile uint8_t *v = p;
    while (n--)
        *v++ = 0;
}

/*
 * Second-order XOR loop between the two generators' interior layers.  The two
 * directions use different cell permutations (transpose; torus shift by
 * (3, 5)) so A and B can never fall into mirrored lock-step.  Each XOR reads a
 * layer the other XOR does not write, so the pair is invertible.
 */
static void couple(clr_dual *d)
{
    for (unsigned i = 0; i < CLR_DIM; i++)
        for (unsigned j = 0; j < CLR_DIM; j++) {
            d->a.cell[CLR_LAYERS - 1][AT(i, j)] ^= d->b.cell[CLR_LAYERS - 2][AT(j, i)];
            d->b.cell[CLR_LAYERS - 1][AT(i, j)] ^= d->a.cell[CLR_LAYERS - 2][AT(i + 3, j + 5)];
        }
}

void clr_dual_permute(clr_dual *d)
{
    for (uint32_t r = 0; r < d->a.rounds; r++) {
        clr_round(&d->a);
        clr_round(&d->b);
        couple(d);
    }
}

void clr_dual_unpermute(clr_dual *d)
{
    for (uint32_t r = 0; r < d->a.rounds; r++) {
        couple(d);                              /* XORs of unchanged layers: self-inverse */
        clr_unround(&d->b);
        clr_unround(&d->a);
    }
}

static void xor_state(clr_state *dst, const clr_state *src)
{
    for (unsigned k = 0; k < CLR_LAYERS; k++)
        for (unsigned c = 0; c < CLR_CELLS; c++)
            dst->cell[k][c] ^= src->cell[k][c];
    dst->gold_u ^= src->gold_u;
    dst->gold_v ^= src->gold_v;
}

void clr_dual_ratchet(clr_dual *d)
{
    clr_state a0 = d->a, b0 = d->b;
    clr_dual_permute(d);
    xor_state(&d->a, &a0);                      /* feed-forward: S <- P(S) xor S */
    xor_state(&d->b, &b0);
    wipe(&a0, sizeof a0);
    wipe(&b0, sizeof b0);
}

void clr_dual_init(clr_dual *d, const void *seed, size_t len, uint32_t rounds)
{
    static const char twin[] = "CLR interior twin";
    memset(d, 0, sizeof *d);
    clr_init(&d->a, seed, len, rounds);
    clr_init(&d->b, seed, len, rounds);
    clr_absorb(&d->b, twin, sizeof twin - 1);   /* domain-separate B from A */
    clr_dual_ratchet(d);
    clr_dual_ratchet(d);
    d->buf_pos = CLR_BLOCK;
}

void clr_dual_absorb(clr_dual *d, const void *data, size_t len)
{
    clr_absorb(&d->b, data, len);
    clr_dual_ratchet(d);
    wipe(d->buf, sizeof d->buf);
    d->buf_pos = CLR_BLOCK;
}

/* OS entropy into the occluded twin only; it reaches the output solely through the loop. */
static void reseed_twin(clr_dual *d)
{
    uint8_t ent[32];
    clr_os_entropy(ent, sizeof ent);
    clr_absorb(&d->b, ent, sizeof ent);
    wipe(ent, sizeof ent);
    d->since_reseed = 0;
}

void clr_dual_init_occluded(clr_dual *d, const void *seed, size_t len, uint32_t rounds, uint32_t reseed_every)
{
    clr_dual_init(d, seed, len, rounds);
    reseed_twin(d);
    clr_dual_ratchet(d);
    clr_dual_ratchet(d);
    d->reseed_every = reseed_every ? reseed_every : CLR_RESEED_BLOCKS;
}

void clr_dual_block(clr_dual *d, uint8_t out[CLR_BLOCK])
{
    if (d->reseed_every && ++d->since_reseed >= d->reseed_every)
        reseed_twin(d);
    clr_dual_ratchet(d);
    for (unsigned b = 0; b < CLR_BLOCK; b++)
        out[b] = d->a.cell[0][b] ^ d->a.cell[0][b + CLR_BLOCK];
}

void clr_dual_fill(clr_dual *d, void *out, size_t len)
{
    uint8_t *p = out;
    while (len) {
        if (d->buf_pos == CLR_BLOCK) {
            clr_dual_block(d, d->buf);
            d->buf_pos = 0;
        }
        size_t n = CLR_BLOCK - d->buf_pos;
        if (n > len)
            n = len;
        memcpy(p, d->buf + d->buf_pos, n);
        wipe(d->buf + d->buf_pos, n);
        d->buf_pos += (uint32_t)n;
        p += n;
        len -= n;
    }
}

uint64_t clr_dual_u64(clr_dual *d)
{
    uint8_t b[8];
    uint64_t v = 0;
    clr_dual_fill(d, b, sizeof b);
    for (int i = 7; i >= 0; i--)
        v = (v << 8) | b[i];
    return v;
}

size_t clr_dual_size(void) { return sizeof(clr_dual); }
