/*
 * CLR -- closed-loop Gold/Galois lattice RNG.
 *
 * State: CLR_LAYERS layers, each an 8x8 torus of GF(2^8) cells, plus a Gold
 * preferred pair of 61-bit Galois LFSRs.  One round closes the loop:
 *
 *   for each layer k:
 *     Gold pair -> 16 bytes of Gold sequence (embedded generator, never output)
 *     layer k   -> two checkerboard half-steps: every cell is rewritten from a
 *                  GF(2^8)-weighted sum of its 4 neighbours through the S-box
 *                  (GF(2^8) inversion), so the layer self-interacts spatially
 *     layer k   -> Reed-Solomon single-error decode of each column; the decoded
 *                  correction (location, magnitude) is applied to layer k+1,
 *                  and the deepest layer's correction is applied to layer 0
 *   layer 0     -> fed back into the Gold pair (the outer layer drives the
 *                  generators that drive it)
 *
 * Every step is invertible (clr_unround undoes clr_round exactly), so the
 * loop never loses entropy no matter how many times it cycles.  Output is
 * taken from layer 0 only, folded (rows 0-3 XOR rows 4-7); layers 1..3 and
 * the Gold pair are never exposed.
 *
 * Status: experimental design, not independently cryptanalysed.
 */
#ifndef CLR_H
#define CLR_H

#include <stddef.h>
#include <stdint.h>

#define CLR_LAYERS 4                    /* layer 0 is the only one ever observed; must be >= 2 */
#define CLR_DIM 8                       /* each layer is a CLR_DIM x CLR_DIM torus */
#define CLR_CELLS (CLR_DIM * CLR_DIM)
#define CLR_BLOCK 32                    /* output bytes per block */
#define CLR_DEFAULT_ROUNDS 8            /* loop passes between output blocks */
#define CLR_SEED_ROUNDS 16              /* loop passes after absorbing a seed */

typedef struct clr_state {
    uint8_t  cell[CLR_LAYERS][CLR_CELLS];
    uint64_t gold_u, gold_v;            /* 61-bit Galois LFSRs, Gold preferred pair */
    uint32_t rounds;                    /* passes of the loop per output block */
    uint32_t buf_pos;
    uint8_t  buf[CLR_BLOCK];
} clr_state;

/* Initialise from a seed of any length; rounds == 0 selects CLR_DEFAULT_ROUNDS. */
void clr_init(clr_state *st, const void *seed, size_t len, uint32_t rounds);

/* Optional: mix external entropy into the outer layer and re-bury it. */
void clr_absorb(clr_state *st, const void *data, size_t len);

/* One pass of the closed loop, and its exact inverse. */
void clr_round(clr_state *st);
void clr_unround(clr_state *st);

/* Output: one raw block, or any number of bytes (buffered). */
void clr_block(clr_state *st, uint8_t out[CLR_BLOCK]);
void clr_fill(clr_state *st, void *out, size_t len);
uint64_t clr_u64(clr_state *st);

size_t clr_state_size(void);

#endif
