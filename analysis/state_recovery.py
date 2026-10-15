#!/usr/bin/env python3
"""
White-box state-recovery attack.

The attacker knows the whole algorithm and sees T output blocks.  They model
the generator as an affine map over GF(2) -- output = M x + c, where x is the
complete 2170-bit hidden state (all four layers and both Gold registers) --
build M by running the algorithm on basis states, and solve for x by Gaussian
elimination.  If the model fits, they recover the buried state and predict
every future block.

This is exact for any design whose dynamics are GF(2)-linear (LFSRs, Gold
codes, multiplication by Galois-field constants), however non-local, layered
or self-coupled.  Against nonlinear dynamics the system of equations has no
solution and the "recovered" state is noise.

Usage: python3 analysis/state_recovery.py [observed_blocks] [predicted_blocks]
"""

import sys

import numpy as np

import clr as clrmod

T = int(sys.argv[1]) if len(sys.argv) > 1 else 16
P = int(sys.argv[2]) if len(sys.argv) > 2 else 16
ROUNDS = 8
n = clrmod.STATE_BITS
rng = np.random.default_rng(1234)


def output_bits(lib, x, nblocks):
    s = clrmod.State.from_bits(x, ROUNDS)
    return np.unpackbits(np.frombuffer(lib.blocks(s, nblocks), dtype=np.uint8), bitorder="little")


def to_int(bits):
    return int.from_bytes(np.packbits(bits, bitorder="little").tobytes(), "little")


def attack(lib, label):
    # 1. Model: affine response of every output bit to every hidden state bit.
    zero = np.zeros(n, dtype=np.uint8)
    y0 = output_bits(lib, zero, T + P)
    cols = np.empty((n, len(y0)), dtype=np.uint8)
    for i in range(n):
        e = zero.copy()
        e[i] = 1
        cols[i] = output_bits(lib, e, T + P) ^ y0
    obs_rows = cols[:, : 128 * 2 * T].T                 # equation r: sum_i M[r,i] x_i = y[r] ^ y0[r]

    # 2. Victim: a secret random state; the attacker sees T blocks.
    secret = rng.integers(0, 2, n, dtype=np.uint8)
    y = output_bits(lib, secret, T + P)

    # 3. Gaussian elimination over GF(2) on rows packed as ints (bit n = right-hand side).
    rows = [to_int(r) | (int(b) << n) for r, b in zip(obs_rows, y[: 256 * T] ^ y0[: 256 * T])]
    pivots, inconsistent = {}, 0
    for r in rows:
        for col, pr in pivots.items():
            if (r >> col) & 1:
                r ^= pr
        low = r & ((1 << n) - 1)
        if not low:
            inconsistent += r >> n                      # 0 = 1: the model does not fit the data
            continue
        col = low.bit_length() - 1
        for c2 in list(pivots):
            if (pivots[c2] >> col) & 1:
                pivots[c2] ^= r
        pivots[col] = r
    xhat = np.zeros(n, dtype=np.uint8)
    for col, pr in pivots.items():
        xhat[col] = (pr >> n) & 1                       # free variables set to 0

    # 4. Predict the P blocks the attacker never saw.
    pred = (y0 ^ (xhat @ cols % 2).astype(np.uint8))[256 * T:]
    future_acc = (pred == y[256 * T:]).mean()
    seen_acc = ((y0 ^ (xhat @ cols % 2).astype(np.uint8))[: 256 * T] == y[: 256 * T]).mean()
    layers = [(xhat[k * 512:(k + 1) * 512] == secret[k * 512:(k + 1) * 512]).mean() for k in range(4)]
    gold = (xhat[2048:] == secret[2048:]).mean()
    print(f"\n{label}")
    print(f"  equations {len(rows)}, rank {len(pivots)} of {n}, contradictory equations {inconsistent}")
    print(f"  model fits the {T} observed blocks:   {seen_acc:7.2%} of bits")
    print(f"  predicts the {P} unseen blocks:       {future_acc:7.2%} of bits")
    print(f"  hidden state bits recovered:   L0 {layers[0]:.0%}  L1 {layers[1]:.0%}  "
          f"L2 {layers[2]:.0%}  L3 {layers[3]:.0%}  Gold {gold:.0%}")


print(f"Attacker sees {T} blocks ({256 * T} bits) and knows the full algorithm; state = {n} bits.")
attack(clrmod.linear_control(), "Same architecture, all-linear control (Gold + Galois constants only):")
attack(clrmod.clr, "CLR as built (Galois inversion S-box + RS decoding in the loop):")
