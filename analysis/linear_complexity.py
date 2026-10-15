#!/usr/bin/env python3
"""
Berlekamp-Massey attack: find the shortest linear recurrence that reproduces
the first N output bits, then use it to predict the next bits.

Any generator built only from GF(2)-linear parts (LFSRs, Gold codes, Galois
field multiplication by constants), however many layers or feedback loops it
has, collapses to one recurrence no longer than its state, so BM recovers it
from ~2x(state bits) of output and predicts everything after.  A generator
with no such structure has linear complexity ~N/2 and BM predicts at chance.

Usage: python3 analysis/linear_complexity.py [N]
"""

import os
import subprocess
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
TEST = 4000
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "tools"))
from gen_tables import P61, derive_pair, galois_bits  # noqa: E402


def berlekamp_massey(bits):
    """Linear complexity and connection polynomial (int, bit i = c_i) of a bit list."""
    C, B, L, m, R = 1, 1, 0, -1, 0
    for n, s in enumerate(bits):
        R = (R << 1) | s                        # bit i of R = s[n - i]
        if (C & R).bit_count() & 1:
            T = C
            C ^= B << (n - m)
            if 2 * L <= n:
                L, m, B = n + 1 - L, n, T
    return L, C


def predict(bits, train, C):
    """Fraction of bits[train:] predicted correctly from the true preceding bits."""
    R, hits = 0, 0
    for s in bits[:train]:
        R = (R << 1) | s
    taps = C >> 1
    for s in bits[train:]:
        hits += ((taps & R).bit_count() & 1) == s
        R = (R << 1) | s
    return hits / (len(bits) - train)


def cli_bits(binary, rounds, nblocks, bit=0, single=True):
    """Output bit number `bit` of each 32-byte block.

    Each block is 256 different readouts of the state, so the right linear
    attack follows one bit position from block to block: for a linear design
    that sequence obeys a recurrence no longer than the state.
    """
    args = [os.path.join(ROOT, "build", binary), "--seed", "linear-complexity",
            "--rounds", str(rounds), "--bytes", str(32 * nblocks)] + (["--single"] if single else [])
    out = subprocess.run(args, check=True, capture_output=True).stdout
    return [(out[32 * t + bit // 8] >> (bit % 8)) & 1 for t in range(nblocks)]


def gold_bits(nbits):
    mu, mv, _ = derive_pair(P61, 61)
    u, _ = galois_bits(mu, 61, 0x123456789ABCDEF, nbits)
    v, _ = galois_bits(mv, 61, 0xFEDCBA987654321, nbits)
    return [a ^ b for a, b in zip(u, v)]


CASES = [
    ("embedded Gold pair on its own (122 state bits)", lambda n: gold_bits(n)),
    ("single loop, all-linear control, output bit 0", lambda n: cli_bits("clr_linear_control", 8, n, 0)),
    ("single loop, all-linear control, output bit 137", lambda n: cli_bits("clr_linear_control", 8, n, 137)),
    ("dual, all-linear control, output bit 0", lambda n: cli_bits("clr_linear_control", 8, n, 0, False)),
    ("single loop, 1 pass per block, output bit 0", lambda n: cli_bits("clr", 1, n, 0)),
    ("single loop, 8 passes per block, output bit 0", lambda n: cli_bits("clr", 8, n, 0)),
    ("single loop, 8 passes per block, output bit 137", lambda n: cli_bits("clr", 8, n, 137)),
    ("dual (default), output bit 0", lambda n: cli_bits("clr", 8, n, 0, False)),
    ("dual (default), output bit 137", lambda n: cli_bits("clr", 8, n, 137, False)),
]

print(f"Berlekamp-Massey on the first {N} bits of each sequence, then predict the next {TEST}")
print(f"(state = 2170 bits single, 4340 dual; a random sequence has linear complexity ~{N // 2} and is predicted ~50%)\n")
for name, gen in CASES:
    bits = gen(N + TEST)
    L, C = berlekamp_massey(bits[:N])
    acc = predict(bits, N, C)
    print(f"  {name:<62} L = {L:>5}   predicted {acc:7.2%}")
