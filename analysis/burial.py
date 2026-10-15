#!/usr/bin/env python3
"""
Burial depth: flip one bit somewhere inside the system and follow the
difference through every layer, pass by pass.  A fully buried difference
leaves every layer (and the output) differing in ~50% of its bits, with no
output bit flipping more or less often than a fair coin.

Usage: python3 analysis/burial.py [trials] [max_rounds]
"""

import sys

import numpy as np

import clr

trials = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
max_rounds = int(sys.argv[2]) if len(sys.argv) > 2 else 10
rng = np.random.default_rng(2026)

ORIGINS = {
    "deepest layer (L3)": lambda s, b: s.cell[3].__setitem__(b % 64, s.cell[3][b % 64] ^ (1 << (b // 64 % 8))),
    "embedded Gold gen.": lambda s, b: setattr(s, "gold_u", s.gold_u ^ (1 << (b % 61))),
    "outer layer (L0)": lambda s, b: s.cell[0].__setitem__(b % 64, s.cell[0][b % 64] ^ (1 << (b // 64 % 8))),
}


def bits(a):
    return np.unpackbits(np.asarray(a, dtype=np.uint8))


cols = ["L0", "L1", "L2", "L3", "Gold", "output"]
for origin, flip in ORIGINS.items():
    frac = np.zeros((max_rounds, len(cols)))
    out_flips = np.zeros((max_rounds, 8 * clr.BLOCK))
    for t in range(trials):
        a = clr.random_state(rng)
        b = a.copy()
        flip(b, int(rng.integers(0, 1 << 30)))
        for r in range(max_rounds):
            clr.step(a)
            clr.step(b)
            la, lb = a.layers(), b.layers()
            for k in range(4):
                frac[r, k] += bits(la[k] ^ lb[k]).mean()
            frac[r, 4] += bits(a.gold_bytes() ^ b.gold_bytes()).mean() * 128 / 122
            d = bits(a.fold() ^ b.fold())
            frac[r, 5] += d.mean()
            out_flips[r] += d
    frac /= trials
    p = out_flips / trials
    sigma = 0.5 / np.sqrt(trials)
    print(f"\n1-bit difference injected into the {origin}  ({trials} trials; ideal = 0.500)")
    print("passes " + "".join(f"{c:>8}" for c in cols) + "   worst output bit  (z)")
    for r in range(max_rounds):
        worst = np.abs(p[r] - 0.5).max()
        print(f"{r + 1:>6} " + "".join(f"{v:>8.3f}" for v in frac[r]) + f"   {0.5 + worst if p[r][np.abs(p[r]-0.5).argmax()] > .5 else 0.5 - worst:>10.3f} ({worst / sigma:4.1f})")
print(f"\nworst-bit z: max over {8 * clr.BLOCK} output bits of |p - 0.5| / sigma; "
      f"for a fair coin this is typically 2.5-3.8 at {trials} trials.")
