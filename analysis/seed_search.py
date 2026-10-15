#!/usr/bin/env python3
"""
Seed search: the attacker knows the whole algorithm, sees the first 32-byte
block of the victim's output, and tries candidate seeds in a natural order
(decimal numbers as text, then every 1- and 2-byte binary seed) until one
reproduces that block.  With the seed found, the attacker regenerates the
victim's entire stream.

Run against both modes of the dual generator:
  --reproducible  output depends on the seed only: a small seed falls at once
  occluded        (default) the interior twin is seeded and reseeded from OS
                  entropy, so the output is not a function of the seed

Usage: python3 analysis/seed_search.py [SEED] [max_decimal]
"""

import itertools
import os
import subprocess
import sys
import time

import clr as clrmod

SEED = sys.argv[1] if len(sys.argv) > 1 else "1"
MAX_DECIMAL = int(sys.argv[2]) if len(sys.argv) > 2 else 1_000_000
STREAM = 1 << 20
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
lib = clrmod.clr


def candidates():
    for k in range(MAX_DECIMAL):
        yield str(k).encode()
    for length in (1, 2):
        for t in itertools.product(range(256), repeat=length):
            yield bytes(t)


def attack(flags, label):
    victim = subprocess.run([os.path.join(ROOT, "build", "clr"), *flags, "--seed", SEED, "--bytes", str(STREAM)],
                            check=True, capture_output=True).stdout
    observed = victim[:clrmod.BLOCK]
    print(f"{label}: `clr {' '.join(flags + ['--seed', SEED])}`; attacker sees {observed.hex()[:32]}...")
    start = time.perf_counter()
    found, tries = None, 0
    for tries, cand in enumerate(candidates(), 1):
        if lib.dual_blocks(lib.dual_new(cand), 1) == observed:
            found = cand
            break
    elapsed = time.perf_counter() - start
    if found is None:
        true_block = lib.dual_blocks(lib.dual_new(SEED.encode()), 1)
        print(f"  no seed reproduces it: {tries:,} candidates tried in {elapsed:.0f} s, "
              f"including the true seed {SEED!r} (its block: {true_block.hex()[:32]}...)\n")
        return
    replay = lib.dual_blocks(lib.dual_new(found), STREAM // clrmod.BLOCK)
    match = sum(a == b for a, b in zip(replay, victim)) / STREAM
    print(f"  seed found after {tries:,} tries ({elapsed * 1e3:.1f} ms): {found!r}")
    print(f"  regenerated the victim's next {STREAM >> 20} MiB: {match:.2%} of bytes identical\n")


attack(["--reproducible"], "Reproducible mode")
attack([], "Occluded mode (default)")

# what the same search costs for larger seeds in reproducible mode, at this machine's rate
t0 = time.perf_counter()
for k in range(200):
    lib.dual_blocks(lib.dual_new(str(k).encode()), 1)
rate = 200 / (time.perf_counter() - t0)
print(f"Reproducible mode, search rate here: {rate:,.0f} seeds/s on one core")
for label, space in (("any 6-digit number", 10 ** 6), ("any 8-byte seed", 2 ** 64), ("32 bytes from the OS", 2 ** 256)):
    secs = space / rate
    when = f"{secs:,.0f} s" if secs < 3600 else f"{secs / 3.15e7:.1e} years"
    print(f"  {label:<22} {float(space):.1e} seeds -> {when} worst case on one core")
