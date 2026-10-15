#!/usr/bin/env python3
"""
Seed search: the attack that works when the seed is small, whatever the mixing.

The victim runs `build/clr --seed SEED`.  The attacker sees only the first
32-byte block, then tries candidate seeds in a natural order (decimal numbers
as text, then every 1- and 2-byte binary seed) until one reproduces that
block.  With the seed found, the attacker regenerates the victim's entire
output stream.

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

victim = subprocess.run([os.path.join(ROOT, "build", "clr"), "--seed", SEED, "--bytes", str(STREAM)],
                        check=True, capture_output=True).stdout
observed = victim[:clrmod.BLOCK]
print(f"Victim: `clr --seed {SEED}` (dual generator). Attacker sees only the first block:")
print(f"  {observed.hex()}\n")


def candidates():
    for k in range(MAX_DECIMAL):
        yield str(k).encode()
    for length in (1, 2):
        for t in itertools.product(range(256), repeat=length):
            yield bytes(t)


lib = clrmod.clr
start = time.perf_counter()
found = None
for tries, cand in enumerate(candidates(), 1):
    if lib.dual_blocks(lib.dual_new(cand), 1) == observed:
        found = cand
        break
elapsed = time.perf_counter() - start

if found is None:
    print(f"no seed found in {tries} tries ({elapsed:.1f} s)")
    sys.exit(1)

replay = lib.dual_blocks(lib.dual_new(found), STREAM // clrmod.BLOCK)
match = sum(a == b for a, b in zip(replay, victim)) / STREAM
print(f"Seed found after {tries} tries ({elapsed * 1e3:.1f} ms): {found!r}")
print(f"Regenerated the victim's next {STREAM >> 20} MiB: {match:.2%} of bytes identical\n")

# what the same search costs for larger seeds, at this machine's measured rate
t0 = time.perf_counter()
for k in range(200):
    lib.dual_blocks(lib.dual_new(str(k).encode()), 1)
rate = 200 / (time.perf_counter() - t0)
print(f"Search rate here: {rate:,.0f} seeds/s on one core (unoptimised Python + C)")
for label, space in (("any 6-digit number", 10 ** 6), ("any 8-byte seed", 2 ** 64), ("32 bytes from the OS", 2 ** 256)):
    secs = space / rate
    if secs < 3600:
        when = f"{secs:,.0f} s"
    else:
        years = secs / 3.15e7
        when = f"{years:.1e} years" if years > 1 else f"{secs / 3600:,.1f} h"
    print(f"  {label:<22} {float(space):.1e} seeds -> {when} worst case on one core")
