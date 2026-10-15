#!/usr/bin/env python3
"""
Backtracking attack: the attacker captures the complete internal state right
after the latest block was delivered and tries to regenerate blocks that were
delivered before it.

  1. Single CLR loop: every step is invertible, so running the loop backwards
     regenerates all past output.
  2. Dual generator, naive: run the invertible part P backwards, ignoring the
     closing feed-forward XOR.
  3. Dual generator, white-box: the attacker knows the algorithm and solves
     P(x) xor x = captured for the previous state x over GF(2), linearising P
     from basis states (exact when the dynamics are linear), and enumerates
     every candidate the solution space allows.

Usage: python3 analysis/backtrack.py [past_blocks]
"""

import itertools
import sys

import numpy as np

import clr as clrmod

PAST = int(sys.argv[1]) if len(sys.argv) > 1 else 64
SEED = b"backtracking target"
ROUNDS = 8
N = 2 * clrmod.STATE_BITS
MAX_KERNEL = 12                                     # enumerate up to 2^12 candidate states


def blocks_to_bits(raw):
    return np.unpackbits(np.frombuffer(raw, dtype=np.uint8), bitorder="little").reshape(-1, 256)


def single_loop():
    lib = clrmod.clr
    s = lib.new(SEED, ROUNDS)
    past = blocks_to_bits(lib.blocks(s, PAST + 1))  # victim delivers PAST + 1 blocks
    captured = s.copy()                             # attacker grabs the state after the last one
    hits = []
    for t in range(PAST - 1, -1, -1):               # regenerate block t from the captured state
        for _ in range(ROUNDS):
            lib.unstep(captured)
        hits.append((np.unpackbits(captured.fold(), bitorder="little") == past[t]).mean())
    return np.array(hits)


def dual_attacks(lib):
    d = lib.dual_new(SEED, ROUNDS)
    lib.dual_blocks(d, 3)
    prev = d.bits()                                 # state that produced the previous block (secret)
    prev_block = clrmod.fold_bits(prev)
    last_block = blocks_to_bits(lib.dual_blocks(d, 1))[0]
    captured = d.copy()
    assert (clrmod.fold_bits(captured.bits()) == last_block).all()

    # naive: run P backwards
    naive = captured.copy()
    lib.dual_unpermute(naive)
    naive_match = (clrmod.fold_bits(naive.bits()) == prev_block).mean()

    # white-box: solve R(x) = P(x) xor x = captured, with R linearised around 0
    def ratchet_bits(x):
        s = clrmod.DualState.from_bits(x, ROUNDS)
        lib.dual_ratchet(s)
        return s.bits()

    zero = np.zeros(N, dtype=np.uint8)
    r0 = ratchet_bits(zero)
    M = np.empty((N, N), dtype=np.uint8)
    for i in range(N):
        e = zero.copy()
        e[i] = 1
        M[:, i] = ratchet_bits(e) ^ r0
    x, rank, contra, basis = clrmod.gf2_solve(M, captured.bits() ^ r0, kernel=True)

    # every candidate previous state consistent with the linear model
    found, tried, best = False, 0, 0.0
    if len(basis) <= MAX_KERNEL:
        for combo in itertools.product((0, 1), repeat=len(basis)):
            cand = x.copy()
            for c, v in zip(combo, basis):
                if c:
                    cand ^= v
            tried += 1
            found |= bool((cand == prev).all())
            best = max(best, (clrmod.fold_bits(cand) == prev_block).mean())
    return naive_match, rank, contra, tried, found, best


print("Attacker captures the complete internal state after the latest block and tries")
print("to regenerate output that was delivered before it (50% = no better than guessing).\n")

h = single_loop()
print("1. Single CLR loop, run backwards:")
print(f"   regenerated {int((h == 1).sum())}/{len(h)} past blocks exactly\n")

for label, lib in (("Dual generator, all-linear control", clrmod.linear_control()),
                   ("Dual generator (CLR)", clrmod.clr)):
    naive, rank, contra, tried, found, best = dual_attacks(lib)
    print(f"{label}:")
    print(f"   2. run P backwards, ignoring the feed-forward: previous block {naive:.1%} of bits")
    print(f"   3. solve P(x) xor x = captured over GF(2): rank {rank}/{N}, contradictions {contra}")
    print(f"      {tried} candidate previous states; true state among them: {'YES' if found else 'no'}; "
          f"best candidate matches the previous block in {best:.1%} of bits\n")
