# closed_loop_feedback_rng

**CLR** is a random number generator that buries itself in its own dynamics.
Gold and Galois-field dynamics run in a closed feedback loop across several
layers. Each layer interacts with itself spatially, and every layer corrects
the next with a Galois code. Only the outer layer is ever seen. The inner
layers and the embedded Gold generators feed the outer layer's inputs and
outputs but never appear in the output.

> Status: experimental design. It passes every attack in `analysis/` (below),
> but it has not been independently cryptanalysed, and passing tests is not a
> proof. Treat it as research until it has been.

## Design

```
  ┌──────────────── Gold feedback (L0 → Gold pair) ────────────────┐
  ▼                                                                │
 Gold pair: u ⊕ v from two 61-bit Galois LFSRs                     │
  │ 16 bytes into every layer, every pass                          │
  ▼                                                                │
 L0 ──RS──▶ L1 ──RS──▶ L2 ──RS──▶ L3 ──RS──▶ back into L0          │
 │                                                                 │
 ├─────────────────────────────────────────────────────────────────┘
 └──▶ output = L0 rows 0–3 ⊕ rows 4–7   (L1–L3 and the Gold pair are never output)
```

| Idea | Mechanism |
|---|---|
| Closed feedback loop | L0 → L1 → L2 → L3 → L0 through Reed–Solomon corrections; the outer layer drives the Gold pair that drives every layer |
| Buried several layers deep | Output is a fold of layer 0 only; layers 1–3 (1536 bits) and the Gold pair (122 bits) are never output |
| Embedded, never-used generators | A Gold preferred pair (degree 61, `x^61+x^5+x^2+x+1` and its 3-decimation) supplies 16 bytes to every layer on every pass |
| Spatial self-interaction | Each layer is an 8×8 torus of GF(2^8) cells, updated in two checkerboard half-steps: `cell ^= S(N ⊕ 2E ⊕ 4S ⊕ 8W ⊕ position ⊕ Gold)`. `S` is GF(2^8) inversion (the AES S-box) |
| Galois-code self-correction | Each column of a layer is decoded as a length-8 Reed–Solomon word; the single-error correction (location `S2/S1`, magnitude `S1²/S2`) is applied to the same column of the next layer. Where it lands depends on the data |
| Re-randomisation | 8 full passes of the loop between output blocks (full burial takes 3) |
| No entropy lost to cycling | Every step is invertible: half-steps are involutions, corrections and feedback are Feistel steps between different layers. `clr_unround` undoes `clr_round` exactly, so the loop never merges states however long it runs |

Total hidden state: 4 × 512 + 122 = 2170 bits.

## What an attacker gets

Each attack runs against CLR and against a **linear control**: the same
architecture, with the S-box replaced by multiplication by a constant and the RS
decode replaced by raw syndrome injection (`make` builds it as
`build/clr_linear_control`). The control isolates what actually defeats pattern
matching.

**Burial depth** (`analysis/burial.py`). Flip one bit anywhere and follow the
difference. Ideal is 0.500 everywhere.

| Passes | Flipped in deepest layer → output | Flipped in Gold generator → output | Flipped in outer layer → output |
|---|---|---|---|
| 1 | 0.109 | 0.478 | 0.290 |
| 2 | 0.498 | 0.500 | 0.501 |
| 3 | 0.501 | 0.500 | 0.500 |

From pass 3 on, every layer and every output bit flip with probability 0.5.
The worst of the 256 output bits stays within the range a fair coin gives.

**Linear-recurrence attack** (`analysis/linear_complexity.py`). Berlekamp–Massey
on 12,000 bits, then predict the next 4,000.

| Sequence | Linear complexity | Predicted |
|---|---|---|
| Embedded Gold pair on its own | 122 | 100% |
| Linear control, per output bit | 1147 | 100% |
| CLR, 1 pass per block | 5999 (random: ~6000) | 48.5% |
| CLR, 8 passes per block | 6000 | 50.3% |

**White-box state recovery** (`analysis/state_recovery.py`). The attacker knows
the whole algorithm, sees 16 blocks, and solves for the entire 2170-bit hidden
state over GF(2).

| | Hidden state recovered (L0 / L1 / L2 / L3 / Gold) | Unseen output predicted |
|---|---|---|
| Linear control | 100% / 100% / 100% / 100% / 100% | 100% |
| CLR | 45% / 47% / 53% / 48% / 47% | 49.9% |

(For CLR, 957 of the 4096 equations contradict each other: no linear model
fits it. The "fits 77% of observed bits" line the script prints is what fitting
2170 free unknowns to 4096 bits gives by construction.)

**PractRand 0.95** (`./build/clr | RNG_test stdin`):

| Build | Result |
|---|---|
| Linear control | **FAIL** at 16 MB (`BRank`, p ≈ 1e-377) |
| CLR, 1 pass per block | no failures through 4 GB (run continuing) |
| CLR, 8 passes per block (default) | no failures through 512 MB (run continuing) |

("Unusual" flags at isolated lengths that never escalate are expected from
any good generator; only escalating flags or FAIL indicate structure.)

**The conclusion.** Spatial non-locality, layering and closed-loop feedback
hide nothing on their own: the linear control has all three and is fully
recovered. What makes the buried state unrecoverable is the Galois-field
nonlinearity (S-box inversion and RS division) inside that structure. With it,
the same attacks return noise, even after a single pass.

## Limits that apply to any deterministic generator

- **The seed is the ceiling.** Mixing never adds entropy, so the output can be
  no less guessable than the seed. Seed with at least 32 bytes from the OS (the
  CLI does this when `--seed` is omitted).
- **No forward secrecy.** The loop is invertible by design, so anyone who
  captures the complete internal state can run it backwards and regenerate past
  output. `clr_absorb` can mix in fresh entropy; a one-way ratchet would be
  needed to protect the past.
- **Not constant-time.** S-box lookups and data-dependent correction locations
  can leak through cache timing on shared hardware.
- **Speed.** About 2.5 MiB/s at 8 passes (portable C, unvectorised).

## Build and run

```sh
make                         # build/clr, build/test_clr, build/libclr.so, linear control
make test                    # unit tests + re-verification of every generated constant
./build/clr --seed abc --hex 64
./build/clr | RNG_test stdin # PractRand, if installed
python3 analysis/burial.py
python3 analysis/linear_complexity.py
python3 analysis/state_recovery.py
```

The C API is in `src/clr.h`: `clr_init`, `clr_fill`, `clr_u64`, `clr_absorb`,
`clr_round`/`clr_unround`.

## Layout

| Path | Contents |
|---|---|
| `src/clr.c`, `src/clr.h` | the generator |
| `src/clr_tables.h` | generated constants (Gold masks, GF(2^8) tables, S-box) |
| `src/clr_cli.c` | command-line streamer |
| `tools/gen_tables.py` | derives the constants. It proves both Gold polynomials primitive and checks the preferred-pair construction's three-valued cross-correlation exhaustively at n = 7 and 13 |
| `tests/test_clr.c` | tables, Gold stepping, word-parallel vs reference steps, exact invertibility over 10,000 random states, seed avalanche, known-answer vector |
| `analysis/` | burial, Berlekamp–Massey and state-recovery attacks (Python, via ctypes) |
