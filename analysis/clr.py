"""ctypes binding to build/libclr.so and the linear control build (run `make` first)."""

import ctypes
import os

import numpy as np

LAYERS, CELLS, BLOCK = 4, 64, 32
STATE_BITS = LAYERS * CELLS * 8 + 2 * 61
M61 = (1 << 61) - 1
_BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build")


class State(ctypes.Structure):
    _fields_ = [
        ("cell", (ctypes.c_uint8 * CELLS) * LAYERS),
        ("gold_u", ctypes.c_uint64),
        ("gold_v", ctypes.c_uint64),
        ("rounds", ctypes.c_uint32),
        ("buf_pos", ctypes.c_uint32),
        ("buf", ctypes.c_uint8 * BLOCK),
    ]

    def copy(self):
        s = State()
        ctypes.memmove(ctypes.byref(s), ctypes.byref(self), ctypes.sizeof(State))
        return s

    def layers(self):
        """Cells as a (LAYERS, CELLS) uint8 array."""
        return np.frombuffer(bytes(self.cell), dtype=np.uint8).reshape(LAYERS, CELLS)

    def gold_bytes(self):
        return np.frombuffer(self.gold_u.to_bytes(8, "little") + self.gold_v.to_bytes(8, "little"), dtype=np.uint8)

    def fold(self):
        """What the next output block would be if it were taken now."""
        c = self.layers()[0]
        return c[:BLOCK] ^ c[BLOCK:]

    def bits(self):
        """All 2170 state bits: cells (LSB first), then gold_u, gold_v."""
        cells = np.unpackbits(self.layers().reshape(-1), bitorder="little")
        gu = np.array([(self.gold_u >> i) & 1 for i in range(61)], dtype=np.uint8)
        gv = np.array([(self.gold_v >> i) & 1 for i in range(61)], dtype=np.uint8)
        return np.concatenate([cells, gu, gv])

    @staticmethod
    def from_bits(x, rounds):
        s = State()
        x = np.asarray(x, dtype=np.uint8)
        ctypes.memmove(s.cell, np.packbits(x[:2048], bitorder="little").tobytes(), LAYERS * CELLS)
        s.gold_u = sum(int(b) << i for i, b in enumerate(x[2048:2109]))
        s.gold_v = sum(int(b) << i for i, b in enumerate(x[2109:2170]))
        s.rounds = rounds
        s.buf_pos = BLOCK
        return s


class DualState(ctypes.Structure):
    """clr_dual: exterior generator a, interior twin b."""
    _fields_ = [
        ("a", State),
        ("b", State),
        ("buf_pos", ctypes.c_uint32),
        ("buf", ctypes.c_uint8 * BLOCK),
    ]

    def copy(self):
        s = DualState()
        ctypes.memmove(ctypes.byref(s), ctypes.byref(self), ctypes.sizeof(DualState))
        return s

    def bits(self):
        return np.concatenate([self.a.bits(), self.b.bits()])

    @staticmethod
    def from_bits(x, rounds):
        s = DualState()
        s.a = State.from_bits(x[:STATE_BITS], rounds)
        s.b = State.from_bits(x[STATE_BITS:], rounds)
        s.buf_pos = BLOCK
        return s


class Lib:
    def __init__(self, name):
        self.lib = ctypes.CDLL(os.path.join(_BUILD, name))
        S, D = ctypes.POINTER(State), ctypes.POINTER(DualState)
        for fn, args in (
            ("clr_init", [S, ctypes.c_char_p, ctypes.c_size_t, ctypes.c_uint32]),
            ("clr_round", [S]),
            ("clr_unround", [S]),
            ("clr_block", [S, ctypes.c_char_p]),
            ("clr_fill", [S, ctypes.c_char_p, ctypes.c_size_t]),
            ("clr_dual_init", [D, ctypes.c_char_p, ctypes.c_size_t, ctypes.c_uint32]),
            ("clr_dual_block", [D, ctypes.c_char_p]),
            ("clr_dual_permute", [D]),
            ("clr_dual_unpermute", [D]),
            ("clr_dual_ratchet", [D]),
        ):
            getattr(self.lib, fn).argtypes = args
            getattr(self.lib, fn).restype = None
        self.lib.clr_state_size.restype = ctypes.c_size_t
        self.lib.clr_dual_size.restype = ctypes.c_size_t
        assert self.lib.clr_state_size() == ctypes.sizeof(State), "State layout does not match clr.h"
        assert self.lib.clr_dual_size() == ctypes.sizeof(DualState), "DualState layout does not match clr.h"

    def dual_new(self, seed: bytes, rounds: int = 0) -> DualState:
        d = DualState()
        self.lib.clr_dual_init(ctypes.byref(d), seed, len(seed), rounds)
        return d

    def dual_blocks(self, d: DualState, n: int) -> bytes:
        buf = ctypes.create_string_buffer(BLOCK)
        out = []
        for _ in range(n):
            self.lib.clr_dual_block(ctypes.byref(d), buf)
            out.append(buf.raw)
        return b"".join(out)

    def dual_permute(self, d: DualState):
        self.lib.clr_dual_permute(ctypes.byref(d))

    def dual_unpermute(self, d: DualState):
        self.lib.clr_dual_unpermute(ctypes.byref(d))

    def dual_ratchet(self, d: DualState):
        self.lib.clr_dual_ratchet(ctypes.byref(d))

    def new(self, seed: bytes, rounds: int = 0) -> State:
        s = State()
        self.lib.clr_init(ctypes.byref(s), seed, len(seed), rounds)
        return s

    def step(self, s: State):
        self.lib.clr_round(ctypes.byref(s))

    def unstep(self, s: State):
        self.lib.clr_unround(ctypes.byref(s))

    def blocks(self, s: State, n: int) -> bytes:
        buf = ctypes.create_string_buffer(BLOCK)
        out = []
        for _ in range(n):
            self.lib.clr_block(ctypes.byref(s), buf)
            out.append(buf.raw)
        return b"".join(out)

    def fill(self, s: State, n: int) -> bytes:
        buf = ctypes.create_string_buffer(n)
        self.lib.clr_fill(ctypes.byref(s), buf, n)
        return buf.raw


clr = Lib("libclr.so")


def linear_control() -> Lib:
    return Lib("libclr_linear_control.so")


def random_state(rng: np.random.Generator, rounds: int = 1) -> State:
    return State.from_bits(rng.integers(0, 2, STATE_BITS, dtype=np.uint8), rounds)


def fold_bits(state_bits):
    """Output block (as bits) that a state's outer layer folds to."""
    cells = np.packbits(state_bits[:512], bitorder="little")
    return np.unpackbits(cells[:BLOCK] ^ cells[BLOCK:], bitorder="little")


def gf2_solve(M, y, kernel=False):
    """Solve M x = y over GF(2) (M: m x n of 0/1).  Free variables are set to 0.

    Returns (x, rank, contradictions): contradictions counts equations that
    reduce to 0 = 1, i.e. that no x satisfies.  With kernel=True, also returns
    a basis of the null space (every solution is x xor a combination of it).
    """
    m, n = M.shape
    aug = np.packbits(np.concatenate([M, y.reshape(-1, 1)], axis=1).astype(np.uint8), axis=1, bitorder="little")
    used = np.zeros(m, dtype=bool)
    pivots = []
    for col in range(n):
        colbits = (aug[:, col >> 3] >> (col & 7)) & 1 == 1
        cand = np.flatnonzero(colbits & ~used)
        if cand.size == 0:
            continue
        p = cand[0]
        used[p] = True
        pivots.append((p, col))
        colbits[p] = False
        aug[colbits] ^= aug[p]
    rhs = (aug[:, n >> 3] >> (n & 7)) & 1
    x = np.zeros(n, dtype=np.uint8)
    for p, col in pivots:
        x[col] = rhs[p]
    result = (x, len(pivots), int(rhs[~used].sum()))
    if not kernel:
        return result
    pivot_cols = {col for _, col in pivots}
    basis = []
    for f in (c for c in range(n) if c not in pivot_cols):
        v = np.zeros(n, dtype=np.uint8)
        v[f] = 1
        for p, col in pivots:
            v[col] = (aug[p, f >> 3] >> (f & 7)) & 1
        basis.append(v)
    return result + (basis,)


# module-level shortcuts for the real generator
new, step, unstep, fill = clr.new, clr.step, clr.unstep, clr.fill
