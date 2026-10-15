#!/usr/bin/env python3
"""
Derive and verify every constant the generator uses, then emit src/clr_tables.h.

  * Gold preferred pair: a primitive degree-61 polynomial p and the minimal
    polynomial q of its 3-decimation (Gold's theorem: n odd, decimation 2^k+1
    with gcd(n, k) = 1 gives a preferred pair).  The same derivation is run at
    n = 7 and n = 13 and the three-valued cross-correlation is checked
    exhaustively, so the code path that produces the n = 61 pair is verified.
  * GF(2^8) tables (AES field, x^8+x^4+x^3+x+1): exp/log for generator 0x03,
    the S-box (multiplicative inverse followed by the AES affine map).
  * Byte-at-a-time stepping tables for the two Galois-configuration LFSRs.

Usage: python3 tools/gen_tables.py [--check]   (--check: verify only, write nothing)
"""

import os
import sys

# ---------------------------------------------------------------- GF(2)[x]

def pmod(a, m):
    dm = m.bit_length() - 1
    while a.bit_length() - 1 >= dm:
        a ^= m << (a.bit_length() - 1 - dm)
    return a


def pmulmod(a, b, m):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a.bit_length() == m.bit_length():
            a ^= m
    return r


def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a


def is_prime(n):
    """Deterministic Miller-Rabin for n < 3.3e24."""
    if n < 2:
        return False
    bases = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for p in bases:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in bases:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def prime_factors(n):
    """Distinct prime factors via Pollard rho (2^61 - 1 is prime; 2^13 - 1 and 2^7 - 1 too)."""
    import math
    if n == 1:
        return []
    if is_prime(n):
        return [n]
    for c in range(1, 100):
        x = y = 2
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n:
            return sorted(set(prime_factors(d) + prime_factors(n // d)))
    raise ValueError(f"could not factor {n}")


def is_primitive(p):
    """Rabin irreducibility test + order check of x modulo p."""
    n = p.bit_length() - 1
    # irreducible iff x^(2^n) = x mod p and gcd(x^(2^(n/r)) - x, p) = 1 for primes r | n
    def x_pow_2k(k):
        t = 0b10
        for _ in range(k):
            t = pmulmod(t, t, p)
        return t
    if x_pow_2k(n) != 0b10:
        return False
    for r in prime_factors(n):
        if pgcd(p, x_pow_2k(n // r) ^ 0b10) != 1:
            return False
    # primitive iff x has order exactly 2^n - 1
    order = (1 << n) - 1
    for r in prime_factors(order):
        e, acc, base = order // r, 1, 0b10
        while e:
            if e & 1:
                acc = pmulmod(acc, base, p)
            base = pmulmod(base, base, p)
            e >>= 1
        if acc == 1:
            return False
    return True


# ---------------------------------------------------------------- LFSRs

def galois_bits(mask, n, state, count):
    """Right-shift Galois LFSR exactly as implemented in src/clr.c."""
    out = []
    for _ in range(count):
        o = state & 1
        state >>= 1
        if o:
            state ^= mask
        out.append(o)
    return out, state


def berlekamp_massey(bits):
    """Returns (L, C) with C as an int: bit i = coefficient of x^i in the connection polynomial."""
    c, b, L, m = 1, 1, 0, -1
    for N, sN in enumerate(bits):
        d = sN
        cc = c >> 1
        i = 1
        while cc:
            if cc & 1:
                d ^= bits[N - i]
            cc >>= 1
            i += 1
        if d:
            t = c
            c ^= b << (N - m)
            if 2 * L <= N:
                L, m, b = N + 1 - L, N, t
    return L, c


def reciprocal(poly, n):
    return int(format(poly, f"0{n + 1}b")[::-1], 2)


def mask_for_poly(poly, n):
    """Galois mask whose output sequence has characteristic polynomial `poly`.

    For the right-shift Galois LFSR the output satisfies the recurrence whose
    connection polynomial is C(x) = 1 + sum_j mask_j x^(j+1), i.e. C = (mask << 1) | 1.
    """
    assert poly & 1 and poly.bit_length() == n + 1
    return poly >> 1


def derive_pair(p, n):
    """Return (mask_u, mask_v): mask_u realises p, mask_v realises its 3-decimation."""
    assert is_primitive(p), "p is not primitive"
    # p as a connection polynomial (constant term 1); find the matching mask and confirm with BM
    mask_u = mask_for_poly(p, n)
    u, _ = galois_bits(mask_u, n, 1, 6 * n)
    L, conn = berlekamp_massey(u)
    assert L == n and conn == p, "mask convention mismatch"
    # 3-decimation u_{3t} and its minimal (connection) polynomial
    long_u, _ = galois_bits(mask_u, n, 1, 3 * 6 * n)
    v_ref = long_u[::3]
    L, q = berlekamp_massey(v_ref)
    assert L == n and is_primitive(q) and q not in (p, reciprocal(p, n))
    mask_v = mask_for_poly(q, n)
    # v-LFSR must reproduce a shift of u_{3t}: same minimal polynomial, n-bit window match
    v, _ = galois_bits(mask_v, n, 1, 6 * n)
    assert berlekamp_massey(v) == (n, q)
    assert (mask_u >> (n - 1)) & 1 and (mask_v >> (n - 1)) & 1, "top tap needed for invertibility"
    return mask_u, mask_v, q


def verify_gold_crosscorrelation(p, n):
    import numpy as np
    mu, mv, _ = derive_pair(p, n)
    N = (1 << n) - 1
    u = np.array(galois_bits(mu, n, 1, N)[0], dtype=np.int64)
    v = np.array(galois_bits(mv, n, 1, N)[0], dtype=np.int64)
    a, b = 1 - 2 * u, 1 - 2 * v
    # periodic cross-correlation at every shift (N is small enough to do directly)
    vals = set()
    for tau in range(N):
        vals.add(int(np.dot(a, np.roll(b, -tau))))
    t = 1 + (1 << ((n + 1) // 2))
    expected = {-1, -t, t - 2}
    assert vals <= expected and len(vals) == 3, f"n={n}: cross-correlation {sorted(vals)} != {sorted(expected)}"
    return sorted(vals)


# ---------------------------------------------------------------- GF(2^8)

def gmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 0x100:
            a ^= 0x11B
    return r


def gf_tables():
    exp, log = [0] * 510, [0] * 256
    x = 1
    for i in range(255):
        exp[i] = exp[i + 255] = x
        log[x] = i
        x = gmul(x, 0x03)
    assert len(set(exp[:255])) == 255, "0x03 must generate GF(2^8)*"
    inv = [0] + [exp[(255 - log[a]) % 255] for a in range(1, 256)]
    sbox = []
    for a in range(256):
        b = inv[a]
        s = b
        for k in range(1, 5):
            s ^= ((b << k) | (b >> (8 - k))) & 0xFF
        sbox.append(s ^ 0x63)
    assert sbox[0x00] == 0x63 and sbox[0x01] == 0x7C and sbox[0x53] == 0xED  # FIPS-197 values
    return exp, log, sbox


def step8_tables(mask, n):
    T, O = [], []
    for b in range(256):
        bits, s = galois_bits(mask, n, b, 8)
        T.append(s)
        O.append(sum(bit << i for i, bit in enumerate(bits)))
    return T, O


# ---------------------------------------------------------------- main

P61 = (1 << 61) | (1 << 5) | (1 << 2) | (1 << 1) | 1     # x^61 + x^5 + x^2 + x + 1
P7 = (1 << 7) | (1 << 1) | 1                             # x^7 + x + 1
P13 = (1 << 13) | (1 << 4) | (1 << 3) | (1 << 1) | 1     # x^13 + x^4 + x^3 + x + 1


def fmt_poly(p):
    terms = [("x^%d" % i if i > 1 else ("x" if i == 1 else "1")) for i in range(p.bit_length() - 1, -1, -1) if (p >> i) & 1]
    return " + ".join(terms)


def c_array(name, ctype, values, per_line, fmt):
    lines = []
    for i in range(0, len(values), per_line):
        lines.append("    " + ", ".join(fmt % v for v in values[i:i + per_line]) + ",")
    return f"static const {ctype} {name}[{len(values)}] = {{\n" + "\n".join(lines) + "\n};\n"


def main():
    check_only = "--check" in sys.argv
    for p, n in ((P7, 7), (P13, 13)):
        vals = verify_gold_crosscorrelation(p, n)
        print(f"n={n:<3} preferred pair verified, cross-correlation values {vals}")
    mask_u, mask_v, q61 = derive_pair(P61, 61)
    print(f"n=61  p(x) = {fmt_poly(P61)}  (primitive)")
    print(f"      q(x) = {fmt_poly(q61)}  (primitive, 3-decimation of p)")
    print(f"      Galois masks u = 0x{mask_u:016X}  v = 0x{mask_v:016X}")

    exp, log, sbox = gf_tables()
    Tu, Ou = step8_tables(mask_u, 61)
    Tv, Ov = step8_tables(mask_v, 61)
    if check_only:
        return

    out = []
    out.append("/* Generated by tools/gen_tables.py -- do not edit. */\n")
    out.append("#ifndef CLR_TABLES_H\n#define CLR_TABLES_H\n\n#include <stdint.h>\n\n")
    out.append(f"/* Gold preferred pair, degree 61 (2^61 - 1 is prime).\n"
               f" *   u: {fmt_poly(P61)}\n"
               f" *   v: {fmt_poly(q61)}  (3-decimation of u)\n */\n")
    out.append(f"#define CLR_GOLD_MASK_U UINT64_C(0x{mask_u:016X})\n")
    out.append(f"#define CLR_GOLD_MASK_V UINT64_C(0x{mask_v:016X})\n\n")
    out.append("/* GF(2^8) mod x^8+x^4+x^3+x+1, generator alpha = 0x03. */\n")
    out.append(c_array("CLR_GF_EXP", "uint8_t", exp, 16, "0x%02X"))
    out.append(c_array("CLR_GF_LOG", "uint8_t", log, 16, "0x%02X"))
    out.append("/* S-box: GF(2^8) inversion followed by the AES affine map. */\n")
    out.append(c_array("CLR_SBOX", "uint8_t", sbox, 16, "0x%02X"))
    out.append("/* Galois LFSRs stepped 8 bits at a time: s' = (s >> 8) ^ T[s & 0xFF], bits out = O[s & 0xFF]. */\n")
    out.append(c_array("CLR_GOLD_TU", "uint64_t", Tu, 4, "UINT64_C(0x%016X)"))
    out.append(c_array("CLR_GOLD_OU", "uint8_t", Ou, 16, "0x%02X"))
    out.append(c_array("CLR_GOLD_TV", "uint64_t", Tv, 4, "UINT64_C(0x%016X)"))
    out.append(c_array("CLR_GOLD_OV", "uint8_t", Ov, 16, "0x%02X"))
    out.append("\n#endif\n")
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "clr_tables.h")
    with open(path, "w") as f:
        f.write("".join(out))
    print(f"wrote {os.path.normpath(path)}")


if __name__ == "__main__":
    main()
