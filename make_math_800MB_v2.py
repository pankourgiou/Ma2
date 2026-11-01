#!/usr/bin/env python3
"""~800 MB HTML of number theory / trigonometry / Fibonacci / combinatorics.
Standard library only. EVERY row is verified by 200 separate checks; the script
aborts on the first failure.

Usage:  python make_math_800MB_v2.py          (800 MB)
        python make_math_800MB_v2.py 30       (custom size in MB, for a quick test)
"""
import sys, time, math
from math import gcd, lcm, isqrt, comb, perm, factorial, isclose
from array import array
from decimal import Decimal, getcontext

SIZE_MB = int(sys.argv[1]) if len(sys.argv) > 1 else 800
TARGET = SIZE_MB * 1_000_000
OUT = f"math_v2_{SIZE_MB}MB.html"
getcontext().prec = 40
MODS = range(3, 28)                      # 25 moduli
Y_VALUES = [0.1 * k for k in range(1, 11)]
MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def build_spf(limit):
    spf = array("I", range(limit + 1))
    for i in range(2, isqrt(limit) + 1):
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def is_prime_mr(n):
    if n < 2:
        return False
    for p in MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in MR_BASES:
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


def fib_pair(k):                          # fast doubling -> (F(k), F(k+1))
    if k == 0:
        return 0, 1
    a, b = fib_pair(k >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if k & 1 else (c, d)


def fib_matrix(k):                        # independent second method
    def mul(x, y):
        return (x[0]*y[0] + x[1]*y[2], x[0]*y[1] + x[1]*y[3],
                x[2]*y[0] + x[3]*y[2], x[2]*y[1] + x[3]*y[3])
    r, m = (1, 0, 0, 1), (1, 1, 1, 0)
    while k:
        if k & 1:
            r = mul(r, m)
        m = mul(m, m); k >>= 1
    return r[1]


def make_row_factory(spf, limit):
    xs = 1.2 / limit

    def factor(m):
        out = []
        while m > 1:
            p = spf[m]; e = 0
            while m % p == 0:
                m //= p; e += 1
            out.append((p, e))
        return out

    def phi_of(m):
        r = m
        for p, _ in factor(m):
            r = r // p * (p - 1)
        return r

    def mu_of(m):
        f = factor(m)
        return 0 if any(e > 1 for _, e in f) else (-1) ** len(f)

    def row(n):
        fs = factor(n)
        prime = spf[n] == n
        phi = n
        tau = 1
        sigma = 1
        for p, e in fs:
            phi = phi // p * (p - 1)
            tau *= e + 1
            sigma *= (p ** (e + 1) - 1) // (p - 1)
        mu = 0 if any(e > 1 for _, e in fs) else (-1) ** len(fs)
        divs = [1]
        for p, e in fs:
            divs = [d * p ** k for d in divs for k in range(e + 1)]
        base = next(a for a in (2, 3, 5, 7, 11, 13, 17, 19, 23) if gcd(a, n) == 1) if n > 2 else 1

        x = 0.1 + n * xs
        s, c = math.sin(x), math.cos(x)
        k = n % 1000 + 1
        m2 = n % 700 + 1
        fk, fk1 = fib_pair(k)
        fkm1 = fk1 - fk
        small = n % 20

        # ---- 20 number-theory checks ----
        A = [
            math.prod(p ** e for p, e in fs) == n,
            all(is_prime_mr(p) for p, _ in fs),
            prime == is_prime_mr(n),
            phi * math.prod(p for p, _ in fs) == n * math.prod(p - 1 for p, _ in fs),
            sum(phi_of(d) for d in divs) == n,
            tau == len(divs),
            sigma == sum(divs),
            all(n % d == 0 for d in divs),
            len(set(divs)) == len(divs),
            (mu != 0) == all(e == 1 for _, e in fs),
            sum(mu_of(d) for d in divs) == 0,
            phi % 2 == 0 or n == 2,
            (phi == n - 1) == prime,
            sigma >= n + 1,
            pow(base, phi, n) == 1 % n,
            lcm(n, n + 1) == n * (n + 1),
            gcd(n, n + 1) == 1,
            isqrt(n) ** 2 <= n < (isqrt(n) + 1) ** 2,
            (spf[n] == 2) == (n % 2 == 0),
            sigma == math.prod(sum(p ** j for j in range(e + 1)) for p, e in fs),
        ]
        # ---- 100 modular checks (4 forms x 25 moduli) ----
        B = []
        for m in MODS:
            B.append(n * n % m == (n % m) * (n % m) % m)
            B.append(n ** 3 % m == pow(n % m, 3, m))
            t = n % (2 * m)
            B.append(n * (n + 1) // 2 % m == t * (t + 1) // 2 % m)
            B.append(pow(n, 5, m) == pow(n % m, 5, m))
        # ---- 40 floating-point / trigonometry checks ----
        ic = lambda a, b: isclose(a, b, rel_tol=1e-9, abs_tol=1e-9)
        C = [
            ic(s * s + c * c, 1.0),
            ic(math.sin(2 * x), 2 * s * c),
            ic(math.cos(2 * x), c * c - s * s),
            ic(math.cos(2 * x), 1 - 2 * s * s),
            ic(math.cos(2 * x), 2 * c * c - 1),
            ic(math.exp(math.log(n)), n),
            ic(math.log(n * n), 2 * math.log(n)),
            ic(math.sqrt(n) ** 2, n),
            ic(math.hypot(s, c), 1.0),
            ic(math.cosh(x) ** 2 - math.sinh(x) ** 2, 1.0),
            ic(math.tanh(x), math.sinh(x) / math.cosh(x)),
            ic(math.tan(x), s / c),
            ic(math.atan(math.tan(x)), x),
            ic(math.asin(s), x),
            ic(math.acos(c), x),
            ic(math.log2(n), math.log(n) / math.log(2)),
            ic(math.log10(n), math.log(n) / math.log(10)),
            ic(math.expm1(x), math.exp(x) - 1),
            ic(math.log1p(x), math.log(1 + x)),
            ic(math.degrees(math.radians(x)), x),
        ]
        C += [ic(math.sin(x + y), s * math.cos(y) + c * math.sin(y)) for y in Y_VALUES]
        C += [ic(math.cos(x + y), c * math.cos(y) - s * math.sin(y)) for y in Y_VALUES]
        # ---- 40 Fibonacci / combinatorics checks ----
        L = fkm1 + fk1
        D = [
            fib_pair(k + 1)[1] == fk1 + fk,
            fk1 * fkm1 - fk * fk == (-1) ** k,
            fib_pair(2 * k)[0] == fk * (2 * fk1 - fk),
            fib_pair(2 * k + 1)[0] == fk * fk + fk1 * fk1,
            gcd(fk, fib_pair(m2)[0]) == fib_pair(gcd(k, m2))[0],
            (fk % 2 == 0) == (k % 3 == 0),
            fib_pair(5 * k)[0] % 5 == 0,
            fib_pair(3 * k)[0] % 2 == 0,
            L * L - 5 * fk * fk == 4 * (-1) ** k,
            fib_matrix(k) == fk,
        ]
        D += [comb(n, j) == comb(n - 1, j - 1) + comb(n - 1, j) for j in range(1, 11)]
        D += [comb(n, j) == comb(n, n - j) for j in range(1, 11)]
        D += [comb(n, j) * (n - j) == comb(n, j + 1) * (j + 1) for j in range(1, 6)]
        D += [
            comb(n, 2) == n * (n - 1) // 2,
            comb(n, 3) * 6 == n * (n - 1) * (n - 2),
            perm(n, 2) == n * (n - 1),
            factorial(small) * comb(20, small) * factorial(20 - small) == factorial(20),
            comb(n, 1) == n,
        ]
        checks = A + B + C + D
        if len(checks) != 200:
            raise SystemExit(f"internal error: {len(checks)} checks instead of 200")
        if not all(checks):
            bad = [i + 1 for i, ok in enumerate(checks) if not ok]
            raise SystemExit(f"CHECK(S) {bad} FAILED for n={n}")

        fstr = "&middot;".join(f"{p}^{e}" if e > 1 else str(p) for p, e in fs)
        return (f"<tr><td>{n}<td>{fstr}<td>{phi}<td>{tau}<td>{sigma}<td>{mu}<td>{int(prime)}"
                f"<td>{Decimal(n).sqrt()}<td>{s!r}<td>{c!r}<td>{math.log(n)!r}"
                f"<td>{comb(n, 10)}<td>{k}<td>{fk % 1_000_000_007}\n")

    return row


def main():
    # rough first guess of the row count, then benchmark with a real sieve
    guess = int(TARGET / 330)
    limit = guess + 1000
    t0 = time.perf_counter()
    spf = build_spf(limit)
    sieve_t = time.perf_counter() - t0
    row = make_row_factory(spf, limit)
    t0 = time.perf_counter()
    mid = max(12, guess // 2)          # sample mid-table rows: they are longer than the first rows
    sample = "".join(row(n) for n in range(mid, mid + 2000))
    per_row_t = (time.perf_counter() - t0) / 2000
    per_row_b = len(sample.encode()) / 2000
    est_rows = int(TARGET / per_row_b)
    if est_rows + 1000 > limit:           # rebuild sieve if the guess was too small
        limit = int(est_rows * 1.1) + 1000
        spf = build_spf(limit); row = make_row_factory(spf, limit)
    est_min = (per_row_t * est_rows + sieve_t) / 60
    print(f"Estimated rows: {est_rows:,} | estimated time: {est_min:.1f} min", flush=True)

    head = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>Math v2 - {SIZE_MB} MB</title>
<style>body{{font-family:system-ui,sans-serif;margin:1rem}}table{{border-collapse:collapse;font:11px monospace}}
td,th{{border:1px solid #8884;padding:2px 6px;text-align:right}}</style></head><body>
<h1>Number theory, trigonometry and combinatorics (~{SIZE_MB} MB)</h1>
<p>About {est_rows:,} rows, one per integer n = 12, 13, 14, ... Columns: prime factorization, Euler phi, number of
divisors, sum of divisors, Moebius mu, is-prime flag, sqrt(n) to 40 digits, sin(x), cos(x), ln(n) with
x = 0.1 + 1.2*n/limit, C(n,10), k = n mod 1000 + 1, and Fibonacci F(k) mod 1,000,000,007.</p>
<p><b>Verification:</b> every row passed <b>200</b> separate checks: 20 number-theory checks (factor product,
Miller-Rabin vs sieve, phi/tau/sigma by two methods, divisor sums, Moebius sums, Euler's theorem),
100 modular checks, 40 floating-point and trigonometric identities, and 40 Fibonacci/Pascal/binomial identities.</p>
<p><b>Estimated generation time on this machine:</b> about {est_min:.1f} minutes.</p>
<p><b>Warning:</b> browsers freeze on files this big. Use a large-file text editor, <code>less</code> or <code>grep</code>.</p>
<table><thead><tr><th>n<th>factors<th>phi<th>tau<th>sigma<th>mu<th>prime<th>sqrt(n)<th>sin x<th>cos x<th>ln n<th>C(n,10)<th>k<th>F(k) mod p</thead><tbody>
"""
    start = time.perf_counter()
    written, n = 0, 12
    with open(OUT, "w", encoding="ascii", newline="\n") as f:
        f.write(head); written += len(head)
        while written < TARGET - 600 and n < limit:
            stop = min(n + 500, limit)
            chunk = "".join(row(i) for i in range(n, stop))
            f.write(chunk); written += len(chunk); n = stop
            print(f"\r{written/1e6:7.1f} MB / {SIZE_MB} MB  {time.perf_counter()-start:6.0f}s", end="", flush=True)
        elapsed = time.perf_counter() - start
        f.write(f"</tbody></table>\n<p><b>Done.</b> {n-12:,} rows, 200/200 checks passed on every row. "
                f"Actual time: {elapsed/60:.1f} min.</p>\n</body></html>\n")
    print(f"\nDone: {OUT}, {n-12:,} rows, {elapsed/60:.1f} min")


main()
