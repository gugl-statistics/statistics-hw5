"""
Homework 5 - Ex. 2: numerical computation of the standard normal CDF Phi(z).
Python 3 + mpmath (only for the 50-digit reference).  Run:  python3 ex2_normal_cdf.py
"""
import math, time, json
import mpmath
mpmath.mp.dps = 50

SQ2PI = math.sqrt(2 * math.pi)
phi = lambda t: math.exp(-0.5 * t * t) / SQ2PI

# ---------------------------------------------------------------- 1. numerical integration of the density
def cdf_trapezoid(z, N=100):
    """Phi(z) = 1/2 + integral_0^z phi(t) dt, trapezoid rule with N panels"""
    h = z / N
    s = 0.5 * (phi(0) + phi(z)) + sum(phi(i * h) for i in range(1, N))
    return 0.5 + h * s

def cdf_simpson(z, N=100):
    if N % 2: N += 1
    h = z / N
    s = phi(0) + phi(z) + 4 * sum(phi((2 * i - 1) * h) for i in range(1, N // 2 + 1)) \
        + 2 * sum(phi(2 * i * h) for i in range(1, N // 2))
    return 0.5 + h * s / 3

# ---------------------------------------------------------------- 2. power series (the one in the Hogg-McKean-Craig notes)
# Phi(z) = 1/2 + 1/sqrt(2 pi) * sum_k (-1)^k z^(2k+1) / (2^k (2k+1) k!)
def cdf_series(z, tol=1e-17, kmax=500):
    term = z                    # (-1)^k z^(2k+1) / (2^k k!)   without the 1/(2k+1)
    s = z
    k = 0
    while k < kmax:
        k += 1
        term *= -z * z / (2 * k)
        add = term / (2 * k + 1)
        s += add
        if abs(add) < tol * abs(s): break
    return 0.5 + s / SQ2PI, k

# ---------------------------------------------------------------- 3. classic polynomial approximation (Abramowitz & Stegun 26.2.17)
def cdf_as26217(z):
    p, b = 0.2316419, (0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429)
    x = abs(z)
    t = 1 / (1 + p * x)
    poly = t * (b[0] + t * (b[1] + t * (b[2] + t * (b[3] + t * b[4]))))
    q = phi(x) * poly           # upper tail Q(x)
    return 1 - q if z >= 0 else q

# ---------------------------------------------------------------- 4. continued fraction for the tail (Laplace)
# Q(x) = 1 - Phi(x) = phi(x) / (x + 1/(x + 2/(x + 3/(x + ...))))   for x > 0
def tail_cf(x, n):
    f = x
    for k in range(n, 0, -1):   # evaluated backwards, n levels deep
        f = x + k / f
    return phi(x) / f

# ---------------------------------------------------------------- 5. my final method: series in the middle, continued fraction in the tails
CROSS = 2.5
def cf_terms(x):                # depth that gives ~1e-16 relative error (tuned below)
    return int(10 + 400 / (x * x))

def cdf_final(z):
    x = abs(z)
    if x <= CROSS:
        return cdf_series(z)[0]
    q = tail_cf(x, cf_terms(x))
    return 1 - q if z > 0 else q

# ---------------------------------------------------------------- 6. library erfc (C math library), for comparison
def cdf_erfc(z):
    return 0.5 * math.erfc(-z / math.sqrt(2))

# ---------------------------------------------------------------- reference
def ref(z):
    return mpmath.ncdf(mpmath.mpf(z))

def relerr(v, z):
    """relative error |v - Phi(z)| / Phi(z); used on the lower tail z <= 0, where Phi(z) is small
    (for z > 0 a double near 1 can't hold the tail anyway, so there we look at the absolute error)"""
    r = ref(z)
    return float(abs(mpmath.mpf(v) - r) / r)

def abserr(v, z): return float(abs(mpmath.mpf(v) - ref(z)))

def main():
    out = {}
    # ---- tune the continued fraction depth
    print("continued fraction: depth needed for rel. error < 2e-16")
    tune = []
    for x in [2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 20.0]:
        r = 1 - ref(x) if x < 8 else mpmath.ncdf(-x)
        n = 1
        while True:
            q = tail_cf(x, n)
            if abs(mpmath.mpf(q) - r) / r < 2e-16 or n > 3000: break
            n += 1
        tune.append(dict(x=x, n=n, formula=cf_terms(x)))
        print(f"x={x:5}: n={n}   (formula {cf_terms(x)})")
    out["cf_tune"] = tune

    # ---- series: terms and accuracy vs z (cancellation)
    print("\npower series: terms needed and error")
    ser = []
    for z in [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]:
        v, k = cdf_series(z)
        biggest = max(z ** (2 * j + 1) / (2 ** j * math.factorial(j) * (2 * j + 1)) for j in range(k + 1))
        ser.append(dict(z=z, terms=k, abs_err=abserr(v, z), rel_tail=relerr(v, z), biggest=biggest / SQ2PI))
        print(f"z={z:4}: {k:3d} terms, biggest term {biggest / SQ2PI:.3e}, abs err {abserr(v, z):.2e}, rel err of tail {relerr(v, z):.2e}")
    out["series"] = ser

    # ---- trapezoid / Simpson convergence at z = 1.96
    print("\ntrapezoid and Simpson at z = 1.96")
    quad = []
    for N in [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        et, es = abserr(cdf_trapezoid(1.96, N), 1.96), abserr(cdf_simpson(1.96, N), 1.96)
        quad.append(dict(N=N, trap=et, simp=es))
        print(f"N={N:5d}: trapezoid err {et:.2e}   Simpson err {es:.2e}")
    out["quad"] = quad

    # ---- accuracy over a grid of z for every method
    print("\nmax abs error over z in [-10, 10], max relative error on the lower tail z in [-5, 0] and [-10, 0]")
    zs = [i / 20 for i in range(-200, 201)]
    meths = {
        "trapezoid N=100": lambda z: cdf_trapezoid(z, 100),
        "Simpson N=100": lambda z: cdf_simpson(z, 100),
        "power series": lambda z: cdf_series(z)[0],
        "A&S 26.2.17": cdf_as26217,
        "series + cont. fraction": cdf_final,
        "math.erfc (C library)": cdf_erfc,
    }
    acc = {}
    for name, f in meths.items():
        ae = [abserr(f(z), z) for z in zs]
        re = [relerr(f(z), z) for z in zs if z <= 0]
        re_mid = max(r for r, z in zip(re, [z for z in zs if z <= 0]) if z >= -5)
        acc[name] = dict(max_abs=max(ae), max_rel_5=re_mid, max_rel_10=max(re))
        print(f"{name:26s} max abs err {max(ae):.2e}   max rel err |z|<=5: {re_mid:.2e}   |z|<=10: {max(re):.2e}")
    out["accuracy"] = acc
    out["curve_z"] = zs
    out["curves_abs"] = {name: [max(abserr(f(z), z), 1e-18) for z in zs] for name, f in meths.items()}
    out["curves_rel"] = {name: [max(relerr(f(z), z), 1e-18) for z in zs if z <= 0] for name, f in meths.items()}

    # ---- speed
    print("\ntime per evaluation (Python, 20,000 random z in [-6, 6])")
    import random
    random.seed(1)
    Z = [random.uniform(-6, 6) for _ in range(20000)]
    speed = {}
    for name, f in meths.items():
        t0 = time.perf_counter()
        for z in Z: f(z)
        us = (time.perf_counter() - t0) / len(Z) * 1e6
        speed[name] = us
        print(f"{name:26s} {us:8.2f} us")
    out["speed"] = speed

    # ---- comparison with web references
    # Wikipedia "Standard normal table" (cumulative table, 5 decimals; tail table, 6 significant digits)
    # Wikipedia "68-95-99.7 rule": P(|Z| < 0.5) = 0.382924922548026, P(|Z| < 1) = 0.682689492137086
    web = [
        (0.0, 0.50000, "z-table"), (0.5, 0.69146, "z-table"), (1.0, 0.84134, "z-table"), (1.5, 0.93319, "z-table"),
        (1.96, 0.97500, "z-table"), (2.0, 0.97725, "z-table"), (2.5, 0.99379, "z-table"), (3.0, 0.99865, "z-table"),
        (3.5, 0.99977, "z-table"),
        (0.5, (1 + 0.382924922548026) / 2, "68-95-99.7 page (15 digits)"),
        (1.0, (1 + 0.682689492137086) / 2, "68-95-99.7 page (15 digits)"),
        # the tail table gives 1 - Phi(z) directly; I keep it as a tail (1 - 1e-16 can't be stored in a double)
        (4.0, 3.16712e-5, "tail table"), (5.0, 2.86652e-7, "tail table"),
        (6.0, 9.86588e-10, "tail table"), (7.0, 1.27981e-12, "tail table"), (8.0, 6.22096e-16, "tail table"),
    ]
    print("\ncomparison with the web values")
    rows = []
    for z, w, src in web:
        mine = cdf_final(z)
        tail_mine = 1 - mine if z > 0 else None
        q_mine = tail_cf(z, cf_terms(z)) if z > CROSS else 1 - mine
        rows.append(dict(z=z, web=w, src=src, mine=mine, tail_mine=q_mine, ref=float(ref(z)), ref_tail=float(1 - ref(z))))
        print(f"z={z:4}  web {w!r:<22} mine {mine:.16f}  tail mine {q_mine:.8e}  ({src})")
    out["web"] = rows
    json.dump(out, open("ex2_results.json", "w"), indent=1)

if __name__ == "__main__":
    main()
