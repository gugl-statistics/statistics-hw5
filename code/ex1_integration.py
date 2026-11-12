"""
Homework 5 - Ex. 1: Riemann vs Lebesgue numerical integration.
Python 3, standard library only.  Run:  python3 ex1_integration.py
"""
import math, time, json
from fractions import Fraction

# ---------------------------------------------------------------- the functions
FUNCS = {
    "sin":     dict(f=math.sin,                                a=0.0, b=math.pi, exact=2.0,               label="sin(x) on [0, pi]"),
    "x2":      dict(f=lambda x: x * x,                         a=0.0, b=1.0,     exact=1 / 3,             label="x^2 on [0, 1]"),
    "step13":  dict(f=lambda x: 1.0 if x < 0.5 else 3.0,       a=0.0, b=1.0,     exact=2.0,               label="1 if x<0.5, 3 otherwise"),
    "abs":     dict(f=lambda x: abs(x - 0.5),                  a=0.0, b=1.0,     exact=0.25,              label="|x - 0.5| on [0, 1]"),
    "invsqrt": dict(f=lambda x: 1 / math.sqrt(x),              a=0.0, b=1.0,     exact=2.0,               label="1/sqrt(x) on (0, 1]", unbounded=True),
    "step10":  dict(f=lambda x: 1.0 if x < 0.5 else 0.0,       a=0.0, b=1.0,     exact=0.5,               label="1 if x<0.5, 0 otherwise"),
    "sinc":    dict(f=lambda x: math.sin(x) / x if x else 1.0, a=0.0, b=math.pi, exact=1.851937051982466, label="sin(x)/x on [0, pi]"),
}

class Counter:
    """wraps f and counts how many times it is evaluated"""
    def __init__(self, f): self.f, self.n = f, 0
    def __call__(self, x): self.n += 1; return self.f(x)

# ---------------------------------------------------------------- Riemann: split the DOMAIN
def riemann(f, a, b, N, rule="mid"):
    h = (b - a) / N
    off = 0.5 if rule == "mid" else 0.0
    return h * sum(f(a + (i + off) * h) for i in range(N))

# ---------------------------------------------------------------- Lebesgue: split the RANGE
# For f >= 0:  integral = sum over bands [t_j, t_j+1] of (t_j+1 - t_j) * mu({f > t})   ("layer cake")
# which is the integral of the simple function  sum_j t_j * 1{t_j < f <= t_j+1}.
# mu({f > t}) is the length of the preimage, found from the sign changes of f - t on a grid
# and refined with bisection.
def level_set_measure(f, t, xs, fs, iters=55):
    total = 0.0
    inside = fs[0] > t
    start = xs[0] if inside else None
    for i in range(1, len(xs)):
        now = fs[i] > t
        if now != inside:                       # the boundary of {f > t} is between xs[i-1] and xs[i]
            lo, hi = xs[i - 1], xs[i]
            for _ in range(iters):
                mid = 0.5 * (lo + hi)
                if (f(mid) > t) == inside: lo = mid
                else: hi = mid
            cross = 0.5 * (lo + hi)
            if inside: total += cross - start
            else: start = cross
            inside = now
    if inside: total += xs[-1] - start
    return total

def _defined(f, x):
    try: f(x); return True
    except (ZeroDivisionError, ValueError): return False

def lebesgue(f, a, b, M, K=1000, levels=None, value="mid"):
    xs = [a + (b - a) * i / K for i in range(K + 1)]
    if not _defined(f, a): xs[0] = a + 1e-300            # e.g. 1/sqrt(x) at 0
    fs = [f(x) for x in xs]
    if levels is None:                                    # M equal bands between 0 and max f
        top = max(fs)
        levels = [top * j / M for j in range(M + 1)]
    total = 0.0
    for j in range(len(levels) - 1):
        t0, t1 = levels[j], levels[j + 1]
        t = 0.5 * (t0 + t1) if value == "mid" else t1     # "lower" = classical lower Lebesgue sum
        total += (t1 - t0) * level_set_measure(f, t, xs, fs)
    return total

def geometric_levels(top, m_lin, m_geo, big=1e12):
    """for an unbounded f: m_lin equal bands up to `top`, then m_geo geometric bands up to `big`"""
    lv = [top * j / m_lin for j in range(m_lin + 1)]
    r = (big / top) ** (1 / m_geo)
    lv += [top * r ** j for j in range(1, m_geo + 1)]
    return lv

# ---------------------------------------------------------------- experiments
def convergence():
    out = {}
    Ns = [11, 31, 101, 301, 1001, 3001, 10001]   # odd N: the jump at 0.5 is NOT on a grid point
    for key, F in FUNCS.items():
        rows = []
        for N in Ns:
            c = Counter(F["f"]); t0 = time.perf_counter()
            R = riemann(c, F["a"], F["b"], N); tr = time.perf_counter() - t0; nr = c.n
            c = Counter(F["f"]); t0 = time.perf_counter()
            if F.get("unbounded"):
                L = lebesgue(c, F["a"], F["b"], N, levels=geometric_levels(10.0, N // 2, N - N // 2))
            else:
                L = lebesgue(c, F["a"], F["b"], N)
            tl = time.perf_counter() - t0; nl = c.n
            rows.append(dict(N=N, riemann=R, lebesgue=L, err_r=abs(R - F["exact"]), err_l=abs(L - F["exact"]),
                             evals_r=nr, evals_l=nl, time_r=tr, time_l=tl))
            print(f"{key:8s} N={N:6d}  R={R:.12f} err {abs(R - F['exact']):.2e} ({nr:>8} evals, {tr * 1e3:8.1f} ms)   "
                  f"L={L:.12f} err {abs(L - F['exact']):.2e} ({nl:>8} evals, {tl * 1e3:8.1f} ms)")
        out[key] = dict(label=F["label"], exact=F["exact"], rows=rows)
    return out

def lower_vs_mid():
    print("\nLebesgue: lower sum vs midpoint-level sum")
    res = {}
    for key in ["x2", "sin", "step13"]:
        F = FUNCS[key]; res[key] = []
        for M in [4, 16, 64, 256, 1024]:
            lo = lebesgue(F["f"], F["a"], F["b"], M, value="lower")
            mi = lebesgue(F["f"], F["a"], F["b"], M, value="mid")
            res[key].append(dict(M=M, lower=lo, mid=mi))
            print(f"{key:7s} M={M:5d}  lower={lo:.10f} (err {abs(lo - F['exact']):.2e})   mid={mi:.10f} (err {abs(mi - F['exact']):.2e})")
    return res

# ---------------------------------------------------------------- pathological functions (exact arithmetic)
def thomae(x: Fraction):
    return Fraction(1, x.denominator)

def smallest_denominator_in(lo: Fraction, hi: Fraction):
    """the fraction with the smallest denominator in [lo, hi] (Stern-Brocot tree)"""
    n = math.ceil(lo)
    if n <= hi: return Fraction(n)
    a, b, c, d = math.floor(lo), 1, math.floor(lo) + 1, 1
    while True:
        m = Fraction(a + c, b + d)
        if m < lo: a, b = a + c, b + d
        elif m > hi: c, d = a + c, b + d
        else: return m

def pathological():
    print("\nThomae: Riemann sums with rational (midpoint) tags, irrational tags, and upper Darboux sums")
    th = []
    for N in [10, 100, 1000, 10000]:
        mid = sum(thomae(Fraction(2 * i + 1, 2 * N)) for i in range(N)) / N
        left = sum(thomae(Fraction(i, N)) for i in range(N)) / N
        upper = sum(Fraction(1, smallest_denominator_in(Fraction(i, N), Fraction(i + 1, N)).denominator) for i in range(N)) / N
        th.append(dict(N=N, mid=float(mid), left=float(left), upper=float(upper), irrational=0.0))
        print(f"N={N:6d}  midpoint tags {float(mid):.6f}  left tags {float(left):.6f}  upper Darboux {float(upper):.6f}  irrational tags 0")
    lev = []
    for t in [0.5, 0.1, 0.01, 0.001]:
        qmax = math.ceil(1 / t) - 1                      # f(p/q) = 1/q > t  <=>  q < 1/t
        count = sum(1 for q in range(1, qmax + 1) for p in range(0, q + 1) if math.gcd(p, q) == 1)
        lev.append(dict(t=t, points=count))
        print(f"t={t}: {{f > t}} has {count} points -> measure 0")
    print("\nDirichlet: rational tags give 1, irrational tags give 0, for every N")
    eps = 1e-3
    print("cover of the first 10^5 rationals with lengths eps/2^(k+1): total", eps * (1 - 0.5 ** 100000), "< eps =", eps)
    return dict(thomae=th, thomae_levels=lev)

if __name__ == "__main__":
    conv = convergence()
    lm = lower_vs_mid()
    path = pathological()
    json.dump(dict(conv=conv, lower_mid=lm, path=path), open("ex1_results.json", "w"), indent=1)
