// ===================================================================
// integrate.js - Riemann (split the domain) vs Lebesgue (split the range)
// ===================================================================

const FUNCS = {
  sin:     { label: "sin(x), [0, π]",            f: Math.sin,                       a: 0, b: Math.PI, exact: 2 },
  x2:      { label: "x², [0, 1]",                f: (x) => x * x,                   a: 0, b: 1, exact: 1 / 3 },
  step13:  { label: "1 if x < 0.5, 3 otherwise",  f: (x) => (x < 0.5 ? 1 : 3),      a: 0, b: 1, exact: 2 },
  abs:     { label: "|x − 0.5|, [0, 1]",         f: (x) => Math.abs(x - 0.5),       a: 0, b: 1, exact: 0.25 },
  invsqrt: { label: "1/√x, (0, 1]",              f: (x) => 1 / Math.sqrt(x),        a: 0, b: 1, exact: 2, unbounded: true },
  step10:  { label: "1 if x < 0.5, 0 otherwise",  f: (x) => (x < 0.5 ? 1 : 0),      a: 0, b: 1, exact: 0.5 },
  sinc:    { label: "sin(x)/x, [0, π]",          f: (x) => (x === 0 ? 1 : Math.sin(x) / x), a: 0, b: Math.PI, exact: 1.851937051982466 },
};

// wrap f to count evaluations
function counted(f) {
  const g = (x) => { g.n++; return f(x); };
  g.n = 0;
  return g;
}

// ---------- Riemann: N subintervals of the domain, value at the midpoint (or left end)
function riemann(f, a, b, N, rule = "mid") {
  const h = (b - a) / N, off = rule === "mid" ? 0.5 : 0;
  let s = 0;
  for (let i = 0; i < N; i++) s += f(a + (i + off) * h);
  return s * h;
}

// ---------- Lebesgue: bands [t_j, t_j+1] of the range, measure of the preimage {f > t}
// integral (f >= 0) = sum_j (t_j+1 - t_j) * mu({f > t_j*})        (layer cake formula)
function levelSetMeasure(f, t, xs, fs, iters = 55, keepIntervals = null) {
  let total = 0, inside = fs[0] > t, start = inside ? xs[0] : null;
  for (let i = 1; i < xs.length; i++) {
    const now = fs[i] > t;
    if (now !== inside) {
      let lo = xs[i - 1], hi = xs[i];
      for (let k = 0; k < iters; k++) {
        const mid = 0.5 * (lo + hi);
        if ((f(mid) > t) === inside) lo = mid; else hi = mid;
      }
      const cross = 0.5 * (lo + hi);
      if (inside) { total += cross - start; if (keepIntervals) keepIntervals.push([start, cross]); }
      else start = cross;
      inside = now;
    }
  }
  if (inside) { total += xs[xs.length - 1] - start; if (keepIntervals) keepIntervals.push([start, xs[xs.length - 1]]); }
  return total;
}

function makeGrid(f, a, b, K) {
  const xs = new Float64Array(K + 1), fs = new Float64Array(K + 1);
  for (let i = 0; i <= K; i++) xs[i] = a + ((b - a) * i) / K;
  if (!isFinite(f(a))) xs[0] = a + 1e-300;
  for (let i = 0; i <= K; i++) fs[i] = f(xs[i]);
  return { xs, fs };
}

function geometricLevels(top, mLin, mGeo, big = 1e12) {
  const lv = [];
  for (let j = 0; j <= mLin; j++) lv.push((top * j) / mLin);
  const r = Math.pow(big / top, 1 / mGeo);
  for (let j = 1; j <= mGeo; j++) lv.push(top * Math.pow(r, j));
  return lv;
}

// returns { value, levels, sets } (sets = preimage intervals of each band, for the picture)
function lebesgue(f, a, b, M, opts = {}) {
  const K = opts.K || 1000, value = opts.value || "mid";
  const { xs, fs } = makeGrid(f, a, b, K);
  let levels = opts.levels;
  if (!levels) {
    let top = 0;
    for (const v of fs) if (v > top) top = v;
    levels = Array.from({ length: M + 1 }, (_, j) => (top * j) / M);
  }
  let total = 0;
  const sets = [];
  for (let j = 0; j + 1 < levels.length; j++) {
    const t0 = levels[j], t1 = levels[j + 1];
    const t = value === "mid" ? 0.5 * (t0 + t1) : t1;
    const iv = opts.keep ? [] : null;
    const mu = levelSetMeasure(f, t, xs, fs, 55, iv);
    total += (t1 - t0) * mu;
    if (opts.keep) sets.push({ t0, t1, mu, iv });
  }
  return { value: total, levels, sets };
}

// ---------- Thomae's function on exact fractions p/q
function gcd(a, b) { while (b) [a, b] = [b, a % b]; return a; }
// fraction with the smallest denominator in [lo, hi], lo = ln/ld, hi = hn/hd (Stern-Brocot)
function smallestDenominator(ln, ld, hn, hd) {
  const n = Math.ceil(ln / ld);
  if (n * hd <= hn) return 1;
  let a = Math.floor(ln / ld), b = 1, c = a + 1, d = 1;
  for (;;) {
    const p = a + c, q = b + d;
    if (p * ld < ln * q) { a = p; b = q; }
    else if (p * hd > hn * q) { c = p; d = q; }
    else return q;
  }
}
function thomaeSums(N) {
  let mid = 0, left = 0, upper = 0;
  for (let i = 0; i < N; i++) {
    const p = 2 * i + 1, q = 2 * N;
    mid += gcd(p, q) / q;                          // f(p/q) = 1 / (q / gcd)
    left += i === 0 ? 1 : gcd(i, N) / N;
    upper += 1 / smallestDenominator(i, N, i + 1, N);
  }
  return { mid: mid / N, left: left / N, upper: upper / N };
}

if (typeof module !== "undefined") module.exports = { FUNCS, counted, riemann, lebesgue, geometricLevels, thomaeSums, smallestDenominator };
