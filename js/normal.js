// ===================================================================
// normal.js - ways to compute Phi(z), the standard normal CDF
// (JavaScript has no Math.erf, so everything is done by hand)
// ===================================================================

const SQ2PI = Math.sqrt(2 * Math.PI);
const pdf = (t) => Math.exp(-0.5 * t * t) / SQ2PI;

// 1. numerical integration of the density on [0, z]
function cdfTrapezoid(z, N = 100) {
  const h = z / N;
  let s = 0.5 * (pdf(0) + pdf(z));
  for (let i = 1; i < N; i++) s += pdf(i * h);
  return 0.5 + h * s;
}
function cdfSimpson(z, N = 100) {
  if (N % 2) N++;
  const h = z / N;
  let s = pdf(0) + pdf(z);
  for (let i = 1; i < N; i++) s += (i % 2 ? 4 : 2) * pdf(i * h);
  return 0.5 + (h * s) / 3;
}

// 2. power series (term-by-term integration of the series of e^(-t^2/2))
function cdfSeries(z, tol = 1e-17, kmax = 500) {
  let term = z, s = z, k = 0;
  while (k < kmax) {
    k++;
    term *= (-z * z) / (2 * k);
    const add = term / (2 * k + 1);
    s += add;
    if (Math.abs(add) < tol * Math.abs(s)) break;
  }
  cdfSeries.terms = k;
  return 0.5 + s / SQ2PI;
}

// 3. Abramowitz & Stegun 26.2.17 (|error| < 7.5e-8)
function cdfAS(z) {
  const p = 0.2316419, b = [0.31938153, -0.356563782, 1.781477937, -1.821255978, 1.330274429];
  const x = Math.abs(z), t = 1 / (1 + p * x);
  const q = pdf(x) * t * (b[0] + t * (b[1] + t * (b[2] + t * (b[3] + t * b[4]))));
  return z >= 0 ? 1 - q : q;
}

// 4. continued fraction for the upper tail Q(x) = 1 - Phi(x), x > 0
//    Q(x) = phi(x) / (x + 1/(x + 2/(x + 3/(x + ...))))
function tailCF(x, n) {
  let f = x;
  for (let k = n; k >= 1; k--) f = x + k / f;
  return pdf(x) / f;
}
const cfTerms = (x) => Math.floor(10 + 400 / (x * x));

// 5. final method: series for |z| <= 2.5, continued fraction outside
const CROSS = 2.5;
function cdfFinal(z) {
  const x = Math.abs(z);
  if (x <= CROSS) return cdfSeries(z);
  const q = tailCF(x, cfTerms(x));
  return z > 0 ? 1 - q : q;
}
// upper tail with full relative precision (for z > 2.5)
function tailFinal(z) {
  return z > CROSS ? tailCF(z, cfTerms(z)) : 1 - cdfFinal(z);
}

const METHODS = {
  trap:   { name: "Trapezoid (N = 100)", f: (z) => cdfTrapezoid(z, 100) },
  simp:   { name: "Simpson (N = 100)",   f: (z) => cdfSimpson(z, 100) },
  series: { name: "Power series",        f: (z) => cdfSeries(z) },
  as:     { name: "A&S 26.2.17",         f: cdfAS },
  final:  { name: "Series + cont. fraction (mine)", f: cdfFinal },
};

if (typeof module !== "undefined") module.exports = { cdfTrapezoid, cdfSimpson, cdfSeries, cdfAS, tailCF, cfTerms, cdfFinal, tailFinal, METHODS };
