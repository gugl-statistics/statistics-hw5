// ===================================================================
// main.js - connects the page to integrate.js, normal.js and charts.js
// ===================================================================
const $ = (id) => document.getElementById(id);
const sci = (x, d = 2) => (x === 0 ? "0" : !isFinite(x) ? "-" : Math.abs(x) >= 1e-3 && Math.abs(x) < 1e4 ? (+x.toPrecision(d + 1)).toString() : x.toExponential(d));
const fmtInt = (x) => Math.round(x).toLocaleString("en-US");

function seg(el, onChange) {
  const btns = [...el.querySelectorAll("button")];
  btns.forEach((b) => (b.onclick = () => { btns.forEach((x) => x.classList.toggle("on", x === b)); onChange && onChange(b.dataset.v); }));
  return { get: () => el.querySelector("button.on").dataset.v };
}

// =================================================== 1.2 demo
const fSel = $("d-f");
for (const [k, F] of Object.entries(FUNCS)) fSel.insertAdjacentHTML("beforeend", `<option value="${k}">${F.label}</option>`);

function drawDemo() {
  const F = FUNCS[fSel.value], N = +$("d-n").value;
  $("d-nv").textContent = N;
  $("d-band").max = N;
  if (+$("d-band").value > N) $("d-band").value = N;
  const band = +$("d-band").value;
  $("d-bv").textContent = band;
  const f = F.f, a = F.a, b = F.b;
  const yTop = F.unbounded ? 6 : Math.max(...Array.from({ length: 400 }, (_, i) => f(a + ((b - a) * (i + 0.5)) / 400))) * 1.08;
  const curve = { type: "fn", f: (x) => Math.min(f(x), yTop * 2), color: COLORS.white, width: 2 };

  // Riemann rectangles
  const h = (b - a) / N, rects = [];
  for (let i = 0; i < N; i++) rects.push([a + i * h, 0, a + (i + 1) * h, f(a + (i + 0.5) * h)]);
  plot($("d-riemann"), { x: [a, b], y: [0, yTop], series: [{ type: "rects", rects, color: COLORS.sky, alpha: 0.45, stroke: "rgba(92,200,245,0.9)" }, curve] });

  // Lebesgue bands (for the picture, 1/sqrt(x) uses equal bands up to the top of the plot)
  const opt = { keep: true };
  if (F.unbounded) opt.levels = Array.from({ length: N + 1 }, (_, j) => (yTop * j) / N);
  const L = lebesgue(f, a, b, N, opt);
  const lr = [], hl = [], marks = [];
  L.sets.forEach((s, j) => {
    for (const [x1, x2] of s.iv) (j === band - 1 ? hl : lr).push([x1, s.t0, x2, s.t1]);
    if (j === band - 1) for (const [x1, x2] of s.iv) marks.push([x1, -yTop * 0.03, x2, 0]);
  });
  plot($("d-lebesgue"), {
    x: [a, b], y: [-yTop * 0.04, yTop],
    series: [
      { type: "rects", rects: lr, color: COLORS.sky, alpha: 0.45, stroke: "rgba(92,200,245,0.9)" },
      { type: "rects", rects: hl, color: COLORS.orange, alpha: 0.75 },
      { type: "rects", rects: marks, color: COLORS.orange, alpha: 1 },
      curve,
    ],
  });

  // numbers (the real Lebesgue value, with the geometric levels for 1/sqrt(x))
  const cr = counted(f), cl = counted(f);
  const R = riemann(cr, a, b, N);
  const lv = F.unbounded ? geometricLevels(10, Math.ceil(N / 2), Math.floor(N / 2) || 1) : null;
  const Lv = lebesgue(cl, a, b, N, { levels: lv }).value;
  const s = L.sets[band - 1];
  $("d-table").innerHTML =
    `<tr><th></th><th class="num">value</th><th class="num">exact</th><th class="num">error</th><th class="num">evaluations of f</th></tr>` +
    `<tr><td>Riemann, ${N} strips</td><td class="num">${R.toFixed(8)}</td><td class="num">${F.exact.toFixed(8)}</td><td class="num">${sci(Math.abs(R - F.exact))}</td><td class="num">${fmtInt(cr.n)}</td></tr>` +
    `<tr><td>Lebesgue, ${N} bands${F.unbounded ? " (half geometric, up to 10¹²)" : ""}</td><td class="num">${Lv.toFixed(8)}</td><td class="num">${F.exact.toFixed(8)}</td><td class="num">${sci(Math.abs(Lv - F.exact))}</td><td class="num">${fmtInt(cl.n)}</td></tr>` +
    `<tr><td colspan="5" style="color:var(--muted)">highlighted band: ${s.t0.toFixed(3)} ≤ t ≤ ${s.t1.toFixed(3)}; preimage {f &gt; ${(0.5 * (s.t0 + s.t1)).toFixed(3)}} = ` +
    (s.iv.length ? s.iv.map(([x1, x2]) => `[${x1.toFixed(4)}, ${x2.toFixed(4)}]`).join(" ∪ ") : "∅") + `, measure ${s.mu.toFixed(5)}</td></tr>`;
}
fSel.value = "sin";
fSel.onchange = drawDemo;
$("d-n").oninput = drawDemo;
$("d-band").oninput = drawDemo;

// =================================================== 1.3 results (Python)
const rSeg = $("r-f");
let rKey = "sin";
Object.entries(EX1.conv).forEach(([k, v]) => rSeg.insertAdjacentHTML("beforeend", `<button data-v="${k}" class="${k === rKey ? "on" : ""}">${FUNCS[k].label}</button>`));
seg(rSeg, (v) => { rKey = v; drawResults(); });

function drawResults() {
  const C = EX1.conv[rKey], rows = C.rows;
  const floor = 1e-17;
  const allE = rows.flatMap((r) => [r.err_r, r.err_l]).map((e) => Math.max(e, floor));
  const ev = rows.flatMap((r) => [r.evals_r, r.evals_l]);
  plot($("r-chart"), {
    x: [Math.min(...ev) / 2, Math.max(...ev) * 2], y: [Math.max(Math.min(...allE) / 3, 1e-17), Math.max(...allE) * 3], xlog: true, ylog: true,
    xlabel: "function evaluations", ylabel: "absolute error",
    series: [
      { type: "line", x: rows.map((r) => r.evals_r), y: rows.map((r) => Math.max(r.err_r, floor)), color: COLORS.sky, width: 2 },
      { type: "points", x: rows.map((r) => r.evals_r), y: rows.map((r) => Math.max(r.err_r, floor)), color: COLORS.sky },
      { type: "line", x: rows.map((r) => r.evals_l), y: rows.map((r) => Math.max(r.err_l, floor)), color: COLORS.orange, width: 2 },
      { type: "points", x: rows.map((r) => r.evals_l), y: rows.map((r) => Math.max(r.err_l, floor)), color: COLORS.orange },
    ],
  });
  let t = `<tr><th class="num">N = M</th><th class="num">Riemann</th><th class="num">error</th><th class="num">evals</th><th class="num">time</th><th class="num">Lebesgue</th><th class="num">error</th><th class="num">evals</th><th class="num">time</th></tr>`;
  for (const r of rows) {
    t += `<tr><td class="num">${fmtInt(r.N)}</td><td class="num">${r.riemann.toFixed(10)}</td><td class="num">${sci(r.err_r)}</td><td class="num">${fmtInt(r.evals_r)}</td><td class="num">${(r.time_r * 1e3).toFixed(2)} ms</td>` +
         `<td class="num">${r.lebesgue.toFixed(10)}</td><td class="num">${sci(r.err_l)}</td><td class="num">${fmtInt(r.evals_l)}</td><td class="num">${(r.time_l * 1e3).toFixed(1)} ms</td></tr>`;
  }
  t += `<tr class="best"><td>exact</td><td class="num" colspan="8" style="text-align:left">${C.exact}</td></tr>`;
  $("r-table").innerHTML = t;
}

(function summary() {
  let t = `<tr><th>f</th><th class="num">exact</th><th class="num">Riemann</th><th class="num">error</th><th class="num">Lebesgue</th><th class="num">error</th><th class="num">evals ratio</th><th>better</th></tr>`;
  for (const [k, C] of Object.entries(EX1.conv)) {
    const r = C.rows[C.rows.length - 1];
    const better = r.err_l < r.err_r / 3 ? "Lebesgue" : r.err_r < r.err_l / 3 ? "Riemann" : "similar";
    t += `<tr><td>${FUNCS[k].label}</td><td class="num">${(+C.exact).toFixed(10)}</td><td class="num">${r.riemann.toFixed(10)}</td><td class="num">${sci(r.err_r)}</td>` +
         `<td class="num">${r.lebesgue.toFixed(10)}</td><td class="num">${sci(r.err_l)}</td><td class="num">×${(r.evals_l / r.evals_r).toFixed(0)}</td><td>${better}</td></tr>`;
  }
  $("r-summary").innerHTML = t;
  const ex = { x2: 1 / 3, sin: 2, step13: 2 };
  let u = `<tr><th>f</th><th class="num">M</th><th class="num">lower sum</th><th class="num">error</th><th class="num">midpoint of band</th><th class="num">error</th></tr>`;
  for (const [k, rows] of Object.entries(EX1.lower_mid))
    for (const r of rows) u += `<tr><td>${FUNCS[k].label}</td><td class="num">${r.M}</td><td class="num">${r.lower.toFixed(8)}</td><td class="num">${sci(Math.abs(r.lower - ex[k]))}</td><td class="num">${r.mid.toFixed(8)}</td><td class="num">${sci(Math.abs(r.mid - ex[k]))}</td></tr>`;
  $("r-lower").innerHTML = u;
})();

// =================================================== 1.4 pathological
const dirTags = seg($("dir-tags"), drawDirichlet);
function drawDirichlet() {
  const N = 10 ** +$("dir-n").value;
  $("dir-nv").textContent = fmtInt(N);
  const rational = dirTags.get() === "rat";
  // (i + 1/2)/N is rational; (i + 1/sqrt 2)/N is irrational for every i, N (a rational multiple of sqrt 2 plus a rational)
  const value = rational ? 1 : 0;
  $("dir-out").innerHTML =
    `Riemann sum with ${fmtInt(N)} strips and ${rational ? "rational" : "irrational"} tags = Σ f(ξ<sub>i</sub>)·(1/N) = <span class="k">${value}</span>\n` +
    `upper Darboux sum = 1, lower Darboux sum = 0 (for every N)\n` +
    `Lebesgue: 1·μ(ℚ∩[0,1]) + 0·μ(irrationals) = 1·0 + 0·1 = <span class="s">0</span>\n` +
    `(note: in floating point, (i + 1/√2)/N would be stored as a rational number m/2^k, so a naive program would give 1 here too)`;
}
$("dir-n").oninput = drawDirichlet;

(function thomaePlot() {
  const xs = [], ys = [];
  for (let q = 1; q <= 60; q++) for (let p = 0; p <= q; p++) if (gcd(p, q) === 1) { xs.push(p / q); ys.push(1 / q); }
  plot($("th-plot"), { x: [0, 1], y: [0, 1.05], series: [{ type: "points", x: xs, y: ys, color: COLORS.sky, r: 2 }] });
  let t = `<tr><th class="num">N</th><th class="num">midpoint tags</th><th class="num">left tags</th><th class="num">upper Darboux</th><th class="num">irrational tags / lower</th></tr>`;
  for (const r of EX1.path.thomae) t += `<tr><td class="num">${fmtInt(r.N)}</td><td class="num">${r.mid.toFixed(6)}</td><td class="num">${r.left.toFixed(6)}</td><td class="num">${r.upper.toFixed(6)}</td><td class="num">0</td></tr>`;
  $("th-table").innerHTML = t;
})();
function drawThomae() {
  const N = Math.round(10 ** +$("th-n").value);
  $("th-nv").textContent = fmtInt(N);
  const s = thomaeSums(N);
  $("th-out").innerHTML = `midpoint tags: <span class="k">${s.mid.toFixed(6)}</span>\nleft tags: ${s.left.toFixed(6)}\nupper Darboux: ${s.upper.toFixed(6)}\nirrational tags / lower Darboux: 0\nLebesgue: <span class="s">0</span> (every {f &gt; t} is finite)`;
}
$("th-n").oninput = drawThomae;

// =================================================== 2.3 tables
(function methodTables() {
  let t = `<tr><th class="num">z</th><th class="num">terms</th><th class="num">biggest term</th><th class="num">absolute error</th><th class="num">rel. error of 1 − Φ(z)</th></tr>`;
  for (const r of EX2.series) t += `<tr><td class="num">${r.z}</td><td class="num">${r.terms}</td><td class="num">${sci(r.biggest)}</td><td class="num">${sci(r.abs_err)}</td><td class="num">${sci(r.rel_tail)}</td></tr>`;
  $("ser-table").innerHTML = t;
  let u = `<tr><th class="num">x</th>${EX2.cf_tune.map((r) => `<th class="num">${r.x}</th>`).join("")}</tr>` +
          `<tr><td>depth needed (rel. error &lt; 2·10⁻¹⁶)</td>${EX2.cf_tune.map((r) => `<td class="num">${r.n}</td>`).join("")}</tr>` +
          `<tr><td>depth used, 10 + 400/x²</td>${EX2.cf_tune.map((r) => `<td class="num">${r.formula}</td>`).join("")}</tr>`;
  $("cf-table").innerHTML = u;
})();

// =================================================== 2.3 calculator
const refAt = (z) => PHI_REF[Math.round(z * 100) + 1000];
function drawCalc() {
  const z = +$("c-z").value / 100;
  $("c-zv").textContent = z.toFixed(2);
  const r = refAt(z);
  const rel = z <= 0;
  let t = `<tr><th>method</th><th class="num">Φ(z)</th><th class="num">${rel ? "relative" : "absolute"} error</th><th class="num">time / eval</th></tr>`;
  for (const [k, M] of Object.entries(METHODS)) {
    const v = M.f(z);
    const e = rel ? Math.abs(v - r) / r : Math.abs(v - r);
    let reps = 0;
    const t0 = performance.now();
    while (performance.now() - t0 < 4) { for (let i = 0; i < 200; i++) M.f(z); reps += 200; }
    const ns = ((performance.now() - t0) / reps) * 1e6;
    t += `<tr${k === "final" ? ' class="best"' : ""}><td>${M.name}</td><td class="num">${v.toPrecision(16)}</td><td class="num">${sci(e)}</td><td class="num">${ns < 1000 ? ns.toFixed(0) + " ns" : (ns / 1000).toFixed(1) + " µs"}</td></tr>`;
  }
  t += `<tr><td>reference (mpmath)</td><td class="num">${r.toPrecision(16)}</td><td></td><td></td></tr>`;
  $("c-table").innerHTML = t;
  cdfSeries(z);
  const tail = z > 0 ? tailFinal(z) : null;
  $("c-note").innerHTML =
    `power series used ${cdfSeries.terms} terms` + (Math.abs(z) > CROSS ? ` · continued fraction depth ${cfTerms(Math.abs(z))}` : "") +
    (z > CROSS ? `\nupper tail 1 − Φ(z) from the continued fraction: <span class="k">${tail.toExponential(10)}</span>, reference Φ(−z) = <span class="s">${refAt(-z).toExponential(10)}</span> (near 1 a double can't store Φ(z) itself with these digits, so the tail is the useful number)` : "") +
    (z <= 0 ? "\nfor z ≤ 0 the error is relative, which is what matters for small probabilities" : "");
  // convergence of trapezoid / Simpson at this z
  const Ns = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048];
  const fl = 1e-17;
  const et = Ns.map((n) => Math.max(Math.abs(cdfTrapezoid(z, n) - r) / (rel ? r : 1), fl));
  const es = Ns.map((n) => Math.max(Math.abs(cdfSimpson(z, n) - r) / (rel ? r : 1), fl));
  plot($("c-conv"), {
    x: [1.5, 3000], y: [Math.max(1e-17, Math.min(...et, ...es) / 3), Math.max(...et, ...es, 1e-15) * 3], xlog: true, ylog: true,
    xlabel: "panels N", ylabel: rel ? "relative error" : "absolute error",
    series: [
      { type: "line", x: Ns, y: et, color: COLORS.sky, width: 2 }, { type: "points", x: Ns, y: et, color: COLORS.sky },
      { type: "line", x: Ns, y: es, color: COLORS.orange, width: 2 }, { type: "points", x: Ns, y: es, color: COLORS.orange },
    ],
  });
}
$("c-z").oninput = drawCalc;
document.querySelectorAll("#methods button[data-z]").forEach((b) => (b.onclick = () => { $("c-z").value = Math.round(+b.dataset.z * 100); drawCalc(); }));

// =================================================== 2.4 accuracy
const ACC_COL = { "trapezoid N=100": COLORS.grey, "Simpson N=100": "#a78bfa", "power series": COLORS.white, "A&S 26.2.17": "#6ee7a8",
                  "series + cont. fraction": COLORS.orange, "math.erfc (C library)": COLORS.sky };
const accMode = seg($("acc-mode"), drawAcc);
function drawAcc() {
  const mode = accMode.get();
  const zs = mode === "rel" ? EX2.curve_z.filter((z) => z <= 0) : EX2.curve_z;
  const curves = mode === "rel" ? EX2.curves_rel : EX2.curves_abs;
  const series = Object.entries(curves).map(([name, ys]) => ({ type: "line", x: zs, y: ys.map((v) => Math.min(Math.max(v, 1e-17), 1e3)), color: ACC_COL[name], width: name.startsWith("series") ? 2.5 : 1.6 }));
  plot($("acc-chart"), { x: [zs[0], zs[zs.length - 1]], y: [1e-17, 1e3], ylog: true, xlabel: "z", ylabel: mode === "rel" ? "relative error (z ≤ 0)" : "absolute error", series });
  $("acc-legend").innerHTML = Object.keys(curves).map((n) => `<span><i style="background:${ACC_COL[n]}"></i>${n}</span>`).join("");
}
(function accTable() {
  let t = `<tr><th>method</th><th class="num">max abs. error, |z| ≤ 10</th><th class="num">max rel. error, −5 ≤ z ≤ 0</th><th class="num">max rel. error, −10 ≤ z ≤ 0</th></tr>`;
  for (const [n, a] of Object.entries(EX2.accuracy))
    t += `<tr${n.startsWith("series") ? ' class="best"' : ""}><td>${n}</td><td class="num">${sci(a.max_abs)}</td><td class="num">${sci(a.max_rel_5)}</td><td class="num">${sci(a.max_rel_10)}</td></tr>`;
  $("acc-table").innerHTML = t;
  let w = `<tr><th class="num">z</th><th>source</th><th class="num">value on the web</th><th class="num">my method</th><th class="num">agrees to the printed digits</th></tr>`;
  for (const r of EX2.web) {
    const isTail = r.src === "tail table";
    const webTxt = isTail ? `1 − Φ = ${r.web.toExponential(5)}` : r.src.startsWith("68") ? r.web.toFixed(15) : r.web.toFixed(5);
    const mineTxt = isTail ? `1 − Φ = ${r.tail_mine.toExponential(8)}` : r.mine.toFixed(16);
    let ok;
    if (isTail) ok = r.tail_mine.toPrecision(6) === r.web.toPrecision(6);
    else if (r.src.startsWith("68")) ok = Math.abs(r.mine - r.web) < 6e-16;
    else ok = r.mine.toFixed(5) === r.web.toFixed(5);
    const src = r.src === "z-table" ? "Wikipedia z-table" : isTail ? "Wikipedia tail table" : "Wikipedia 68-95-99.7";
    w += `<tr><td class="num">${r.z}</td><td>${src}</td><td class="num">${webTxt}</td><td class="num">${mineTxt}</td><td class="num">${ok ? '<span class="ok">✓</span>' : '<span class="bad">✗</span>'}</td></tr>`;
  }
  $("web-table").innerHTML = w;
})();

// =================================================== 2.5 speed
(function speedTable() {
  let t = `<tr><th>method (Python)</th><th class="num">µs per evaluation</th></tr>`;
  for (const [n, us] of Object.entries(EX2.speed)) t += `<tr${n.startsWith("series") ? ' class="best"' : ""}><td>${n}</td><td class="num">${us.toFixed(2)}</td></tr>`;
  $("sp-table").innerHTML = t;
})();
$("sp-run").onclick = () => {
  $("sp-run").disabled = true;
  $("sp-js").innerHTML = `<tr><td>running…</td></tr>`;
  setTimeout(() => {
    const Z = Float64Array.from({ length: 200000 }, () => -6 + 12 * Math.random());
    let t = `<tr><th>method (JavaScript, this browser)</th><th class="num">ns per evaluation</th></tr>`;
    let sink = 0;
    for (const [k, M] of Object.entries(METHODS)) {
      const n = k === "trap" || k === "simp" ? 20000 : Z.length;
      const t0 = performance.now();
      for (let i = 0; i < n; i++) sink += M.f(Z[i]);
      const ns = ((performance.now() - t0) / n) * 1e6;
      t += `<tr${k === "final" ? ' class="best"' : ""}><td>${M.name}</td><td class="num">${ns.toFixed(0)}</td></tr>`;
    }
    $("sp-js").innerHTML = t + (sink === 42 ? "" : "");
    $("sp-run").disabled = false;
  }, 30);
};

// =================================================== init
drawDemo();
drawResults();
drawDirichlet();
drawThomae();
drawCalc();
drawAcc();

let rt;
window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(() => { drawDemo(); drawResults(); drawCalc(); drawAcc(); }, 150); });

const links = [...document.querySelectorAll("nav.toc a")];
const obs = new IntersectionObserver((entries) => {
  entries.forEach((e) => { if (e.isIntersecting) links.forEach((a) => a.classList.toggle("active", a.getAttribute("href") === "#" + e.target.id)); });
}, { rootMargin: "-40% 0px -55% 0px" });
document.querySelectorAll("main section").forEach((s) => obs.observe(s));
