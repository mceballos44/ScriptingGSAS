# Builds a self-contained HTML report from the pipeline outputs:
# output/seq_results.csv, output/sample_info.csv and the two stats_*.csv files.
#
# Run from the repo root after main.py:
#     python -m scripts.report
# It runs the statistics first (scripts/stats.py) so the report always matches
# the current results, then writes output/report.html. Open it in a browser.

### Written by Mauricio Ceballos, Joester Group, Northwestern University

import json
from datetime import datetime
import numpy as np
import pandas as pd
from scripts.config import OUTPUT_DIR, CTE_SETTINGS, STATS_SETTINGS
from scripts import analysis, stats

REPORT_FILE = OUTPUT_DIR / "report.html"


def _scan_status(results, info):
    """Label each scan: used in the analysis, freezing scan, excluded scan, not refined"""
    freeze = dict(zip(info['Sample'], info['T_freeze']))
    excluded = CTE_SETTINGS.get('excluded_scans') or {}
    status = []
    for row in results.itertuples():
        if getattr(row, 'Refined', True) is False or pd.isna(row.A):
            status.append('not refined')
        elif CTE_SETTINGS.get('exclude_freeze_scan') and row.T == freeze.get(row.Sample):
            status.append('freezing scan')
        elif row.T in excluded.get(row.Sample, []):
            status.append('excluded scan')
        else:
            status.append('used')
    return status


def _records(df, columns, digits=6):
    """DataFrame -> list of dicts with NaN as None, rounded floats"""
    out = []
    for rec in df[columns].to_dict('records'):
        clean = {}
        for k, v in rec.items():
            if isinstance(v, (float, np.floating)):
                clean[k] = None if not np.isfinite(v) else round(float(v), digits)
            elif isinstance(v, (np.integer,)):
                clean[k] = int(v)
            elif isinstance(v, (np.bool_,)):
                clean[k] = bool(v)
            else:
                clean[k] = v
        out.append(clean)
    return out


def build_report(results_file, sample_info_file, per_sample, tests, out_file, run_label):
    results = pd.read_csv(results_file)
    info = pd.read_csv(sample_info_file)
    results['Status'] = _scan_status(results, info)

    data = {
        'run': run_label,
        'generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
        't_ref': STATS_SETTINGS['t_ref'],
        'excluded_scans': CTE_SETTINGS.get('excluded_scans') or {},
        'exclude_freeze_scan': bool(CTE_SETTINGS.get('exclude_freeze_scan')),
        'scans': _records(results, ['Sample', 'T', 'A', 'sigmaA', 'C', 'sigmaC',
                                    'Rwp', 'Status']),
        'samples': _records(per_sample, [
            'Sample', 'Group', 'T_start', 'T_freeze', 'N_points',
            'CTE_A', 'sigma_CTE_A', 'CTE_C', 'sigma_CTE_C',
            'A_Tref', 'sigma_A_Tref', 'C_Tref', 'sigma_C_Tref',
            'C_over_A_Tref', 'sigma_C_over_A_Tref', 'Tref_extrapolated']
            + [c for c in ['Spottiness'] if c in per_sample]),
        'tests': _records(tests, [
            'Quantity', 'Comparison', 'N_A', 'N_B', 'Mean_A', 'Mean_B',
            'Diff_A_minus_B', 'CI95_low', 'CI95_high', 'p_permutation',
            'p_mann_whitney', 'p_holm'], digits=8),
    }
    html = TEMPLATE.replace('__DATA__', json.dumps(data))
    out_file.write_text(html, encoding='utf-8')
    return out_file


TEMPLATE = r'''<title>Larval Ice Lattice Report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Sans+Condensed:wght@500;600&display=swap">
<style>
/* Layout: single reading column (~68ch) for text, figures run to 1080px; sections stacked by question */
:root {
  --bg: #f7f8f8;
  --surface: #fcfcfb;
  --ink: #101418;
  --ink-2: #4b5157;
  --muted: #7e858b;
  --rule: #dfe3e5;
  --axis: #c2c8cc;
  --accent: #1c5cab;
  --afp: #2a78d6;
  --wt: #eb6834;
  --good: #006300;
  --chip-bg: #e9eef2;
  --font-body: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --font-display: "IBM Plex Sans Condensed", "IBM Plex Sans", system-ui, sans-serif;
  --font-data: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #0f1113; --surface: #171a1c; --ink: #eef0f1; --ink-2: #c1c6ca; --muted: #8a9196;
    --rule: #2a2f33; --axis: #3a4045; --accent: #86b6ef; --afp: #3987e5; --wt: #d95926;
    --good: #0ca30c; --chip-bg: #22282c; color-scheme: dark;
  }
}
:root[data-theme="dark"] {
  --bg: #0f1113; --surface: #171a1c; --ink: #eef0f1; --ink-2: #c1c6ca; --muted: #8a9196;
  --rule: #2a2f33; --axis: #3a4045; --accent: #86b6ef; --afp: #3987e5; --wt: #d95926;
  --good: #0ca30c; --chip-bg: #22282c; color-scheme: dark;
}
body { background: var(--bg); color: var(--ink); font-family: var(--font-body); font-size: 15px; line-height: 1.55; }
.wrap { max-width: 1080px; margin: 0 auto; padding-inline: 20px; padding-block: 40px 64px; display: grid; gap: 48px; }
header { display: grid; gap: 10px; }
.eyebrow { font-family: var(--font-data); font-size: 12px; letter-spacing: .06em; text-transform: uppercase; color: var(--muted); }
h1 { font-family: var(--font-display); font-weight: 600; font-size: clamp(28px, 4.2vw, 40px); line-height: 1.1; margin: 0; text-wrap: balance; }
h2 { font-family: var(--font-display); font-weight: 600; font-size: 22px; margin: 0; text-wrap: balance; }
h3 { font-size: 15px; font-weight: 600; margin: 0; }
p { margin: 0; max-width: 68ch; }
.lede { color: var(--ink-2); font-size: 17px; }
section { display: grid; gap: 14px; min-width: 0; }
.note { color: var(--ink-2); font-size: 14px; max-width: 68ch; }
.figure { background: var(--surface); border: 1px solid var(--rule); border-radius: 6px; padding: 16px; min-width: 0; }
.legend { display: flex; flex-wrap: wrap; gap: 6px 18px; font-size: 13px; color: var(--ink-2); align-items: center; }
.legend span { display: inline-flex; align-items: center; gap: 6px; }
.legend svg { overflow: visible; }
.chart svg { display: block; width: 100%; height: auto; overflow: visible; }
.chart text { fill: var(--muted); font-family: var(--font-data); font-size: 11px; }
.chart .label { fill: var(--ink-2); font-family: var(--font-body); font-size: 12px; }
.chart .title { fill: var(--ink); font-family: var(--font-body); font-size: 13px; font-weight: 600; }
.grid line { stroke: var(--rule); }
.axis-line { stroke: var(--axis); }
.findings { display: grid; gap: 0; border-top: 1px solid var(--rule); }
.finding { display: grid; grid-template-columns: minmax(0, 1.4fr) minmax(0, 1fr) auto; gap: 6px 20px; padding: 12px 0; border-bottom: 1px solid var(--rule); align-items: baseline; }
.finding .what { font-weight: 500; }
.finding .what small { display: block; color: var(--muted); font-weight: 400; font-size: 13px; }
.finding .num { font-family: var(--font-data); font-size: 13px; color: var(--ink-2); font-variant-numeric: tabular-nums; }
.chip { font-family: var(--font-data); font-size: 11px; letter-spacing: .04em; text-transform: uppercase; padding: 3px 8px; border-radius: 999px; background: var(--chip-bg); color: var(--ink-2); white-space: nowrap; }
.chip.sig { background: var(--good); color: var(--surface); }
.chip.sugg { outline: 1px solid var(--accent); color: var(--accent); background: transparent; }
@media (max-width: 640px) { .finding { grid-template-columns: 1fr; } }
.table-wrap { overflow-x: auto; border: 1px solid var(--rule); border-radius: 6px; background: var(--surface); }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { padding: 7px 10px; text-align: right; white-space: nowrap; border-bottom: 1px solid var(--rule); }
th { font-weight: 600; color: var(--ink-2); background: var(--surface); position: sticky; top: 0; }
th:first-child, td:first-child, th.l, td.l { text-align: left; }
td { font-family: var(--font-data); font-variant-numeric: tabular-nums; }
td.l { font-family: var(--font-body); }
tr:last-child td { border-bottom: none; }
.swatch { width: 10px; height: 10px; border-radius: 50%; display: inline-block; }
.methods { display: grid; gap: 10px; }
.methods ul { margin: 0; padding-left: 18px; max-width: 72ch; display: grid; gap: 6px; color: var(--ink-2); }
code { font-family: var(--font-data); font-size: 13px; }
#tip { position: fixed; pointer-events: none; background: var(--surface); color: var(--ink); border: 1px solid var(--axis); border-radius: 6px; padding: 8px 10px; font-size: 12px; line-height: 1.4; box-shadow: 0 4px 16px rgba(0,0,0,.12); max-width: 260px; z-index: 10; }
#tip b { font-weight: 600; }
#tip .m { font-family: var(--font-data); }
.two { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 16px; }
:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
</style>

<div class="wrap">
  <header>
    <div class="eyebrow" id="eyebrow"></div>
    <h1>Ice in AFP and wild-type larvae</h1>
    <p class="lede" id="lede"></p>
  </header>

  <section id="findings-sec">
    <h2>Findings</h2>
    <p class="note">One value per larva. Differences are first group minus second, with a bootstrap 95% interval. "Holm" is the p-value corrected for running all of these tests together.</p>
    <div class="findings" id="findings"></div>
  </section>

  <section>
    <h2>Freezing temperature</h2>
    <p class="note">Warmest scan with ice for each larva, by group and start temperature. Each run cools in 1 K steps, so values are known to within 1 K. Bars show the mean of each column.</p>
    <div class="figure"><div class="legend" id="legend-freeze"></div><div class="chart" id="chart-freeze"></div></div>
  </section>

  <section>
    <h2>Expansion and lattice per larva</h2>
    <p class="note" id="note-cte"></p>
    <div class="figure"><div class="legend" id="legend-cte"></div><div class="chart" id="chart-cte"></div></div>
  </section>

  <section>
    <h2>Lattice parameters against temperature</h2>
    <p class="note">Every scan of every larva. Filled points were used in the fits; hollow points are freezing scans or scans excluded after inspection. Absolute a and c carry a per-larva offset from sample position (see Methods); slopes and c/a do not.</p>
    <div class="figure"><div class="legend" id="legend-lattice"></div>
      <div class="two"><div class="chart" id="chart-a"></div><div class="chart" id="chart-c"></div></div></div>
  </section>

  <section>
    <h2>Fit quality</h2>
    <p class="note">Weighted profile R-factor for each scan. Compare within a larva rather than as an absolute value. The excluded scans stand out at 17&ndash;25%.</p>
    <div class="figure"><div class="legend" id="legend-rwp"></div><div class="chart" id="chart-rwp"></div></div>
  </section>

  <section>
    <h2>All tests</h2>
    <div class="table-wrap"><table id="table-tests"></table></div>
  </section>

  <section>
    <h2>Per-larva values</h2>
    <div class="table-wrap"><table id="table-samples"></table></div>
  </section>

  <section class="methods">
    <h2>Methods and caveats</h2>
    <ul id="methods"></ul>
  </section>
</div>
<div id="tip" hidden></div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
<script>
const DATA = __DATA__;
const TREF = DATA.t_ref;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const groupColor = g => css(g === 'AFP' ? '--afp' : '--wt');
const fmt = (v, d) => v == null ? '–' : Number(v).toFixed(d);
const fmtP = p => p == null ? '–' : (p < 0.001 ? '<0.001' : Number(p).toFixed(3));
const sampleNum = s => parseInt(s.replace(/\D/g, ''), 10);
const samples = DATA.samples.slice().sort((a, b) =>
  a.Group.localeCompare(b.Group) || sampleNum(a.Sample) - sampleNum(b.Sample));
const startOf = Object.fromEntries(samples.map(s => [s.Sample, s.T_start]));
const groupOf = Object.fromEntries(samples.map(s => [s.Sample, s.Group]));

// ---------- text ----------
document.getElementById('eyebrow').textContent = `Pipeline run ${DATA.run} · report ${DATA.generated}`;
const nA = samples.filter(s => s.Group === 'AFP').length, nW = samples.filter(s => s.Group === 'WT').length;
document.getElementById('lede').textContent =
  `Sequential Rietveld refinement of ice Ih in ${samples.length} live larvae (${nA} AFP, ${nW} WT), ` +
  `each cooled in 1 K steps from 260 K or 280 K. Coefficients of thermal expansion (CTE) of a and c, ` +
  `the lattice at ${TREF} K and freezing temperatures, compared between groups and start temperatures.`;
document.getElementById('note-cte').textContent =
  `CTE from a weighted straight-line fit of a and c against temperature (units 10⁻⁶ K⁻¹), and c/a read off that fit at ${TREF} K. ` +
  `Bars are ±1 standard uncertainty. Vertical lines mark each group's mean. c/a cancels any error that scales a and c together, so it is the most robust structural comparison.`;

function niceName(q) {
  return q.replace('T_freeze (K)', 'Freezing temperature').replace(/^CTE a \(1e-6\/K\)/, 'CTE of a (10⁻⁶ K⁻¹)')
    .replace(/^CTE c \(1e-6\/K\)/, 'CTE of c (10⁻⁶ K⁻¹)').replace(/\(A\)$/, '(Å)');
}
function verdict(t) {
  if (t.p_holm != null && t.p_holm < 0.05) return ['sig', 'Significant'];
  if (t.p_permutation < 0.1) return ['sugg', 'Suggestive'];
  return ['', 'No difference'];
}
function diffText(t) {
  const isT = t.Quantity.startsWith('T_freeze');
  const isCTE = t.Quantity.startsWith('CTE');
  const isSpot = t.Quantity.startsWith('Ring');
  const d = isT || isCTE ? 1 : isSpot ? 2 : (t.Quantity.startsWith('c/a') ? 5 : 4);
  const unit = isT ? ' K' : (isCTE || isSpot || t.Quantity.startsWith('c/a')) ? '' : ' Å';
  return `${t.Diff_A_minus_B >= 0 ? '+' : ''}${fmt(t.Diff_A_minus_B, d)}${unit} [${fmt(t.CI95_low, d)}, ${fmt(t.CI95_high, d)}]`;
}
const order = t => (verdict(t)[0] === 'sig' ? 0 : verdict(t)[0] === 'sugg' ? 1 : 2);
const fbox = document.getElementById('findings');
DATA.tests.slice().sort((a, b) => order(a) - order(b)).forEach(t => {
  const [cls, label] = verdict(t);
  const div = document.createElement('div');
  div.className = 'finding';
  div.innerHTML = `<div class="what">${niceName(t.Quantity)}<small>${t.Comparison}</small></div>
    <div class="num">${diffText(t)}<br>p ${fmtP(t.p_permutation)} · Holm ${fmtP(t.p_holm)}</div>
    <div><span class="chip ${cls}">${label}</span></div>`;
  fbox.appendChild(div);
});

// ---------- legends ----------
function marker(svg, x, y, group, start, r = 5) {
  return svg.append('circle').attr('cx', x).attr('cy', y).attr('r', r)
    .attr('fill', start === 280 ? css('--surface') : groupColor(group))
    .attr('stroke', groupColor(group)).attr('stroke-width', 2);
}
function legend(id, items) {
  const el = d3.select('#' + id).html('');
  items.forEach(it => {
    const span = el.append('span');
    const svg = span.append('svg').attr('width', 14).attr('height', 14);
    it.draw(svg);
    span.append('span').text(it.text);
  });
}
function standardLegend(id, withShapes = true, withStatus = false) {
  const items = [
    { text: 'AFP', draw: s => s.append('circle').attr('cx', 7).attr('cy', 7).attr('r', 5).attr('fill', groupColor('AFP')) },
    { text: 'WT', draw: s => s.append('circle').attr('cx', 7).attr('cy', 7).attr('r', 5).attr('fill', groupColor('WT')) },
  ];
  if (withShapes) items.push(
    { text: 'filled: start 260 K', draw: s => marker(s, 7, 7, 'AFP', 260).attr('fill', css('--muted')).attr('stroke', css('--muted')) },
    { text: 'hollow: start 280 K', draw: s => marker(s, 7, 7, 'AFP', 280).attr('stroke', css('--muted')) });
  if (withStatus) items.push(
    { text: 'used in fits', draw: s => s.append('circle').attr('cx', 7).attr('cy', 7).attr('r', 4).attr('fill', css('--muted')) },
    { text: 'freezing / excluded scan', draw: s => s.append('circle').attr('cx', 7).attr('cy', 7).attr('r', 4).attr('fill', css('--surface')).attr('stroke', css('--muted')).attr('stroke-width', 2) });
  legend(id, items);
}

// ---------- tooltip ----------
const tip = document.getElementById('tip');
function showTip(ev, html) {
  tip.innerHTML = html; tip.hidden = false;
  const pad = 14, w = tip.offsetWidth, h = tip.offsetHeight;
  let x = ev.clientX + pad, y = ev.clientY + pad;
  if (x + w > window.innerWidth - 8) x = ev.clientX - w - pad;
  if (y + h > window.innerHeight - 8) y = ev.clientY - h - pad;
  tip.style.left = x + 'px'; tip.style.top = y + 'px';
}
const hideTip = () => { tip.hidden = true; };

function frame(id, height, margin) {
  const el = document.getElementById(id);
  el.innerHTML = '';
  const width = Math.max(300, el.clientWidth);
  const svg = d3.select(el).append('svg').attr('viewBox', `0 0 ${width} ${height}`)
    .attr('role', 'img');
  return { svg, width, height, m: margin, iw: width - margin.l - margin.r, ih: height - margin.t - margin.b };
}
function yGrid(g, y, iw, ticks, f) {
  const gg = g.append('g').attr('class', 'grid');
  y.ticks(ticks).forEach(v => {
    gg.append('line').attr('x1', 0).attr('x2', iw).attr('y1', y(v)).attr('y2', y(v));
    g.append('text').attr('x', -8).attr('y', y(v)).attr('dy', '0.32em').attr('text-anchor', 'end').text(f(v));
  });
}

// ---------- freezing temperature strip ----------
function drawFreeze() {
  const F = frame('chart-freeze', 300, { t: 16, r: 16, b: 44, l: 48 });
  const g = F.svg.append('g').attr('transform', `translate(${F.m.l},${F.m.t})`);
  const cats = [['AFP', 260], ['AFP', 280], ['WT', 260], ['WT', 280]];
  const x = d3.scaleBand().domain(cats.map(c => c.join(' '))).range([0, F.iw]).padding(0.25);
  const vals = samples.map(s => s.T_freeze);
  const y = d3.scaleLinear().domain([d3.min(vals) - 1, d3.max(vals) + 1]).range([F.ih, 0]);
  yGrid(g, y, F.iw, 6, v => v + ' K');
  cats.forEach(([grp, st]) => {
    const key = grp + ' ' + st;
    const pts = samples.filter(s => s.Group === grp && s.T_start === st);
    const cx = x(key) + x.bandwidth() / 2;
    const counts = {};
    pts.forEach(p => {
      const k = p.T_freeze; const i = counts[k] = (counts[k] || 0) + 1;
      p._dx = ((i % 2 ? 1 : -1) * Math.floor(i / 2)) * 13;
    });
    const mean = d3.mean(pts, p => p.T_freeze);
    g.append('line').attr('x1', cx - x.bandwidth() / 2).attr('x2', cx + x.bandwidth() / 2)
      .attr('y1', y(mean)).attr('y2', y(mean)).attr('stroke', groupColor(grp)).attr('stroke-width', 2).attr('opacity', .55);
    pts.forEach(p => {
      marker(g, cx + p._dx, y(p.T_freeze), grp, st, 5.5)
        .on('mousemove', ev => showTip(ev, `<b>${p.Sample}</b><br>froze at <span class="m">${p.T_freeze} K</span><br>start ${p.T_start} K`))
        .on('mouseleave', hideTip);
    });
    g.append('text').attr('class', 'label').attr('x', cx).attr('y', F.ih + 20).attr('text-anchor', 'middle').text(`${grp}, start ${st} K`);
    g.append('text').attr('x', cx).attr('y', F.ih + 36).attr('text-anchor', 'middle').text(`n = ${pts.length}, mean ${fmt(mean, 1)} K`);
  });
  g.append('line').attr('class', 'axis-line').attr('x1', 0).attr('x2', F.iw).attr('y1', F.ih).attr('y2', F.ih);
}

// ---------- per-larva dot rows ----------
function drawPerSample() {
  const cols = [
    { key: 'CTE_A', sig: 'sigma_CTE_A', title: 'CTE of a (10⁻⁶ K⁻¹)', d: 1 },
    { key: 'CTE_C', sig: 'sigma_CTE_C', title: 'CTE of c (10⁻⁶ K⁻¹)', d: 1 },
    { key: 'C_over_A_Tref', sig: 'sigma_C_over_A_Tref', title: `c/a at ${TREF} K`, d: 5 },
  ];
  const el = document.getElementById('chart-cte');
  const narrow = el.clientWidth < 620;
  const rowH = 20, gapG = 14;
  const labelW = 64;
  const nCols = narrow ? 1 : 3;
  const colGap = 44;
  const rowsH = samples.length * rowH + gapG;
  const blockH = rowsH + 52;
  const height = narrow ? blockH * 3 : blockH;
  const F = frame('chart-cte', height, { t: 8, r: 12, b: 4, l: labelW });
  const colW = (F.iw - colGap * (nCols - 1)) / nCols;
  const yRow = i => i * rowH + (samples[i].Group === 'WT' ? gapG : 0) + rowH / 2 + 22;
  cols.forEach((c, ci) => {
    const ox = narrow ? 0 : ci * (colW + colGap);
    const oy = narrow ? ci * blockH : 0;
    const g = F.svg.append('g').attr('transform', `translate(${F.m.l + ox},${F.m.t + oy})`);
    const lo = d3.min(samples, s => s[c.key] - (s[c.sig] || 0));
    const hi = d3.max(samples, s => s[c.key] + (s[c.sig] || 0));
    const pad = (hi - lo) * 0.06;
    const x = d3.scaleLinear().domain([lo - pad, hi + pad]).nice(c.d > 2 ? 3 : 4).range([0, colW]);
    g.append('text').attr('class', 'title').attr('x', 0).attr('y', 8).text(c.title);
    const gg = g.append('g').attr('class', 'grid');
    const nt = c.d > 2 ? 3 : 4;
    x.ticks(nt).forEach(v => {
      gg.append('line').attr('x1', x(v)).attr('x2', x(v)).attr('y1', 18).attr('y2', rowsH + 22);
      g.append('text').attr('x', x(v)).attr('y', rowsH + 38).attr('text-anchor', 'middle').text(d3.format(c.d > 2 ? '.4f' : 'd')(v));
    });
    if (ci === 0 || narrow) samples.forEach((s, i) => {
      g.append('text').attr('class', 'label').attr('x', -10).attr('y', yRow(i)).attr('dy', '0.32em').attr('text-anchor', 'end').text(s.Sample);
    });
    ['AFP', 'WT'].forEach(grp => {
      const idx = samples.map((s, i) => [s, i]).filter(([s]) => s.Group === grp);
      const mean = d3.mean(idx, ([s]) => s[c.key]);
      g.append('line').attr('x1', x(mean)).attr('x2', x(mean))
        .attr('y1', yRow(idx[0][1]) - rowH / 2).attr('y2', yRow(idx[idx.length - 1][1]) + rowH / 2)
        .attr('stroke', groupColor(grp)).attr('stroke-width', 2).attr('stroke-dasharray', '4 3').attr('opacity', .7);
    });
    samples.forEach((s, i) => {
      const yy = yRow(i), col = groupColor(s.Group);
      if (s[c.sig] != null) g.append('line').attr('x1', x(s[c.key] - s[c.sig])).attr('x2', x(s[c.key] + s[c.sig]))
        .attr('y1', yy).attr('y2', yy).attr('stroke', col).attr('stroke-width', 2).attr('stroke-linecap', 'round');
      marker(g, x(s[c.key]), yy, s.Group, s.T_start, 4.5);
      g.append('rect').attr('x', 0).attr('y', yy - rowH / 2).attr('width', colW).attr('height', rowH).attr('fill', 'transparent')
        .on('mousemove', ev => showTip(ev, `<b>${s.Sample}</b> (${s.Group}, start ${s.T_start} K)<br>${c.title}: <span class="m">${fmt(s[c.key], c.d)} ± ${fmt(s[c.sig], c.d)}</span><br>${s.N_points} scans in fit${s.Tref_extrapolated && c.key === 'C_over_A_Tref' ? `<br>${TREF} K is outside this larva's range (extrapolated)` : ''}`))
        .on('mouseleave', hideTip);
    });
  });
}

// ---------- lattice / Rwp against temperature ----------
function drawVsT(id, key, sigKey, yLabel, digits) {
  const F = frame(id, 300, { t: 26, r: 14, b: 34, l: 58 });
  const g = F.svg.append('g').attr('transform', `translate(${F.m.l},${F.m.t})`);
  const pts = DATA.scans.filter(d => d[key] != null && groupOf[d.Sample]);
  const x = d3.scaleLinear().domain(d3.extent(pts, d => d.T)).nice().range([0, F.iw]);
  const y = d3.scaleLinear().domain(d3.extent(pts, d => d[key])).nice().range([F.ih, 0]);
  g.append('text').attr('class', 'title').attr('x', -F.m.l + 4).attr('y', -10).text(yLabel);
  yGrid(g, y, F.iw, 5, v => d3.format(`.${digits}f`)(v));
  x.ticks(6).forEach(v => g.append('text').attr('x', x(v)).attr('y', F.ih + 18).attr('text-anchor', 'middle').text(v + ' K'));
  g.append('line').attr('class', 'axis-line').attr('x1', 0).attr('x2', F.iw).attr('y1', F.ih).attr('y2', F.ih);
  const line = d3.line().x(d => x(d.T)).y(d => y(d[key]));
  d3.groups(pts, d => d.Sample).forEach(([sample, rows]) => {
    rows.sort((a, b) => a.T - b.T);
    const col = groupColor(groupOf[sample]);
    const used = rows.filter(r => r.Status === 'used');
    g.append('path').attr('d', line(used)).attr('fill', 'none').attr('stroke', col).attr('stroke-width', 1.5).attr('opacity', .55);
    rows.forEach(r => {
      const isUsed = r.Status === 'used';
      g.append('circle').attr('cx', x(r.T)).attr('cy', y(r[key])).attr('r', isUsed ? 3 : 4)
        .attr('fill', isUsed ? col : css('--surface')).attr('stroke', col).attr('stroke-width', isUsed ? 0 : 2);
      g.append('circle').attr('cx', x(r.T)).attr('cy', y(r[key])).attr('r', 8).attr('fill', 'transparent')
        .on('mousemove', ev => showTip(ev, `<b>${sample}</b> at ${r.T} K<br>${yLabel}: <span class="m">${fmt(r[key], digits + 1)}${sigKey && r[sigKey] != null ? ' ± ' + fmt(r[sigKey], digits + 1) : ''}</span><br>${r.Status}`))
        .on('mouseleave', hideTip);
    });
  });
}

// ---------- tables ----------
function table(id, headers, rows) {
  const t = d3.select('#' + id).html('');
  t.append('thead').append('tr').selectAll('th').data(headers).join('th').attr('class', h => h[2] || null).text(h => h[0]);
  const tb = t.append('tbody');
  rows.forEach(r => {
    const tr = tb.append('tr');
    headers.forEach(h => {
      const td = tr.append('td').attr('class', h[2] || null);
      const v = h[1](r);
      if (v instanceof Node) td.node().appendChild(v); else td.text(v);
    });
  });
}
function groupCell(g) {
  const span = document.createElement('span');
  span.innerHTML = `<span class="swatch" style="background:${groupColor(g)}"></span> ${g}`;
  return span;
}
function drawTables() {
  table('table-tests', [
    ['Quantity', t => niceName(t.Quantity), 'l'], ['Comparison', t => t.Comparison, 'l'],
    ['n', t => `${t.N_A} / ${t.N_B}`], ['Mean', t => `${fmt(t.Mean_A, 4)} / ${fmt(t.Mean_B, 4)}`],
    ['Difference [95% CI]', diffText], ['p (perm.)', t => fmtP(t.p_permutation)],
    ['p (Mann-Whitney)', t => fmtP(t.p_mann_whitney)], ['p (Holm)', t => fmtP(t.p_holm)],
  ], DATA.tests);
  table('table-samples', [
    ['Larva', s => s.Sample, 'l'], ['Group', s => groupCell(s.Group), 'l'], ['Start', s => s.T_start + ' K'],
    ['Froze', s => s.T_freeze + ' K'], ['Scans', s => s.N_points],
    ['CTE a', s => `${fmt(s.CTE_A, 1)} ± ${fmt(s.sigma_CTE_A, 1)}`], ['CTE c', s => `${fmt(s.CTE_C, 1)} ± ${fmt(s.sigma_CTE_C, 1)}`],
    [`a at ${TREF} K (Å)`, s => fmt(s.A_Tref, 4)], [`c at ${TREF} K (Å)`, s => fmt(s.C_Tref, 4)],
    [`c/a at ${TREF} K`, s => fmt(s.C_over_A_Tref, 5) + (s.Tref_extrapolated ? ' *' : '')],
  ].concat(samples.some(s => s.Spottiness != null) ? [['Spottiness', s => fmt(s.Spottiness, 2)]] : []), samples);
}

// ---------- methods ----------
const ex = Object.entries(DATA.excluded_scans).map(([s, ts]) => `${s} ${ts.join(', ')} K`).join('; ');
[
  'GSAS-II sequential Rietveld refinement per larva: the coldest scan is refined first (background, sample displacement, cell, uniaxial microstrain), then every scan in turn, each starting from the previous result. Lattice parameters come from the Dij terms.',
  'Sample displacement is refined once per larva and then held fixed; refining it at every temperature traded off against a and c and doubled their scatter. Zero shift is fixed at 0.',
  (DATA.exclude_freeze_scan ? 'The freezing scan of each larva (partly liquid) is left out of the fits. ' : '') + (ex ? `Scans excluded after inspection: ${ex}.` : ''),
  'Statistics use one value per larva. Permutation tests on the difference in means; AFP vs WT labels are shuffled only within the same start temperature. Bootstrap 95% intervals; Mann-Whitney as a rank-based check; Holm correction across all tests.',
  'Start temperature coincides with collection batch (larvae 2–8 started at 260 K, 15–24 at 280 K), so a batch or cohort difference would look the same as a start-temperature effect.',
  'The LaB6 calibration was collected separately and each larva sits at a slightly different position, so absolute a and c carry a per-larva scale offset. CTE (a relative slope) and c/a are not affected by it.',
  `* ${TREF} K lies outside the measured range of this larva and is extrapolated from its fit.`,
].forEach(t => { if (t) { const li = document.createElement('li'); li.textContent = t; document.getElementById('methods').appendChild(li); } });

// ---------- render ----------
function render() {
  standardLegend('legend-freeze');
  standardLegend('legend-cte');
  standardLegend('legend-lattice', false, true);
  standardLegend('legend-rwp', false, true);
  drawFreeze();
  drawPerSample();
  drawVsT('chart-a', 'A', 'sigmaA', 'a (Å)', 3);
  drawVsT('chart-c', 'C', 'sigmaC', 'c (Å)', 3);
  drawVsT('chart-rwp', 'Rwp', null, 'Rwp (%)', 0);
  drawTables();
}
render();
let rt; window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(render, 150); });
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', render);
new MutationObserver(render).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
</script>
'''


if __name__ == "__main__":
    per_sample, tests = stats.run_stats(**STATS_SETTINGS)
    per_sample.to_csv(stats.PER_SAMPLE_FILE, index=False)
    tests.to_csv(stats.TESTS_FILE, index=False)
    run_time = datetime.fromtimestamp(analysis.RESULTS_FILE.stat().st_mtime)
    out = build_report(
        analysis.RESULTS_FILE, analysis.SAMPLE_INFO_FILE, per_sample, tests,
        REPORT_FILE, run_time.strftime('%Y-%m-%d %H:%M')
    )
    print(f"Report written to {out}")
