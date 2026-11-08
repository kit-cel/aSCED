"""Builds the interactive HTML page for the beam-search IR-VN-splitting
waterfall sweep (sweeps/beam_search_ir_vn_waterfall.py's output), following
sweeps/PLOTTING_STYLE.md's data format / colour tokens / trace conventions.
"""
import json
from pathlib import Path

HERE = Path(__file__).parent
DATA_PATH = HERE / "beam_search_ir_vn_data.json"
OUT_PATH = HERE.parent / "beam_search_ir_vn_waterfall.html"

data = json.loads(DATA_PATH.read_text())
setting = data["beam_search_ir_vn_zc11_r4n132"]
split_histograms = data["split_histograms"]
meta = data["meta"]

DATA_JS = json.dumps({"setting": setting, "splitHistograms": split_histograms, "meta": meta})

HTML = """<title>Beam-search IR-VN Splitting</title>
<style>
:root{--bg:#f6f7f9;--panel:#fff;--ink:#14171c;--ink-2:#4a505a;--ink-3:#737a85;--rule:#dfe2e7;--grid:#e7e9ed;
      --c-a:#2a78d6;--c-b:#eb6834;--c-base:#1baf7a;--c-extra:#4a3aa7;--c-ref:#14171c}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){color-scheme:dark;--bg:#121417;--panel:#1a1d21;
      --ink:#eef0f3;--ink-2:#b4bac3;--ink-3:#8a909a;--rule:#2c3036;--grid:#262a30;
      --c-a:#3987e5;--c-b:#d95926;--c-base:#199e70;--c-extra:#9085e9;--c-ref:#eef0f3}}
:root[data-theme="dark"]{color-scheme:dark;--bg:#121417;--panel:#1a1d21;
      --ink:#eef0f3;--ink-2:#b4bac3;--ink-3:#8a909a;--rule:#2c3036;--grid:#262a30;
      --c-a:#3987e5;--c-b:#d95926;--c-base:#199e70;--c-extra:#9085e9;--c-ref:#eef0f3}

*{box-sizing:border-box}
body{background:var(--bg);color:var(--ink);margin:0;font-family:system-ui,sans-serif;padding:0 16px 48px;line-height:1.5}
.wrap{max-width:1080px;margin:0 auto}
h1{font-size:1.5rem;font-weight:600;margin:28px 0 4px;letter-spacing:-.01em}
h2{font-size:1.1rem;font-weight:600;margin:0 0 6px}
.lede{color:var(--ink-2);max-width:72ch;margin:0 0 20px;font-size:.95rem}
.mono{font-family:ui-monospace,'SF Mono',Consolas,monospace}
.card{background:var(--panel);border:1px solid var(--rule);border-radius:10px;padding:20px 22px;margin-bottom:24px}
.note{color:var(--ink-3);font-size:.85rem;margin-top:10px}
.plot{width:100%;height:560px}
.hist-plot{width:100%;height:420px}
table{border-collapse:collapse;width:100%;font-size:.88rem;margin-top:8px}
th,td{text-align:right;padding:4px 8px;border-bottom:1px solid var(--rule)}
th:first-child,td:first-child{text-align:left}
th{color:var(--ink-2);font-weight:600}
.small{font-size:.8rem;color:var(--ink-3)}
</style>
<div class="wrap">
  <h1>Beam-search decoder: IR-VN splitting, first waterfall comparison</h1>
  <p class="lede">
    aSCED-48 TCOM-paper code family (5G BG2, Z<sub>c</sub>=11, k=66, base n=154).
    Beam_Search_Decoder on an IR-VN-extended PCM (base H + 2 appended HARQ
    candidate-row blocks, each with its own new degree-1 IR-VN column,
    n_ext=176) with <b>all</b> columns eligible for splitting, vs. a fresh
    NMSA-32 baseline (base H, no extension) and the project's existing
    aSCED-48 reference curve (reused, not rerun). Click a legend entry to
    hide/show a curve.
  </p>
  <div class="card">
    <h2>FER vs. E<sub>b</sub>/N<sub>0</sub></h2>
    <div class="plot" id="err-plot"></div>
    <p class="note">Hollow markers = fewer than 300 frame errors at that point (preliminary budget: target_errors=200,
      max_transmissions capped per SNR point -- see page note below). aSCED-48 reference (black, dash-dot) is the
      project's existing target_errors=1000 run, NOT rerun here.</p>
  </div>
  <div class="card">
    <h2>Split-location histogram (aggregate across all SNR points)</h2>
    <div class="hist-plot" id="hist-plot"></div>
    <table id="hist-table"></table>
    <p class="note">Bucketed by VN column index: punctured state VNs = the leading 2&times;Z<sub>c</sub>=22 standard 5G
      systematic-puncturing columns; other state VNs = the remaining, normally-transmitted base-code columns;
      IR VNs = the 22 newly appended degree-1 columns.</p>
  </div>
  <div class="card">
    <p class="small mono" id="meta-note"></p>
  </div>
</div>
<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>
<script>
const DATA = __DATA_JS__;
const tok = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

function errTraces(d) {
  const col = {"Beam search": tok('--c-a'), "Baseline": tok('--c-base')};
  const dashes = ["solid", "dash", "dot", "dashdot"], syms = ["circle", "square", "diamond", "triangle-up"], cnt = {};
  const tr = d.curves.map(c => {
    const i = (cnt[c.group] = (cnt[c.group] || 0) + 1) - 1;
    return {
      name: c.name, x: c.pts.map(p => p[0]), y: c.pts.map(p => p[1]), mode: 'lines+markers',
      line: {color: col[c.group], width: 2, dash: dashes[i % 4]},
      marker: {size: 7, symbol: syms[i % 4],
               color: c.pts.map(p => p[2] < 300 ? 'rgba(0,0,0,0)' : col[c.group]),
               line: {width: 1.6, color: col[c.group]}},
      customdata: c.pts.map(p => p[2]),
      hovertemplate: '%{y:.3e} · %{customdata} err<extra>' + c.name + '</extra>',
    };
  });
  const R = d.reference;
  tr.push({name: R.name, x: R.pts.map(p => p[0]), y: R.pts.map(p => p[1]), mode: 'lines+markers',
           line: {color: tok('--c-ref'), dash: 'dashdot', width: 2}, marker: {size: 5, color: tok('--c-ref')},
           hovertemplate: '%{y:.3e}<extra>' + R.name + '</extra>'});
  return tr;
}

function layout(d) {
  const ax = {gridcolor: tok('--grid'), linecolor: tok('--rule'), zerolinecolor: tok('--rule'),
              tickfont: {color: tok('--ink-2'), size: 11}};
  return {paper_bgcolor: tok('--panel'), plot_bgcolor: tok('--panel'),
          font: {family: 'system-ui, sans-serif', color: tok('--ink')}, margin: {l: 60, r: 10, t: 8, b: 44},
          xaxis: Object.assign({title: {text: d.xlabel, font: {color: tok('--ink-2'), size: 12}}}, ax),
          yaxis: Object.assign({title: {text: d.ylabel, font: {color: tok('--ink-2'), size: 12}}, type: 'log',
                                 exponentformat: 'power'}, ax),
          legend: {orientation: 'h', y: -0.2, font: {size: 11, color: tok('--ink-2')}, bgcolor: 'rgba(0,0,0,0)'},
          hovermode: 'closest',
          hoverlabel: {bgcolor: tok('--panel'), bordercolor: tok('--rule'), font: {color: tok('--ink'), size: 12}}};
}
const cfg = {responsive: true, displaylogo: false, modeBarButtonsToRemove: ['lasso2d', 'select2d']};

function histTraces(sh) {
  const names = Object.keys(sh);
  const cats = ['punctured_state_vns', 'other_state_vns', 'ir_vns'];
  const labels = {'punctured_state_vns': 'Punctured state VNs (2·Zc)', 'other_state_vns': 'Other state VNs',
                  'ir_vns': 'IR VNs'};
  const colors = {'punctured_state_vns': tok('--c-b'), 'other_state_vns': tok('--c-base'), 'ir_vns': tok('--c-extra')};
  return cats.map(cat => ({
    name: labels[cat], type: 'bar', x: names, y: names.map(n => 100 * sh[n].aggregate[cat] / sh[n].aggregate.total),
    marker: {color: colors[cat]},
    hovertemplate: '%{y:.2f}%<extra>' + labels[cat] + '</extra>',
  }));
}

function histLayout() {
  const ax = {gridcolor: tok('--grid'), linecolor: tok('--rule'), tickfont: {color: tok('--ink-2'), size: 11}};
  return {barmode: 'stack', paper_bgcolor: tok('--panel'), plot_bgcolor: tok('--panel'),
          font: {family: 'system-ui, sans-serif', color: tok('--ink')}, margin: {l: 50, r: 10, t: 8, b: 60},
          xaxis: ax, yaxis: Object.assign({title: {text: '% of all splits', font: {color: tok('--ink-2'), size: 12}}}, ax),
          legend: {orientation: 'h', y: -0.25, font: {size: 11, color: tok('--ink-2')}, bgcolor: 'rgba(0,0,0,0)'}};
}

function drawAll() {
  Plotly.react('err-plot', errTraces(DATA.setting), layout(DATA.setting), cfg);
  Plotly.react('hist-plot', histTraces(DATA.splitHistograms), histLayout(), cfg);

  const rows = ['<tr><th>Config</th><th>Total splits</th><th>Punctured state</th><th>Other state</th><th>IR VNs</th></tr>'];
  Object.entries(DATA.splitHistograms).forEach(([name, h]) => {
    const a = h.aggregate, t = a.total;
    const pct = x => (100 * x / t).toFixed(1) + '%';
    rows.push(`<tr><td>${name}</td><td>${t}</td><td>${a.punctured_state_vns} (${pct(a.punctured_state_vns)})</td>` +
              `<td>${a.other_state_vns} (${pct(a.other_state_vns)})</td><td>${a.ir_vns} (${pct(a.ir_vns)})</td></tr>`);
  });
  document.getElementById('hist-table').innerHTML = rows.join('');

  const m = DATA.meta;
  document.getElementById('meta-note').textContent =
    `${m.note} Zc=${m.Zc}, base n=${m.number_vns}, n_ext=${m.n_ext}, Delta=${m.Delta}, k=${m.k}, ` +
    `num_ad_blocks=${m.num_ad_blocks}, SNR grid=${m.snr_grid.join(', ')} dB. ` +
    `Sweep wall-clock: ${m.wall_clock_total_seconds.toFixed(1)}s.`;
}
drawAll();
try { matchMedia('(prefers-color-scheme: dark)').addEventListener('change', drawAll) } catch (e) {}
new MutationObserver(drawAll).observe(document.documentElement, {attributes: true, attributeFilter: ['data-theme']});
</script>
"""

OUT_PATH.write_text(HTML.replace("__DATA_JS__", DATA_JS))
print(f"Wrote {OUT_PATH}")
