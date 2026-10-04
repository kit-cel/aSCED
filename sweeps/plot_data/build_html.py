import json, csv
from pathlib import Path

OUT = Path(__file__).parent

def csv_rows(name):
    return list(csv.DictReader(open(OUT / name)))

sweepA = json.load(open(OUT / "sweepA_data.json"))
groupsize = json.load(open(OUT / "sweepB_groupsize.json"))
ordering = json.load(open(OUT / "sweepB_ordering.json"))
raw = json.load(open(OUT / "raw_numbers.json"))
greedy = json.load(open(OUT / "greedy_data.json"))
greedy_raw = csv_rows("greedy_raw.csv")
pcmfirst = json.load(open(OUT / "pcmfirst_data.json"))
pcmfirst_raw = csv_rows("pcmfirst_raw.csv")
scheduling = json.load(open(OUT / "scheduling_data.json"))
scheduling_raw = csv_rows("scheduling_raw.csv")
alpha = json.load(open(OUT / "alpha_data.json"))
alpha_raw = csv_rows("alpha_raw.csv")

DATA_JS = json.dumps({
    "sweepA": sweepA,
    "groupsize": groupsize,
    "ordering": ordering,
    "raw": raw,
    "greedy": greedy,
    "greedyRaw": greedy_raw,
    "pcmfirst": pcmfirst,
    "pcmfirstRaw": pcmfirst_raw,
    "scheduling": scheduling,
    "schedulingRaw": scheduling_raw,
    "alpha": alpha,
    "alphaRaw": alpha_raw,
}, indent=None)

HTML = """<title>Sequential aSCED Decoding</title>
<style>
:root{--bg:#f6f7f9;--panel:#fff;--ink:#14171c;--ink-2:#4a505a;--ink-3:#737a85;--rule:#dfe2e7;--grid:#e7e9ed;
      --c-a:#2a78d6;--c-b:#eb6834;--c-base:#1baf7a;--c-extra:#4a3aa7;--c-ref:#14171c;
      --good:#1baf7a;--warn:#c98a1a}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){color-scheme:dark;--bg:#121417;--panel:#1a1d21;
      --ink:#eef0f3;--ink-2:#b4bac3;--ink-3:#8a909a;--rule:#2c3036;--grid:#262a30;
      --c-a:#3987e5;--c-b:#d95926;--c-base:#199e70;--c-extra:#9085e9;--c-ref:#eef0f3;
      --good:#2bd696;--warn:#e3a53a}}
:root[data-theme="dark"]{color-scheme:dark;--bg:#121417;--panel:#1a1d21;
      --ink:#eef0f3;--ink-2:#b4bac3;--ink-3:#8a909a;--rule:#2c3036;--grid:#262a30;
      --c-a:#3987e5;--c-b:#d95926;--c-base:#199e70;--c-extra:#9085e9;--c-ref:#eef0f3;
      --good:#2bd696;--warn:#e3a53a}

*{box-sizing:border-box}
body{background:var(--bg);color:var(--ink);margin:0;font-family:'IBM Plex Sans',system-ui,sans-serif;
     padding:0 16px 48px;line-height:1.5}
.wrap{max-width:1180px;margin:0 auto}
h1{font-size:1.5rem;font-weight:600;margin:28px 0 4px;letter-spacing:-.01em}
h2{font-size:1.15rem;font-weight:600;margin:0}
.lede{color:var(--ink-2);max-width:72ch;margin:0 0 28px;font-size:.95rem}
.mono{font-family:'IBM Plex Mono',monospace}
.card{background:var(--panel);border:1px solid var(--rule);border-radius:10px;padding:20px 22px;margin-bottom:28px}
.card > header{display:flex;align-items:baseline;gap:10px;flex-wrap:wrap;margin-bottom:6px}
.sub{color:var(--ink-3);font-size:.82rem;font-weight:400}
.pill{font-size:.72rem;padding:2px 9px;border-radius:999px;font-weight:600;letter-spacing:.02em}
.pill.good{background:color-mix(in srgb, var(--good) 18%, transparent);color:var(--good)}
.pill.warn{background:color-mix(in srgb, var(--warn) 18%, transparent);color:var(--warn)}
.verdict{color:var(--ink-2);font-size:.88rem;margin:4px 0 16px;max-width:85ch}
.pair{display:grid;grid-template-columns:1.35fr 1fr;gap:14px}
@media (max-width:860px){.pair{grid-template-columns:1fr}}
.plot{height:480px}
.quad{display:grid;grid-template-columns:1fr 1fr;gap:14px}
@media (max-width:760px){.quad{grid-template-columns:1fr}}
.bars{height:220px}
details{margin-top:14px}
summary{cursor:pointer;color:var(--ink-2);font-size:.85rem;font-weight:500}
summary:hover{color:var(--ink)}
table{border-collapse:collapse;width:100%;font-size:.78rem;margin-top:10px}
.tscroll{overflow-x:auto}
th,td{padding:4px 9px;text-align:right;border-bottom:1px solid var(--rule);white-space:nowrap;
      font-variant-numeric:tabular-nums}
th{color:var(--ink-3);font-weight:600;text-align:right;position:sticky;top:0;background:var(--panel)}
td:first-child,th:first-child{text-align:left}
td:nth-child(2),th:nth-child(2){text-align:left}
.toc{display:flex;gap:8px;flex-wrap:wrap;margin:10px 0 24px;position:sticky;top:8px;z-index:5}
.toc a{font-size:.78rem;color:var(--ink-2);background:var(--panel);border:1px solid var(--rule);
       border-radius:999px;padding:4px 11px;text-decoration:none}
.toc a:hover{color:var(--ink);border-color:var(--ink-3)}
.note{font-size:.78rem;color:var(--ink-3);max-width:90ch;margin-top:10px}
.legend-key{display:flex;gap:16px;flex-wrap:wrap;font-size:.78rem;color:var(--ink-2);margin:2px 0 14px}
.legend-key span{display:inline-flex;align-items:center;gap:5px}
.dot{width:9px;height:9px;border-radius:50%;display:inline-block}
</style>

<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>

<div class="wrap">
  <h1>Sequential &amp; m-converged aSCED decoding</h1>
  <p class="lede">
    aSCED-48 and aSCED-384 (RL / QC-block 5G LDPC, C<sub>5G</sub>(132,66), Z<sub>c</sub>=11) decoded with the new
    <span class="mono">DecoderSelector</span> + <span class="mono">MConvergedPolicy</span> group-sequential stopping
    mechanism, compared against full-parallel ensemble decoding and against Fig. 3 of the paper. Six studies in
    total: the original selector/group-size sweeps (A, B), plus four follow-ups &mdash; greedy RL/QC-block
    selection, PCM-first ensembles, BP scheduling/iteration count, and the MSA normalization constant (&alpha;).
    All curves below are interactive &mdash; click a legend entry to hide/show, double-click to isolate. Hollow
    markers mark points with &lt;300 collected frame errors.
  </p>

  <div class="toc">
    <a href="#sweepA">Sweep A &middot; selector comparison</a>
    <a href="#sweepB1">Sweep B &middot; latency/complexity vs. group size</a>
    <a href="#sweepB2">Sweep B &middot; ordering rule</a>
    <a href="#greedy">Greedy block selection</a>
    <a href="#pcmfirst">PCM-first (aSCED-49/385)</a>
    <a href="#scheduling">Scheduling &amp; iterations</a>
    <a href="#alpha">&alpha; (norm_factor) sweep</a>
    <a href="#raw">Raw numbers</a>
  </div>

  <section id="sweepA">
    <h2 style="margin-bottom:2px">Sweep A &mdash; selector comparison</h2>
    <p class="note" style="margin-top:2px;margin-bottom:18px">
      Job 532495, 70/70 complete, n_simul=132 (C<sub>5G</sub>(132,66) &mdash; the earlier n_simul=110 sweep was
      a different, higher-rate code and has been discarded). members_per_group=8 fixed; target_num_converged &isin; {2, 6}; full parallel = no
      grouping (one group, all paths). Complexity = <span class="mono">average_ensemble_effort</span> (total BP
      iterations per decoded frame, summed over the ensemble &mdash; the complexity axis). Literature curves from
      Fig. 3 are shown for context, hidden by default.
    </p>
    <div id="settings-sweepA"></div>
  </section>

  <section id="sweepB1">
    <h2 style="margin-bottom:2px">Sweep B, part 1 &mdash; best latency vs. best complexity</h2>
    <p class="note" style="margin-top:2px;margin-bottom:10px">
      Job 532496, 56/56 complete, n_simul=132. aSCED-48, syndrome_sequential, target_num_converged=6, members_per_group swept
      {1,2,4,8,16,24}; the mpg=48 endpoint reuses full-parallel (one group spanning the whole ensemble &mdash; no
      early stop is possible). Complexity = <span class="mono">average_ensemble_effort</span> (total BP iterations
      per frame, summed over the ensemble). Latency = <span class="mono">average_ensemble_latency</span> (sum of
      per-executed-group max effort &mdash; assumes full hardware parallelism <em>within</em> a group). Each curve
      point is the group size's own SNR interpolated to the stated FER (log-FER, log-metric interpolation between
      the two bracketing simulated SNR points &mdash; never extrapolated).
    </p>
    <div class="pair">
      <div class="plot" id="groupsize-plot"></div>
      <div id="groupsize-legend" style="padding-top:8px">
        <div class="legend-key">
          <span><span class="dot" style="background:var(--c-a)"></span>FER = 10<sup>&minus;1</sup></span>
          <span><span class="dot" style="background:var(--c-b)"></span>FER = 10<sup>&minus;3</sup></span>
        </div>
        <p class="note">Point color encodes group size (darker = larger groups = more hardware parallelism, less
        stopping granularity). Lower-left is better on both axes; neither endpoint dominates &mdash; full sequential
        (mpg=1) minimizes complexity, full parallel (mpg=48) minimizes latency.</p>
      </div>
    </div>
  </section>

  <section id="sweepB2">
    <h2 style="margin-bottom:2px">Sweep B, part 2 &mdash; ordering rule at fixed group size</h2>
    <p class="note" style="margin-top:2px;margin-bottom:10px">
      Same job, members_per_group=4 fixed, target_num_converged=6, selector &isin;
      {fixed_sequential, random_sequential, syndrome_sequential}. Same interpolation procedure as above, read off at
      the two FER operating points.
    </p>
    <div class="quad">
      <div><div class="sub" style="margin-bottom:4px">Complexity @ FER=10&#8315;&#185; &middot; BP it./frame</div><div class="bars" id="ord-eff-01"></div></div>
      <div><div class="sub" style="margin-bottom:4px">Latency @ FER=10&#8315;&#185; &middot; BP it./frame</div><div class="bars" id="ord-lat-01"></div></div>
      <div><div class="sub" style="margin-bottom:4px">Complexity @ FER=10&#8315;&#179; &middot; BP it./frame</div><div class="bars" id="ord-eff-001"></div></div>
      <div><div class="sub" style="margin-bottom:4px">Latency @ FER=10&#8315;&#179; &middot; BP it./frame</div><div class="bars" id="ord-lat-001"></div></div>
    </div>
    <p class="note">syndrome_sequential orders groups by estimated closeness to convergence (lowest syndrome weight
    first), so the paths most likely to converge run early &mdash; the lower bars above are the direct effect.</p>
  </section>

  <section id="greedy">
    <h2 style="margin-bottom:2px">Greedy RL/QC-block selection</h2>
    <p class="note" style="margin-top:2px;margin-bottom:18px">
      Job 535455, 70/70 complete, n_simul=132, members_per_group=8, target_num_converged &isin; {2, 6} (full parallel
      ignores target). <b>Baseline</b> = the reference script's arbitrary choice of "the next 2 (resp. 4) subsequent
      RL/QC blocks from offset 0". <b>Greedy</b> = blocks picked by round-by-round greedy forward selection over all
      34 candidates, ranked by combined FER rank at 1.8/3.3dB (see
      <span class="mono">sweeps/greedy_block_search_results.md</span>) &mdash; winning order aSCED-48: blocks
      [21, 16]; aSCED-384: blocks [31, 14, 17, 3]. Only block 16 appears in both lists. Complexity =
      <span class="mono">average_ensemble_effort</span>.
    </p>
    <p class="note" style="margin-bottom:18px"><b>Caveat (documented, not re-verified here):</b> 2 of the 34
      candidate blocks (indices 10 and 19) reproducibly segfault channel-code-lib2's BP decoder &mdash; both are
      degenerate check-node-degree-1 blocks (every row weight exactly 1). Both search orchestrators automatically
      excluded them from ranking, so they never corrupted the search; they were simply permanently out of
      contention. This is a C++-side bug, out of scope for this page (tracked separately).</p>
    <div id="settings-greedy"></div>
  </section>

  <section id="pcmfirst">
    <h2 style="margin-bottom:2px">PCM-first ensembles (aSCED-49/385)</h2>
    <p class="note" style="margin-top:2px;margin-bottom:10px">
      aSCED-49/385 = aSCED-48/384 plus one extra decoding path on the plain, unmodified PCM <span class="mono">H</span>
      (no auxiliary rows &mdash; a single path, not a batch), prepended at ensemble index 0 so it is scheduled first
      under <span class="mono">fixed_sequential</span>. Job 534062, 42/42 complete. Scoped to
      <span class="mono">full_parallel</span> + <span class="mono">fixed_sequential</span> only (pinning index 0
      first needs <span class="mono">members_per_group=7</span> so 49=7&times;7 / 385=7&times;55 divide evenly;
      <span class="mono">syndrome_sequential</span> reorders dynamically per word and was out of scope).
    </p>
    <p class="verdict" style="background:color-mix(in srgb, var(--warn) 12%, transparent);border:1px solid var(--warn);border-radius:8px;padding:10px 14px;margin-bottom:16px">
      <b>Confounded comparison &mdash; not yet corrected.</b> The PCM-first sweep necessarily used
      <span class="mono">members_per_group=7</span>, while the baseline below (plain aSCED-48/384) uses
      <span class="mono">members_per_group=8</span> &mdash; so the gap plotted here mixes <i>two</i> effects: (1)
      adding the pinned PCM-first path, and (2) shrinking the group size 8&rarr;7 (finer grouping alone tends to
      lower latency/effort). A de-confounding control (plain aSCED-48/384 at mpg=7, no PCM path &mdash; job 535646)
      was running at the time this page was built (16/28 tasks complete, 10 running, 1 pending; check
      <span class="mono">sacct -j 535646</span>). Once it lands, this section should be rebuilt to show the clean,
      isolated PCM-first effect (pcmfirst mpg=7 vs. mpg=7 control) separately from the pure group-size effect
      (mpg=7 control vs. mpg=8 baseline). Until then, treat the gap below as an upper bound on the PCM-first effect,
      not a clean attribution.
    </p>
    <div id="settings-pcmfirst"></div>
  </section>

  <section id="scheduling">
    <h2 style="margin-bottom:2px">BP scheduling &amp; iteration count</h2>
    <p class="note" style="margin-top:2px;margin-bottom:18px">
      Job 534114, 126/126 complete, n_simul=132, full_parallel only (one clean variable at a time, decoupled from the
      selector work above). <span class="mono">scheduling_type</span> &isin; {flooding, row_layered (natural order),
      row_layered (appended-rows block first)} &times; <span class="mono">max_iterations</span> &isin; {32, 16}, for
      3 variants: standalone NMSA (BP on <span class="mono">H</span> alone), aSCED-48, aSCED-384 (full parallel).
      Solid = 32 iterations, dashed = 16 iterations. <b>Headline:</b> row-layered scheduling at 16 iterations
      roughly matches or beats flooding at 32 iterations on both FER and effort &mdash; i.e. the same error-rate
      floor at about half the iteration budget.
    </p>
    <div id="settings-scheduling"></div>
  </section>

  <section id="alpha">
    <h2 style="margin-bottom:2px">MSA normalization constant (&alpha;) sweep</h2>
    <p class="note" style="margin-top:2px;margin-bottom:10px">
      Job 534481, 189/189 complete, n_simul=132. <span class="mono">BP_config.norm_factor</span> (&alpha;, used when
      <span class="mono">cn_update_type="msa"</span>) swept over 9 dyadic values in [0.5, 1.0], all else held at the
      established defaults (flooding, max_iterations=32). One line per SNR point (darker = higher SNR / lower FER);
      the dashed vertical guide marks &alpha;=0.75, the long-standing default. <b>Headline:</b> &alpha;=0.75 is at or
      very near the FER-minimizing value for all three variants &mdash; it wins outright at every SNR &ge;2.5dB for
      all three, and is within noise of the (slightly different) winner at the lower, noisier SNR points.
    </p>
    <div class="quad" id="settings-alpha"></div>
  </section>

  <section id="raw">
    <h2>Raw numbers</h2>
    <p class="note" style="margin-top:2px">CSV/JSON of everything below were sent alongside this page.</p>
    <details open>
      <summary>Sweep B operating points (interpolated) &mdash; group-size trade-off</summary>
      <div class="tscroll" id="tbl-groupsize"></div>
    </details>
    <details>
      <summary>Sweep B operating points (interpolated) &mdash; ordering comparison</summary>
      <div class="tscroll" id="tbl-ordering"></div>
    </details>
    <details>
      <summary>All 119 simulated (config &times; SNR) data points &mdash; Sweeps A/B</summary>
      <div class="tscroll" id="tbl-raw"></div>
    </details>
    <details>
      <summary>Greedy block selection &mdash; baseline vs. greedy, raw (config &times; SNR) points</summary>
      <div class="tscroll" id="tbl-greedy"></div>
    </details>
    <details>
      <summary>PCM-first &mdash; raw (config &times; SNR) points</summary>
      <div class="tscroll" id="tbl-pcmfirst"></div>
    </details>
    <details>
      <summary>Scheduling study &mdash; raw (config &times; SNR) points</summary>
      <div class="tscroll" id="tbl-scheduling"></div>
    </details>
    <details>
      <summary>&alpha; sweep &mdash; raw (config &times; SNR) points</summary>
      <div class="tscroll" id="tbl-alpha"></div>
    </details>
  </section>

  <p class="note">CPU-time-per-frame was not instrumented separately from decoding effort in these runs, in any of
  the six studies above; the "complexity"/"effort" axis is BP-iteration effort (the paper's own cost proxy), not
  wall-clock CPU seconds. All runs: n_simul=132 (C<sub>5G</sub>(132,66)), target_errors=200, max_transmissions=2e6
  per point, unless a section's own note states otherwise.</p>
</div>

<script>
const DATA = __DATA_JS__;
const tok = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const cfg = {responsive: true, displaylogo: false, modeBarButtonsToRemove: ['lasso2d', 'select2d']};

function axisStyle() {
  return {gridcolor: tok('--grid'), linecolor: tok('--rule'), zerolinecolor: tok('--rule'),
          tickfont: {color: tok('--ink-2'), family: 'IBM Plex Mono, monospace', size: 11}};
}
function baseLayout(xlabel, ylabel, xlog, ylog) {
  const ax = axisStyle();
  return {paper_bgcolor: tok('--panel'), plot_bgcolor: tok('--panel'),
          font: {family: 'IBM Plex Sans, system-ui, sans-serif', color: tok('--ink')},
          margin: {l: 60, r: 10, t: 8, b: 44},
          xaxis: Object.assign({title: {text: xlabel, font: {color: tok('--ink-2'), size: 12}}, type: xlog ? 'log' : 'linear'}, ax),
          yaxis: Object.assign({title: {text: ylabel, font: {color: tok('--ink-2'), size: 12}}, type: ylog ? 'log' : 'linear', exponentformat: 'power'}, ax),
          legend: {font: {size: 11, color: tok('--ink-2')}, bgcolor: 'rgba(0,0,0,0)'}, hovermode: 'closest',
          hoverlabel: {bgcolor: tok('--panel'), bordercolor: tok('--rule'), font: {color: tok('--ink'), family: 'IBM Plex Mono, monospace', size: 12}}};
}

const GROUP_COLOR = {"Baseline": () => tok('--c-base'), "Fixed order": () => tok('--c-a'), "Syndrome order": () => tok('--c-b'),
  "Baseline blocks": () => tok('--c-base'), "Greedy blocks": () => tok('--c-a'),
  "No PCM-first (mpg=8, baseline)": () => tok('--c-base'), "PCM-first (mpg=7)": () => tok('--c-a'),
  "flooding": () => tok('--c-base'), "row-layered (natural)": () => tok('--c-a'), "row-layered (appended-first)": () => tok('--c-b')};
const DASHES = ["solid", "dash", "dot", "dashdot", "longdash"];
const SYMS = ["circle", "square", "diamond", "triangle-up", "x"];

function curveTraces(d, cpu) {
  const cnt = {};
  return d.curves.map(c => {
    const i = (cnt[c.group] = (cnt[c.group] || 0) + 1) - 1, few = c.pts.length <= 2;
    const col = GROUP_COLOR[c.group] ? GROUP_COLOR[c.group]() : tok('--c-extra');
    const yv = cpu ? c.pts.map(p => p[3]) : c.pts.map(p => p[1]);
    return {
      name: c.name, legendgroup: c.name, showlegend: !cpu,
      x: c.pts.map(p => p[0]), y: yv,
      customdata: c.pts.map(p => p[2]),
      hovertemplate: (cpu ? '%{y:.4g} it/frame' : '%{y:.3e}') + ' &middot; %{customdata} err<extra>' + c.name + '</extra>',
      mode: few ? 'markers' : 'lines+markers',
      line: {color: col, width: 2, dash: DASHES[i % 5]},
      marker: {size: few ? 10 : 6, symbol: SYMS[i % 5],
               color: c.pts.map(p => p[2] < 300 ? 'rgba(0,0,0,0)' : col),
               line: {width: 1.6, color: col}},
      visible: c.visible ? true : 'legendonly'
    };
  });
}

function settingTraces(d, cpu) {
  const tr = curveTraces(d, cpu);
  if (!cpu) {
    const R = d.reference;
    if (R) {
      tr.push({name: R.name, x: R.pts.map(p => p[0]), y: R.pts.map(p => p[1]), mode: 'lines+markers',
               line: {color: tok('--c-ref'), dash: 'dashdot', width: 2}, marker: {size: 5, color: tok('--c-ref')},
               hovertemplate: '%{y:.3e}<extra>' + R.name + '</extra>'});
    }
    Object.entries(d.paper || {}).forEach(([n, s]) => tr.push({name: n + ' (paper)', x: s.map(p => p[0]), y: s.map(p => p[1]),
             mode: 'lines+markers', line: {color: tok('--ink-3'), dash: 'dot', width: 1.3}, marker: {size: 4, color: tok('--ink-3')},
             visible: 'legendonly'}));
  }
  return tr;
}

function drawSetting(id, d) {
  const errLayout = Object.assign(baseLayout(d.xlabel, d.ylabel, d.xlog, true), {height: 480,
    legend: {orientation: 'h', y: -0.22, font: {size: 10.5, color: tok('--ink-2')}, bgcolor: 'rgba(0,0,0,0)'}});
  const cpuLayout = Object.assign(baseLayout(d.xlabel, 'BP iterations / frame (effort)', d.xlog, true),
    {height: 480, showlegend: false});
  Plotly.react('err-' + id, settingTraces(d, false), errLayout, cfg);
  Plotly.react('cpu-' + id, settingTraces(d, true), cpuLayout, cfg);
  const el = document.getElementById('err-' + id);
  if (!el._synced) {
    el._synced = true;
    el.on('plotly_restyle', () => {
      const vis = {}; el.data.forEach(t => vis[t.name] = t.visible);
      const c = document.getElementById('cpu-' + id);
      Plotly.restyle(c, {visible: c.data.map(t => vis[t.name] === undefined ? true : vis[t.name])});
    });
  }
}

function buildSweepA() {
  const host = document.getElementById('settings-sweepA');
  Object.entries(DATA.sweepA).forEach(([id, d]) => {
    const K = id.replace('asced', '');
    const sec = document.createElement('div');
    sec.className = 'card';
    sec.innerHTML = `
      <header><h2>aSCED-${K} <span class="sub">RL 5G LDPC &middot; C<sub>5G</sub>(132,66), Z<sub>c</sub>=11</span></h2></header>
      <div class="pair">
        <div class="plot" id="err-${id}"></div>
        <div class="plot" id="cpu-${id}"></div>
      </div>`;
    host.appendChild(sec);
    drawSetting(id, d);
  });
}

// Generic FER+effort "pair" card renderer, reused by greedy / pcmfirst / scheduling
// (same shape as buildSweepA, parameterized over host id, data object and a title function).
function buildPairedSettings(hostId, dataObj, titleFor) {
  const host = document.getElementById(hostId);
  host.innerHTML = '';
  Object.entries(dataObj).forEach(([id, d]) => {
    const uid = hostId + '-' + id;
    const sec = document.createElement('div');
    sec.className = 'card';
    sec.innerHTML = `
      <header><h2>${titleFor(id)}</h2></header>
      <div class="pair">
        <div class="plot" id="err-${uid}"></div>
        <div class="plot" id="cpu-${uid}"></div>
      </div>`;
    host.appendChild(sec);
    drawSetting(uid, d);
  });
}

function buildGreedy() {
  buildPairedSettings('settings-greedy', DATA.greedy, id =>
    `aSCED-${id.replace('asced', '')} <span class="sub">baseline (subsequent blocks) vs. greedy-selected blocks</span>`);
}

function buildPcmfirst() {
  buildPairedSettings('settings-pcmfirst', DATA.pcmfirst, id =>
    `${id.replace('K', 'K=')} family <span class="sub">aSCED-${id === 'K48' ? '48 vs. 49' : '384 vs. 385'}</span>`);
}

function buildScheduling() {
  buildPairedSettings('settings-scheduling', DATA.scheduling, id => {
    const label = id === 'nmsa' ? 'Standalone NMSA' : `aSCED-${id.replace('asced', '')}`;
    return `${label} <span class="sub">full parallel</span>`;
  });
}

function buildGroupsize() {
  const mpgs = DATA.groupsize.mpg_values;
  const lo = Math.log2(Math.min(...mpgs)), hi = Math.log2(Math.max(...mpgs));
  const colorsFor = pts => pts.map(p => (Math.log2(p.mpg) - lo) / (hi - lo));
  const mk = (key, name, color, symbol) => {
    const pts = DATA.groupsize.operating_points[key].slice().sort((a, b) => a.mpg - b.mpg);
    return {
      name, x: pts.map(p => p.latency), y: pts.map(p => p.effort),
      mode: 'lines+markers', line: {color, width: 1.6, dash: 'dot'},
      marker: {size: 13, symbol, color: colorsFor(pts), colorscale: [[0, '#cfe0f5'], [1, color]],
               cmin: 0, cmax: 1, line: {width: 1.4, color}},
      text: pts.map(p => 'mpg=' + p.mpg),
      customdata: pts.map(p => [p.mpg, p.snr.toFixed(2), p.frame_errors]),
      hovertemplate: 'mpg=%{customdata[0]} &middot; SNR&approx;%{customdata[1]}dB<br>latency %{x:.4g} &middot; effort %{y:.4g}<extra>' + name + '</extra>'
    };
  };
  const traces = [mk('0.1', 'FER = 1e-1', tok('--c-a'), 'circle'), mk('0.001', 'FER = 1e-3', tok('--c-b'), 'diamond')];
  const layout = Object.assign(baseLayout('average ensemble latency [BP it./frame]', 'average ensemble effort [BP it./frame]', true, true),
    {height: 480, legend: {orientation: 'h', y: -0.18}});
  Plotly.react('groupsize-plot', traces, layout, cfg);
}

function buildOrderingBars(divId, key, metric, color0) {
  const pts = DATA.ordering.operating_points[key];
  const names = {"fixed_sequential": "fixed", "random_sequential": "random", "syndrome_sequential": "syndrome"};
  const colors = {"fixed_sequential": tok('--c-a'), "random_sequential": tok('--c-extra'), "syndrome_sequential": tok('--c-b')};
  const sorted = pts.slice().sort((a, b) => a[metric] - b[metric]);
  const layout = Object.assign(baseLayout(null, null, false, false), {height: 220, showlegend: false,
    margin: {l: 90, r: 30, t: 6, b: 30},
    xaxis: Object.assign({title: null}, axisStyle()),
    yaxis: Object.assign({title: null, automargin: true}, axisStyle())});
  const trace = {
    x: sorted.map(p => p[metric]), y: sorted.map(p => names[p.selector]), type: 'bar', orientation: 'h',
    marker: {color: sorted.map(p => colors[p.selector])},
    hovertemplate: '%{x:.4g} it/frame<extra></extra>'
  };
  Plotly.react(divId, [trace], layout, cfg);
}

function buildOrdering() {
  buildOrderingBars('ord-eff-01', '0.1', 'effort');
  buildOrderingBars('ord-lat-01', '0.1', 'latency');
  buildOrderingBars('ord-eff-001', '0.001', 'effort');
  buildOrderingBars('ord-lat-001', '0.001', 'latency');
}

// Parameter-sweep plot (PLOTTING_STYLE.md section 8): x = alpha (continuous decoder
// parameter), y = absolute FER (log), one line per channel point (SNR), sequential
// single-hue ramp (darker = higher SNR = lower FER = "better channel").
function buildAlpha() {
  const host = document.getElementById('settings-alpha');
  host.innerHTML = '';
  const names = {nmsa: 'Standalone NMSA', asced48: 'aSCED-48', asced384: 'aSCED-384'};
  Object.entries(DATA.alpha.variants).forEach(([variant, lines]) => {
    const div = document.createElement('div');
    div.innerHTML = `<div class="sub" style="margin-bottom:4px">${names[variant] || variant}</div><div class="plot" id="alpha-${variant}" style="height:400px"></div>`;
    host.appendChild(div);
    const snrVals = lines.map(l => l.snr);
    const lo = Math.min(...snrVals), hi = Math.max(...snrVals);
    const traces = lines.filter(l => l.pts.length).map(l => {
      const t = hi === lo ? 1 : (l.snr - lo) / (hi - lo);
      // light -> dark single-hue blue ramp, darker = higher SNR (lower FER, "better channel")
      const light = 82 - 54 * t; // lightness 82% (low SNR) -> 28% (high SNR)
      const color = `hsl(215, 70%, ${light}%)`;
      const few = l.pts.length <= 2;
      return {
        name: `SNR=${l.snr.toFixed(1)}dB`, x: l.pts.map(p => p[0]), y: l.pts.map(p => p[1]),
        customdata: l.pts.map(p => p[2]),
        hovertemplate: '%{y:.3e} &middot; %{customdata} err<extra>SNR=' + l.snr.toFixed(1) + 'dB</extra>',
        mode: few ? 'markers' : 'lines+markers',
        line: {color, width: 2},
        marker: {size: few ? 9 : 6, color: l.pts.map(p => p[2] < 300 ? 'rgba(0,0,0,0)' : color), line: {width: 1.6, color}}
      };
    });
    const layout = Object.assign(baseLayout('norm_factor (α)', 'FER', false, true), {
      height: 400, legend: {orientation: 'h', y: -0.24, font: {size: 10, color: tok('--ink-2')}, bgcolor: 'rgba(0,0,0,0)'},
      shapes: [{type: 'line', x0: 0.75, x1: 0.75, yref: 'paper', y0: 0, y1: 1,
                line: {color: tok('--ink-3'), dash: 'dash', width: 1.3}}],
      annotations: [{x: 0.75, y: 1, yref: 'paper', yanchor: 'bottom', showarrow: false, text: 'α=0.75 (default)',
                      font: {size: 10, color: tok('--ink-3')}}]
    });
    Plotly.react('alpha-' + variant, traces, layout, cfg);
  });
}

function table(rows, cols) {
  const head = '<tr>' + cols.map(c => `<th>${c.label}</th>`).join('') + '</tr>';
  const body = rows.map(r => '<tr>' + cols.map(c => `<td>${c.fmt ? c.fmt(r[c.key]) : (r[c.key] == null ? '&ndash;' : r[c.key])}</td>`).join('') + '</tr>').join('');
  return `<table><thead>${head}</thead><tbody>${body}</tbody></table>`;
}
const f3 = x => x == null ? '&ndash;' : (+x).toPrecision(4);
const fsnr = x => (+x).toFixed(2);

function buildTables() {
  document.getElementById('tbl-groupsize').innerHTML = table(
    [].concat(...Object.entries(DATA.groupsize.operating_points).map(([fer, pts]) => pts.map(p => Object.assign({fer}, p)))),
    [{key: 'fer', label: 'target FER'}, {key: 'mpg', label: 'members/group'}, {key: 'snr', label: 'SNR* [dB]', fmt: fsnr},
     {key: 'latency', label: 'latency', fmt: f3}, {key: 'effort', label: 'effort', fmt: f3}, {key: 'frame_errors', label: 'min FE (bracket)'}]);
  document.getElementById('tbl-ordering').innerHTML = table(
    [].concat(...Object.entries(DATA.ordering.operating_points).map(([fer, pts]) => pts.map(p => Object.assign({fer}, p)))),
    [{key: 'fer', label: 'target FER'}, {key: 'selector', label: 'selector'}, {key: 'snr', label: 'SNR* [dB]', fmt: fsnr},
     {key: 'latency', label: 'latency', fmt: f3}, {key: 'effort', label: 'effort', fmt: f3}, {key: 'frame_errors', label: 'min FE (bracket)'}]);
  document.getElementById('tbl-raw').innerHTML = table(DATA.raw,
    [{key: 'variant', label: 'variant'}, {key: 'selector', label: 'selector'}, {key: 'members_per_group', label: 'mpg'},
     {key: 'target_num_converged', label: 'target'}, {key: 'snr_db', label: 'SNR [dB]'}, {key: 'fer', label: 'FER', fmt: f3},
     {key: 'frame_errors', label: 'FE'}, {key: 'trials', label: 'trials'},
     {key: 'average_ensemble_effort', label: 'effort', fmt: f3}, {key: 'average_ensemble_latency', label: 'latency', fmt: f3},
     {key: 'average_number_converged_path', label: 'avg #converged', fmt: f3}]);
  document.getElementById('tbl-greedy').innerHTML = table(DATA.greedyRaw,
    [{key: 'source', label: 'source'}, {key: 'variant', label: 'variant'}, {key: 'selector', label: 'selector'},
     {key: 'target_num_converged', label: 'target'}, {key: 'snr_db', label: 'SNR [dB]'}, {key: 'fer', label: 'FER', fmt: f3},
     {key: 'frame_errors', label: 'FE'}, {key: 'trials', label: 'trials'},
     {key: 'average_ensemble_effort', label: 'effort', fmt: f3}, {key: 'average_ensemble_latency', label: 'latency', fmt: f3}]);
  document.getElementById('tbl-pcmfirst').innerHTML = table(DATA.pcmfirstRaw,
    [{key: 'source', label: 'source'}, {key: 'variant', label: 'variant'}, {key: 'selector', label: 'selector'},
     {key: 'members_per_group', label: 'mpg'}, {key: 'target_num_converged', label: 'target'},
     {key: 'snr_db', label: 'SNR [dB]'}, {key: 'fer', label: 'FER', fmt: f3}, {key: 'frame_errors', label: 'FE'},
     {key: 'trials', label: 'trials'}, {key: 'average_ensemble_effort', label: 'effort', fmt: f3},
     {key: 'average_ensemble_latency', label: 'latency', fmt: f3}]);
  document.getElementById('tbl-scheduling').innerHTML = table(DATA.schedulingRaw,
    [{key: 'variant', label: 'variant'}, {key: 'scheduling_config', label: 'schedule'}, {key: 'max_iterations', label: 'maxiter'},
     {key: 'snr_db', label: 'SNR [dB]'}, {key: 'fer', label: 'FER', fmt: f3}, {key: 'frame_errors', label: 'FE'},
     {key: 'trials', label: 'trials'}, {key: 'average_ensemble_effort', label: 'effort', fmt: f3},
     {key: 'average_ensemble_latency', label: 'latency', fmt: f3}]);
  document.getElementById('tbl-alpha').innerHTML = table(DATA.alphaRaw,
    [{key: 'variant', label: 'variant'}, {key: 'norm_factor', label: 'α'}, {key: 'snr_db', label: 'SNR [dB]'},
     {key: 'fer', label: 'FER', fmt: f3}, {key: 'frame_errors', label: 'FE'}, {key: 'trials', label: 'trials'},
     {key: 'average_ensemble_effort', label: 'effort', fmt: f3}]);
}

function buildAll() {
  buildSweepA();
  buildGroupsize();
  buildOrdering();
  buildGreedy();
  buildPcmfirst();
  buildScheduling();
  buildAlpha();
}

buildAll();
buildTables();

try {
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', buildAll);
} catch (e) {}
new MutationObserver(buildAll)
  .observe(document.documentElement, {attributes: true, attributeFilter: ['data-theme']});
</script>
"""

HTML = HTML.replace("__DATA_JS__", DATA_JS)
dest = Path("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/sweeps/plot_sequential_results.html")
dest.write_text(HTML)
print("wrote", dest, len(HTML), "bytes")
