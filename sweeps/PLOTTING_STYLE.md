# Plotting style: interactive error-rate and complexity plots for decoder studies

A style guide for result pages that compare decoders (e.g. ensemble decoders against a near-optimal reference) by
error rate and complexity. It was developed for aSCED on quantum codes and is written generally so it applies to
classical / wireless settings, where the reference is an ML decoder or OSD.

Working reference implementation (QEC version, same structure): `~/qec_claude_logs/build_plots.py` (per-setting
error-rate + complexity plots, parameter-sweep plots) and `~/qec_claude_logs/build_sweep_data.py`.

## 0. Vocabulary (map to your setting)

| Term in this guide | QEC setting | Wireless / classical setting |
|---|---|---|
| **setting** | a code | a code + channel (+ modulation, block length) |
| **channel parameter** (x axis) | physical error probability ε | Eb/N0 or SNR in dB (or BSC crossover probability) |
| **error rate** (y axis) | logical error rate (LER) | BLER / FER (and BER if needed) |
| **reference decoder** | Tesseract long beam (≈ most-likely-error) | ML decoder if feasible, otherwise OSD (state the order, e.g. OSD-2) |
| **baselines** | MWPM, BP, paper curves | BP / min-sum, SCL (polar), OSD of lower order, paper curves |
| **ensemble size K** | number of aSCED paths | number of ensemble decoders / list size |
| **complexity** | CPU time per trial | CPU time per decoded frame (optionally iterations or operations per frame) |
| **match** | reference within or above the 95 % CI of our curve at the chosen operating points | same |

Choose the operating points for the match criterion up front (e.g. the two highest SNRs or lowest ε where the
reference has reliable statistics) and keep them fixed across all plots of a setting.

## 1. What every plot must show

| Rule | How |
|---|---|
| Interactive | Plotly (`plotly.js-dist-min@2.35.2` from `cdn.jsdelivr.net/npm/`), one self-contained HTML page |
| Reference always on | black dash-dot line, never hidden. Say what it is in the name: `ML`, `OSD-2`, or e.g. `Tesseract (long beam)`. If the reference is only an approximation of ML and that was not verified, mark it (e.g. `OSD-2, non-ML`); if it was checked against ML, say so (`≈ ML (checked)`) |
| Literature curves available but hidden | `visible: 'legendonly'`, grey dotted, name suffix ` (paper)` |
| Baselines | standard decoders of the field (e.g. BP / min-sum, SCL), green family; further special baselines violet |
| Don't overload | show the key curves; secondary curves start hidden (`legendonly`), and list the left-out curves with the reason under the plot |
| Complexity | next to every error-rate plot, a complexity plot (CPU per frame) of the same curves; the usual goal is reference-level error rate at lower cost |
| Few errors | points with < 300 errors (frame errors) get hollow markers; capped points are allowed but flagged |
| Recent additions | suffix ` (new)` in the curve name |
| Changes of method | a highlighted note at the top of the page when a method changes (e.g. a new schedule that changes CPU numbers) |

## 2. Data format

One JSON object per setting:

```json
{
  "LDPC_648_rate_1/2_AWGN": {
    "xlabel": "Eb/N0 [dB]", "xlog": false, "ylabel": "BLER",
    "curves": [
      {"name": "BP, 50 it.", "group": "Baseline", "visible": true, "pts": [[x, err_rate, errors, cpu_seconds_total], ...]},
      {"name": "ensemble K=16 · layered", "group": "Ensemble", "visible": true, "pts": [...]},
      {"name": "ensemble K=64 · layered (new)", "group": "Ensemble", "visible": false, "pts": [...]}
    ],
    "reference": {"name": "ML", "pts": [[x, err_rate], ...]},
    "paper": {"Author et al. 2024": [[x, err_rate], ...]}
  }
}
```

- `pts`: `[x, error_rate, errors, cpu_seconds_total]`. CPU per frame in ms = `1000 * cpu * error_rate / errors`
  (frames = errors / error_rate); store the frame count explicitly if error_rate can be 0.
- `cpu_seconds_total` should be thread-seconds (wall time × threads), so multi-threaded runs are comparable.
- Merge several runs of the same configuration (e.g. extension runs) by summing errors, frames and CPU per x before
  building `pts`, but only if the runs draw independent noise (independent RNG seeds).
- Use as many `group`s as there are decoder families (e.g. `Ensemble`, `Hard ensemble`, `Soft ensemble`, `Baseline`).

## 3. Colours and theme

Define colours as CSS tokens and read them in JS, so light and dark mode both work:

```css
:root{--bg:#f6f7f9;--panel:#fff;--ink:#14171c;--ink-2:#4a505a;--ink-3:#737a85;--rule:#dfe2e7;--grid:#e7e9ed;
      --c-a:#2a78d6;--c-b:#eb6834;--c-base:#1baf7a;--c-extra:#4a3aa7;--c-ref:#14171c}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){color-scheme:dark;--bg:#121417;--panel:#1a1d21;
      --ink:#eef0f3;--ink-2:#b4bac3;--ink-3:#8a909a;--rule:#2c3036;--grid:#262a30;
      --c-a:#3987e5;--c-b:#d95926;--c-base:#199e70;--c-extra:#9085e9;--c-ref:#eef0f3}}
:root[data-theme="dark"]{ /* same dark values again, so an explicit toggle wins */ }
body{background:var(--bg);color:var(--ink)}
```

```js
const tok = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
try { matchMedia('(prefers-color-scheme: dark)').addEventListener('change', drawAll) } catch (e) {}
new MutationObserver(drawAll).observe(document.documentElement, {attributes: true, attributeFilter: ['data-theme']});
```

| Series | Colour token | Line | Marker |
|---|---|---|---|
| Decoder family A (e.g. hard / main method) | `--c-a` (blue) | solid / dash / dot / dashdot / longdash, cycling within the family | circle, square, diamond, triangle-up, x |
| Decoder family B (e.g. soft / variant) | `--c-b` (orange) | same cycling | same |
| Standard baselines | `--c-base` (green) | dash / dot | small circles |
| Special baseline | `--c-extra` (violet) | solid | star |
| Reference (ML / OSD) | `--c-ref` (ink) | dashdot, width 2 | small circles |
| Literature | `--ink-3` (grey) | dot, width 1.3, `legendonly` | small circles |

Colour encodes the family; line style and marker distinguish curves within a family. Curves with only 1–2
points are drawn as markers only (`mode: 'markers'`, larger marker). At most ~8 categorical colours; beyond that,
fold curves into hidden legend entries or separate plots.

## 4. Error-rate trace (per curve)

```js
function curveTraces(d, cpu) {   // cpu = false: error-rate plot, true: complexity plot
  const col = {"Ensemble": tok('--c-a'), "Soft ensemble": tok('--c-b'), "Baseline": tok('--c-base')};
  const dashes = ["solid", "dash", "dot", "dashdot", "longdash"], syms = ["circle", "square", "diamond", "triangle-up", "x"], cnt = {};
  return d.curves.map(c => {
    const i = (cnt[c.group] = (cnt[c.group] || 0) + 1) - 1, few = c.pts.length <= 2;
    const ms = c.pts.map(p => p[3] > 0 && p[2] > 0 ? 1000 * p[3] * p[1] / p[2] : null);
    return {
      name: c.name, legendgroup: c.name, showlegend: !cpu,
      x: c.pts.map(p => p[0]), y: cpu ? ms : c.pts.map(p => p[1]),
      customdata: c.pts.map((p, k) => [p[2], ms[k] == null ? 'n/a' : ms[k].toFixed(ms[k] < 10 ? 2 : 0)]),
      hovertemplate: (cpu ? '%{y:.3g} ms/frame' : '%{y:.3e}') + ' · %{customdata[0]} err · %{customdata[1]} ms/frame<extra>' + c.name + '</extra>',
      mode: few ? 'markers' : 'lines+markers',
      line: {color: col[c.group], width: 2, dash: dashes[i % 5]},
      marker: {size: few ? 10 : 6, symbol: syms[i % 5], color: c.pts.map(p => p[2] < 300 ? 'rgba(0,0,0,0)' : col[c.group]),
               line: {width: 1.6, color: col[c.group]}},          // hollow marker = fewer than 300 errors
      visible: c.visible ? true : 'legendonly'                    // secondary curves start hidden
    };
  });
}
```

Append the reference (always visible) and the hidden literature curves:

```js
const R = d.reference;
tr.push({name: R.name, x: R.pts.map(p => p[0]), y: R.pts.map(p => p[1]), mode: 'lines+markers',
         line: {color: tok('--c-ref'), dash: 'dashdot', width: 2}, marker: {size: 5, color: tok('--c-ref')},
         hovertemplate: '%{y:.3e}<extra>' + R.name + '</extra>'});
Object.entries(d.paper || {}).forEach(([n, s]) => tr.push({name: n + ' (paper)', x: s.map(p => p[0]), y: s.map(p => p[1]),
         mode: 'lines+markers', line: {color: tok('--ink-3'), dash: 'dot', width: 1.3}, marker: {size: 4, color: tok('--ink-3')},
         visible: 'legendonly'}));
```

If a reference decoder is only available at some points (ML is often infeasible at low error rates), plot what exists
and state the range in the text; do not extrapolate it.

## 5. Layout

```js
function layout(d, ytitle, ylog) {
  const ax = {gridcolor: tok('--grid'), linecolor: tok('--rule'), zerolinecolor: tok('--rule'),
              tickfont: {color: tok('--ink-2'), family: 'IBM Plex Mono, monospace', size: 11}};
  return {paper_bgcolor: tok('--panel'), plot_bgcolor: tok('--panel'),
          font: {family: 'IBM Plex Sans, system-ui, sans-serif', color: tok('--ink')}, margin: {l: 60, r: 10, t: 8, b: 44},
          xaxis: Object.assign({title: {text: d.xlabel, font: {color: tok('--ink-2'), size: 12}}, type: d.xlog ? 'log' : 'linear'}, ax),
          yaxis: Object.assign({title: {text: ytitle, font: {color: tok('--ink-2'), size: 12}}, type: ylog ? 'log' : 'linear', exponentformat: 'power'}, ax),
          legend: {font: {size: 11, color: tok('--ink-2')}, bgcolor: 'rgba(0,0,0,0)'}, hovermode: 'closest',
          hoverlabel: {bgcolor: tok('--panel'), bordercolor: tok('--rule'), font: {color: tok('--ink'), family: 'IBM Plex Mono, monospace', size: 12}}};
}
const cfg = {responsive: true, displaylogo: false, modeBarButtonsToRemove: ['lasso2d', 'select2d']};
```

- Error rate: log y axis. x: SNR / Eb/N0 in dB on a linear axis (error rate falls to the right); error probability ε
  on a linear or log axis (error rate falls to the left). Keep the orientation consistent across all plots of a page.
- Many curves: legend below the plot (`legend: {orientation: 'h', y: -0.2}`) and a plot height of ≈ 560 px.
- Complexity plot: log y, `'CPU per frame [ms]'`, `showlegend: false` (the legend lives on the error-rate plot).
- Use `Plotly.react` (not `newPlot`) so redraws after a theme change are cheap.

## 6. Activating / deactivating curves (legend toggles, synced plots)

- Click a legend entry to hide / show a curve; double-click to isolate it (Plotly default).
- `visible: 'legendonly'` makes a curve start hidden but keeps it in the legend.
- The complexity plot has no legend of its own; keep it in sync with the error-rate plot:

```js
Plotly.react('err-' + id, errTraces, errLayout, cfg);
Plotly.react('cpu-' + id, cpuTraces, cpuLayout, cfg);
const el = document.getElementById('err-' + id);
if (!el._synced) {                              // attach once; react() keeps the handler
  el._synced = true;
  el.on('plotly_restyle', () => {
    const vis = {}; el.data.forEach(t => vis[t.name] = t.visible);
    const c = document.getElementById('cpu-' + id);
    Plotly.restyle(c, {visible: c.data.map(t => vis[t.name] === undefined ? true : vis[t.name])});
  });
}
```

Traces are matched by `name`, so every curve needs a unique, readable name, e.g.
`ensemble K=64 · layered · clustering (new)`.

## 7. Page layout per setting

```html
<section class="card" id="LDPC_648">
  <header><h2>LDPC (648, 324) <span class="sub">rate 1/2, AWGN, BPSK</span></h2><span class="pill good">matched</span></header>
  <p class="verdict">One-paragraph verdict with the key numbers (error rate / reference at the operating points, CPU per frame).</p>
  <div class="pair">                              <!-- grid: 1.35fr 1fr; one column below 860 px -->
    <div class="plot" id="err-LDPC_648"></div>
    <div class="plot" id="cpu-LDPC_648"></div>
  </div>
  <details><summary>Left out of this plot, and why</summary><ul><li>…</li></ul></details>
</section>
```

- A sticky table of contents with one chip per setting (coloured dot = verdict).
- Status pills: `good` = matches the reference at the operating points, `warn` = partly / tentative, `bad` = not reached.
- Phone width: plots stack to one column, heights ≈ 480 px, 16 px side gutter, no horizontal page scroll.

## 8. Parameter-sweep plots (separate from the error-rate-vs-channel plots)

| Sweep | x axis | y axis | Traces |
|---|---|---|---|
| Ensemble size K | K, log axis with `dtick: Math.log10(2)` (powers of two) | error rate / reference at one operating point, log, 95 % CI bars, dashed line at 1 | one per configuration (K removed from the label); family by colour, variant by dash; hollow < 300 errors. Companion plot: CPU per frame vs K |
| A continuous decoder parameter (bias, damping, scaling factor, list size, …) | parameter (categories if the values are irregular) | absolute error rate, log | one line per channel point, sequential single-hue ramp (e.g. `hsl(215, 70→50 %, 28→73 %)`, darker = better channel) |
| A structural parameter (e.g. number of splitters, sub-decoders per ensemble) | integer ticks | error rate / reference at one operating point, reference line at 1 | one per variant; companion CPU plot |
| Algorithm variants (schedule, selection rule, …) | horizontal bars, sorted by value | error rate / reference at one operating point (log x), CI bars, dashed line at 1 | bar colour = family; label `K=256 · layered · clustering`; div height `90 + 26 × bars` px |

Reference line at 1 in a ratio plot:

```js
tr.push({x: [x0, x1], y: [1, 1], mode: 'lines', line: {color: tok('--ink-3'), dash: 'dash', width: 1}, hoverinfo: 'skip', showlegend: false});
// for horizontal bar charts use a shape instead:
L.shapes = [{type: 'line', x0: 1, x1: 1, yref: 'paper', y0: 0, y1: 1, line: {color: tok('--ink-3'), dash: 'dash', width: 1}}];
```

CI bars from the Wilson interval, divided by the reference value at that point:

```python
def ci95(k, n):                      # Wilson score interval for k errors in n frames
    import math
    p, z = k / n, 1.96
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)
```

```js
error_y: {type: 'data', symmetric: false, array: hi - ratio, arrayminus: ratio - lo, thickness: 1}
```

If the reference has no value at the chosen point, interpolate it in log(error rate) between its neighbours
(never extrapolate).

## 9. Comparing against the reference on shared samples

When a decoder appears to beat a near-ML reference, check it on the **same** noise realisations before believing it:
decode one set of sampled frames with both decoders, then compare

- failures of each and the paired counts (both fail / only ours / only reference), with a McNemar z-score
  `(only_ref - only_ours) / sqrt(only_ref + only_ours)`;
- on frames where only the reference fails, the likelihood (or distance / weight) of each decoder's output against
  the transmitted word: if the reference returns a *less* likely word than ours, the reference is not ML there.

Reproduce the reference's published numbers on these shared samples first (same decoder configuration and code
construction as its authors); a hand-built configuration can be much weaker than the real one.

## 10. Text around the plots

- A short lede under the title: what is plotted, which reference and baselines, how to toggle.
- A small-print note: hollow markers = < 300 errors; CPU measured on the stated hardware (mixed node types → compare
  only large ratios); what the reference label means (ML / OSD order, verified or not).
- Per setting: the verdict paragraph and the "left out, and why" list (every omitted configuration with its reason).
- Numbers in the text: error rate / reference with two decimals, errors in parentheses where fewer than 300, CPU per
  frame in ms.

## 11. Checklist before publishing

1. Reference visible on every error-rate plot and named precisely (ML, OSD-order, …); literature curves hidden but present.
2. Every curve name unique (the legend sync matches by name).
3. Hollow markers where errors < 300; the left-out list is filled.
4. Axis orientation consistent across the page (SNR in dB vs error probability).
5. Parse-check the page's JS (e.g. `uv run --with esprima python -c "import esprima; esprima.parseScript(js)"`)
   and look at the rendered page in light and dark mode.
6. Page title: short name (e.g. `Ensemble Decoding Plots`), not a sentence.
