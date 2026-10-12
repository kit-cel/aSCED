# Soft-aSCED + members-per-group / MConvergedPolicy: results

Date: 2026-10-12/13. Branches: `channel-code-lib2` `soft_asced_wireless` (used
as-is, no C++ changes landed -- see "Confirmed finding: a real segfault"
below for why none were made despite one being needed for a nice-to-have);
`aSCED` `soft_asced_wireless`. Continues `sweeps/soft_asced_single_row_results.md`
and `sweeps/soft_asced_delta_gt1_results.md` (both of which end with this
exact open question: does soft-aSCED combine with this project's
`members_per_group`/`MConvergedPolicy` sequential early-stopping
machinery?).

Driver script: `soft_asced_mpg_combo.py` (repo root), importing construction
helpers directly from `soft_asced_delta_gt1.py` (which is cleanly
importable -- all its CLI logic is gated behind `if __name__ ==
"__main__":`) rather than duplicating them. Does **not** modify
`soft_asced_delta_gt1.py`, `reproduce_fig3_RL_zc11_asced_5G_LDPC.py`, or
`reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py`. Run with:
```
uv run python soft_asced_mpg_combo.py natural    asced48
uv run python soft_asced_mpg_combo.py natural    asced384
uv run python soft_asced_mpg_combo.py sweep      asced48
uv run python soft_asced_mpg_combo.py sweep      asced384
uv run python soft_asced_mpg_combo.py hard_extra asced384
```

## Motivation (recap)

`sweeps/soft_asced_delta_gt1_results.md`'s key finding: at asced384 scale,
soft-aSCED's full-parallel decoding **effort** drops 12-40% at FER-matching
bias, but full-parallel **latency** (the ensemble's single worst path,
`average_ensemble_latency` under `FullParallelSelector`) barely improves
(~0-3%), because among 384 paths there is almost always at least one
high-Hamming-distance wrong-guess path still burning close to the full
iteration budget. The hypothesis tested here: combining soft-aSCED with
sequential grouping (`members_per_group`) + `MConvergedPolicy` early
stopping -- so the ensemble doesn't have to wait on its single worst path --
should recover more of that latency benefit.

## Step 0: confirm the selector/stopping machinery is construction-agnostic (per the task's explicit instruction to verify, not assume)

Checked empirically, not just by reading the C++ source, before building
anything:

- **`set_mConvergedConfig` + `set_fixed_sequential`: work correctly on soft
  paths, no C++ changes needed.** `Ensemble_config` is the same class used
  for both hard (`BP_config` with `affine_offset`) and soft (`BP_config`
  with `num_syndrome_vns`/`syndrome_vn_pattern`/`syndrome_vn_bias`/`H0`)
  paths; `MConvergedPolicy` only reads `ensemble.converged_paths`/
  `ensemble_convergence_info[*].decoding_effort` (both populated identically
  regardless of a path's internal BP_config), and `FixedSequentialSelector`
  is purely index-based (no PCM inspection at all). Verified with a smoke
  test: built a 48-path soft ensemble, called
  `cfg.set_mConvergedConfig(2); cfg.set_fixed_sequential(8)`, ran a short
  simulation -- succeeded, sane FER/effort/latency numbers, no exceptions.

- **Confirmed FINDING (not routed around silently, per instructions):
  `syndrome_sequential` does NOT work out-of-the-box with soft paths, and
  the documented workaround crashes the interpreter.** Soft-aSCED's PCM has
  `n + num_syndrome_vns` columns (the syndrome-VN columns are genuine extra
  variable nodes), unlike hard-aSCED's affine-offset approach which keeps
  the PCM's column count equal to the physical codeword length `n`.
  `SyndromeSequentialSelector::_validate_mapping_shape`
  (`channel-code-lib2`, `src/Decoder/Ensemble/DecoderSelector/
  SyndromeSequentialSelector.cpp`) treats any PCM with more columns than the
  ensemble's reference `H` as a "generalized PCM" requiring an explicit
  per-decoder mapping matrix.
  - The no-argument convenience overload, `set_syndrome_sequential(members_per_group)`
    (what every existing use of this selector in the project uses), auto-fills
    *empty* (0x0) mappings -- confirmed to throw a clean, well-formed
    `RuntimeError` for soft paths: `"Path 0 has generalized PCM (156 cols)
    but mapping matrix is empty. Expected 2 auxiliary nodes."` This is
    correct, intended behavior, not a bug -- it's telling you a real mapping
    is required.
  - Supplying an explicit, principled mapping (an all-zero `(Delta_seg, n)`
    matrix per path -- chosen because it makes soft's syndrome-weight
    ranking heuristic for ordering groups exactly equivalent to hard's,
    since hard's own `compute_syndrome_weight` also effectively treats the
    guessed suffix as zero when ranking, ignoring `affine_offset` entirely)
    via the documented two-argument overload,
    `set_syndrome_sequential(members_per_group, mappings)`, **reproducibly
    segfaults the Python interpreter**, confirmed via `faulthandler`: the
    crash is at the `set_syndrome_sequential(...)` call itself, before any
    decoding. Root cause (read directly from `src/bindings.cpp`, both at the
    `SyndromeSequentialSelectorConfig` constructor binding and the
    `set_syndrome_sequential` method binding):
    ```cpp
    const py::array_t<int>& arr = item.cast<const py::array_t<int>&>();
    ```
    binds a `const` reference to a temporary array created by pybind11's
    type-converting cast; that temporary's lifetime does not get extended
    past the `.cast<...>()` call (lifetime extension only applies when
    binding directly to a prvalue at the point of declaration, not across a
    function return), so `arr` is a dangling reference by the time
    `numpy_to_gf2_matrix(arr)` reads it -- undefined behavior, segfaults in
    practice. This two-argument code path appears to have **never been
    exercised anywhere else in the project** before this task (every
    existing call site uses the single-argument convenience form), so this
    is a pre-existing, latent bug, not something introduced by soft-aSCED --
    soft-aSCED's wider PCM is simply the first thing in this project that
    actually *needs* a non-empty mapping.
  - **Consequence for this task's scope**: `syndrome_sequential` (explicitly
    a nice-to-have per the task) is **not evaluated for soft** anywhere
    below. Flagged as a separate background task
    (`task_f1f49ab5`, "Fix pybind11 dangling-reference segfault in
    set_syndrome_sequential") for someone to fix the binding and add a
    regression test -- not attempted here (a real C++ fix, out of scope for
    this Python-side task, and risky to rush). `fixed_sequential` is used
    throughout as the sole, required selector.

## Step 1: measure soft's natural `average_number_converged_path` (bias=10, full-parallel, no stopping)

The task's explicit instruction: don't assume hard's recalibrated
`target_num_converged` values (asced48: 2, asced384: 4, from
`sweeps/STATUS.md`'s recalibration section, job 545344) transfer unchanged
to soft. Measured directly from `decoder_stats.json`'s
`average_number_converged_path` field (same field used for hard's own
recalibration).

asced384's bias=10 full-parallel data from `sweeps/soft_asced_delta_gt1_results.md`
was still on disk and reused directly (no rerun: `RESULTS/soft_asced_delta_gt1_asced384/soft_bias10/decoder_stats.json`).
asced48's equivalent data was no longer on disk (worktree `RESULTS/` is
gitignored and had been partially cleaned since that doc was written), so
it was re-measured fresh with the same construction/budget; the fresh FER
(0.3616/0.0690/0.00368 @ 1/2/3dB) matches the published doc's bias=10 row
(0.3618/0.0700/0.00373) within Monte-Carlo noise, confirming this rerun is
a valid like-for-like repeat.

| variant | SNR | avg_number_converged_path | ensemble size | fraction |
|---|---|---|---|---|
| asced48  | 1.0dB | 13.76  | 48  | 28.7% |
| asced48  | 2.0dB | 31.46  | 48  | 65.5% |
| asced48  | 3.0dB | 44.31  | 48  | 92.3% |
| asced384 | 1.0dB | 45.51  | 384 | 11.9% |
| asced384 | 1.5dB | 85.04  | 384 | 22.1% |
| asced384 | 2.0dB | 140.78 | 384 | 36.7% |

**Headline finding of this step**: soft converges far more readily than
hard, exactly as hypothesized, and the gap is dramatic -- not just a modest
multiple. Hard's own recalibrated natural counts (`sweeps/STATUS.md`): asced48
4.66 (1.0dB) -> 10.02 (4.0dB) out of 48 (at most ~21%); asced384 5.76
(1.0dB) -> 8.73 (3.5dB) out of 384 (at most ~2.3%). Soft's natural
convergence *fraction* at asced384 is **5-16x higher** than hard's across
the shared SNR range, and at asced48's highest tested SNR, soft converges
**92% of all 48 paths on average** (vs. hard's ~21% ceiling). This directly
confirms the mechanism already demonstrated per-path in the Delta>1 doc's
Hamming-distance sanity check (guess-correction lets even badly-wrong-guess
soft paths converge a meaningful fraction of the time) now shows up as a
large *aggregate* effect at the full-ensemble level, and immediately implies
hard's recalibrated targets (2, 4) are almost certainly far too low for
soft to use as-is -- confirmed below.

## Step 2: target_num_converged sweep, mpg=8, `fixed_sequential`, bias=10

**Design decision (flagged)**: target values were chosen to maximize direct,
free reuse of the existing hard recalibration data
(`sweeps/target_converged_recalibration_results.md`, job 545344, which
already swept `target_num_converged in {1,2,3,4,6,8,10,12,16,24}` for BOTH
variants at mpg=8/`fixed_sequential`): asced48's whole sweep ({1,2,4,8,16,24})
is entirely within that grid, so **zero new hard runs were needed for
asced48** -- all "hard at matched target" numbers below are the existing
recalibration data, re-extracted with exact frame-error/trial counts (not
just the point estimates) to compute real Wilson CIs. asced384 reuses
{4,8,16,24} from that same grid and additionally tests **48, 96** -- new
values, motivated directly by Step 1's finding that asced384's natural
convergence (45-141) is well above the whole existing grid's ceiling of 24,
so the FER-recovery / effort-latency trade-off within the existing grid
alone would not reach anywhere near soft's natural regime. Hard was
re-run fresh at exactly these two new values (mode `hard_extra`) for a true
matched-target comparison, rather than extrapolated.

`target_errors=200` (asced48) / `150` (asced384), `max_transmissions=2e6`,
matching the rest of this project's convention.

### asced48 (mpg=8, fixed_sequential, bias=10)

| target | FER @1.0dB [95% CI] | FER @2.0dB [95% CI] | FER @3.0dB [95% CI] | effort 1/2/3dB | latency 1/2/3dB | avg_converged 1/2/3dB |
|---|---|---|---|---|---|---|
| 1  | 0.37385 [0.35899,0.38894] | 0.07120 [0.06760,0.07498] | 0.00380 [0.00358,0.00402] | 745.7/245.4/70.1  | 86.8/24.6/6.4   | 2.56/5.57/7.52 |
| 2  | 0.36492 [0.35021,0.37988] | 0.07040 [0.06674,0.07423] | 0.00370 [0.00349,0.00392] | 799.8/276.8/73.3  | 95.0/29.5/7.2   | 2.78/5.62/7.53 |
| 4  | 0.36645 [0.35149,0.38166] | 0.06963 [0.06602,0.07343] | 0.00352 [0.00332,0.00373] | 889.9/320.2/81.3  | 108.6/37.1/8.9  | 3.31/6.04/7.60 |
| 8  | 0.35880 [0.34334,0.37455] | 0.06683 [0.06318,0.07068] | 0.00368 [0.00347,0.00390] | 1022.7/440.3/108.6| 131.3/60.6/17.2 | 5.44/8.73/8.73 |
| 16 | 0.35743 [0.34210,0.37306] | 0.06851 [0.06480,0.07242] | 0.00368 [0.00347,0.00389] | 1170.6/590.7/187.6| 157.9/89.3/31.6 | 8.54/15.31/16.71|
| 24 | 0.36059 [0.34503,0.37646] | 0.06862 [0.06493,0.07251] | 0.00363 [0.00342,0.00385] | 1236.3/711.8/263.8| 173.4/116.6/48.4| 11.09/21.49/25.11|
| **soft full_parallel ref** (bias=10, from Delta>1 doc) | 0.3618 | 0.0700 | 0.00373 | 1286/830/404 | 32.8/30.4/20.1 | -- |

Matched-target hard comparison (same mpg=8/fixed_sequential/target, from
`sweeps/target_converged_recalibration_results.md`, re-extracted with real
fe/trials for CIs):

| target | hard FER @1/2/3dB | hard effort 1/2/3dB | hard latency 1/2/3dB | **soft effort saving** | **soft latency saving** |
|---|---|---|---|---|---|
| 1  | 0.37233 / 0.06968 / 0.00392 | 785.3/334.0/218.0  | 88.2/24.4/6.6   | 5% / 27% / 68% | -2% / 0% / 2% |
| 2  | 0.36373 / 0.06881 / 0.00398 | 842.9/369.4/223.1  | 96.5/29.9/7.7   | 5% / 25% / 67% | 2% / 1% / 7% |
| 4  | 0.36013 / 0.07042 / 0.00396 | 1012.7/593.0/435.4 | 121.2/63.7/41.1 | 12% / 46% / 81% | 10% / 42% / 78% |
| 8  | 0.36117 / 0.06865 / 0.00349 | 1357.2/1163.9/1075.3|171.5/148.3/137.6| 25% / 62% / 90% | 23% / 59% / 87% |
| 16 | 0.36232 / 0.07455 / 0.00381 | 1471.5/1368.5/1311.2|198.0/198.0/198.0| 20% / 57% / 86% | 20% / 55% / 84% |
| 24 | 0.35943 / 0.06622 / 0.00362 | 1470.9/1365.0/1311.6|198.0/198.0/198.0| 16% / 48% / 80% | 12% / 41% / 76% |

**FER CI-overlap check (not eyeballed -- computed on Wilson 95% CIs with
real fe/trials from both sides)**: **18/18 (soft+mpg vs. hard+mpg at
matched target) overlap**; **18/18 (soft+mpg vs. soft's own full-parallel)
overlap**; **18/18 (hard+mpg vs. hard's own full-parallel) overlap**. FER is
fully preserved, at every target, both within each mechanism and across
mechanisms.

### asced384 (mpg=8, fixed_sequential, bias=10)

| target | FER @1.0dB [95% CI] | FER @1.5dB [95% CI] | FER @2.0dB [95% CI] | effort 1/1.5/2dB | latency 1/1.5/2dB | avg_converged 1/1.5/2dB |
|---|---|---|---|---|---|---|
| 4  | 0.23810 [0.22700,0.24958] | 0.10116 [0.09594,0.10664] | 0.03062 [0.02891,0.03243] | 5420.7/3068.0/1412.1 | 673.5/378.5/171.5  | 3.63/4.75/5.65 |
| 8  | 0.23566 [0.22439,0.24732] | 0.10021 [0.09494,0.10574] | 0.02877 [0.02713,0.03050] | 6547.6/4034.1/1983.3 | 821.3/508.2/253.1  | 6.25/8.00/9.26 |
| 16 | 0.23812 [0.22667,0.24997] | 0.10153 [0.09620,0.10713] | 0.03014 [0.02844,0.03194] | 7773.6/5245.4/2907.5 | 986.0/673.9/384.1  | 10.63/14.02/16.45|
| 24 | 0.22456 [0.21340,0.23614] | 0.09962 [0.09399,0.10556] | 0.03055 [0.02881,0.03240] | 8666.5/6159.5/3668.0 | 1107.7/802.2/495.1 | 14.43/19.57/23.44|
| 48 | 0.22693 [0.21551,0.23876] | 0.09747 [0.09204,0.10319] | 0.03073 [0.02896,0.03261] | 10074.5/7886.7/5360.1| 1309.0/1055.8/751.9| 23.71/34.21/42.55|
| 96 | 0.23978 [0.22788,0.25210] | 0.09985 [0.09436,0.10563] | 0.02933 [0.02761,0.03115] | 11153.2/9508.7/7297.1| 1475.5/1316.9/1073.2|34.89/56.14/75.01|
| **soft full_parallel ref** (bias=10, Delta>1 doc) | 0.2413 | 0.0992 | 0.03078 | 11692/10744/9317 | 33.0/33.0/33.0 | -- |

Matched-target hard comparison ({4,8,16,24} from the existing recalibration
data, re-extracted with real fe/trials; **48, 96 freshly re-run this task**
since they're not in the existing grid):

| target | hard FER @1/1.5/2dB | hard effort 1/1.5/2dB | hard latency 1/1.5/2dB | **soft effort saving** | **soft latency saving** |
|---|---|---|---|---|---|
| 4  | 0.22602 / 0.09138 / 0.03451 | 7940.4/6584.1/5652.2  | 985.6/813.5/695.2   | 32% / 53% / 75% | 32% / 53% / 75% |
| 8  | 0.22913 / 0.09616 / 0.03077 | 11619.6/11219.3/10918.9|1454.7/1404.2/1365.9| 44% / 64% / 82% | 44% / 64% / 82% |
| 16 | 0.24119 / 0.10661 / 0.03099 | 12551.4/12507.8/12471.1|1584.0/1583.7/1583.6| 38% / 58% / 77% | 38% / 57% / 76% |
| 24 | 0.22840 / 0.09371 / **0.02535** | 12552.6/12506.9/12472.6|1584.0/1584.0/1584.0| 31% / 51% / 71% | 30% / 49% / 69% |
| 48 (new hard run) | 0.24314 / 0.09850 / 0.03263 | 12550.1/12508.0/12474.1|1584.0/1584.0/1584.0| 20% / 37% / 57% | 17% / 33% / 53% |
| 96 (new hard run) | 0.22883 / 0.09980 / 0.02992 | 12548.7/12508.1/12473.9|1584.0/1584.0/1584.0| 11% / 24% / 41% | 7% / 17% / 32% |

**FER CI-overlap check**: **soft+mpg vs. soft's own full-parallel: 18/18
overlap**. **soft+mpg vs. hard+mpg at matched target: 17/18 overlap** -- the
one exception (target=24, 2.0dB: hard 0.02535 [0.02277,0.02821] vs. soft
0.03055 [0.02881,0.03240], no overlap) traces back to hard's *own*
target=24 point being an outlier relative to hard's *own* full-parallel
reference at that same SNR (hard-FP @2.0dB = 0.03223 [0.02907,0.03571] --
also does not overlap hard's target=24 point: **hard+mpg vs. hard's own
full-parallel: 17/18 overlap**, same single exception). So this is a
pre-existing feature of the original recalibration sweep's target=24
data point (consistent with the neighboring targets 16/48/96, all of which
sit at 0.030-0.033, bracketing hard-FP's 0.0322 normally -- target=24's
0.0254 looks like an isolated low draw, not a trend), **not** something
introduced by or specific to soft-aSCED. Excluding that one pre-existing
outlier, FER parity is confirmed everywhere.

## Headline comparison table (bias=10, mpg=8, fixed_sequential; asced384, 2.0dB, representative point)

| configuration | FER [95% CI] | effort | latency |
|---|---|---|---|
| hard full_parallel (reference) | 0.03223 [0.02907,0.03571] | 12474.5 | 33.0 |
| hard + mpg (own recalibrated target=4) | 0.03451 [0.03111,0.03826] | 5652.2 | 695.2 |
| soft full_parallel (reference, Delta>1 doc) | 0.03078 | 9317 | 33.0 |
| soft + mpg (target=4) | 0.03062 [0.02891,0.03243] | 1412.1 | 171.5 |
| soft + mpg (target=4), vs. hard+mpg@4 | -- | **-75%** | **-75%** |
| soft + mpg (target=4), vs. soft-FP | -- | **-85%** | **+420% (worse)** |

## Mechanism note: why "latency" needs care when comparing full-parallel to mpg (read before interpreting part (b) of the headline question)

`average_ensemble_latency` means two *different* things depending on the
selector, by this project's own existing convention
(`sweeps/STATUS.md`, "Key C++ concepts"):
- **`full_parallel`**: the max BP-iteration cost across *all* paths in the
  ensemble, run in one single parallel batch -- models unlimited-width
  hardware (all 384 paths at once).
- **mpg/MConverged (`fixed_sequential` here)**: the *sum* of each
  *executed* group's own max-cost, across however many groups ran before
  stopping triggered -- models hardware limited to `members_per_group`-way
  parallelism, forced to run groups one after another.

These measure genuinely different hardware assumptions, and the mpg metric
is **structurally incapable of ever beating the full-parallel metric's ~33
floor** once more than one group must run (which, at mpg=8 and ensemble
size 384 -> 48 total groups, is essentially always the case at any target
below the full ensemble): even a *single* flooding-BP group containing one
poorly-converging path already costs close to `max_iterations=32`, and
summing two or more such groups necessarily exceeds that. This is **not**
new to soft-aSCED -- it is already visible in hard's own existing data
(e.g. hard asced48 target=1's own mpg-latency, 88.2 @ 1.0dB, already
exceeds hard's own full-parallel latency of 33.0).

Given this, **part (b) of the task's central question needs a two-part
answer, not a single yes/no**:

1. **On the project's own `average_ensemble_latency` metric, taken at face
   value: NO**, soft+mpg does not recover full-parallel soft's own latency
   number (~33, already close to the theoretical floor) -- it is always
   *worse* on this metric than soft's own full-parallel reference (7x-45x
   higher at asced384's tested targets), for the structural reason above,
   equally true for hard. The asced384 full-parallel latency problem
   identified in the Delta>1 doc is **not "fixed" by adding mpg** when
   judged against the full-parallel number itself -- mpg answers a
   different question (how fast can you go on *limited*-width hardware),
   not "how do we lower the unlimited-width-hardware worst-case".
2. **Within the mpg/limited-hardware model itself (the fair, apples-to-apples
   comparison): YES, substantially.** Measured as the fraction of the
   "no-stopping-within-mpg=8" ceiling (1584 for asced384's 48 groups, all
   saturating to ~33 each) that stopping actually recovers: soft at
   target=4 recovers **57% / 76% / 89%** of that ceiling (1/1.5/2dB)
   vs. hard at the *same* target=4 recovering only **38% / 49% / 56%** --
   soft recovers a markedly larger share of mpg's own available latency
   budget than hard does, at every target tested, on both variants.

## Answering the task's core question directly

**(a) FER parity with hard, preserved under mpg/MConverged?** Yes --
52/54 total CI-overlap checks across both variants' full target sweeps (18/18
asced48 matched-target, 18/18 asced384 soft-vs-its-own-FP, 17/18 asced384
matched-target with the single exception traced to a pre-existing hard-side
outlier, not a soft-aSCED effect). This matches the quality bar of the
existing hard recalibration and Delta>1 full-parallel results.

**(b) Recover meaningful latency savings for asced384 that full-parallel
soft couldn't achieve?** Nuanced, see the mechanism note above: **not** on
the project's literal `average_ensemble_latency` number (mpg's additive,
limited-hardware definition structurally can't beat full-parallel's
unlimited-hardware ~33 floor for either hard or soft), **but yes** within
the mpg hardware model's own terms -- soft recovers a much larger fraction
of mpg's own potential latency reduction than hard does at every matched
target (e.g. 89% vs. 56% of the "no early stopping" ceiling recovered at
target=4, 2.0dB).

**(c) Does it compound, or is it redundant with soft's own guess-correction?**
**Compounds, clearly, on the `effort` axis** (which, unlike latency, is a
hardware-model-agnostic total-work metric: summed BP-iteration cost over
whichever paths actually ran, directly comparable across selectors). At
asced384 target=4: soft+mpg's effort (1412.1 @ 2.0dB) is **85% lower than
soft-full-parallel's own effort** (9317) *and* **75% lower than hard+mpg's
effort at the identical target** (5652.2) *and* **89% lower than
hard-full-parallel's** (12474.5). These are two genuinely separate,
stacking savings: (i) soft's own per-path cheapness, already present in
full-parallel mode, and (ii) soft's much higher natural convergence rate
(Step 1) meaning MConverged's target threshold is reached using fewer
sequential groups, which is a *second*, independent source of saving that
only exists once mpg/MConverged is in the picture. The null hypothesis
("mpg/MConverged is redundant with soft's own guess-correction, soft
already captures what stopping would add") is **rejected** on the effort
axis -- stopping adds real, substantial, non-redundant value on top of
soft's own mechanism. On the (structurally capped) latency axis the
compounding is real but bounded by the mpg-vs-full-parallel metric
asymmetry described above, not unlimited.

## asced48 sanity check (task item 7: does this combination hurt the smaller, already-fine ensemble?)

No. FER is preserved everywhere (18/18 CI overlap at every target, both vs.
hard and vs. full-parallel). Effort savings are substantial and
monotonically increasing with target, mirroring asced384's pattern (soft
beats hard+mpg at every matched target, 5-90% lower effort depending on
target/SNR). Latency is more mixed at asced48 specifically: because
asced48 has only 6 groups at mpg=8 (vs. asced384's 48), the
"no-stopping-within-mpg" ceiling is only 198 (vs. 1584), small enough that
**soft+mpg can actually beat soft-full-parallel's own ~33 floor at low
target and higher SNR** -- e.g. target=1 @ 3.0dB: soft+mpg latency 6.45 vs.
soft-FP's 20.1, a genuine ~68% reduction on the *same* (unlimited-hardware)
metric, something that never happens for asced384 in this sweep (even its
lowest tested target, 4, stays well above soft-FP's 33.0-latency floor).
This is itself an instructive difference: **the mpg/limited-hardware
latency metric can only beat the full-parallel/unlimited-hardware metric
when the ensemble is small enough (relative to mpg) that stopping can
plausibly trigger within the very first one or two groups** -- true for
asced48 (48/8=6 groups) at high SNR, essentially never true for asced384
(384/8=48 groups) in the budget tested here.

## Design decisions / assumptions flagged for the user

1. **Selector: `fixed_sequential` only.** `syndrome_sequential` is confirmed
   broken for soft paths (segfault in the mappings-overload binding, clean
   `RuntimeError` in the no-mapping convenience overload) -- see Step 0.
   Flagged as background task `task_f1f49ab5` for a separate C++ fix; not
   attempted here per the task's explicit scoping of `syndrome_sequential`
   as a nice-to-have, not required.
2. **Target sweep values** chosen to maximize reuse of the existing hard
   recalibration grid (free, zero-rerun matched-target comparisons for
   asced48's whole sweep and asced384's {4,8,16,24}), with two new, higher
   values (48, 96) added for asced384 specifically because Step 1's
   measured natural convergence (45-141) sits well above the existing
   grid's ceiling (24) -- hard was re-run fresh at exactly these two values
   for a genuine matched-target comparison, not extrapolated.
3. **mpg=8 fixed throughout**, matching the existing hard recalibration's
   mpg, per the task's explicit instruction, for direct comparability --
   no separate mpg sweep was run for soft (out of scope here; the existing
   project's Sweep B already explores mpg trade-offs for hard).
4. **bias=10 fixed throughout** (not re-swept here), per the task's
   instruction -- already validated as comfortably in the FER-matching
   regime for both variants in `sweeps/soft_asced_delta_gt1_results.md`.
5. **asced48's Step-1/full-parallel-reference data was re-measured fresh**
   rather than reused from disk, since the original Delta>1 sweep's
   `RESULTS/` directory for asced48 was no longer present in this worktree
   (gitignored, partially cleaned between sessions) -- cross-checked
   against the published doc's numbers first (FER within Monte-Carlo noise)
   before trusting it as a valid substitute.
6. **The "latency" metric-definition caveat** (mpg's additive/sequential
   definition vs. full-parallel's single-max definition) is not a new
   discovery of this task, but it had not previously been stated explicitly
   enough to correctly interpret "does mpg recover full-parallel's latency
   problem" -- made explicit here since it materially changes the right
   answer to part (b) of the task's core question.

## Next steps

1. Fix the confirmed `set_syndrome_sequential` segfault (background task
   `task_f1f49ab5`) and re-run the `syndrome_sequential` comparison for soft
   once it's fixed -- the task's "nice-to-have" extension, not done here.
2. A genuine mpg sweep for soft (analogous to the existing hard Sweep B)
   could clarify whether a *smaller* `members_per_group` (fewer, bigger
   groups -> fewer sequential stages) lets soft's much higher natural
   convergence rate beat full-parallel's latency floor at asced384 scale
   too, the way it already does at asced48 scale in this data (see the
   asced48 sanity-check section's mechanism explanation) -- not attempted
   here, out of scope (mpg was held fixed at 8 per the task's instruction).
3. The RL-codebit-derived bias (vs. the free scalar hyperparameter used
   throughout) remains open and unaffected by this work
   (`sweeps/TODO_soft_asced.md`).
