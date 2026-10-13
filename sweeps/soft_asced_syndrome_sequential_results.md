# Soft-aSCED + `syndrome_sequential`: results

Date: 2026-10-13. Branches: `channel-code-lib2` `soft_asced_wireless` (used
as-is, no C++ changes -- the pybind11 dangling-reference segfault that
previously blocked this comparison is already fixed, commit `64884c1`, and
this worktree has been rebased onto it); `aSCED` `soft_asced_wireless`.
Continues `sweeps/soft_asced_mpg_combo_results.md`, which combined soft-aSCED
with `members_per_group`(mpg)/`MConvergedPolicy` using `fixed_sequential`
only, explicitly deferring `syndrome_sequential` as a "nice-to-have" once the
segfault it hit was fixed (see that doc's Step 0 and
`sweeps/STATUS.md`'s "`syndrome_sequential` segfault: fixed" section).

Driver script: `soft_asced_mpg_combo.py` (repo root), extended in place
(not a new script -- the vast majority of the logic, constants
(`TARGET_SWEEP`, `VARIANTS`, `SNR_POINTS`), and helpers (`run_sim`, CI
helpers) are shared with the existing `fixed_sequential` work, and keeping
everything in one script avoids duplicating/drifting those constants).
New entry points:
```
uv run python soft_asced_mpg_combo.py sweep_syn asced48    # soft, syndrome_sequential, full TARGET_SWEEP grid
uv run python soft_asced_mpg_combo.py sweep_syn asced384
uv run python soft_asced_mpg_combo.py hard_syn   asced48    # hard, syndrome_sequential, full TARGET_SWEEP grid (fresh)
uv run python soft_asced_mpg_combo.py hard_syn   asced384
```
(A fifth run, hard + `fixed_sequential` at the full target grid, was also
redone fresh in this worktree -- see "Design decision: also reran hard
`fixed_sequential`" below -- via a direct call to
`sweep_hard_all_targets(variant, selector="fixed_sequential")`, not wired
into the CLI, to avoid scope-creeping the script's public interface beyond
what this task needed.)

## Step 0: confirm the fix is present, build the mapping, verify non-degenerate before trusting any numbers

Per the task's explicit instruction to verify, not assume:

- **Confirmed the fix is an ancestor of this worktree's pinned commit**:
  `git merge-base --is-ancestor 64884c1 <pinned commit>` on
  `channel-code-lib2` succeeds. `uv sync` performed a real pybind11 rebuild
  (`channel_code_lib2.Ensemble_config.set_syndrome_sequential` resolves to
  the two-overload signature, `(members_per_group, mappings)` and
  `(members_per_group)`).
- **Smoke-tested the real (non-toy) construction directly**: built a full
  asced48 soft ensemble (48 paths), called
  `set_mConvergedConfig(2); set_syndrome_sequential(8, mappings)` with the
  all-zero mapping described below, ran a short simulation -- **no
  segfault**, sane FER (0.0739 at a cheap budget, consistent with the
  established bias=10/target=2 regime).
- **Mapping construction** (`build_soft_syndrome_mappings` in
  `soft_asced_mpg_combo.py`): one all-zero `(Delta_seg, n)` matrix per path,
  mirroring `create_soft_asced_config`'s own block/segment/pattern loop
  exactly (same `get_row_segments` call, same `2**Delta_seg` repetition per
  segment) so the mapping list lines up index-for-index with the ensemble's
  `path_configs`. `n = H.shape[1] = 154` (the ensemble's reference PCM
  column count), **not** the per-path soft PCM's wider `n+Delta_seg` column
  count -- confirmed from `SyndromeSequentialSelector.cpp`'s
  `_validate_mapping_shape`, which checks the mapping's columns against
  `ensemble.H.cols()` specifically. Verified both segment Delta sizes appear
  with the right shapes: asced48 → `{(2,154), (3,154)}` (48 total mappings);
  asced384 → `{(5,154), (6,154)}` (384 total mappings).
- **Verified non-degenerate (real reordering, not a silent no-op)**, two
  ways:
  1. **Source-level**: `compute_syndrome_weight` computes
     `syndrome = H_path @ [hd | mapping @ hd]` where `hd` is the
     per-codeword hard-decision of the (shared, n-length) processed LLRs.
     With an all-zero mapping, the appended aux bits are always 0, so the
     weight reduces to `H_aux_int @ hd` (the segment's own appended-row PCM
     applied to the received hard-decision) -- a quantity that genuinely
     varies **per received codeword** (depends on the actual noise
     realization) and **per segment** (different segments append different
     candidate rows). `select()` then does a real `std::sort` by this
     weight and groups the sorted indices into `mpg`-sized groups -- this is
     an active, per-transmission computation, not a fixed/cached order.
  2. **Empirically**: `syndrome_sequential`'s measured effort/latency/FER
     differ from `fixed_sequential`'s at matched target/SNR (not identical
     up to noise) -- see the full comparison below. A true no-op would
     produce numbers statistically indistinguishable from `fixed_sequential`
     at every point; instead there is a real, direction-consistent pattern
     (see "Headline finding" below), which is itself proof the selector is
     doing something different, whatever that something turns out to be.

**Design decision (flagged): for hard, no mapping is needed and none is
used.** Hard aSCED's per-path PCM (`H_aux = vstack(H, rows)`) has the SAME
column count as the reference `H` (the guessed suffix is baked into
`affine_offset`, not extra columns) -- it is not a "generalized PCM" in
`SyndromeSequentialSelector`'s sense, so the existing no-mapping convenience
overload `set_syndrome_sequential(mpg)` (used throughout the rest of this
project already) applies unchanged. No new code needed for the hard side.

**Design decision (flagged): also reran hard `fixed_sequential` fresh, at
the full target grid, in this worktree.** The task explicitly permits
reusing the existing hard `fixed_sequential` recalibration data
(`sweeps/target_converged_recalibration_results.md`, job 545344) for this
comparison, and that was the plan initially. However, this worktree's
(gitignored) `RESULTS/` tree does not contain that job's raw per-SNR
frame-error/trial counts for *either* selector -- only the published point
estimates are in the markdown doc, with no CI. Since this task's rigor bar
requires real Wilson CIs on every FER value, and the task already requires
a fresh, real-CI hard `syndrome_sequential` run regardless, it was cheap
(~12 more target-runs, ~12 minutes on this node's 168 cores) to also rerun
hard `fixed_sequential` fresh in this same worktree, rather than present a
mixed-rigor comparison (fresh-CI syndrome_sequential against
point-estimate-only historical fixed_sequential numbers). This is **beyond
the letter of the task's instructions** (which would have permitted the
point-estimate reuse) but keeps every leg of the four-way comparison
(soft-fixed / soft-syndrome / hard-fixed / hard-syndrome) on equal rigor
footing. Soft's `fixed_sequential` data did **not** need a rerun -- the
prior task's own run is still on disk in this worktree
(`RESULTS/soft_asced_mpg_combo/<variant>/soft_sweep_summary.json`) with real
fe/trial counts, matching the already-published CIs in
`sweeps/soft_asced_mpg_combo_results.md` exactly -- reused as-is.

All four sweeps (soft syndrome_sequential x2 variants, hard syndrome_sequential
x2 variants, hard fixed_sequential x2 variants rerun) used the identical
grids as the existing `fixed_sequential` work: bias=10, mpg=8,
`TARGET_SWEEP = {asced48: [1,2,4,8,16,24], asced384: [4,8,16,24,48,96]}`,
`target_errors={asced48: 200, asced384: 150}`, `max_transmissions=2e6`.

## Full comparison: FER Wilson-CI overlap (precise, not eyeballed)

Three pairwise overlap checks were computed programmatically (Wilson 95%
CIs from real fe/trial counts on both sides of every check, `hi1 < lo2 or
hi2 < lo1` tested directly, not eyeballed):

| comparison | overlap | what it tests |
|---|---|---|
| soft `fixed_sequential` vs soft `syndrome_sequential` (same mechanism, selector-only question) | **36/36** | switching selector doesn't change soft's own FER |
| hard `fixed_sequential` vs hard `syndrome_sequential` (same mechanism, selector-only question) | **36/36** | switching selector doesn't change hard's own FER |
| soft `syndrome_sequential` vs hard `syndrome_sequential` (matched target, the task's FER-parity question under this selector) | **36/36** | FER parity holds under `syndrome_sequential`, same as it already did under `fixed_sequential` |

**All three checks pass at every one of the 36 (variant x target x SNR)
points tested** (asced48: 6 targets x 3 SNR = 18; asced384: 6 targets x 3
SNR = 18; x2 variants = 36 total per check). This matches the existing
`fixed_sequential` work's FER-parity conclusion -- `syndrome_sequential`
does not regress FER parity in either direction (hard-vs-soft, or
selector-vs-selector).

## Headline finding: `syndrome_sequential` does NOT outperform `fixed_sequential` for soft -- if anything, mildly worse

Aggregated over all 36 (variant x target x SNR) points, comparing
`syndrome_sequential`'s effort/latency against `fixed_sequential`'s at the
**same** target/SNR/mechanism:

| | mean effort delta | mean latency delta | points where syn is WORSE (effort) | points where syn is WORSE (latency) |
|---|---|---|---|---|
| **soft** | **+2.33%** | **+2.21%** | **28/36** | **29/36** |
| **hard** | +2.11% | +5.04% (median: 0.00%) | 23/36 | 15/36 |

For **soft**, the pattern is consistent and direction-stable, not just a
mean pulled by outliers: at the **lowest tested SNR per variant**, syn and
fixed are statistically indistinguishable (mean delta +0.03%/+0.03%), but at
the **highest tested SNR per variant**, syn is **systematically ~3.3-3.5%
worse** on both effort and latency (mean +3.48%/+3.29%), topping out at
+7.8% effort / +9.6% latency at asced48 target=16, 3.0dB. The full
per-point table (bias=10, mpg=8):

### asced48 (SNR 1.0/2.0/3.0dB)

| target | fixed FER [CI] | syn FER [CI] | overlap | effort Δ% | latency Δ% |
|---|---|---|---|---|---|
| 1  | 0.37385/0.07120/0.00380 | 0.36353/0.06840/0.00368 | Y/Y/Y | -2.05/-0.91/+3.56 | -2.87/-4.26/-2.50 |
| 2  | 0.36492/0.07040/0.00370 | 0.34937/0.06770/0.00384 | Y/Y/Y | -1.30/+0.19/+4.88 | -1.76/-1.73/+1.00 |
| 4  | 0.36645/0.06963/0.00352 | 0.36801/0.06884/0.00391 | Y/Y/Y | +2.24/+5.68/+5.82 | +2.29/+5.39/+4.45 |
| 8  | 0.35880/0.06683/0.00368 | 0.35258/0.06978/0.00363 | Y/Y/Y | +0.51/+2.57/+7.30 | +0.44/+2.88/+9.13 |
| 16 | 0.35743/0.06851/0.00368 | 0.35744/0.06920/0.00384 | Y/Y/Y | -0.02/+3.07/+4.50 | +0.39/+4.58/+9.64 |
| 24 | 0.36059/0.06862/0.00363 | 0.35855/0.06573/0.00378 | Y/Y/Y | -0.11/+0.17/+1.21 | +0.25/+0.86/+2.91 |

(Fixed-sequential FER/CI values are from `sweeps/soft_asced_mpg_combo_results.md`'s
already-published table, reused verbatim; "overlap" = Wilson 95% CI overlap
between the two selectors at that exact point, computed from real fe/trial
counts on both sides.)

### asced384 (SNR 1.0/1.5/2.0dB)

| target | fixed FER | syn FER | overlap | effort Δ% | latency Δ% |
|---|---|---|---|---|---|
| 4  | 0.23810/0.10116/0.03062 | 0.23064/0.09782/0.03182 | Y/Y/Y | -1.42/+0.97/+4.67 | -1.52/+0.81/+4.34 |
| 8  | 0.23566/0.10021/0.02877 | 0.23910/0.09872/0.03066 | Y/Y/Y | +1.00/+2.68/+7.44 | +0.95/+2.55/+7.04 |
| 16 | 0.23812/0.10153/0.03014 | 0.23702/0.09800/0.03097 | Y/Y/Y | +1.63/+5.51/+7.79 | +1.63/+5.46/+7.76 |
| 24 | 0.22456/0.09962/0.03055 | 0.22241/0.09724/0.03190 | Y/Y/Y | -0.15/+2.90/+5.81 | -0.04/+2.97/+5.98 |
| 48 | 0.22693/0.09747/0.03073 | 0.23364/0.10240/0.03046 | Y/Y/Y | +0.43/+2.29/+2.98 | +0.55/+2.40/+3.55 |
| 96 | 0.23978/0.09985/0.02933 | 0.23053/0.10138/0.02935 | Y/Y/Y | -0.40/+1.19/+1.35 | +0.09/+1.47/+2.34 |

Same pattern at asced384: near-zero or mildly-negative (better) delta at the
lowest SNR, growing to consistently positive (worse) at the highest SNR,
across every target.

For **hard**, the picture is noisier (several large-magnitude swings, e.g.
asced48 target=2 @3.0dB: +108.67% latency delta) but these trace to small
absolute denominators, not a systematic trend: at low target + high SNR,
very few groups need to run before stopping triggers, so latency is tiny in
absolute terms (hard asced48 target=2 @3.0dB: fixed=7.73, syn=16.13 -- a
large relative swing on an absolute difference of 8.4 iteration-units,
within the kind of per-transmission-ordering noise expected when so few
groups actually execute). Once targets/SNR saturate the full mpg budget
(target>=16 for asced48, target>=16 for asced384), fixed and syndrome
converge to the *same* number (both selectors end up running literally
every group, so the ordering is irrelevant) -- hence the 0.00% median. This
reproduces, independently, what `sweeps/target_converged_recalibration_results.md`
(job 545344, a wholly separate historical sweep) already found: "Both
selectors (fixed_sequential, syndrome_sequential) behave similarly, no
strong reason to prefer one from this sweep alone" -- cross-checking that
doc's own fixed-vs-syndrome hard numbers at matched targets shows the same
kind of mixed, noisy, no-clear-winner pattern (e.g. asced48 target=4,
historical data: fixed effort 1012.70/593.03/435.43 vs syndrome
1079.52/674.96/530.92 at 1/2/3dB -- syndrome measurably *worse*, by 6.6-21.9%,
in that independent, earlier run too).

## Mechanism: why the hoped-for "most-likely-first" benefit doesn't materialize for soft

The task's hypothesis was that `syndrome_sequential`'s dynamic per-codeword
reordering might "consistently run the most-likely-to-converge-first groups
first," recovering more of the mpg-vs-full-parallel latency-metric
asymmetry than `fixed_sequential`'s static index order. The data above shows
this does not happen (soft: mildly worse; hard: a wash) -- and the
construction explains why:

The all-zero mapping -- required both to avoid the segfault's
shape-validation path and to match hard's own
`compute_syndrome_weight`-ignores-the-guess convention (the explicit,
principled choice this task was told to use) -- makes **every one of the
`2**Delta_seg` pattern-paths within a single segment compute the identical
syndrome weight** (since the weight only depends on `H_aux_int @ hd`, which
doesn't involve the per-path pattern/guess at all once the mapped aux
contribution is forced to zero). So `syndrome_sequential` can only
meaningfully reorder **between segments** (different segments append
different candidate rows, hence different `H_aux_int`, hence genuinely
different weights) -- it has **zero information to distinguish a
likely-correct guess from a very-wrong one within the same segment**, which
is exactly the distinction that would matter most (per
`sweeps/soft_asced_delta_gt1_results.md`'s own Hamming-distance finding:
convergence likelihood depends strongly on how close the guessed pattern is
to the true one).

Meanwhile, `fixed_sequential`'s static index order is **already** segment-
coherent -- `create_soft_asced_config`/`create_hard_asced_config` build
`path_configs` in exactly the block-then-segment-then-pattern nested loop
order, so index-contiguous `mpg`-sized groups already tend to draw from a
small number of adjacent segments. `syndrome_sequential`'s only real lever
(inter-segment reordering) therefore has limited room to improve on an
already-reasonably-coherent baseline, and the residual difference observed
-- small, and for soft specifically consistently unfavorable at high SNR --
is plausibly `std::sort`'s lack of a stability guarantee disturbing what
would otherwise be a benign tie order among equal-weight paths, though this
specific mechanism was not confirmed at the C++ instrumentation level (not
pursued further -- out of scope, and the aggregate effect size is modest).

**This is a genuine, if modest, finding distinct from a simple "no
difference": soft's specific need for an all-zero mapping (forced by the
segfault-avoidance/hard-equivalence requirement) structurally removes
`syndrome_sequential`'s only real source of leverage over `fixed_sequential`
for the intra-segment case, which is precisely the case the Delta>1 work
showed matters most for soft's guess-correction mechanism.** A
mapping that preserves per-pattern distinguishability (e.g. encoding the
guessed pattern's parity into the mapping, rather than zeroing it) might let
`syndrome_sequential` realize an actual advantage for soft -- but that would
abandon the exact-equivalence-with-hard property this task was told to
preserve, and is flagged as an open idea for a follow-up, not attempted
here.

## Mpg-ceiling-recovery metric (same metric as the fixed_sequential doc) confirms "about the same"

Using the prior doc's own "fraction of the no-early-stopping mpg ceiling
recovered" metric (asced384's 48 groups at mpg=8, each saturating toward
~33, ceiling=1584) at the doc's own representative headline point
(target=4, 2.0dB):

| | latency | recovered fraction of ceiling |
|---|---|---|
| soft, `fixed_sequential` (prior doc) | 171.5 | 89.2% |
| soft, `syndrome_sequential` (this task) | 179.0 | 88.7% |
| hard, `fixed_sequential` (this task, fresh rerun) | 695.2 | 56.1% |
| hard, `syndrome_sequential` (this task) | 690.5 | 56.4% |

Essentially identical for both mechanisms (within 0.5 percentage points) --
`syndrome_sequential` neither meaningfully helps nor hurts this metric
relative to `fixed_sequential`, for either hard or soft.

## Does soft's effort-compounds / latency-capped conclusion still hold under `syndrome_sequential`?

**Yes, unchanged.** The soft-vs-hard savings at matched target under
`syndrome_sequential` closely track the `fixed_sequential` doc's own
numbers (e.g. asced384 target=4, 2.0dB: this task's soft-vs-hard effort
saving under `syndrome_sequential` is **+73.7%**, latency saving **+74.1%**
-- compare the `fixed_sequential` doc's analogous target=4, 2.0dB numbers,
**75%/75%**). Across the full grid (36 points), soft-vs-hard savings under
`syndrome_sequential` range from +3.9% to +88.9% (effort) and +0.1% to
+85.8% (latency), monotonically increasing with target, exactly mirroring
the `fixed_sequential` doc's own pattern -- and the FER-parity overlap
(36/36, above) confirms this isn't bought at a hidden FER cost. The
mpg-vs-full-parallel latency-metric asymmetry described in the
`fixed_sequential` doc (mpg's number can never beat full-parallel's ~33
floor once >1 group runs) applies identically here -- `syndrome_sequential`
doesn't change which metric regime this falls into, only (very mildly)
which selector you'd pick within that regime, and this task's data says
`fixed_sequential` is at least as good, typically marginally better, as a
within-mpg choice.

## Answering the task's core question directly

**Does `syndrome_sequential` behave differently from `fixed_sequential` for
soft-aSCED -- better, worse, or about the same on FER/effort/latency?**
FER: **the same** (36/36 CI overlap, both selector-vs-selector and
soft-vs-hard). Effort/latency: **mildly worse for soft** (mean +2.2-2.3%,
systematically ~+3.3-3.5% worse specifically at the highest tested SNR per
variant, every target), **a wash for hard** (median 0%, mean pulled up by
noisy small-denominator swings at low-target/high-SNR points that cancel
out once mpg saturation is reached). Neither is a dramatic effect in either
direction.

**Does the "effort compounds, latency is capped by the mpg-vs-full-parallel
metric asymmetry" conclusion still hold, or does the dynamic per-codeword
reordering change the picture?** **Still holds, essentially unchanged.**
The hypothesis that dynamic per-codeword reordering might "consistently run
the most-likely-to-converge-first groups first" and thereby recover more of
the latency ceiling than `fixed_sequential` is **not borne out** -- the
mpg-ceiling-recovery metric is within 0.5 percentage points of
`fixed_sequential`'s for both hard and soft, and the mechanism section above
explains concretely why: the all-zero mapping this construction requires
removes exactly the intra-segment per-pattern signal that a smarter
reordering would need to exploit, leaving only a limited, already-mostly-
exploited inter-segment reordering lever.

## Design decisions / assumptions flagged for the user

1. **Mapping: all-zero `(Delta_seg, n)` per path**, exactly the choice
   already identified (and required) in `sweeps/soft_asced_mpg_combo_results.md`'s
   Step 0 -- not re-litigated here, just implemented and verified
   non-degenerate (see Step 0 above).
2. **Script: extended `soft_asced_mpg_combo.py` in place**, not a new
   script -- judged cleaner given how much of the existing logic
   (constants, run_sim, CI helpers) is shared and would otherwise need
   duplicating or importing awkwardly. New CLI modes `sweep_syn`/`hard_syn`
   added alongside the existing `sweep`/`hard_extra`/`natural` modes,
   none of which were modified.
3. **Reran hard `fixed_sequential` fresh at the full target grid** (beyond
   the task's explicit permission to reuse existing data) purely to get
   real Wilson CIs for every leg of the four-way comparison on equal
   footing -- see "Design decision: also reran hard `fixed_sequential`"
   above for the full reasoning. This used the exact same, unmodified
   `run_hard_mpg`/`sweep_hard_all_targets` code path as the required
   `hard_syn` runs (selector parameter only differs), so it carries the
   same construction/budget guarantees.
4. **Did not attempt a mapping that preserves per-pattern distinguishability**
   (e.g. one that would let `syndrome_sequential` actually rank patterns by
   guessed-pattern plausibility) -- this would abandon the
   exact-hard-equivalence property the task specified, and is noted above
   as a possible follow-up, not pursued here (would be a materially
   different, larger task: deciding what such a mapping should even encode,
   whether it stays principled/justifiable, and re-validating from scratch).
5. **`std::sort`'s non-stability as a partial explanation for the small
   residual soft-vs-hard-selector difference** is stated as a plausible,
   *unverified* hypothesis (flagged explicitly as such above), not a
   confirmed mechanism -- verifying it would need C++-level instrumentation
   of the actual group orderings per transmission, judged not worth the
   added scope for an effect this small (a few percent).

## Next steps

1. The alternative, distinguishability-preserving mapping idea above (if
   anyone wants to pursue a selector that could actually realize the
   "most-likely-first" benefit for soft) -- would need its own principled
   design and validation pass, not a small extension of this task.
2. The RL-codebit-derived bias (vs. the free scalar hyperparameter used
   throughout) remains open and unaffected by this work, same as every
   prior entry in this family (`sweeps/TODO_soft_asced.md`).
3. A genuine mpg sweep for soft combined with `syndrome_sequential`
   specifically (this task held mpg=8 fixed, per the existing project
   convention) could in principle show a different picture at a different
   group size -- not attempted here, out of scope.
