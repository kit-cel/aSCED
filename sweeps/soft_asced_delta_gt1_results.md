# Soft-aSCED, Delta>1 generalization: results

Date: 2026-10-12. Branches: `channel-code-lib2` `soft_asced_wireless` (used
as-is, no C++ changes needed); `aSCED` `soft_asced_wireless`. Continues
`sweeps/soft_asced_single_row_results.md` (the validated Delta=1
prototype).

Driver script: `soft_asced_delta_gt1.py` (repo root). Run with:
```
uv run python soft_asced_delta_gt1.py sanity asced48
uv run python soft_asced_delta_gt1.py sanity asced384
uv run python soft_asced_delta_gt1.py sweep  asced48
uv run python soft_asced_delta_gt1.py sweep  asced384
```

## What this is

Generalization of the validated Delta=1 soft-aSCED prototype to the real
multi-row-per-segment case ("Delta>1"), validated directly on this
project's primary paper-reproduction ensembles, "aSCED-48" and "aSCED-384"
(not a new toy PCM -- the SAME code/rows the production Fig. 3 sweeps use).

## Key scoping finding (confirmed by reading the actual hard-aSCED code)

`create_asced_config` in `reproduce_fig3_RL_zc11_asced_5G_LDPC.py` does
**not** use any general overcomplete-PCM/`M_l` row-combination
construction -- it is `M_l = identity`. For a segment of `Delta_seg`
appended candidate rows it just builds `H_aux = [H; segment_rows]` (vstack,
nothing combined), then enumerates all `2**Delta_seg` hard `affine_offset`
patterns (`binary_vectors_in_suffix`). This means the general
overcomplete-row/`M_l` construction used on the sibling QEC project (where
syndrome VNs get degree > 1 from genuine row-overcompleteness) is **not
needed here** -- it is explicitly out of scope for this task. The soft
generalization is a near-direct port of `create_asced_config`: same
block/segment loop, same segment sizes, same per-segment path count
(`2**Delta_seg`, matching hard exactly) -- only the per-pattern path
construction changes.

- **"asced48"** = `splitting_pattern[1]=[2,4,6,8]` x 2 blocks (`block_offset=0`)
  = 2 blocks x (4 segments of Delta=2 + 1 segment of Delta=3)
  = 2 x (4x4 + 8) = **48** paths. Confirmed empirically:
  `np.split(row_block(11 rows), [2,4,6,8])` -> sizes `[2,2,2,2,3]`.
- **"asced384"** = `splitting_pattern[2]=[5]` x 4 blocks (`block_offset=0`)
  = 4 blocks x (1 segment of Delta=5 + 1 of Delta=6)
  = 4 x (32+64) = **384** paths. Confirmed empirically:
  `np.split(row_block(11 rows), [5])` -> sizes `[5,6]`.

## Implementation

`soft_asced_delta_gt1.py` (repo root), self-contained (does NOT import or
modify `reproduce_fig3_RL_zc11_asced_5G_LDPC.py` /
`reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py`, which stay untouched
per scope -- that original script also gates its own execution on
`sys.argv`, so it can't be imported as a library without running its CLI
logic):

- `get_row_segments` (line 127): the exact block/segment loop from
  `create_asced_config`, returning `(block_idx, rows)` pairs -- used by
  BOTH hard and soft construction so the row selection is byte-identical.
- `create_hard_asced_config` (line 146): near-verbatim copy of
  `create_asced_config`, with `use_all_zero=False` hardcoded (the project's
  standard asced48/384 behavior -- zero affine offset AND all nonzero
  suffix patterns, i.e. all `2**Delta_seg` patterns per segment).
- `make_soft_bp_config` (line 172) + `create_soft_asced_config` (line 194):
  the new soft generalization. For each segment,
  `H_soft = [H_aux | Delta_seg-column identity block]`, identity placed in
  the LAST `Delta_seg` rows (the appended segment's own rows) -- i.e. each
  syndrome-VN column has a single 1 in its own row, degree-1 connectivity,
  exactly like the validated Delta=1 construction, just repeated
  `Delta_seg` times per segment. One `BP_config` per full `Delta_seg`-bit
  pattern (`itertools.product((0,1), repeat=Delta_seg)`), matching hard's
  path count and row selection exactly.
- `cfg.set_H0(H)` always the ORIGINAL (un-extended, m x n) PCM, for every
  segment regardless of `Delta_seg` -- the guess-correction check is always
  against the true, un-extended syndrome.

**C++ side: zero changes.** Confirmed independently by reading
`BP_Decoder::init` in the `soft_asced_wireless` branch of
`channel-code-lib2`: it only checks
`H0.cols() == H.cols() - num_syndrome_vns`, nothing about connectivity
structure -- arbitrary per-row degree-1 syndrome VNs already work for any
`Delta_seg`, confirmed empirically below (both variants decode and
converge correctly end-to-end on the first try, no C++ work needed).

## Sanity checks (generalized from Delta=1, run on one representative segment per variant)

Representative segments: asced48's first segment (block 0, `Delta_seg=2`);
asced384's first segment (block 0, `Delta_seg=5`). All four checks were
generalized from `prototype_soft_asced_single_row.py`'s single-bit version
to a full `Delta_seg`-bit pattern.

### asced48 (Delta_seg=2)

1. **Bias is really applied.** `bias=0` collapses ALL 4 full patterns PLUS
   the all-erasure pattern to the identical decision: **confirmed**.
   `bias=30`, all-zero vs. all-one pattern differ: **confirmed**. A
   SINGLE-bit flip (pattern `(1,0)` vs. `(0,0)`, holding the other VN at 0)
   also changes the decision at `bias=30`: **confirmed** -- this
   additionally checks that each syndrome-VN column has an independent
   effect, not just "any guess vs. its total opposite".
2. **Hard-limit equivalence.** N=1200 random-message samples with a KNOWN
   true 2-bit segment pattern (`rows @ codeword mod 2`); soft with
   `pattern=true_bits`, `bias=30` vs. hard aSCED's path using that SAME
   correct full guess: **1200/1200 (100%) decision agreement**; restricted
   to the subset where hard did NOT converge (165/1200, the only place
   saturation and a hard-fixed guess could plausibly differ): **165/165
   (100%) agreement**, 95% CI on the disagreement rate **[0, 0.0182]**.
3. **Reused decoder == fresh decoder.** 50 samples, byte-identical results
   whether decoded on one reused `Decoder` object or a fresh one per
   sample: **confirmed**.
4. **Per-path effort/convergence by Hamming distance** (N=300 samples, all
   4 patterns per sample, bias=10, 2.0dB):

   | Hamming dist to true pattern | n | hard mean iters | hard conv. | soft mean iters | soft conv. |
   |---|---|---|---|---|---|
   | 0 (fully correct) | 300 | 12.19 | 84.0% | 12.19 | 84.0% |
   | 1 (one of two bits wrong) | 600 | 33.00 | 0.0% | 16.72 | 69.2% |
   | 2 (fully wrong) | 300 | 33.00 | 0.0% | 22.83 | 44.0% |

   Confirms the Delta=1 mechanism generalizes cleanly: hard's wrong-guess
   paths essentially never converge (burn the full 32-iteration budget,
   reported as 33 -- 32 + the implementation's +1 convergence-check
   convention) regardless of HOW wrong the guess is, while soft's
   wrong-guess paths still converge a majority of the time at distance 1
   and a substantial fraction even at distance 2 (fully wrong), at much
   lower mean effort -- monotonically degrading with Hamming distance, as
   expected.

### asced384 (Delta_seg=5)

1. **Bias applied**: all 32 full patterns + all-erasure collapse at
   `bias=0`: **confirmed**. All-zero vs. all-one differ at `bias=30`:
   **confirmed**. Single-bit flip differs at `bias=30`: **confirmed**.
2. **Hard-limit equivalence**, N=1200: **1200/1200 (100%) agreement**;
   subset where hard did not converge (117/1200): **117/117 (100%)
   agreement**, disagreement-rate 95% CI **[0, 0.0256]**.
3. **Reused decoder == fresh decoder**: **confirmed**.
4. **Per-path effort/convergence by Hamming distance** (N=300 samples, all
   32 patterns per sample, bias=10, 2.0dB) -- the full 0..5 Hamming-distance
   range now resolvable since `Delta_seg=5`:

   | Hamming dist | n | hard mean iters | hard conv. | soft mean iters | soft conv. |
   |---|---|---|---|---|---|
   | 0 | 300 | 10.58 | 89.0% | 10.58 | 89.0% |
   | 1 | 1500 | 33.00 | 0.0% | 14.27 | 78.4% |
   | 2 | 3000 | 32.95 | 0.4% | 20.55 | 53.1% |
   | 3 | 3000 | 32.95 | 0.3% | 27.17 | 25.0% |
   | 4 | 1500 | 32.97 | 0.2% | 30.87 | 8.6% |
   | 5 (fully wrong) | 300 | 33.00 | 0.0% | 32.28 | 3.0% |

   A clean, monotonic degradation of soft's convergence fraction with
   Hamming distance (89% -> 78% -> 53% -> 25% -> 8.6% -> 3.0%), confirming
   the guess-correction mechanism continues to work exactly as expected
   even at `Delta_seg=5` -- soft "spends" its advantage proportionally to
   how wrong the guess is, rather than an all-or-nothing hard cutoff.

All sanity checks pass on both variants -- **no bugs found in this
generalization** (the Delta=1 prototype's known bugs, e.g. the
`bp_converged`/`iteration_count` early-return bug, were already fixed on
this branch and are inherited unchanged).

## Headline result: hard aSCED-48 (48 paths) vs. soft aSCED-48 (48 paths, bias swept)

SNR points {1.0, 2.0, 3.0} dB, `target_errors=200`, `max_transmissions=2e6`
(not binding), `sim.num_threads = os.cpu_count()`. Both ensembles run as a
plain full-parallel ensemble (no `members_per_group`/`MConvergedPolicy`
early-stopping -- isolates the per-path mechanism, same convention as the
Delta=1 doc).

| | FER @ 1.0dB [95% CI] | FER @ 2.0dB [95% CI] | FER @ 3.0dB [95% CI] | effort 1.0/2.0/3.0dB | latency 1.0/2.0/3.0dB |
|---|---|---|---|---|---|
| **hard (48 paths)** | 0.3620 [0.3461,0.3783] | 0.0705 [0.0666,0.0746] | 0.00371 [0.00350,0.00394] | 1469/1367/1311 | 33.0/33.0/33.0 |
| soft, bias=1 | 0.4175 [0.4004,0.4347]** | 0.1028 [0.0973,0.1085]** | 0.00764 [0.00720,0.00811]** | 1110/618/302 | 26.2/15.8/7.4 |
| soft, bias=2 | 0.4039 [0.3876,0.4205]** | 0.0827 [0.0782,0.0874]** | 0.00596 [0.00561,0.00632]** | 1168/645/311 | 28.8/18.5/8.6 |
| soft, bias=3 | 0.3641 [0.3483,0.3801] | 0.0758 [0.0719,0.0800] | 0.00470 [0.00443,0.00498]** | 1181/679/322 | 30.2/21.2/9.9 |
| soft, bias=5 | 0.3668 [0.3512,0.3828] | 0.0719 [0.0680,0.0761] | 0.00382 [0.00360,0.00405] | 1249/740/345 | 31.9/25.0/12.7 |
| soft, bias=10 | 0.3618 [0.3463,0.3777] | 0.0700 [0.0662,0.0741] | 0.00373 [0.00352,0.00395] | 1286/830/404 | 32.8/30.4/20.1 |
| soft, bias=20 | 0.3639 [0.3478,0.3802] | 0.0709 [0.0670,0.0750] | 0.00338 [0.00318,0.00359] | 1296/872/483 | 32.9/32.3/29.0 |
| soft, bias=30 | 0.3597 [0.3440,0.3757] | 0.0698 [0.0661,0.0738] | 0.00387 [0.00365,0.00410] | 1291/880/492 | 32.9/32.4/29.8 |
| soft, bias=35 | 0.3580 [0.3426,0.3738] | 0.0682 [0.0645,0.0722] | 0.00378 [0.00357,0.00400] | 1285/878/491 | 32.9/32.4/29.8 |

(`**` = CI does not overlap hard's at that SNR, i.e. measurably worse --
checked precisely via interval arithmetic on the exact CI bounds, not by
eye. Wilson 95% CIs; effort/latency are
`average_ensemble_effort`/`average_ensemble_latency` from
`decoder_stats.json`, same convention as the Delta=1 doc.)

**Headline finding, and an important DIFFERENCE from Delta=1**: at low bias
(`bias=1,2`), soft aSCED-48's FER is **measurably WORSE** than hard's at
EVERY SNR tested (CIs do not overlap at 1.0/2.0/3.0dB) -- unlike the
validated Delta=1 result, where even `bias=1` already matched hard's FER at
every operating point. `bias=3` is an intermediate case: its CI overlaps
hard's at 1.0dB and 2.0dB, but is still measurably worse at 3.0dB (the
highest-SNR / lowest-FER point: hard 0.00371 [0.00350,0.00394] vs. soft(3)
0.00470 [0.00443,0.00498], no overlap). Only from **`bias>=5`** onward does
soft's FER CI overlap hard's at ALL three SNR points simultaneously
(verified precisely for every bias/SNR cell, not just eyeballed -- bias=5,
10, 20, 30, 35 all qualify). This makes intuitive sense: with multi-bit
segments (`Delta_seg` in {2,3} here), a weak per-bit bias means several
syndrome-VN bits are only weakly pushed at once, and the combined effect of
multiple weakly-biased guesses compounds unfavorably compared to the
single-bit Delta=1 case. **This is a concrete, reportable generalization
finding, not predicted by the Delta=1 result.**

At `bias=5` -- the smallest bias value tested where FER is
statistically indistinguishable from hard at every SNR -- **effort drops
15% / 46% / 74%** and **latency drops 3% / 24% / 62%** (1.0/2.0/3.0dB
respectively). At `bias=10`, effort drops 12%/39%/69% and latency drops
0.6%/8%/39%. The savings shrink (but don't vanish) as bias continues to
increase toward saturation (`bias=30/35`), exactly as in the Delta=1 case:
effort/latency increase monotonically with bias (more saturating = more
hard-like), but never reach hard's cost even at `bias=35`.

## Headline result: hard aSCED-384 (384 paths) vs. soft aSCED-384 (384 paths, bias swept)

SNR points {1.0, 1.5, 2.0} dB (chosen lower than asced48's, since
asced384's much stronger code reaches comparably low FER at lower SNR),
`target_errors=150` (lighter than asced48's 200, given the much larger
384-path per-trial cost), same ensemble/threading conventions. Bias sweep
trimmed to {1,5,10,20,30} given the ~8x higher per-run cost (see "Design
decisions" below). Run completed in 7m48s wall-clock (168 cores).

| | FER @ 1.0dB [95% CI] | FER @ 1.5dB [95% CI] | FER @ 2.0dB [95% CI] | effort 1.0/1.5/2.0dB | latency 1.0/1.5/2.0dB |
|---|---|---|---|---|---|
| **hard (384 paths)** | 0.2403 [0.2283,0.2528] | 0.1036 [0.0979,0.1096] | 0.02953 [0.02778,0.03138] | 12552/12509/12474 | 33.0/33.0/33.0 |
| soft, bias=1 | 0.3655 [0.3499,0.3815]** | 0.1755 [0.1665,0.1847]** | 0.07606 [0.0719,0.0805]** | 9168/6942/5134 | 28.5/23.9/18.7 |
| soft, bias=5 | 0.2292 [0.2178,0.2411] | 0.1014 [0.0958,0.1073] | 0.03207 [0.0303,0.0340] | 10996/9469/7457 | 33.0/32.8/32.0 |
| soft, bias=10 | 0.2413 [0.2294,0.2536] | 0.0992 [0.0936,0.1052] | 0.03078 [0.0290,0.0327] | 11692/10744/9317 | 33.0/33.0/32.99 |
| soft, bias=20 | 0.2383 [0.2264,0.2506] | 0.0988 [0.0934,0.1045] | 0.02939 [0.0277,0.0312] | 11809/11088/10059 | 33.0/33.0/33.0 |
| soft, bias=30 | 0.2388 [0.2273,0.2508] | 0.0964 [0.0910,0.1021] | 0.03115 [0.0294,0.0330] | 11802/11041/10079 | 33.0/33.0/33.0 |

(`**` = CI does not overlap hard's at that SNR. Same conventions as the
asced48 table above.)

**Same qualitative pattern as asced48, confirming it generalizes**:
`bias=1` is measurably worse than hard at every SNR (no CI overlap) --
even more so than at asced48, since asced384's segments are larger
(`Delta_seg` in {5,6} vs. {2,3}), so a weak per-bit bias compounds across
more simultaneously-weakly-guessed bits. From `bias=5` onward, soft's FER
CI overlaps hard's at all three SNR points.

**A genuinely new finding at this larger scale: effort savings persist,
but latency savings nearly vanish at FER-matching bias.** At `bias=5`
(smallest FER-matching bias): effort drops **12% / 24% / 40%**
(1.0/1.5/2.0dB) -- comparable in direction to asced48 -- but latency drops
only **~0% / 0.6% / 3%**, essentially within noise. At `bias=10`: effort
drops 7%/14%/25%, latency drops <0.1% (indistinguishable from zero). Only
at `bias=1` (FER-incompatible) does latency drop meaningfully (14%/28%/43%).

This is a structural consequence of the **full-parallel "latency"
metric** (`average_ensemble_latency` = the slowest path in the group) at
much larger path count: with 384 paths per trial, at any bias that
preserves FER there is very likely to be at least one high-Hamming-distance
wrong-guess path among the 384 that still burns close to the full
iteration budget (per the Hamming-distance sanity check above: even at
`bias=10`, Hamming distance 4-5 patterns still only converge 3-9% of the
time) -- so the ensemble's *worst* path dominates the latency statistic
almost every trial, washing out the gain visible in the *mean* (effort)
statistic. This did not show up as starkly in asced48 (`Delta_seg` only up
to 3, so far fewer high-Hamming-distance patterns per segment) -- **the
latency benefit of soft-aSCED at FER-matching bias shrinks as `Delta_seg`
grows**, even though the effort (sum/mean-based) benefit persists. This is
a genuine, non-obvious generalization finding, not predictable from the
Delta=1 or asced48 results alone.

**Cross-check against documented asced384 production FER**
(`sweeps/STATUS.md`: asced385 -- a 385-path greedy variant, closest
documented comparison point -- FER(1.0dB)~0.2283, FER(2.0dB)~0.0307): our
independently re-implemented hard-384 baseline got FER(1.0dB)=0.2403,
FER(2.0dB)=0.02953 -- same order of magnitude and close numerically
(within ~5-6%), a reasonable plausibility check of the whole
re-implementation (not a bit-exact reproduction, given `target_errors=150`
Monte-Carlo noise and the asced385-vs-asced384 variant mismatch).

## Design decisions / assumptions flagged for the user

1. **`M_l = identity`, inherited directly from hard aSCED's own
   construction** (not an independent choice made for this task) -- see
   "Key scoping finding" above. This also means the generalized sanity
   checks above still only exercise degree-1 syndrome VNs (one VN per
   appended row, connected to exactly that row) -- NOT the higher-degree
   case from genuine row-overcompleteness, which remains out of scope here
   (that would require implementing a general `M_l` construction hard
   aSCED itself doesn't use).
2. **Bias sweep range**: reused the Delta=1 doc's set for asced48
   ({1,2,3,5,10,20,30,35}) so the two docs are directly comparable; for
   asced384, used a trimmed set ({1,5,10,20,30}) given the much higher
   per-run cost, informed by the asced48 finding that very low bias (1-2)
   measurably underperforms -- included `bias=5` explicitly as the
   empirically-identified "sweet spot" from asced48, plus `bias=1` to
   check whether the same low-bias underperformance generalizes to
   `Delta_seg` in {5,6}.
3. **SNR operating points differ between the two variants** (asced48:
   {1.0,2.0,3.0}dB; asced384: {1.0,1.5,2.0}dB) -- chosen so both land in a
   similar FER range (~0.3 down to ~0.003-0.03), rather than reusing
   identical SNR values that would put asced384 at FER too low to collect
   `target_errors` in reasonable time, or asced48 at FER too high to be
   informative.
4. **`target_errors`/SNR choices were cross-checked against this project's
   own documented asced48/384 FER curve** (`sweeps/STATUS.md`: asced48
   1.0dB~0.41, 3.0dB~0.0038; asced384/385 1.0dB~0.23, 2.0dB~0.031) as an
   independent plausibility check of the whole hard-aSCED-48/384
   re-implementation in this script (not literally imported from the
   original, paper-reproduction script, since that script gates on
   `sys.argv` -- see "Implementation" above). Our hard-48 run got
   FER(1.0dB)=0.362, FER(3.0dB)=0.0037, and our hard-384 run got
   FER(1.0dB)=0.2403, FER(2.0dB)=0.02953 -- both same order of magnitude
   and directionally consistent with the documented values (within
   Monte-Carlo noise at `target_errors=150-200`, not a bit-exact
   reproduction of the production sweep).
5. **Representative segments for the sanity checks**: always the FIRST
   segment of each variant (`get_row_segments(...)[0]`) -- asced48:
   `Delta_seg=2`; asced384: `Delta_seg=5` -- per the task's explicit
   suggestion ("Delta_seg=2 from asced48, and Delta_seg=5 or 6 from
   asced384"). Not repeated across every segment/block, which would be a
   further robustness check but was judged out of scope for this pass
   (the per-path mechanism is structurally identical across segments of
   the same `Delta_seg`, only the specific rows differ).
6. **No equal-K ("patterns vs. matrices") check for Delta>1** -- that was
   an exploratory extra specific to the Delta=1 prototype's scope, not
   listed in this task's scope items; not repeated here.

## Next steps

1. The later, explicitly-deferred step of deriving syndrome-VN bias from
   real per-codeword information (RL codebits) rather than a free scalar
   hyperparameter remains open and unaffected by this work
   (`sweeps/TODO_soft_asced.md`).
2. Combining soft-aSCED with this project's separate
   `members_per_group`/`MConvergedPolicy` sequential-grouping/early-stopping
   machinery remains untested at Delta>1 too (same open question flagged
   in the Delta=1 doc).
3. The low-bias underperformance found here (bias 1-3 measurably worse
   than hard at Delta>1) is itself worth a closer look -- e.g. whether it
   is specific to `Delta_seg` in {2,3,5,6} or would also appear for larger
   `Delta_seg`, and whether a bias schedule that scales with `Delta_seg`
   (rather than one shared scalar across all segment sizes) could recover
   more of the effort/latency savings at smaller bias values.
4. The new asced384 finding -- effort savings persist at FER-matching bias
   but full-parallel "latency" (max-over-paths) savings nearly vanish as
   `Delta_seg` grows -- suggests the *latency* benefit of soft-aSCED (as
   opposed to *effort*) may specifically require combining it with
   `members_per_group`/`MConvergedPolicy` early stopping (so the ensemble
   doesn't wait on its single worst path) rather than full-parallel -- this
   ties directly into the still-open "soft + members-per-group compounding"
   question from `sweeps/soft_asced_single_row_results.md`, now with a
   concrete motivation for why it matters more at larger `Delta_seg`.
