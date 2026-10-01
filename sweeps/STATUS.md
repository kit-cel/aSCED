# Sequential/mconverged aSCED experiment — status

Working context: aSCED (Python orchestration, this repo) + channel-code-lib2
(C++ decoder, cloned at `/home/pj9034/channel-code-lib2`). Both on git branch
`claude_sequential` (paired branches, aSCED's `pyproject.toml` pins
channel-code-lib2 to that branch). Both pushed to their remotes.
This worktree: `/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478`.

## What this is

Extending aSCED (paper: "Affine Subcode Ensemble Decoding of Linear Block
Codes") with a `DecoderSelector` (groups ensemble paths into sequentially-
executed batches) + `MConvergedPolicy` stopping (stop once M paths converge),
to measure real decoding-effort/latency savings vs. full-parallel decoding,
for the RL/QC-block 5G LDPC ensembles aSCED-48 and aSCED-384 (code
C_5G(132,66), from `reproduce_fig3_RL_zc11_asced_5G_LDPC.py`, originally on
aSCED's `v1.0.0` branch, reused here).

## Key C++ concepts (channel-code-lib2)

- `Ensemble_Decoder::decode()`: executes `DecoderSelector::select()`'s groups
  in order; `break`s out of the loop once `StoppingPolicy::should_stop()` —
  **this is a real stop**, skipped groups' decoders never run.
- Within the *triggering* group, all members still run to their own natural
  BP stopping point; `MConvergedPolicy` only clips *reported* effort
  (`clipped_decoding_effort`) to an estimate — see `Ensemble_Decoder.cpp`'s
  header comment. Real per-path early truncation was NOT requested/built.
- Stats to read per run: `decoder_stats.json`'s `average_ensemble_effort`
  (total BP-iteration cost, the "complexity" axis) and
  `average_ensemble_latency` (sum of per-group max-effort, the "latency"
  axis, assuming full hardware parallelism within a group).
- `target_num_converged` MUST be a reachable absolute count, not a fraction
  of ensemble size — see Calibration below.

## Bugs found and fixed this session (channel-code-lib2, claude_sequential)

1. Missing virtual destructors on `Encoder`, `BP_Scheduling`, `BP_CN_Processing`
   (UB via unique_ptr<Base> deletion) — fixed, valgrind-verified.
2. `Simulation_Env.cpp` hardcoded `num_threads=1` (debug leftover, reverted
   real multithreading) — removed; `num_threads` also newly exposed to Python.
3. `set_syndrome_sequential(members_per_group)` (no-mapping overload) always
   threw — fixed (auto-fills decoder_set_size empty mappings).
4. `MConvergedConfig(0)` now throws instead of silent UB via assert (compiled
   out in Release).
5. `channel_coding_lib` (shared) was missing nlohmann_json's include dir,
   silently breaking every `-Dtesting=ON` test target — fixed.
6. Added `tests/test_ensemble_stopping.cpp` (MConvergedPolicy boundary logic,
   selector partition correctness, syndrome-fix regression). All pass.

## Calibration finding (important, don't repeat this mistake)

First sweep used `target_fraction_converged=0.5` → target_num_converged=24
(aSCED-48) / 192 (aSCED-384). **Never reachable**: measured
`average_number_converged_path` over 1.0-4.0dB stayed at 3.5-10 (aSCED-48)
and 4.8-8.7 (aSCED-384) regardless of ensemble size (these are overcomplete-
PCM subcode paths; individual-path exact convergence is intrinsically rare
and does NOT scale with ensemble size). Stopping never triggered; all three
selectors matched full_parallel's effort within noise.
Fixed: switched to absolute `target_num_converged`, calibrated from that
data to `{2, 6}`. Verified target=2 triggers real savings (asced48@2.5dB:
effort 1350→383, ~3.5x).

## Race condition found and fixed

Tasks sharing one `save_dir` across SNR points (all 7 SNR tasks of one
config launched in the same SLURM array concurrency wave) raced on the
C++ side's non-atomic JSON read-modify-write (`update_json`), silently
dropping entries — confirmed via `sacct` showing identical start timestamps
for all 7 tasks, 3 of them (fastest-finishing) lost. Fixed: `sim.save_dir`
now includes the SNR value, one directory per run, no shared file. Lost
`asced48/full_parallel` points (1.0/1.5/2.0dB) were backfilled directly.
**If building result-merging code, always merge across per-SNR subdirectories
rather than trusting a single shared FER.json to be complete for runs from
early sweeps** (pre-fix data lives at `<config>/FER.json` directly with no
`snr_*` subdir; post-fix data lives at `<config>/snr_<X>/FER.json`).

## CRITICAL BUG — n_simul=110 was wrong, all n110 data is invalid (found 2026-10-01/02)

Built a results page (`sweeps/plot_sequential_results.html`) comparing our
full_parallel aSCED-48/384 against the paper's Fig. 3 curves and they didn't
match at all (our FER systematically 1.4-6x worse, growing with SNR —
not noise, >300 FE at every point). Root-cause process (all verified, not
guessed):
1. Ran the **original unmodified** `reproduce_fig3_RL_zc11_asced_5G_LDPC.py`
   (from `v1.0.0`) standalone — same bad numbers. Not something broken by
   the sequential-script port.
2. Ran **standalone NMSA** (bypasses aSCED ensemble logic entirely, just BP
   on the base PCM) — also mismatched the paper's NMSA curve. So the bug
   isn't in aSCED/selector/stopping logic at all, it's upstream of all of
   that.
3. Checked out aSCED's actual `v1.0.0` working tree (`/home/pj9034/aSCED`,
   which has its own `.venv` built against channel-code-lib2 pinned at the
   **exact** commit `18faaa12886fb32e5f1cb2daeb6507080370d54a` that
   `v1.0.0`'s own `uv.lock` requests) and ran NMSA/asced48 there — **same
   mismatch**. Rules out any channel-code-lib2 version regression between
   `v1.0.0` and `sequential`/`claude_sequential`.
4. User found it: `n_simul` is **not an offset**, it's the final transmitted
   codeword length `n` directly (`number_vn_simul = n_simul + 2*Zc`, then
   puncturing removes exactly those `2*Zc` bits back off, so printed
   `n == n_simul`). Running with `n_simul=132` (not 110) reproduces the
   paper's NMSA curve closely (e.g. FER@2.0dB: ours 0.1801 vs paper 0.1798).
   **`n_simul=110` was actually simulating the higher-rate `C_5G(110,66)`
   code (rate 0.6) instead of the paper's `C_5G(132,66)` (rate 0.5)** — a
   genuinely different, weaker-protection code, fully explaining the
   systematically worse FER.

Verified the fix on `claude_sequential` too: `asced48 full_parallel` with
`n_simul=132` gives FER 0.371/0.187/0.071 @ 1.0/1.5/2.0dB vs paper's
0.357/0.185/0.069 — matches well.

**Impact: every number in both completed sweeps (jobs 531557, 532206, both
under `RESULTS/.../..._n110/...`) was generated with the wrong code and is
not comparable to the paper.** The *relative* comparisons within those
sweeps (full_parallel vs. sequential selectors, across group sizes/
orderings) are still internally self-consistent since everything in them
used the same (wrong) code — but don't reuse the absolute FER/effort/
latency numbers, and the published plot (`sweeps/plot_sequential_results.html`,
built from `..._n110` data) is now known-stale and needs to be rebuilt once
new data lands. Do not delete the `_n110` RESULTS directories without
asking — kept for now as a record, but treat them as invalid for anything
paper-comparable.

Fixed: `N_SIMUL = 132` in both `sweeps/generate_sequential_rl_manifest.py`
and `sweeps/generate_parallelism_ordering_manifest.py`, manifests
regenerated (same 70/56 row counts, grid unchanged — the SNR range 1.0-4.0dB
still brackets both FER=1e-1 and FER=1e-3 for the corrected code, confirmed
from the paper's own aSCED-48 numbers).

## Relaunched sweeps (2026-10-02, n_simul=132 — supersedes the n110 jobs below)

### Sweep A — main sweep (job 532495, `sweeps/sequential_rl_manifest.csv`, 70 tasks) — SUBMITTED, pending/running
variants {asced48, asced384} x selectors {full_parallel, fixed_sequential,
syndrome_sequential} x SNR {1.0..4.0dB step 0.5} x target_num_converged
{2,6} (full_parallel ignores target, 1 run/SNR only).
Output root: `RESULTS/fig_x_zc11_r4_seq_<selector>_mpg8_n132/<variant>_<selector>[_mpg8_target<N>]/`
`target_num_converged` values {2,6} reused as-is from the old (wrong-code)
calibration — not re-derived yet; if mconverged doesn't trigger on the
corrected code, recalibrate the same way as last time (read measured
`average_number_converged_path` from this run, pick new targets, re-run).

### Sweep B — parallelism + ordering (job 532496, `sweeps/parallelism_ordering_manifest.csv`, 56 tasks) — SUBMITTED, pending/running
asced48 only, `n_simul=132`. Same structure as before:
- Parallelism/group-size trade-off: `members_per_group in {1,2,4,8,16,24}`,
  selector fixed to `syndrome_sequential`, `target_num_converged=6` fixed,
  SNR 1.0-4.0dB. (mpg=48 endpoint = reuse asced48 full_parallel from Sweep A,
  not rerun.)
- Ordering comparison: `members_per_group=4` fixed, selector in
  {fixed_sequential, random_sequential} (syndrome_sequential@mpg=4 already
  in the parallelism sweep above), target=6, same SNR grid.

Once both complete: rebuild `sweeps/plot_data/build_data*.py` +
`sweeps/plot_data/build_html.py` pointing at the `..._n132` RESULTS
directories (same extraction logic — still always prefer per-SNR
subdirectories over flat top-level files, see race-condition note below,
which is independent of the n_simul bug and still applies), republish
`sweeps/plot_sequential_results.html`, and re-verify the "vs. literature"
overlay actually lines up this time. Re-derive the FER=1e-1/1e-3 SNR
operating points from the new full_parallel baseline rather than reusing
the old ≈2.15dB/≈3.7dB estimates (they were based on the wrong code, though
likely close since rate-0.5 vs rate-0.6 shifts the waterfall but not
drastically at these FERs — confirm, don't assume).

## Infrastructure

- `reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py`: CLI
  `<variant> <n_simul> <snr_start> <snr_end> <selector> <members_per_group> <target_num_converged>`.
  `n_simul=132` reproduces the exact C_5G(132,66)/Zc=11 code used in Fig. 3
  (n_simul IS the final transmitted length n directly, not an offset — see
  the n_simul bug section above; 110 was wrong).
  `target_errors=200`, `max_transmissions=2e6` (bounded on purpose — see
  commit messages for why the paper's own 1000/3e8 defaults are impractical
  at high SNR).
- `sweeps/generate_sequential_rl_manifest.py`, `sweeps/generate_parallelism_ordering_manifest.py`:
  regenerate the respective manifest CSVs.
- `sweeps/run_manifest_array.sbatch`: generic runner, `sbatch --array=1-N%10 sweeps/run_manifest_array.sbatch <manifest.csv>`.
  (`sweeps/sequential_rl_array.sbatch` is the original, sweep-A-specific
  version, kept as historical record of exactly what was submitted.)
- `smoke_test_sequential.py` (repo root): quick 4-thread functional smoke
  test, not part of the sweeps.

## Final verification pass (before compacting)

- Both repos clean, fully committed and pushed; aSCED's `uv.lock` pin
  (`c0fd4c1...`) exactly matches channel-code-lib2's `claude_sequential` HEAD.
- Found and fixed one more bug: an earlier `echo "sweeps/logs/" >> .gitignore`
  landed on a file with no trailing newline and silently merged into the
  existing `RESULTS/` line, producing the single broken pattern
  `RESULTS/sweeps/logs/` — which **unignored `RESULTS/`** (the actual
  sweep output, ~1.5MB and growing) from that commit onward. Fixed to two
  separate lines; verified both `RESULTS/` and `sweeps/logs/` are ignored
  again via `git check-ignore`. Confirmed via `git ls-files` that nothing
  under either path was ever accidentally committed in the meantime.
- Cross-checked every config's data completeness by unioning each config's
  flat `FER.json` keys with any `snr_*/FER.json` subdirectory keys (handles
  configs that have both pre- and post-race-fix data, e.g. the asced48
  backfill): confirmed exactly 0/7 missing for all Sweep B configs (56/56),
  and 0/7 missing for all Sweep A configs **except** the 3 still-running
  asced384@4.0dB tasks (6/7, missing only 4.0dB, which matches `sacct`
  showing those 3 as the only non-COMPLETED tasks) — i.e. no further silent
  data loss beyond what was already known and documented above.

## Not yet done / next steps

1. **Build the interactive plot.** Style guide verbatim in
   `sweeps/PLOTTING_STYLE.md` (Plotly, light/dark theme, hollow markers
   <300 errors, reference always visible, etc.) — load this before writing
   any plot code. Literature/baseline curves (NMSA, NMSA-352, AED-11,
   SCED-11/43, aSCED-11/48/384 from the paper, OSD-4 as ~ML reference)
   verbatim tikz in `sweeps/fig3_literature_data.tex` — parse this into the
   style guide's JSON `pts` format (note: most points lack per-point
   error/CPU counts; only a handful have inline comments like "200FE",
   "1133FE" — use those where present, leave others without hollow-marker
   data).
2. Plot deliverables wanted:
   - Sweep A: FER-vs-SNR comparing full_parallel/fixed_sequential/
     syndrome_sequential for both variants, overlaid with the literature
     curves (esp. aSCED-48/384 "(paper)" curves as a sanity cross-check,
     though exact reproduction was explicitly not required).
   - Sweep B part 1: latency-vs-complexity trade-off (average_ensemble_latency
     vs average_ensemble_effort) across members_per_group, evaluated at
     FER=1e-1 and FER=1e-3 (interpolate SNR per curve, then read off
     latency/complexity at that SNR) — this is the main deliverable the user
     asked for ("best latency vs best complexity trade-off").
   - Sweep B part 2: same latency/complexity comparison across orderings
     {fixed, random, syndrome} at members_per_group=4 fixed.
3. Once the plot's built and reviewed, decide whether to extend any of this
   to asced384 (deferred for cost reasons — ~8x slower per earlier timing).
4. Known TODOs not yet acted on (from earlier code review, still valid):
   - RL/QC-block-aligned (variable-size) grouping instead of uniform
     `members_per_group` — natural aSCED "batches" have variable size
     (2^Δ paths, Δ varies 2/3 for aSCED-48, 5/6 for aSCED-384); current
     selectors only support fixed group size.
   - Real per-path early truncation (vs. current effort-estimate-only
     stopping within the triggering group) — scope question raised earlier,
     not resolved either way; revisit if the latency numbers from Sweep B
     suggest it matters.
   - `enc_gf2`/`gf4_ops`/`gf4_vector` C++ tests still fail to build (stale
     `#include` paths from an old refactor) — unrelated pre-existing rot,
     not touched.
