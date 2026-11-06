# Shortening: scoping (2026-11-06)

Read the actual code first (`Codes/generate_5G_LDPC.py`, every `reproduce_*.py`
/ `*_example.py` caller, and the paired `channel-code-lib2` C++ repo's git
history) rather than assuming shortening needs new machinery. Key finding:
**the existing construction-time shortening is already mathematically correct
and is the intended design — but it has literally never been exercised with a
non-empty shortening set anywhere in this repo.** The "properly implementing
shortening" request most likely needs a concrete target scenario, not a
redesign of the mechanism itself.

## What currently exists

`generate_5G_LDPC(BG_no, K_target, N_target)` (`Codes/generate_5G_LDPC.py:10`)
computes, among other things:

- `num_shortened_bits = K_b*Z - K_target` (line 86) — the gap between the
  lifted base graph's natural message capacity (`K_b*Z`, fixed by the chosen
  lifting size `Z`) and the requested `K_target`.
- `S = K_target + np.arange(num_shortened_bits)` (line 114) — the column
  indices of the message bits to shorten, i.e. the standard 5G NR convention:
  the *last* message-bit positions in the base-graph message block are the
  ones fixed to a known value (conventionally 0) and not transmitted.

`S` (returned as `s`) is then consumed by
`get_final_matrices_and_message_bit_pucturing(H, s, p)` (line 125):

- `keep_cols = np.setdiff1d(np.arange(H.shape[1]), s)` then `H = H[:, keep_cols]`
  (lines 160-161) — the shortened columns are **deleted from H entirely**,
  and `G = gf2(H).null_space()` is recomputed from the reduced `H`, giving a
  code with `k == K_target` (not `K_b*Z`).

This is mathematically sound: a message bit that is known to be 0 to both
encoder and decoder contributes nothing to any parity-check equation, so
deleting its column from `H` (and re-deriving `G`) is exactly equivalent to
keeping the column and assigning it infinite-reliability LLR / fixing it to 0
during decoding, just without the wasted variable node. This is a standard
technique in LDPC simulation codebases, and it correctly reproduces the 5G
NR semantic distinction between *shortening* (value known a priori, bit not
transmitted, **removed from the code entirely** here) and *puncturing*
(value unknown a priori, bit not transmitted, **still decoded** — these stay
as live columns and are instead handled at simulation time via
`sim.puncturing(message_bit_pucturing)`, which erases/marks those bits at the
channel but still asks the decoder to recover them).

So: nothing is "missing" in the sense of an unhandled bit-treatment case.
Shortened bits are not present in the simulated `H`/`G` at all, so there is no
decoder-side step to add (no infinite-LLR injection, no bit-comparison
exclusion needed — they're simply not part of `k`/`n` for the simulated
code).

## What's actually missing: shortening is dead code in practice

Grepped every caller of `generate_5G_LDPC(...)` repo-wide:

- `bp_example.py`, `AED_example.py`, `MBBP_and_stopping_example.py`,
  `reproduce_fig3_asced_5G_LDPC.py`, `sequential_asced_5G_LDPC.py`,
  `smoke_test_sequential.py` all call it with the **same literal**
  `(BG_no=2, K_target=66, N_target=132)`.
- For `BG_no=2`, `K_target=66`: `K_b=6` (falls in the `K_target <= 192`
  bracket), and the minimal valid lifting size satisfying `K_b*Z >= 66` is
  exactly `Z=11` (from the `[11, 22, 44, ...]` set) — so `K_b*Z = 66 ==
  K_target` exactly, giving `num_shortened_bits = 0` and `S = []` every time.
- `reproduce_fig4_asced_CCSDS_LDPC.py` and `toy_example.py` import
  `generate_5G_LDPC` but never call it (CCSDS AList file / hand-built
  Hamming code instead) — dead imports.
- The whole "RL zc11 / zc6" family (`reproduce_fig3_RL_zc11_asced_5G_LDPC*.py`,
  `reproduce_fig5_scatter_plot_zc6_5G_LDPC*.py`, and every `sweeps/*.py`
  script — `greedy_block_search*.py`, `alpha_sweep.py`, `scheduling_study.py`,
  `scout_zc6_fer_targets.py` — i.e. essentially all of the project's actual
  sweep/paper-reproduction work) **doesn't call `generate_5G_LDPC` at all**.
  It loads a precomputed max-rank lifted PCM
  (`Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz`) and builds `H`/`G`
  directly via a custom HARQ-style incremental-redundancy row/column slice
  (`H = gf2(H_full[:m, :number_vns])`), with no shortening concept in that
  path whatsoever — only `message_bit_pucturing = np.arange(2*Zc)` (the fixed
  first-2Z puncturing block).

So, across the entire repo as it stands today, `num_shortened_bits > 0` has
never actually been hit at runtime — the shortening branch in
`generate_5G_LDPC.py` (lines 112-117) and the column-deletion branch in
`get_final_matrices_and_message_bit_pucturing` (lines 143-161) are exercised
with `s = []` in every real invocation. The logic is plausible by inspection
but **untested with real data** in this codebase.

## A confusing but informative dead end: `channel-code-lib2`'s stale examples

`channel-code-lib2/examples/demo*.py` (`demo.py`, `demo_FerVsSnr.py`,
`demo_DecoderSweep.py`, `demo_decoder_comparison.py`,
`demo_encdec_awgn.py`) call `sim.shortening(s)` as if `Simulation_Env` had a
dedicated runtime shortening method, analogous to `sim.puncturing(p)`. It
doesn't — there is no `shortening` method in the current
`include/channel_coding_lib/Simulation_Env.h` or `src/bindings.cpp`, so these
calls would raise `AttributeError` if actually run today.

`git log` on `channel-code-lib2` explains why: a `Simulation_Env::shortening`
method existed historically, but commit `d367315` ("BIGGER CHANGES: removed
message bit shortening and parity bit shortening from c++-> should now be
handled in e.g. python; message bit puncturing still supported; ...",
authored by Jonathan Mandelbaum, 2026-06-12, an ancestor of the current
`claude_sequential` HEAD in both repos) deliberately removed it, with the
explicit intent that shortening be handled in Python instead. `src/Simulation_Env.cpp:51`
still has a leftover comment `// random_data[shortened bits] = 0` marking
where the old C++ handling used to live.

**This confirms the aSCED-side column-deletion approach is the intended
design, not a workaround** — but the `examples/demo*.py` scripts in
`channel-code-lib2` were never updated after that refactor and are stale/
broken as written. (Flagging this for awareness; fixing those examples is a
`channel-code-lib2` concern, not in scope here, and that repo's own
`claude_sequential`/`main` branches are off-limits per this task's
constraints regardless.)

One documentation mismatch worth noting: this repo's own `README.md` (line
39) lists "puncturing, shortening, bit comparison scope" together as
`Simulation_Env` "evaluation settings," implying parallel runtime treatment —
but only puncturing is actually a `Simulation_Env`-level call
(`sim.puncturing(...)`); shortening is entirely a pre-`Simulation_Env`,
code-construction-time step (`k`/`n` are already shrunk before
`Simulation_Env(k, n, ...)` is constructed). This may be part of what
prompted the "properly implement shortening" request — the README's phrasing
overstates symmetry with puncturing that doesn't exist in the implementation,
even though the implementation itself looks correct.

## Open questions (blocking for Stage 2)

1. **Which script(s)/code variant is this actually for?** The only path that
   has shortening wired up at all (`generate_5G_LDPC` +
   `get_final_matrices_and_message_bit_pucturing`) is the `k=66,n=132`
   examples family (`bp_example.py`, `AED_example.py`,
   `sequential_asced_5G_LDPC.py`, etc.), and that family never actually
   triggers shortening (`Z=11` divides `K_target` exactly). The "RL
   zc11/zc6" reproduce scripts and all `sweeps/*.py` work — i.e. essentially
   all of the project's real experiments — use a completely different,
   hand-built PCM-slicing construction with no shortening concept at all. Is
   the ask to (a) exercise/verify the existing `generate_5G_LDPC` path with a
   `K_target` that actually needs shortening, (b) add an equivalent
   shortening step to the RL zc11/zc6 HARQ-style construction, or (c)
   something else?
2. **What concrete `(BG_no, K_target, N_target)` (or equivalent) needs to
   work that doesn't today?** Since shortening has never run with non-empty
   `S`, is there a specific target code size driving this (e.g. a real 3GPP
   TBS/code-rate combination that doesn't land on a lifting-size boundary),
   or is this purely about validating/hardening the existing untested
   branch?
3. **Does this interact with the aSCED/soft-aSCED ensemble-decoding work**
   (`STATUS.md`, `TODO_soft_asced.md` on `soft_asced_wireless`), or is it a
   standalone code-construction concern (e.g. just making
   `generate_5G_LDPC`/`bp_example.py`-style scripts correct for arbitrary
   `K_target`)? Neither `sweeps/STATUS.md` nor `sweeps/TODO_soft_asced.md`
   mentions shortening at all — it doesn't appear to be a previously
   flagged/deferred gap in this project's history, so this seems to be a new
   ask rather than picking up a documented TODO.
4. **Is any `channel-code-lib2` (C++) change anticipated?** Based on the
   history above, the project's own prior decision was explicitly to keep
   shortening Python-only and out of the C++ `Simulation_Env`. Should Stage 2
   stay within that constraint (Python-only, construction-time column
   deletion, no C++ changes), or is there a reason to reconsider that
   decision (e.g. performance of decoding the full lifted PCM, needing
   shortened-bit info downstream for some other purpose)?
5. **Should the stale `channel-code-lib2/examples/demo*.py` `sim.shortening(s)`
   calls be fixed/removed as a side effect, or is that explicitly out of
   scope** (noting again that `channel-code-lib2`'s `claude_sequential`/
   `main`/`sequential` branches are off-limits without explicit instruction)?
6. **If a non-trivial `K_target` is picked for testing, is there a reference
   FER curve or existing result to validate against** (so "properly
   implemented" has a concrete correctness check beyond "doesn't crash" —
   e.g. comparing FER of the shortened code against an equivalent
   differently-constructed code of the same effective `(k,n)`)?

## Bottom line

By inspection, the mechanism in `Codes/generate_5G_LDPC.py` +
`get_final_matrices_and_message_bit_pucturing` already does the right thing
for shortening, and matches the project's own prior architectural decision
(C++ repo commit `d367315`) to keep it Python-side. The real gap is that it
is **completely untested** (zero real usages with non-empty `S`) and
**disconnected from the project's actual current experiment family** (RL
zc11/zc6), which has no shortening concept at all. Stage 2 should be scoped
around answering the open questions above, not around redesigning the
column-deletion mechanism.

## Owner's answers to the open questions (2026-11-06)

1. **Which script/variant** — owner found the question confusing as posed.
   Inferring the narrowest reasonable scope from answer 2: this is about the
   `generate_5G_LDPC.py` / `get_final_matrices_and_message_bit_pucturing`
   construction path itself (the mechanism that already has shortening
   wired up), made to work for **arbitrary** `(BG_no, K_target, N_target)`
   — not about retrofitting the RL zc11/zc6 HARQ-style scripts, which use an
   entirely separate precomputed-PCM construction with no shortening concept
   at all. That family is **out of scope for this Stage 2** unless the owner
   asks for it explicitly later (flagged as a possible follow-up at the end
   of this roadmap).
2. **Concrete target** — "I want to run for arbitrary 5G LDPC codes."
   Confirms: generalize across the parameter space, don't hardcode or
   validate only a single new target.
3. **Interaction with aSCED/soft-aSCED ensemble work** — "Once you specify
   shortening, all paths in the ensemble use the same knowledge, so no real
   interaction required but rather already handled via passing correct LLRs
   and adjusting the encoder." Confirmed: standalone code-construction
   concern, no ensemble-specific handling needed.
4. **C++ changes anticipated?** — owner asked back whether they might be
   mistaken in expecting none. Structurally, shortened columns are deleted
   from `H`/`G` before anything reaches `channel_code_lib2` (per this doc's
   earlier finding), so the C++ side never sees a shortened bit — same as
   just using a smaller code, and consistent with the owner's own prior
   decision in `channel-code-lib2` commit `d367315`. But since this path has
   never run with non-empty `S`, this must be **verified empirically, not
   assumed** — see Step 3 below.
5. **Stale `channel-code-lib2/examples/demo*.py` `sim.shortening(s)`
   calls** — "Yes," fix/remove them as part of this work.
6. **Reference FER curve to validate against** — owner will supply one
   later; don't block on it. Tracked as a TODO (Step 7 below).

## Stage 2 implementation roadmap (drafted 2026-11-06, not yet executed)

Ordered so that each step is cheap/fast to run and fails fast before
spending time on the next, expensive step. Nothing in this roadmap has been
run yet — this is the plan for review, per the coordinator's instruction, not
a report of results.

**Step 1 — Python-only construction sweep (no decoder, no C++).**
Call `generate_5G_LDPC(BG_no, K_target, N_target)` across a grid chosen to
force non-trivial, varying `num_shortened_bits`, covering both base graphs:
- `BG_no=2`: `K_target` values picked strictly *between* successive
  `K_b*Z` boundaries (e.g. ~50, 70, 100, 150, 300, 600 — crossing all four
  `K_b` brackets: 6, 8, 9, 10), so `num_shortened_bits` ranges from small to
  close to `Z-1` (its structural maximum, since a larger gap would have
  picked a smaller/different `Z`).
- `BG_no=1`: same idea at the larger `K_b=22` scale (`K_target` in the
  low-to-mid thousands, up to `K_max=8448`).
- For each `K_target`, vary `N_target` to also force `num_punctured_bits`
  to be zero in some combos and nonzero in others, so shortening and
  puncturing interact both separately and simultaneously.
- Edge cases to include explicitly: `num_shortened_bits=0` (today's only
  exercised case, must still pass — regression), and the largest
  `num_shortened_bits` achievable for a given `BG_no`/bracket.
**Validates:** the generator itself (`S`/`P` index computation, `N_b`
adjustment for excess bits) doesn't produce malformed or out-of-range
indices for arbitrary inputs. Pure numpy/galois — fast, no decoder.

**Step 2 — Validate `get_final_matrices_and_message_bit_pucturing` on every
combo from Step 1.**
For each: assert `k == K_target` exactly; assert resulting `n` matches
`N_target` after puncturing; assert `H @ G.T == 0` (mod 2) i.e. G is a valid
generator for the reduced H; check no unexpected rank deficiency (reduced H
should still have full row rank `n-k`, or the code should explicitly handle
redundant/all-zero rows if a check happens to touch only shortened
variables); and confirm the existing assertions in that function
(`parity_bit_puncturing[0] > s[-1]`, parity-puncturing-must-be-trailing-
columns) hold generically rather than only by coincidence at `K_target=66`.
**Validates:** the shortening+puncturing *interaction* logic for arbitrary
inputs — this is the step most likely to surface an edge case that happens
to be masked by the one combo exercised today.

**Step 3 — End-to-end decode through `channel_code_lib2`, the empirical
verification the owner asked for (answer 4).**
For a representative subset of Step 1/2's combos (not the full grid — pick
~4-6 that span small/large `num_shortened_bits`, BG1 and BG2, and
shortening-with-puncturing vs. shortening-alone), run a real decode:
- Construct `Simulation_Env(k, n, ...)` as in `bp_example.py`, with both
  `PCM_Encoder_config` and `G_Encoder_config`, plus the all-zero-codeword
  path, matching `bp_example.py`'s own cross-check pattern.
- Call `sim.puncturing(message_bit_pucturing)` and run `get_error_rates`
  at a couple of SNR points with a small frame count (e.g. 100-500 frames
  — enough to confirm the pipeline runs and produces a sane, non-degenerate
  FER, not a full-fidelity cluster-scale run).
- Check: no crash/exception in either encoder path; `FER_AZ`, `FER_PCM`,
  `FER_G` are all finite, non-degenerate, and consistent with each other
  (large divergence between the three would indicate a hidden construction
  bug, same cross-check `bp_example.py` already does for the non-shortened
  case); no implicit minimum-code-length assumption in `Simulation_Env`/BP
  config is violated for small shortened `k`.
**Validates:** whether a C++ change is actually needed — concretely
confirming or refuting the owner's (and this doc's) expectation of "no,"
rather than assuming it.

**Step 4 — Plausibility cross-check pending the real reference curve.**
Until the owner supplies a reference FER curve (answer 6, tracked as Step 7
TODO below), sanity-check each Step 3 result against a same-length/similar-
rate code built *without* going through shortening (e.g. a `K_target`
that happens to land exactly on a `K_b*Z` boundary with comparable `(k,n)`).
This is not a strict equality check (the codes differ), just a check that
FER magnitude/slope is in a plausible range for the rate/length — a
stand-in sanity gate, not the final correctness bar.

**Step 5 — Only if Step 3 surfaces a real break.**
If something in Step 3 can't be resolved Python-side, stop, write up exactly
what broke and why a C++ change would be needed, and take that back to the
owner/coordinator before touching `channel-code-lib2` at all — whose
`claude_sequential`/`main`/`sequential` branches remain off-limits without
explicit instruction regardless of what Step 3 finds. Per the owner's own
prior architectural decision (commit `d367315`) and this doc's analysis,
the expectation is that this step is a no-op.

**Step 6 — Fix the stale `channel-code-lib2/examples/demo*.py` files**
(`demo.py`, `demo_FerVsSnr.py`, `demo_DecoderSweep.py`,
`demo_decoder_comparison.py`, `demo_encdec_awgn.py`): replace each
`sim.shortening(s)` call with the correct current pattern — i.e. construct
the already-shortened `H`/`G` via `generate_5G_LDPC` +
`get_final_matrices_and_message_bit_pucturing` *before* constructing
`Simulation_Env`, matching how aSCED's own example scripts
(`bp_example.py` et al.) already do it, and drop the dead
`sim.shortening(...)` call entirely. Note: this lives in a different repo
(`channel-code-lib2`) and its `claude_sequential`/`main`/`sequential`
branches are off-limits per the hard constraint on this task — needs a
non-default branch there (coordinate with the coordinating session on
which one) before editing.

**Step 7 — TODO (non-blocking).**
Add an explicit TODO (here or wherever this project tracks open items) that
the owner will supply a reference FER curve later; once available, wire it
in as a concrete regression check (e.g. an assertion in a test/reproduce
script comparing simulated FER against the reference within tolerance),
rather than relying solely on the Step 4 plausibility check.

**Step 8 — Deferred, not part of this Stage 2.**
Once Steps 1-7 land cleanly, consider (separately, later, only if asked)
whether the `generate_5G_LDPC`-based construction should also be made
available to the RL zc11/zc6 family, which currently bypasses it entirely.
Not pursued now per the narrowed scope in answer 1/this doc's Step-1
inference.

**Explicitly not in this roadmap:** no implementation has started. This is
the plan for review before any code is touched, per the coordinator's
instruction.
