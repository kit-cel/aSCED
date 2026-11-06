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
