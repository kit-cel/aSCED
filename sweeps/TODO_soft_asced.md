# TODO: soft-aSCED for Wireless

## Status

Not yet scoped / not started. The user will provide implementation details.
A companion effort is already active on the quantum side: repo
`/home/pj9034/quantum-error-correction`, branch/worktree
`claude/soft-asced-qec-decoder-460f88`
(`.claude/worktrees/soft-asced-qec-decoder-460f88`) — that session is the
"Soft-aSCED expert for the wireless port" the user can be @-mentioned
(session_id `local_2e4ec81e-4615-4304-8fdb-a09d450a4271`). Ask it, or the
user, for implementation details before starting.

## Core idea (as given by the user, 2026-10-09)

- Standard (hard) aSCED: appending a linearly independent row to the PCM
  hard-splits the solution space into two subcodes — a linear subcode
  (guessed syndrome bit = 0) and an affine subcode (guessed syndrome bit =
  1). Covering both guesses costs one decoding path per subcode per
  appended row (2x path blow-up per row).
- Soft-aSCED: treat the guessed syndrome bit as a binary RV with its own
  LLR instead of hard-fixing it.
  - Fixing that LLR to +/-infty recovers the existing hard-aSCED behavior
    (fully pinned to one subcode).
  - A finite LLR (e.g. a bias of +1) softly biases decoding toward one
    subcode without hard-committing — even when the true codeword isn't
    actually in that subcode, the soft push can still help rather than
    locking in a wrong guess.
  - Goal: cut the path-count blow-up (currently 2x per appended row) while
    hopefully preserving most of the error-correction gain, by representing
    syndrome-guess uncertainty continuously instead of discretely.

## New TODO items

- [x] Scope soft-aSCED for Wireless (channel-code-lib2 / aSCED) — scoped
      2026-10-09, see "Design spec" and "Scope decisions" below.
- [x] Implement single-splitter (Delta=1) soft-aSCED prototype on aSCED-48
      — done 2026-10-09. `channel-code-lib2` branch `soft_asced_wireless`
      (off `claude_sequential` @ e5f69b5, pushed @ 21d31b7); `aSCED` branch
      `soft_asced_wireless` (off `claude_sequential`). FER matches hard
      aSCED within noise; effort/latency drop 35-68%/47-81%. Independently
      verified (reran tests + sanity checks + a spot-check sweep). Full
      results: `sweeps/soft_asced_single_row_results.md`. Not merged into
      `claude_sequential` — stays a separate exploratory branch.
- [ ] Generalize to Delta>1 splitters.
- [ ] Once soft-aSCED is implemented: interpret the RL (raptor-like)
      codebits as syndrome VNs (variable nodes), per the soft idea above —
      i.e. derive the syndrome-VN bias from real per-codeword information
      instead of a free hyperparameter. Explicitly deferred past the first
      prototype.

## Scope decisions (user, 2026-10-09)

- Soft-aSCED is a new opt-in mode alongside existing hard aSCED selectors —
  not a replacement. Lets the two be benchmarked directly.
- First pass: syndrome-VN bias is a free scalar hyperparameter, swept
  experimentally (like the MSA-normalization sweep already done elsewhere
  in this project). The RL-codebit connection is explicitly a later step,
  not part of this first implementation.
- First validation target: a single appended splitter row (Delta=1), on
  aSCED-48 (not aSCED-384 or the Zc=6 Fig.5 codes).

## Delta=1 specifics / correction (expert session, 2026-10-09)

The expert initially said Delta=1 was untested on their side; they then
corrected this. Actual state: Delta=1 data exists (Reproduce/sweeps.py:
884-935, 1048-1090; Results/delta_sweep on the quantum repo), but always
with BOTH patterns (0 and 1) run per matrix (L = K/2 matrices, 2 soft paths
per matrix — one biased toward each pattern), never with an all-zero/
single-pattern setup. So "all-zero pattern + bias~10 at Delta=1" (the
original toy-case default) is specifically untested territory.

Results (LER / Tesseract-reference ratio, ~330 errors/point, CI ~+/-11%):
BB_144 (K=32, bias=10): Delta1=1.32/1.38 vs Delta2-all-patterns=1.28/1.26
vs Delta4-rand4=1.20/1.17. Toric_128 (K=128, bias=3): Delta1=1.85/1.26 vs
Delta2-all=1.56/1.19 vs Delta3-rand4=1.44/1.12. Delta=1 works but is the
weakest of Delta 1-4 at equal K (up to ~15% worse on BB_144, clearly worse
on Toric_128).

**Implication for our toy case**: at Delta=1, run BOTH pattern-0-biased
and pattern-1-biased soft paths per matrix (2 paths total — same path
count as hard aSCED's 2 hard paths). The win at Delta=1 isn't fewer paths,
it's that a wrong-guess soft path can self-correct ("guess correction")
instead of being hard-locked to the wrong subcode. All-zero/single-pattern
numbers at Delta=1 can still be reported as an exploratory extra, just not
as the primary hard-vs-soft comparison.

## Design spec (from the "Soft-aSCED expert for the wireless port" session, 2026-10-09)

Mechanism (ground truth, from the sibling quantum/QEC project that already
implements this): the guessed syndrome bit becomes one extra **binary
variable node** z_j per splitter row, with its own prior/LLR — it is NOT
injected into the syndrome or check targets.

- `H_soft = [H_l | M_l[:, m:]]`: Delta extra columns, one z_j per splitter.
  `H_l` is the overcomplete m_oc x n matrix (`H_l = M_l . H0_l mod 2`,
  `H0_l = [H; S_l]`); z_j connects to every overcomplete row k with
  `M_l[k, m+j] = 1` (usually many rows, not one).
- Check k reads: `H_l[k].x XOR sum_j M_l[k,m+j].z_j = M_l[k,:m].s`.
- Wireless case (vs. the quantum/coset case): for a codeword, `H.c = 0`, so
  every check target is 0 — only the z_j carry the subcode choice.
- The guess enters ONLY as z_j's prior: pattern bit 0 -> +bias, bit 1 ->
  -bias, empty pattern -> 0 (erasure); sign convention LLR>0 = bit 0.
- Output = physical part only (first n columns of the estimate).
- Validity is checked against the ORIGINAL (un-extended) H and the true
  syndrome, NOT against the guess -- this "guess correction" (a path can
  converge to a valid codeword that violates its own guess) is why soft is
  cheaper than hard: wrong-guess paths tend to stop early.
- Ensemble handling is unchanged: just run N independent BP decodes (one
  per matrix/pattern) and pick the first/best valid one.
- Binary simplification: no GF4/BP4 wrapper needed on our side (we have no
  GF4 nodes at all) — z_j is an ordinary binary VN, already handled by the
  existing generic `variable_node_update`.
- Hard/soft equivalence (bias -> +/-infinity recovers hard aSCED) is
  believed correct but was **not** verified bit-exactly even on the
  quantum side — worth actually testing this on our side.
- Recommended defaults (from their experience, not gospel): all-zero init
  pattern, fixed bias ~10 (3-10 roughly tied), prefer layered-CN BP
  schedule over flooding (ties into this project's own earlier
  flooding-vs-layered-CN TODO, see `sweeps/scheduling_study.py` /
  `sweeps/STATUS.md`).
- Known bugs on their side, avoid repeating: bias/mode not copied from
  config in `init()`; Serial-VN iteration counter not reset between
  trials; resumable `decode(k)` skipping an iteration per chunk boundary.

## Soft + members-per-group/MConverged compounding (expert session, 2026-10-10)

For when we test combining soft-aSCED with this project's existing
members-per-group/`MConvergedPolicy` ensemble-level stopping (currently
untested, see STATUS.md "open question"): on the QEC side, stopping ran on
top of already-soft configurations and found (CPU vs. the soft baseline,
at the benchmark LER):
- GB_126 soft K64 MW, groups of 8, M=4: 6-9x.
- GB_254 soft K16 MW, groups of 2, M=1: 4-6x.
- BB_288 soft K16 + ELC-B, groups of 2, M=1: ~3.5x.
- BB_144 soft K16 + ELC-B, M=8=K/2: only ~2x (ELC needs many converged paths).
- GB_46 soft K64 MW: no free saving (many paths converge in the same
  iteration; the minimum-weight one is often in a later group).

Caveats from the expert: this is "stopping still saved a lot on top of
soft", not a measured product of two independent factors (no hard-vs-soft
comparison WITH stopping at equal K was run on the same code); the two
mechanisms overlap in principle (soft already makes wrong-guess paths
cheap, so stopping's extra saving is smaller after soft than after hard).
**At K=2 (our current Delta=1 prototype) there is nothing to gain**: the
only grouping option is groups of 1/M=1, i.e. "run path 0, run path 1 only
if path 0 fails" -- saves at most the second path, costs FER whenever path
1 would've been the better codeword. Worth testing only once K>=8 (i.e.
once we've generalized past Delta=1 or run multiple matrices).

## Delta>1 generalization: scoping (2026-10-10)

Read the actual existing hard-aSCED Delta>1 machinery directly
(`create_asced_config` in `reproduce_fig3_RL_zc11_asced_5G_LDPC.py`)
rather than assuming the QEC side's general overcomplete-PCM construction
was needed. Key finding: **it isn't**. Hard aSCED on this project's
codebase already uses the simplest possible construction — `M_l =
identity` (`H_aux = [H; rows]`, the literal appended candidate rows
vstacked onto `H`, no row-combination at all) — exactly the same
simplification the Delta=1 prototype already used. Concretely:

- `candidate_rows` is organized in Zc-row (=11) "blocks" (the RL/QC-block
  pool from the greedy-search work). `num_used_blocks` picks how many
  blocks get appended; within each used block, `split_pattern` (e.g.
  `[2,4,6,8]`) cuts the 11 rows into `row_segments` via `np.split`.
- Each segment of `Delta_seg` rows becomes its OWN independent
  `H_aux = [H; segment_rows]`, run as a full `2**Delta_seg`-path hard
  sub-ensemble (`binary_vectors_in_suffix` enumerates every non-zero
  suffix pattern + the separately-added all-zero one).
- "asced48" = `splitting_pattern[1]=[2,4,6,8]` x 2 blocks = 2 x (4 segments
  of Delta=2 + 1 segment of Delta=3) = 2 x (4x4 + 8) = 2x24 = **48**.
  "asced384" = `splitting_pattern[2]=[5]` x 4 blocks = 4 x (1 segment of
  Delta=5 + 1 of Delta=6) = 4 x (32+64) = **384**. (Confirms these labels
  are total path COUNTS across multiple independent small-Delta segments,
  not one single flat Delta value.)

**Implication**: the C++ side already needs zero changes (confirmed
independently, not just taking the Delta=1 agent's word for it — 
`BP_Decoder::init` only checks `H0.cols() == H.cols() - num_syndrome_vns`,
nothing about connectivity structure, so arbitrary degree-1-per-row
syndrome VNs already work). The Delta>1 generalization is a **near-direct
port of `create_asced_config`**: write `create_soft_asced_config(...)`
with the identical block/segment loop, swapping the per-pattern path
construction from hard (`affine_offset` on `H_aux`) to soft
(`H_soft = [H_aux | Delta_seg-column identity block]`,
`num_syndrome_vns=Delta_seg`, `syndrome_vn_pattern=pattern bits`,
`syndrome_vn_bias=bias`, `H0=H`) — same segment sizes, same per-segment
`2**Delta_seg` path count as hard, for direct comparability. The general
overcomplete-row/`M_l` construction (raising syndrome-VN degree above 1,
as on the QEC side) is NOT required for this and is deferred as a later,
separate enhancement, not a precondition.

**Recommended plan** (to confirm with the user before implementing):
1. Validate first on **asced48** (Delta_seg in {2,3}, 2**Delta_seg in
   {4,8} -- small, tractable, matches the project's primary already-
   characterized baseline), then **asced384** (Delta_seg in {5,6},
   2**Delta_seg in {32,64} -- still tractable, no combinatorial blowup).
   Match hard's exact `2**Delta_seg` path count per segment (apples-to-
   apples with existing hard baselines) rather than a sampled subset —
   tractable at this scale; a sampled-subset option only becomes relevant
   for much larger segments (e.g. the unused "no split" Delta=10 variant,
   2**10=1024 — out of scope here).
2. Bias: keep as one shared swept free hyperparameter across all segments
   in this first pass (same scope as Delta=1), not yet per-segment-size-
   tuned, even though the expert noted best bias drifts with Delta on
   their side — a refinement to consider later if results suggest it
   matters, not a blocker now.
3. Re-run the Delta=1 sanity-check suite, generalized: bias-applied,
   hard-limit equivalence (now pattern = the full Delta_seg-bit true
   vector, bias saturated), reused-decoder determinism, and the per-path
   effort-mechanism check (now split by Hamming distance of the guess from
   the true pattern, not just correct/wrong, since a partially-wrong
   multi-bit guess is a new regime Delta=1 couldn't exercise).
4. RL-codebit-derived-bias and soft+mpg/MConverged compounding remain
   separate, later, deferred items (unaffected by this).

## Reference

- Quantum analog repo granted read access 2026-10-09:
  `/home/pj9034/quantum-error-correction` — similar codebase structure
  (decoder base class, config objects) to the wireless
  `channel-code-lib2`/`aSCED` pair. Core reference files (read at branch
  `claude/soft-asced-qec-decoder-460f88` HEAD, not its diff/history):
  `include/Decoder/BP/BP_Decoder.h:30-77`,
  `src/Decoder/BP/BP_Decoder.cpp` (init/init_syndrome/fill_priors/
  check_syndrome_match), `PythonScripts/run_soft_asced.py:190-261`.
- Related session: "Soft-aSCED expert for the wireless port" (session_id
  `local_2e4ec81e-4615-4304-8fdb-a09d450a4271`), working in the quantum
  repo's `claude/soft-asced-qec-decoder-460f88` branch/worktree.
