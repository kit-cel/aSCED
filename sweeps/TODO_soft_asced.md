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
- [ ] Implement single-splitter (Delta=1) soft-aSCED prototype on aSCED-48
      — in progress, background agent launched 2026-10-09 on
      `channel-code-lib2` branch `soft_asced_wireless` (worktree
      `/home/pj9034/channel-code-lib2/.claude/worktrees/soft_asced_wireless`,
      off `claude_sequential` @ e5f69b5) + a matching aSCED agent worktree.
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
