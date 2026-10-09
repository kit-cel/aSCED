# Soft-aSCED, single-splitter (Delta=1) prototype: results

Date: 2026-10-09. Branches: `channel-code-lib2` `soft_asced_wireless` (off
`claude_sequential` @ `e5f69b5`, pushed @ `21d31b7`); `aSCED`
`soft_asced_wireless` (off `claude_sequential`), `pyproject.toml` pinned to
the above via a `git`+`branch` source (not a local path).

Driver script: `prototype_soft_asced_single_row.py` (repo root). Run with
`uv run python prototype_soft_asced_single_row.py`.

**Independent verification (orchestrating session, 2026-10-09)**: the
implementing agent's worktree turned out to have branched from a stale
`main` ref instead of `claude_sequential` (an infrastructure mixup, not a
correctness issue) -- its C++ branch was correctly based, but its aSCED-side
commit was re-extracted onto a properly-based `soft_asced_wireless` branch
rather than merged as-is. Before doing so: reran the full C++ test suite
(`test_soft_asced_syndrome_vn` + all pre-existing tests) directly, all
green; independently reran all three sanity checks fresh (different RNG
draws) -- identical qualitative results (bias=0 collapses patterns, bias=30
makes them diverge, 100/100 hard-limit agreement at `num_samples=100`,
reused-decoder determinism confirmed); and independently reran a reduced
spot-check of the main sweep (`target_errors=50`, bias in {3, 20}) -- FER
and the effort/latency reduction both reproduced within Monte-Carlo noise of
the table below. The findings reported here are corroborated, not just
transcribed.

## What this is

A first, minimal (single appended row, Delta=1) prototype of **soft-aSCED**:
instead of hard-baking an aSCED subcode guess into an extended syndrome
(today's "hard aSCED", needing one full BP path per guess), the guess
becomes the prior LLR of one extra ordinary binary variable node (the
"syndrome VN"), appended as a new column of the PCM. A path can still
converge to a codeword that violates its own guess ("guess correction"),
which is why soft aSCED can be cheaper than hard aSCED. Full mechanism
description lives in `channel-code-lib2`'s
`include/channel_coding_lib/Decoder/BP/BP_Decoder.h` (soft-aSCED comment
block on `BP_config`) and `BP_Decoder::syndrome_check` in
`src/Decoder/BP/BP_Decoder.cpp`.

Code under test: the "aSCED-48" family's underlying code, `C_5G(132,66)`
(n=154 before puncturing / n_transmitted=132 after puncturing the first
2*Zc=22 bits, k=66), loaded with the exact same parameters
(`n_simul=132, Zc=11, bg_vn=8, bg_cn=4`) as
`reproduce_fig3_RL_zc11_asced_5G_LDPC.py` uses for its own `asced48`
variant.

## C++ changes (channel-code-lib2, branch `soft_asced_wireless`)

- `BP_config` / `BP_Decoder`: new fields `num_syndrome_vns`,
  `syndrome_vn_pattern`, `syndrome_vn_bias`, `H0` (+ `set_H0()`). When
  `num_syndrome_vns == 0` (the default), all existing behavior (plain BP,
  hard aSCED) is completely unchanged -- verified by a dedicated regression
  test (`tests/test_soft_asced_syndrome_vn.cpp`,
  `"num_syndrome_vns=0 leaves plain BP config/behavior unchanged"`) and by
  the full existing test suite still passing (`test_bp_degree1_check`,
  `test_ensemble_stopping`, `test_gfq`, `test_gf2_vector`, `test_gf2_matrix`
  -- `test_enc_gf2`/`test_gf4_*` remain pre-existing, unrelated build rot per
  `sweeps/HANDOFF_decode_from_llrs.md`, not touched).
- `BP_Decoder::init`: when `num_syndrome_vns > 0`, fills the trailing
  `num_syndrome_vns` entries of the AVN scratch buffer (`llr_augmented`)
  once, from `syndrome_vn_pattern`/`syndrome_vn_bias` (sign convention:
  pattern 0 -> `+bias`, pattern 1 -> `-bias`, anything else -> `0`,
  i.e. erasure). This reuses the existing AVN ("additional variable node")
  machinery that `use_avns` already provides for ssPCM-style decoders, so no
  new padding/truncation logic was needed -- `get_codeword()` already
  returns only the physical part (`transmitted_bits` columns), matching the
  "output = physical part only" requirement.
- `BP_Decoder::syndrome_check`: when `num_syndrome_vns > 0`, checks the
  physical part of the codeword estimate against `H0` (the ORIGINAL,
  un-extended PCM) and an all-zero true syndrome, instead of the full
  extended `H`/`affine_offset`. This is the "guess correction" property and
  is a pure addition gated on the new field -- the `num_syndrome_vns == 0`
  branch is untouched.
- Exposed via `src/bindings.cpp`: the four new `BP_config` fields/method,
  plus (new, additive) a generic `Decoder` binding (`decode`,
  `get_codeword`, `.converged`, `.decoding_effort`) and a `build_decoder(config)`
  factory function, so a Python driver can run one decoder (BP or Ensemble)
  on an explicit LLR vector -- needed for the sanity checks below, which
  require feeding the *same* noisy sample through two different decoder
  configs, something `Simulation_Env`'s own Monte-Carlo loop does not expose.
  This happens to also be exactly the "expose `decode(LLRs)->codeword`"
  task scoped (but not implemented) in `sweeps/HANDOFF_decode_from_llrs.md`.

### Bug found and fixed (all three BP scheduling classes)

`BP_Flooding::decode`, `BP_Row_Layered::decode`, `BP_Column_Layered::decode`
all have an early-return path: if the *initial* hard-decision (straight from
the channel LLRs, before any BP iteration) already satisfies
`syndrome_check()`, the method returns immediately -- but it did so
**without ever setting `bp_converged = true` or resetting `iteration_count`**,
leaving both at whatever stale value they held (`false`/a previous call's
iteration count on a reused decoder, zero-initialized-but-wrong on a fresh
one). This is a genuine, pre-existing bug, not introduced by this work, but
it was *latent* in practice: hitting the "already valid before any
iteration" path is rare when checking the full extended syndrome (hard
aSCED, plain BP). Soft-aSCED's guess-correction `syndrome_check` (checking
only the physical part against `H0`, ignoring the syndrome VN) hits this
path far more often -- exactly the scenario where soft aSCED is cheap,
since the physical bits can already be a valid codeword with zero BP
iterations needed, independent of whatever the syndrome VN guesses. Found
via the `test_soft_asced_syndrome_vn.cpp` "guess correction" test case,
which initially failed with `converged == false` despite an obviously
correct decoded codeword. Fixed in all three scheduling classes
(set `bp_converged = true; iteration_count = 0;` before the early return).
All existing tests still pass after the fix.

## Python driver: `prototype_soft_asced_single_row.py`

Builds `H` (88x154), picks the single splitter row `s = candidate_rows[0:1, :]`
(first candidate row of the first QC-block; verified linearly independent of
`H` via `np.linalg.matrix_rank`), forms `H_l = [H; s]` (89x154), and:

- **Hard aSCED** (baseline, unchanged mechanism): 2 `BP_config`s on `H_l`,
  `affine_offset = 0...0` and `0...01` respectively (exactly
  `reproduce_fig3_RL_zc11_asced_5G_LDPC.py`'s own `create_asced_config`
  pattern for a single, unsplit row).
- **Soft aSCED**: `H_soft = [H_l | col]` (89x155) where `col` is zero except
  a 1 in the splitter's own row -- one trailing syndrome-VN column. 2
  `BP_config`s on `H_soft` (`use_avns=True`, `num_syndrome_vns=1`, `H0=H`),
  one with `syndrome_vn_pattern=[0]`, one with `[1]`, both swept over
  `syndrome_vn_bias`.

Both use identical BP settings otherwise (`cn_update_type="msa"`,
`scheduling_type="flooding"`, `norm_factor=0.75`, `max_iterations=32`,
`early_stopping=True` -- the project's existing aSCED-48 defaults), so the
comparison isolates the hard-vs-soft splitter mechanism itself.

### Design decision (flagged, not fully specified in the task)

**M_l = identity for this Delta=1 prototype** (no row-overcompleteness):
`H_l = [H; s]` directly, so the syndrome VN connects to **exactly one**
check row (the splitter's own), i.e. **syndrome-VN degree = 1**. The
general construction (`H_l = M_l . H0_l`, `M_l` combining original rows with
the splitter to reach an overcomplete `m_oc > m+1` row count) would give the
syndrome VN higher degree, as on the sibling QEC project (e.g. degree 23 on
one of their codes) -- but building that general overcomplete-row machinery
was out of scope for "smallest possible toy case" (explicit scope item 3).
Per the coordinator's relayed caveat, this low degree means the `bias=0`
(erasure) data point should not be over-interpreted as representative of
the *general* soft-aSCED construction; it is reported here only as a
mechanism sanity check (see below), not as a sweep data point.

## Sanity checks (run before trusting the main sweep)

All three explicitly-requested checks pass:

1. **Bias is really applied.** At `bias=0`, patterns 0, 1 and erasure (-1)
   produce bit-identical decoded output (prior literally 0 all three ways:
   `+0 == -0 == 0`) -- confirmed. At `bias=30`, pattern 0 vs. pattern 1
   produce *different* decisions on the same (weak-evidence) input --
   confirmed. (Matches the quantum project's own R1 bug check: "every path
   silently ran at a hardcoded default bias regardless of config" did not
   reproduce here.)
2. **Hard-limit equivalence** (new evidence -- never bit-exactly verified
   even on the sibling QEC project): on 200 random-message samples with a
   *known* true splitter bit (`s . c`), running the soft path with
   `pattern = true bit`, `bias = 30` (saturating, not hard-fixing) against
   the existing hard-aSCED path using that same correct guess:
   **200/200 (100%) decision agreement** overall, and **170/170 (100%)**
   agreement restricted to samples where both paths converged. This is
   strong empirical support for the claimed bias-\>infinity equivalence to
   hard aSCED, on real wireless-code data, for the first time.
3. **Reused decoder gives identical results to fresh-per-sample decoding.**
   Decoding the same 30 samples once on one reused `BP_Decoder` object (with
   an unrelated sample decoded in between) and once with a fresh decoder
   per sample: byte-for-byte identical results -- confirmed, no stale-state
   leak.

## Main result: hard vs. soft aSCED, bias swept

2 SNR operating points (2.0 dB, 3.0 dB), `target_errors=200`,
`max_transmissions=1e5` per run (cheap/prototype budget, not a production
sweep). Both hard and soft run as a genuine 2-path ensemble
(`FullParallelSelector`, the default) -- per the coordinator's corrected
scoping: at Delta=1, soft aSCED's apples-to-apples comparison against hard
aSCED runs **both** patterns (one path biased toward 0, one toward 1), not
a single unbiased path, so the path *count* is identical (2 vs. 2) and the
only thing varying is the mechanism (hard-fixed guess vs. soft-biased guess
+ guess correction).

| | FER @ 2.0dB | FER @ 3.0dB | effort @ 2.0dB | effort @ 3.0dB | latency @ 2.0dB | latency @ 3.0dB |
|---|---|---|---|---|---|---|
| **hard aSCED (2 paths)** | 0.1564 | 0.01585 | 45.16 | 39.04 | 32.98 | 32.99 |
| soft, bias=1 | 0.1599 | 0.01666 | 25.68 | 12.48 | 13.44 | 6.41 |
| soft, bias=2 | 0.1594 | 0.01529 | 26.34 | 12.55 | 14.13 | 6.56 |
| soft, bias=3 | 0.1539 | 0.01652 | 26.53 | 12.84 | 14.52 | 6.81 |
| soft, bias=5 | 0.1566 | 0.01648 | 27.50 | 13.19 | 15.42 | 7.15 |
| soft, bias=10 | 0.1428 | 0.01567 | 28.51 | 13.80 | 16.66 | 7.79 |
| soft, bias=20 | 0.1498 | 0.01470 | 29.23 | 14.64 | 17.35 | 8.67 |
| soft, bias=30 | 0.1531 | 0.01588 | 29.41 | 14.86 | 17.40 | 8.86 |
| soft, bias=35 (~PHI_UPPER_LIMIT) | 0.1551 | 0.01569 | 29.61 | 14.93 | 17.53 | 8.91 |

("effort"/"latency" = `average_ensemble_effort`/`average_ensemble_latency`
from `decoder_stats.json`, i.e. total BP-iteration cost and max-per-group
iteration cost, matching this project's existing convention, see
`sweeps/STATUS.md`.)

**Headline finding**: across the whole bias sweep, soft aSCED's FER is
statistically indistinguishable from hard aSCED's (all values within
Monte-Carlo noise for `target_errors=200`, ~7% relative standard error;
no monotonic trend vs. bias), while **decoding effort drops 35-43% at
2.0dB and 62-68% at 3.0dB**, and **latency drops 47-59% at 2.0dB and
73-81% at 3.0dB**, for every bias value tested. Effort/latency increase
monotonically with bias (as expected -- a more saturating guess behaves
more like a hard-fixed one) but **never reach hard aSCED's cost, even at
bias=35**. This is because the saving isn't only "fewer iterations before
giving up on a wrong guess" -- it's that soft aSCED's guess-correction
`syndrome_check` (against `H0` only) is a strictly *easier* condition to
satisfy than hard aSCED's (against the full extended system including the
guessed row), so soft paths can early-stop as soon as the physical bits
alone form a valid codeword, regardless of what the syndrome VN itself
converged to. That is a structural advantage of the mechanism, not just a
property of "weak bias," and is itself a new, concrete finding from this
prototype (not stated in the original design spec).

### Exploratory (unvalidated elsewhere, reported separately per instruction): single-pattern, 1 path only

A single soft path (pattern=0 only, no second pattern-1 path) was also run
as a cheap extra data point, explicitly **not** the primary comparison
(the quantum project has no Delta=1 data for this all-zero-pattern-only
setup either):

| bias | FER @ 2.0dB | FER @ 3.0dB | effort @ 2.0dB | effort @ 3.0dB |
|---|---|---|---|---|
| 1 | 0.1767 | 0.0199 | 12.67 | 6.23 |
| 5 | 0.2195 | 0.0270 | 13.89 | 6.55 |
| 10 | 0.2429 | 0.0360 | 14.62 | 6.90 |
| 30 | 0.2466 | 0.0503 | 14.69 | 7.42 |

As expected, a single unbiased/weakly-biased soft path is **cheaper but
meaningfully worse in FER** than both the hard baseline and the 2-pattern
soft ensemble (FER degrades further as bias increases, since the single
path commits more strongly to a guess it has only a 50% chance of being
right about, with no second path to cover the other case). This matches
the sibling project's finding that Delta=1 with a single path is the
weakest configuration among Delta=1..4, and confirms the 2-pattern
comparison above (not this one) is the right one to draw conclusions from.

## Design decisions / assumptions flagged for the user

1. **M_l = identity (syndrome-VN degree 1)** for this Delta=1 prototype --
   see "Design decision" above. Generalizing to Delta>1 (see below) is
   exactly where the general `M_l`/row-overcompleteness construction would
   need to be implemented, raising syndrome-VN degree.
2. **"aSCED-48" interpretation**: the task description's "aSCED-48
   (C_5G(48,24)-derived...)" is read as *the underlying code of the
   existing `asced48` decoder variant* in
   `reproduce_fig3_RL_zc11_asced_5G_LDPC.py` (`C_5G(132,66)`, n_simul=132)
   -- not a literally different (48,24) code. That script's own `asced48`
   variant natively uses a Delta=4 split (`splitting_pattern[1] = [2,4,6,8]`,
   batch2); for this prototype we instead build our *own* single appended
   row (Delta=1) from the same base PCM / candidate rows, per the explicit
   "smallest possible toy case" scope, rather than reusing asced48's native
   split.
3. **2-path soft comparison, not 1-path**, per the coordinator's corrected
   scoping message (relayed mid-task): primary comparison runs both
   syndrome-VN patterns (0 and 1) as 2 soft paths, matching hard aSCED's 2
   paths exactly; the win is FER-preserving cheaper convergence per path,
   not fewer paths. The single-pattern variant is reported only as a
   clearly-labeled exploratory extra.
4. **Erasure/bias=0 not included in the main sweep**, given the degree-1
   syndrome VN caveat above -- it is exercised only inside sanity check 1
   as a mechanism check, not treated as a meaningful FER/effort data point.
5. **Bug fix to `bp_converged`/`iteration_count` on the "already valid"
   early-return path** in all three BP scheduling classes (see above) --
   judged in-scope since it directly affects whether soft-aSCED's FER/effort
   numbers are reported correctly, and is a strict correctness fix (verified
   not to change any existing passing test's result).

## Next steps for generalizing to Delta>1 splitters

1. Implement the general overcomplete-PCM construction (`H0_l = [H; S_l]`
   with `Delta` splitter rows, `M_l` expressing an `m_oc x (m+Delta)`
   overcomplete `H_l = M_l . H0_l`), matching
   `examples/Codes/overcomplete.py` in `channel-code-lib2` and/or the
   existing hard-aSCED row-combination logic in
   `reproduce_fig3_RL_zc11_asced_5G_LDPC.py`'s `create_asced_config`
   (`np.split`/`row_segments`) as a starting point for which row
   combinations to generate.
2. `H_soft = [H_l | M_l[:, m:]]`: this part of the C++ side needs **no
   change** -- `num_syndrome_vns` already generalizes to `Delta > 1` as
   written (it is already a `size_t` count, and `syndrome_vn_pattern`
   already a per-VN vector); only the *Python-side matrix construction*
   needs generalizing, not `BP_config`/`BP_Decoder`.
3. Ensemble size at `Delta > 1`: decide whether to run all `2^Delta`
   pattern combinations (matching hard aSCED's full blow-up, maximally
   apples-to-apples) or a sampled subset (matching the sibling project's
   "rand4" conventions for larger Delta) -- this is a path-count/compute
   trade-off the user should weigh in, not something this prototype's
   Delta=1 scope needed to resolve.
4. Re-run the three sanity checks at `Delta > 1` (especially hard-limit
   equivalence, which becomes `2^Delta` hard paths vs. `2^Delta` (or fewer)
   soft paths) and re-check the erasure/degree relationship now that
   syndrome VNs will have degree > 1 from genuine row-overcompleteness.
5. Per `sweeps/TODO_soft_asced.md`: the later, explicitly-deferred step of
   deriving syndrome-VN bias from real per-codeword information (RL
   codebits) rather than a free scalar hyperparameter remains open and
   unaffected by this work.
