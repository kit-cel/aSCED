# pylint: disable=invalid-name
"""Soft-aSCED, single-splitter (Delta=1) prototype on aSCED-48.

Compares the existing (hard) aSCED splitter mechanism -- the guessed
subcode bit baked directly into an extended syndrome / ``affine_offset``,
needing one full BP path per guess -- against the new *soft* mechanism,
where the guess instead becomes the prior LLR of one extra ordinary binary
variable node (the "syndrome VN"), appended as a new column of the PCM.
See channel-code-lib2's ``BP_Decoder.h`` (soft-aSCED comment block) for the
full mechanism and ``BP_Decoder::syndrome_check`` for the "guess
correction" property that makes soft aSCED's wrong-guess paths cheaper
than hard aSCED's (they can still converge to a valid codeword instead of
needing a second, oppositely-guessed path).

Code / setup: reuses the exact code-loading block of
``reproduce_fig3_RL_zc11_asced_5G_LDPC.py`` (n_simul=132, Zc=11, bg_vn=8,
bg_cn=4 -- the "aSCED-48"/"aSCED-384" code family, C_5G(132,66) after HARQ
extension, k=66 message bits, n=132 transmitted bits after puncturing the
first 2*Zc=22 bits). Design decision: that script's own "asced48" decoder
variant uses Delta=4 splitter rows (natively); for this Delta=1 prototype
we build our OWN single appended row from the same base PCM / candidate
rows instead of reusing asced48's native split, per the task's explicit
"smallest possible toy case" scope. See sweeps/soft_asced_single_row_results.md
for the full write-up of what this script finds.

Usage:
    uv run python prototype_soft_asced_single_row.py
"""

import json
import os
import shutil

import galois
import numpy as np

gf2 = galois.GF2

import channel_code_lib2 as ccl

# --------------------------------------------------------------------------
# Code setup (verbatim parameters from reproduce_fig3_RL_zc11_asced_5G_LDPC.py,
# n_simul=132 -- the "aSCED-48" code family)
# --------------------------------------------------------------------------

n_simul = 132
remove = 4
bg_vn = 12 - remove  # 8
bg_cn = 4
Zc = 11

code = np.load("Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz")
H_full = code["h"].astype(int)
H_full = np.delete(H_full, np.arange(6 * Zc, 10 * Zc), axis=1)

number_vn_simul = n_simul + 2 * Zc
number_vns_start = bg_vn * Zc + 2 * Zc
harq_part = (number_vn_simul - number_vns_start) // Zc
m = Zc * bg_cn + harq_part * Zc
number_vns = number_vns_start + harq_part * Zc
assert number_vn_simul == number_vns

H = np.array(gf2(H_full[:m, :number_vns])).astype(int)  # (88, 154)
G = np.array(gf2(H).null_space()).astype(int)
k, n = G.shape  # k=66, n=154
assert n == number_vns

candidate_rows = H_full[m:, :number_vns]  # future possible splitter rows
message_bit_pucturing = np.arange(2 * Zc, dtype=int)  # 22 punctured bits
n_transmitted = n - len(message_bit_pucturing)
assert n_transmitted == n_simul

print(f"Code: n={n} (n_transmitted={n_transmitted}), k={k}, m={m} "
      f"(aSCED-48 family, n_simul={n_simul})")

enc_cfg = ccl.PCM_Encoder_config(H, k, n)

# --------------------------------------------------------------------------
# Single splitter row (Delta=1): the first candidate row of the first
# QC-block. H_l = [H; s] (m+1 x n). H0 (for soft's guess-correction check)
# is the ORIGINAL, un-extended H (m x n) -- never H_l.
# --------------------------------------------------------------------------

s = candidate_rows[0:1, :]  # (1, n)
H_l = np.vstack([H, s]).astype(int)  # (m+1, n)
assert np.linalg.matrix_rank(gf2(H_l)) == m + 1, "splitter row must be linearly independent of H"

# H_soft = [H_l | M_l[:, m:]] with M_l = identity (Delta=1, no row-
# overcompleteness in this prototype -- see design-decision note in the
# results doc): one trailing column, 1 only in the splitter's own row.
col = np.zeros((m + 1, 1), dtype=int)
col[-1, 0] = 1
H_soft = np.hstack([H_l, col]).astype(int)  # (m+1, n+1)


def make_hard_bp_config(affine_offset):
    cfg = ccl.BP_config(H_l)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    cfg.affine_offset = list(affine_offset)
    return cfg


def make_soft_bp_config(pattern, bias):
    cfg = ccl.BP_config(H_soft)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    cfg.use_avns = True
    cfg.num_syndrome_vns = 1
    cfg.syndrome_vn_pattern = [pattern]
    cfg.syndrome_vn_bias = float(bias)
    cfg.set_H0(H)
    return cfg


def make_hard_ensemble_config():
    guess0 = make_hard_bp_config(np.zeros(m + 1, dtype=int))
    guess1_offset = np.zeros(m + 1, dtype=int)
    guess1_offset[-1] = 1
    guess1 = make_hard_bp_config(guess1_offset)
    return ccl.Ensemble_config(H, [guess0, guess1])


def make_soft_ensemble_config(bias):
    # Per the coordinator's corrected scoping message: at Delta=1, run BOTH
    # patterns (one soft path biased toward 0, one toward 1) for an
    # apples-to-apples 2-path comparison against hard aSCED's 2 paths --
    # NOT a single all-zero/unbiased path.
    path0 = make_soft_bp_config(0, bias)
    path1 = make_soft_bp_config(1, bias)
    return ccl.Ensemble_config(H, [path0, path1])


# --------------------------------------------------------------------------
# Sanity checks (per the coordinator's follow-up, run BEFORE trusting the
# main sweep -- these are the actual validation the soft-aSCED expert
# flagged as previously untested / bug-prone on the sibling QEC project).
# --------------------------------------------------------------------------

NEAR_ZERO = 1e-10  # matches channel_coding_lib's safe::NEAR_ZERO erasure-LLR magnitude


def awgn_llrs(codeword_bits, snr_db, rng):
    """BPSK + AWGN channel LLRs for the FULL (un-punctured) n-length
    codeword, with the punctured positions overwritten by a near-zero
    (erasure) LLR of random sign -- this matches
    Simulation_Env::simulate_transmission exactly: decode() is always
    called with n (not n_transmitted) LLRs, puncturing only erases some of
    them rather than shrinking the vector. Sign convention: LLR < 0 =>
    hard bit 1."""
    r = k / n_transmitted  # rate is computed from the TRANSMITTED length
    snr_lin = 10 ** (snr_db / 10.0)
    L_c = 4.0 * r * snr_lin  # standard BPSK/AWGN channel LLR scaling (Es/N0 -> L_c)
    sigma = 1.0 / np.sqrt(2.0 * r * snr_lin)
    bpsk = 1.0 - 2.0 * np.asarray(codeword_bits, dtype=float)
    noisy = bpsk + rng.normal(0.0, sigma, size=len(codeword_bits))
    llrs = noisy * L_c
    for i in message_bit_pucturing:
        llrs[i] = NEAR_ZERO if rng.integers(0, 2) == 0 else -NEAR_ZERO
    return llrs.tolist()


def encode_random_message(rng):
    msg = rng.integers(0, 2, size=k)
    codeword = np.array(gf2(msg) @ gf2(G)).astype(np.int64)  # (n,) full codeword incl. punctured bits
    return msg, codeword


def sanity_check_1_bias_applied():
    """bias=0 must give exactly the erasure-path result, and changing the
    bias must actually change outcomes (their R1 bug: bias silently
    ignored)."""
    print("\n--- Sanity check 1: bias is really applied ---")
    rng = np.random.default_rng(1)
    _, codeword = encode_random_message(rng)
    llrs = awgn_llrs(codeword, snr_db=0.5, rng=rng)  # low SNR: ambiguous evidence

    results = {}
    for pattern, bias in [(0, 0.0), (1, 0.0), (-1, 0.0), (0, 30.0), (1, 30.0)]:
        cfg = make_soft_bp_config(pattern, bias)
        dec = ccl.build_decoder(cfg)
        dec.decode(llrs)
        results[(pattern, bias)] = tuple(dec.get_codeword())

    bias0_agree = (
        results[(0, 0.0)] == results[(1, 0.0)] == results[(-1, 0.0)]
    )
    bias_changes_outcome = results[(0, 30.0)] != results[(1, 30.0)]
    print(f"  bias=0 collapses pattern 0/1/erasure to the same decision: {bias0_agree}")
    print(f"  bias=30 pattern 0 vs pattern 1 differ (bias has an effect): {bias_changes_outcome}")
    return bias0_agree, bias_changes_outcome


def sanity_check_2_hard_limit_equivalence(num_samples=3000, snr_db=2.0):
    """On samples with a KNOWN true splitter bit, pattern = true bit with
    bias saturated at +/-30 (via syndrome_vn_bias=30) should (nearly) agree
    with the existing hard-aSCED path using that same correct guess.

    NOTE (per the soft-aSCED expert's methodological correction): agreement
    on N trials only bounds the DISAGREEMENT rate, it does not establish
    exact equivalence -- 0 disagreements in N trials bounds the true
    disagreement rate at < 3/N (95%, "rule of three"). num_samples=3000
    (vs. the original 200) tightens that bound to <0.1%. Non-converged-hard
    trials are intentionally included in the unconditional count (not just
    both-converged), since that's exactly where saturation (+/-30) and a
    hard-fixed bit could plausibly differ."""
    print("\n--- Sanity check 2: hard-limit equivalence (bias=30 saturation vs. hard aSCED) ---")
    rng = np.random.default_rng(2)
    agreements = 0
    both_converged_agreements = 0
    both_converged_count = 0
    hard_not_converged_count = 0
    hard_not_converged_agreements = 0
    for _ in range(num_samples):
        _, codeword = encode_random_message(rng)
        true_bit = int((s[0] @ codeword) % 2)

        llrs = awgn_llrs(codeword, snr_db=snr_db, rng=rng)

        hard_offset = np.zeros(m + 1, dtype=int)
        hard_offset[-1] = true_bit
        hard_cfg = make_hard_bp_config(hard_offset)
        hard_dec = ccl.build_decoder(hard_cfg)
        hard_dec.decode(llrs)
        hard_cw = hard_dec.get_codeword()

        soft_cfg = make_soft_bp_config(true_bit, 30.0)
        soft_dec = ccl.build_decoder(soft_cfg)
        soft_dec.decode(llrs)
        soft_cw = soft_dec.get_codeword()

        agree = hard_cw == soft_cw
        if agree:
            agreements += 1
        if hard_dec.converged and soft_dec.converged:
            both_converged_count += 1
            if agree:
                both_converged_agreements += 1
        if not hard_dec.converged:
            # Per the soft-aSCED expert: this is the ONLY subset where
            # saturation (bias=30) and a hard-fixed bit could plausibly
            # differ, since a converged hard path already pins the decision.
            hard_not_converged_count += 1
            if agree:
                hard_not_converged_agreements += 1

    disagreements = num_samples - agreements
    agreement_rate = agreements / num_samples
    disagreement_ci = wilson_interval(disagreements, num_samples)
    conditional_rate = (
        both_converged_agreements / both_converged_count if both_converged_count else float("nan")
    )
    hard_not_converged_agreement_rate = (
        hard_not_converged_agreements / hard_not_converged_count
        if hard_not_converged_count else float("nan")
    )
    hard_not_converged_disagree_ci = wilson_interval(
        hard_not_converged_count - hard_not_converged_agreements, hard_not_converged_count
    )
    print(f"  overall decision agreement (correct guess, bias=30 vs. hard): "
          f"{agreement_rate:.4f} ({agreements}/{num_samples})")
    print(f"  => disagreement rate 95% CI: [{disagreement_ci[0]:.5f}, {disagreement_ci[1]:.5f}] "
          f"({disagreements} observed disagreements)")
    print(f"  agreement given BOTH converged: {conditional_rate:.4f} "
          f"({both_converged_agreements}/{both_converged_count})")
    print(f"  agreement given HARD DID NOT CONVERGE (the subset where saturation vs. "
          f"hard-fixed CAN differ): {hard_not_converged_agreement_rate:.4f} "
          f"({hard_not_converged_agreements}/{hard_not_converged_count}); "
          f"disagreement 95% CI [{hard_not_converged_disagree_ci[0]:.5f}, "
          f"{hard_not_converged_disagree_ci[1]:.5f}]")
    return (agreement_rate, conditional_rate, disagreement_ci, num_samples,
            hard_not_converged_agreement_rate, hard_not_converged_count,
            hard_not_converged_disagree_ci)


def sanity_check_3_reused_decoder(num_samples=50, snr_db=2.0):
    """Running the same set of samples once as a batch on ONE reused
    decoder object, and once with a fresh decoder per sample, must give
    identical results (catches stale-per-trial-state bugs)."""
    print("\n--- Sanity check 3: reused decoder vs. fresh-per-sample decoder ---")
    rng = np.random.default_rng(3)
    samples = []
    for _ in range(num_samples):
        _, codeword = encode_random_message(rng)
        llrs = awgn_llrs(codeword, snr_db=snr_db, rng=rng)
        samples.append(llrs)

    cfg = make_soft_bp_config(0, 10.0)

    reused_decoder = ccl.build_decoder(cfg)
    reused_results = []
    for llrs in samples:
        reused_decoder.decode(llrs)
        reused_results.append((tuple(reused_decoder.get_codeword()), reused_decoder.converged))

    fresh_results = []
    for llrs in samples:
        fresh_cfg = make_soft_bp_config(0, 10.0)
        fresh_decoder = ccl.build_decoder(fresh_cfg)
        fresh_decoder.decode(llrs)
        fresh_results.append((tuple(fresh_decoder.get_codeword()), fresh_decoder.converged))

    identical = reused_results == fresh_results
    print(f"  reused-decoder results identical to fresh-decoder-per-sample: {identical}")
    if not identical:
        for i, (r, f) in enumerate(zip(reused_results, fresh_results)):
            if r != f:
                print(f"    mismatch at sample {i}: reused={r} fresh={f}")
    return identical


def sanity_check_4_effort_mechanism(num_samples=500, snr_db=2.0, bias=10.0):
    """Direct, per-path evidence for WHY soft aSCED's effort/latency drops
    (per the coordinator's follow-up question -- confirm this is really the
    'wrong' path converging earlier, not some other confound). Splits each
    path's own BP iteration count (`decoding_effort`) and convergence flag
    by whether that path's guess matches the true splitter bit or not,
    separately for hard and soft (same bias/SNR as used elsewhere)."""
    print("\n--- Sanity check 4: per-path effort/convergence, correct vs. wrong guess ---")
    rng = np.random.default_rng(4)
    buckets = {
        "hard_correct": ([], []), "hard_wrong": ([], []),
        "soft_correct": ([], []), "soft_wrong": ([], []),
    }
    for _ in range(num_samples):
        _, codeword = encode_random_message(rng)
        true_bit = int((s[0] @ codeword) % 2)
        wrong_bit = 1 - true_bit
        llrs = awgn_llrs(codeword, snr_db=snr_db, rng=rng)

        for guess, key in [(true_bit, "hard_correct"), (wrong_bit, "hard_wrong")]:
            offset = np.zeros(m + 1, dtype=int)
            offset[-1] = guess
            dec = ccl.build_decoder(make_hard_bp_config(offset))
            dec.decode(llrs)
            buckets[key][0].append(dec.decoding_effort)
            buckets[key][1].append(bool(dec.converged))

        for guess, key in [(true_bit, "soft_correct"), (wrong_bit, "soft_wrong")]:
            dec = ccl.build_decoder(make_soft_bp_config(guess, bias))
            dec.decode(llrs)
            buckets[key][0].append(dec.decoding_effort)
            buckets[key][1].append(bool(dec.converged))

    summary = {}
    for key, (effort, conv) in buckets.items():
        mean_effort = sum(effort) / len(effort)
        converged_frac = sum(conv) / len(conv)
        summary[key] = {"mean_effort": mean_effort, "converged_frac": converged_frac}
        print(f"  {key:14s} mean_effort={mean_effort:6.2f}  converged_frac={converged_frac:.3f}")
    return summary


# --------------------------------------------------------------------------
# Equal-K comparison: same path BUDGET (K=2), spent either as 2 patterns on
# 1 matrix (our headline comparison) or 1 pattern (0 only) on 2 DIFFERENT
# matrices (2 independently-built splitter rows). Per the soft-aSCED
# expert's correction: this -- not our original 1-path-vs-2-path comparison
# -- is the actual like-for-like test of "spend K on patterns vs. matrices",
# and on their side it was code-dependent (helped some codes, hurt others).
# --------------------------------------------------------------------------

s2 = candidate_rows[1:2, :]  # second splitter row candidate, for the equal-K check
H_l_2 = np.vstack([H, s2]).astype(int)
assert np.linalg.matrix_rank(gf2(H_l_2)) == m + 1, "second splitter row must be linearly independent of H"
col_2 = np.zeros((m + 1, 1), dtype=int)
col_2[-1, 0] = 1
H_soft_2 = np.hstack([H_l_2, col_2]).astype(int)


def make_soft_bp_config_matrix2(pattern, bias):
    cfg = ccl.BP_config(H_soft_2)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    cfg.use_avns = True
    cfg.num_syndrome_vns = 1
    cfg.syndrome_vn_pattern = [pattern]
    cfg.syndrome_vn_bias = float(bias)
    cfg.set_H0(H)
    return cfg


def make_equal_k_two_matrices_config(bias):
    """K=2: pattern 0 only, on 2 DIFFERENT splitter-row matrices."""
    path_m1 = make_soft_bp_config(0, bias)
    path_m2 = make_soft_bp_config_matrix2(0, bias)
    return ccl.Ensemble_config(H, [path_m1, path_m2])


# --- Generalized version, for repeating the equal-K check over several
# independent splitter-row draws (per the expert's point that a single
# matrix pair's result can be decided by matrix-to-matrix variance rather
# than the patterns-vs-matrices trade itself). ---

def make_soft_bp_config_for_row(row, pattern, bias):
    H_l_row = np.vstack([H, row]).astype(int)
    col_row = np.zeros((H_l_row.shape[0], 1), dtype=int)
    col_row[-1, 0] = 1
    H_soft_row = np.hstack([H_l_row, col_row]).astype(int)
    cfg = ccl.BP_config(H_soft_row)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    cfg.use_avns = True
    cfg.num_syndrome_vns = 1
    cfg.syndrome_vn_pattern = [pattern]
    cfg.syndrome_vn_bias = float(bias)
    cfg.set_H0(H)
    return cfg


def make_two_patterns_config_for_row(row, bias):
    return ccl.Ensemble_config(H, [
        make_soft_bp_config_for_row(row, 0, bias),
        make_soft_bp_config_for_row(row, 1, bias),
    ])


def make_pattern0_two_matrices_config(row_a, row_b, bias):
    return ccl.Ensemble_config(H, [
        make_soft_bp_config_for_row(row_a, 0, bias),
        make_soft_bp_config_for_row(row_b, 0, bias),
    ])


def equal_k_multi_pair_check(num_pairs=6, bias=10.0, snr_points=(2.0, 3.0),
                              target_errors=100, results_root="RESULTS/soft_asced_prototype_n132"):
    """Repeat the K=2 'patterns vs. matrices' comparison over several
    independent splitter-row draws (disjoint across pairs, drawn from the
    374-row candidate pool). Addresses the expert's point #3: a single pair
    can be decided by matrix-to-matrix variance, not the trade itself."""
    print("\n" + "=" * 70)
    print(f"Equal-K multi-pair check: {num_pairs} independent splitter-row draws, bias={bias}")
    print("=" * 70)
    rng = np.random.default_rng(5)
    idx_pool = rng.permutation(candidate_rows.shape[0])
    results = []
    for pair_i in range(num_pairs):
        r1_idx, r2_idx = int(idx_pool[2 * pair_i]), int(idx_pool[2 * pair_i + 1])
        row1 = candidate_rows[r1_idx:r1_idx + 1, :]
        row2 = candidate_rows[r2_idx:r2_idx + 1, :]

        tp_cfg = make_two_patterns_config_for_row(row1, bias)
        tp_fer, _, _, tp_counts = run_sim(
            tp_cfg, snr_points, os.path.join(results_root, f"mp_{pair_i}_2pat"),
            target_errors=target_errors,
        )
        tm_cfg = make_pattern0_two_matrices_config(row1, row2, bias)
        tm_fer, _, _, tm_counts = run_sim(
            tm_cfg, snr_points, os.path.join(results_root, f"mp_{pair_i}_2mat"),
            target_errors=target_errors,
        )

        row_result = {"row_indices": (r1_idx, r2_idx), "two_patterns": {"fer": tp_fer, "counts": tp_counts},
                      "pattern0_two_matrices": {"fer": tm_fer, "counts": tm_counts}}
        results.append(row_result)
        for snr in snr_points:
            print(f"  pair {pair_i} (rows {r1_idx},{r2_idx}) @ {snr}dB: "
                  f"2pat/1mat={_fer_ci_str(tp_fer, tp_counts, snr)}  |  "
                  f"pat0/2mat={_fer_ci_str(tm_fer, tm_counts, snr)}")

    # Per the soft-aSCED expert's correction: the SNR points within one pair
    # reuse the same matrix pair, so they are NOT independent samples. Count
    # independent PAIRS (does 2-patterns/1-matrix win at every SNR in this
    # pair?), not pair x SNR combinations, and report a sign test over the
    # pairs rather than treating per-pair FER margins as individually
    # significant (at this target_errors budget each point's own CI is wide).
    pair_winners = []  # +1 = 2pat/1mat wins this pair (all SNRs agree), -1 = 2mat wins, 0 = mixed/tie
    for r in results:
        signs = []
        for snr in snr_points:
            tp = r["two_patterns"]["fer"].get(snr) or r["two_patterns"]["fer"].get(str(snr))
            tm = r["pattern0_two_matrices"]["fer"].get(snr) or r["pattern0_two_matrices"]["fer"].get(str(snr))
            if tp is None or tm is None:
                continue
            signs.append(1 if tp < tm else (-1 if tm < tp else 0))
        if signs and all(sgn == signs[0] for sgn in signs):
            pair_winners.append(signs[0])
        else:
            pair_winners.append(0)  # mixed across SNRs within this pair -> ambiguous

    n_decided = sum(1 for w in pair_winners if w != 0)
    n_2pat = sum(1 for w in pair_winners if w == 1)
    n_2mat = sum(1 for w in pair_winners if w == -1)
    # Exact two-sided binomial sign-test p-value under H0: each decided pair
    # independently favors either side with probability 0.5.
    import math
    k = max(n_2pat, n_2mat)
    if n_decided > 0:
        tail = sum(math.comb(n_decided, i) for i in range(k, n_decided + 1)) / (2 ** n_decided)
        sign_test_p = min(1.0, 2 * tail)
    else:
        sign_test_p = float("nan")
    print(f"\n  Across {num_pairs} INDEPENDENT pairs (not pair x SNR): "
          f"2-patterns/1-matrix favored in {n_2pat}/{n_decided} decided pairs, "
          f"pattern0/2-matrices in {n_2mat}/{n_decided} "
          f"({num_pairs - n_decided} pair(s) mixed across SNRs)")
    print(f"  Two-sided sign-test p-value: {sign_test_p:.4f}")
    return {"pair_results": results, "pair_winners": pair_winners, "sign_test_p": sign_test_p}


# --------------------------------------------------------------------------
# Main bias sweep / hard-vs-soft FER comparison
# --------------------------------------------------------------------------

def _lookup_by_snr(d, snr):
    # decoder_stats.json / stats.json keys are produced by the C++ side's own
    # double_to_string (e.g. "2" for 2.0, "2.5" for 2.5) -- match by parsed
    # float value instead of assuming a particular string format.
    for k_str, v in d.items():
        if abs(float(k_str) - snr) < 1e-9:
            return v
    return None


def wilson_interval(num_events, num_trials, z=1.96):
    """Wilson score interval for a binomial proportion (95% by default).
    Falls back to the "rule of three" upper bound when num_events == 0,
    since Wilson degenerates to (0, ~0) there and understates the true
    uncertainty for small/zero counts."""
    if num_trials == 0:
        return (float("nan"), float("nan"))
    if num_events == 0:
        return (0.0, 3.0 / num_trials)
    phat = num_events / num_trials
    denom = 1 + z**2 / num_trials
    center = (phat + z**2 / (2 * num_trials)) / denom
    half = z * np.sqrt(phat * (1 - phat) / num_trials + z**2 / (4 * num_trials**2)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def run_sim(ensemble_cfg, snr_points, save_dir, target_errors=200, max_transmissions=int(1e5)):
    sim = ccl.Simulation_Env(k, n, "all")
    sim.target_errors = target_errors
    sim.max_transmissions = max_transmissions
    try:
        sim.num_threads = os.cpu_count() or 1
    except AttributeError:
        pass
    sim.auto_save = True
    sim.save_dir = save_dir
    sim.puncturing(message_bit_pucturing)
    sim.init(enc_cfg, ensemble_cfg, False)
    sim.get_error_rates(np.array(snr_points, dtype=float))

    fer = dict(sim.error_rates["FER-SNR"])
    effort = {}
    latency = {}
    counts = {}  # snr -> (frame_errors, trials), for CIs

    decoder_stats_path = os.path.join(save_dir, "decoder_stats.json")
    if os.path.exists(decoder_stats_path):
        with open(decoder_stats_path) as f:
            stats = json.load(f)
        for snr in snr_points:
            effort[snr] = _lookup_by_snr(stats.get("average_ensemble_effort", {}), snr)
            latency[snr] = _lookup_by_snr(stats.get("average_ensemble_latency", {}), snr)

    stats_path = os.path.join(save_dir, "stats.json")
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            stats = json.load(f)
        for snr in snr_points:
            fe = _lookup_by_snr(stats.get("frame_errors", {}), snr)
            tr = _lookup_by_snr(stats.get("trials", {}), snr)
            counts[snr] = (fe, tr)

    return fer, effort, latency, counts


def _fer_ci_str(fer_by_snr, counts_by_snr, snr):
    fe_tr = counts_by_snr.get(snr)
    if not fe_tr or fe_tr[0] is None:
        return f"{fer_by_snr.get(snr)}"
    fe, tr = fe_tr
    lo, hi = wilson_interval(fe, tr)
    return f"{fer_by_snr.get(snr):.5f} [{lo:.5f}, {hi:.5f}] ({fe} err / {tr} trials)"


def main():
    print("=" * 70)
    print("Soft-aSCED single-row (Delta=1) prototype -- sanity checks")
    print("=" * 70)
    b0_ok, bchange_ok = sanity_check_1_bias_applied()
    (agree_rate, cond_agree_rate, disagree_ci, hard_limit_n,
     hnc_agree_rate, hnc_count, hnc_disagree_ci) = sanity_check_2_hard_limit_equivalence()
    reuse_ok = sanity_check_3_reused_decoder()
    effort_mechanism = sanity_check_4_effort_mechanism()

    print("\n" + "=" * 70)
    print("Main sweep: hard aSCED (2 paths) vs. soft aSCED (2 paths, bias swept)")
    print("=" * 70)

    snr_points = [2.0, 3.0]
    results_root = "RESULTS/soft_asced_prototype_n132"
    shutil.rmtree(results_root, ignore_errors=True)

    print("\nRunning hard aSCED baseline...")
    hard_cfg = make_hard_ensemble_config()
    hard_fer, hard_effort, hard_latency, hard_counts = run_sim(
        hard_cfg, snr_points, os.path.join(results_root, "hard")
    )
    for snr in snr_points:
        print(f"  hard FER @ {snr}dB: {_fer_ci_str(hard_fer, hard_counts, snr)}")
    print(f"  hard avg ensemble effort: {hard_effort}")
    print(f"  hard avg ensemble latency: {hard_latency}")

    bias_values = [1, 2, 3, 5, 10, 20, 30, 35]
    soft_results = {}
    for bias in bias_values:
        print(f"\nRunning soft aSCED, bias={bias}...")
        soft_cfg = make_soft_ensemble_config(bias)
        fer, effort, latency, counts = run_sim(
            soft_cfg, snr_points, os.path.join(results_root, f"soft_bias{bias}")
        )
        soft_results[bias] = {"fer": fer, "effort": effort, "latency": latency, "counts": counts}
        for snr in snr_points:
            print(f"  soft(bias={bias}) FER @ {snr}dB: {_fer_ci_str(fer, counts, snr)}")
        print(f"  soft(bias={bias}) avg ensemble effort: {effort}")
        print(f"  soft(bias={bias}) avg ensemble latency: {latency}")

    print("\n" + "=" * 70)
    print("Equal-K check (K=2): 2 patterns on 1 matrix vs. pattern-0-only on 2 matrices")
    print("(per the soft-aSCED expert's correction -- this, not 1-path-vs-2-path, is")
    print(" the actual like-for-like 'spend K on patterns vs. matrices' comparison)")
    print("=" * 70)
    equal_k_bias = 10  # the expert's recommended starting point
    print(f"\nRunning 2-patterns-1-matrix (our headline config) at bias={equal_k_bias}...")
    two_patterns_cfg = make_soft_ensemble_config(equal_k_bias)
    tp_fer, tp_effort, tp_latency, tp_counts = run_sim(
        two_patterns_cfg, snr_points, os.path.join(results_root, "equalk_2patterns_1matrix")
    )
    print(f"\nRunning pattern-0-only-2-matrices at bias={equal_k_bias}...")
    two_matrices_cfg = make_equal_k_two_matrices_config(equal_k_bias)
    tm_fer, tm_effort, tm_latency, tm_counts = run_sim(
        two_matrices_cfg, snr_points, os.path.join(results_root, "equalk_pattern0_2matrices")
    )
    for snr in snr_points:
        print(f"  2-patterns/1-matrix  FER @ {snr}dB: {_fer_ci_str(tp_fer, tp_counts, snr)}")
        print(f"  pattern0/2-matrices  FER @ {snr}dB: {_fer_ci_str(tm_fer, tm_counts, snr)}")

    multi_pair_results = equal_k_multi_pair_check(
        num_pairs=6, bias=equal_k_bias, snr_points=snr_points, results_root=results_root,
    )

    summary = {
        "sanity_checks": {
            "bias0_collapses_patterns": bool(b0_ok),
            "bias_changes_outcome": bool(bchange_ok),
            "hard_limit_agreement_rate": agree_rate,
            "hard_limit_conditional_agreement_rate": cond_agree_rate,
            "hard_limit_disagreement_rate_95ci": list(disagree_ci),
            "hard_limit_num_samples": hard_limit_n,
            "hard_limit_hard_not_converged_count": hnc_count,
            "hard_limit_hard_not_converged_agreement_rate": hnc_agree_rate,
            "hard_limit_hard_not_converged_disagreement_rate_95ci": list(hnc_disagree_ci),
            "reused_decoder_matches_fresh": bool(reuse_ok),
            "effort_mechanism_by_guess_correctness": effort_mechanism,
        },
        "snr_points": snr_points,
        "hard": {"fer": hard_fer, "effort": hard_effort, "latency": hard_latency, "counts": hard_counts},
        "soft": soft_results,
        "equal_k_check": {
            "bias": equal_k_bias,
            "two_patterns_one_matrix": {"fer": tp_fer, "effort": tp_effort, "latency": tp_latency, "counts": tp_counts},
            "pattern0_two_matrices": {"fer": tm_fer, "effort": tm_effort, "latency": tm_latency, "counts": tm_counts},
            "multi_pair_check": multi_pair_results,
        },
        "splitter_row_weight": int(s.sum()),
        "syndrome_vn_degree": 1,  # see design-decision note: M_l = identity for this Delta=1 prototype
    }
    out_path = os.path.join(results_root, "summary.json")
    os.makedirs(results_root, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nWrote summary to {out_path}")


if __name__ == "__main__":
    main()
