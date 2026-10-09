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


def sanity_check_2_hard_limit_equivalence(num_samples=200, snr_db=2.0):
    """On samples with a KNOWN true splitter bit, pattern = true bit with
    bias saturated at +/-30 (via syndrome_vn_bias=30) should (nearly) agree
    with the existing hard-aSCED path using that same correct guess."""
    print("\n--- Sanity check 2: hard-limit equivalence (bias=30 saturation vs. hard aSCED) ---")
    rng = np.random.default_rng(2)
    agreements = 0
    both_converged_agreements = 0
    both_converged_count = 0
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

        if hard_cw == soft_cw:
            agreements += 1
        if hard_dec.converged and soft_dec.converged:
            both_converged_count += 1
            if hard_cw == soft_cw:
                both_converged_agreements += 1

    agreement_rate = agreements / num_samples
    conditional_rate = (
        both_converged_agreements / both_converged_count if both_converged_count else float("nan")
    )
    print(f"  overall decision agreement (correct guess, bias=30 vs. hard): "
          f"{agreement_rate:.4f} ({agreements}/{num_samples})")
    print(f"  agreement given BOTH converged: {conditional_rate:.4f} "
          f"({both_converged_agreements}/{both_converged_count})")
    return agreement_rate, conditional_rate


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


# --------------------------------------------------------------------------
# Main bias sweep / hard-vs-soft FER comparison
# --------------------------------------------------------------------------

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
    stats_path = os.path.join(save_dir, "decoder_stats.json")
    effort = {}
    latency = {}
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            stats = json.load(f)

        def _lookup(d, snr):
            # decoder_stats.json keys are produced by the C++ side's own
            # double_to_string (e.g. "2" for 2.0, "2.5" for 2.5) -- match by
            # parsed float value instead of assuming a particular string
            # format.
            for k_str, v in d.items():
                if abs(float(k_str) - snr) < 1e-9:
                    return v
            return None

        for snr in snr_points:
            effort[snr] = _lookup(stats.get("average_ensemble_effort", {}), snr)
            latency[snr] = _lookup(stats.get("average_ensemble_latency", {}), snr)
    return fer, effort, latency


def main():
    print("=" * 70)
    print("Soft-aSCED single-row (Delta=1) prototype -- sanity checks")
    print("=" * 70)
    b0_ok, bchange_ok = sanity_check_1_bias_applied()
    agree_rate, cond_agree_rate = sanity_check_2_hard_limit_equivalence()
    reuse_ok = sanity_check_3_reused_decoder()

    print("\n" + "=" * 70)
    print("Main sweep: hard aSCED (2 paths) vs. soft aSCED (2 paths, bias swept)")
    print("=" * 70)

    snr_points = [2.0, 3.0]
    results_root = "RESULTS/soft_asced_prototype_n132"
    shutil.rmtree(results_root, ignore_errors=True)

    print("\nRunning hard aSCED baseline...")
    hard_cfg = make_hard_ensemble_config()
    hard_fer, hard_effort, hard_latency = run_sim(
        hard_cfg, snr_points, os.path.join(results_root, "hard")
    )
    print(f"  hard FER: {hard_fer}")
    print(f"  hard avg ensemble effort: {hard_effort}")
    print(f"  hard avg ensemble latency: {hard_latency}")

    bias_values = [1, 2, 3, 5, 10, 20, 30, 35]
    soft_results = {}
    for bias in bias_values:
        print(f"\nRunning soft aSCED, bias={bias}...")
        soft_cfg = make_soft_ensemble_config(bias)
        fer, effort, latency = run_sim(
            soft_cfg, snr_points, os.path.join(results_root, f"soft_bias{bias}")
        )
        soft_results[bias] = {"fer": fer, "effort": effort, "latency": latency}
        print(f"  soft(bias={bias}) FER: {fer}")
        print(f"  soft(bias={bias}) avg ensemble effort: {effort}")
        print(f"  soft(bias={bias}) avg ensemble latency: {latency}")

    summary = {
        "sanity_checks": {
            "bias0_collapses_patterns": bool(b0_ok),
            "bias_changes_outcome": bool(bchange_ok),
            "hard_limit_agreement_rate": agree_rate,
            "hard_limit_conditional_agreement_rate": cond_agree_rate,
            "reused_decoder_matches_fresh": bool(reuse_ok),
        },
        "snr_points": snr_points,
        "hard": {"fer": hard_fer, "effort": hard_effort, "latency": hard_latency},
        "soft": soft_results,
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
