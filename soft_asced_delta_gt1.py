# pylint: disable=invalid-name
"""Soft-aSCED, Delta>1 (multi-row-per-segment) generalization, validated
directly on this project's primary paper-reproduction ensembles ("aSCED-48"
and "aSCED-384").

This generalizes ``prototype_soft_asced_single_row.py``'s Delta=1 result
(see ``sweeps/soft_asced_single_row_results.md``) to the real multi-row
splitter segments used by ``reproduce_fig3_RL_zc11_asced_5G_LDPC.py``'s own
``create_asced_config``.

Key scoping finding (see ``sweeps/soft_asced_delta_gt1_results.md`` for the
full writeup): hard aSCED's existing Delta>1 machinery does NOT use any
general overcomplete-PCM/``M_l`` row-combination construction -- it is
``M_l = identity``. For a segment of ``Delta_seg`` appended candidate rows
it just builds ``H_aux = [H; segment_rows]`` (vstack, nothing combined),
then enumerates ALL ``2**Delta_seg`` hard ``affine_offset`` patterns. This
script mirrors that exact block/segment loop
(``create_asced_config``'s structure, duplicated here rather than imported,
since the original script gates on ``sys.argv`` and is otherwise left
untouched per scope) and only changes the per-pattern path construction:
hard bakes the guess into ``affine_offset``; soft instead adds
``Delta_seg`` ordinary binary "syndrome VN" columns
(``H_soft = [H_aux | Delta_seg-column identity block]``, degree-1
connectivity per column, exactly like the validated Delta=1 case, just
repeated ``Delta_seg`` times per segment) with a shared scalar
``syndrome_vn_bias``.

"asced48" = ``splitting_pattern[1] = [2,4,6,8]``, ``num_used_blocks=2``:
2 blocks x (4 segments of Delta=2 + 1 segment of Delta=3) = 2*(4*4+8) = 48
paths.
"asced384" = ``splitting_pattern[2] = [5]``, ``num_used_blocks=4``:
4 blocks x (1 segment of Delta=5 + 1 of Delta=6) = 4*(32+64) = 384 paths.

Usage:
    uv run python soft_asced_delta_gt1.py sanity asced48
    uv run python soft_asced_delta_gt1.py sweep asced48
    uv run python soft_asced_delta_gt1.py sweep asced384
"""

import itertools
import json
import os
import shutil
import sys

import galois
import numpy as np

gf2 = galois.GF2

import channel_code_lib2 as ccl

# --------------------------------------------------------------------------
# Code setup (verbatim parameters from reproduce_fig3_RL_zc11_asced_5G_LDPC.py
# / prototype_soft_asced_single_row.py, n_simul=132 -- the "aSCED-48"/
# "aSCED-384" code family, both built from the SAME underlying code, only
# differing in which/how many candidate rows get appended).
# --------------------------------------------------------------------------

n_simul = 132
remove = 4
bg_vn = 12 - remove  # 8
bg_cn = 4
Zc = 11
block_offset = 0

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

number_additional_cyclic_blocks = candidate_rows.shape[0] // Zc
assert number_additional_cyclic_blocks * Zc == candidate_rows.shape[0]

print(f"Code: n={n} (n_transmitted={n_transmitted}), k={k}, m={m} "
      f"(aSCED-48/384 family, n_simul={n_simul})")

enc_cfg = ccl.PCM_Encoder_config(H, k, n)

# splitting_pattern[1] = [2,4,6,8] -> asced48; splitting_pattern[2] = [5] -> asced384
SPLIT_PATTERN_48 = [2, 4, 6, 8]
SPLIT_PATTERN_384 = [5]

BP_DEFAULTS = dict(
    early_stopping=True,
    max_iterations=32,
    cn_update_type="msa",
    scheduling_type="flooding",
    norm_factor=0.75,
)


def _apply_bp_defaults(cfg):
    for key, val in BP_DEFAULTS.items():
        setattr(cfg, key, val)
    return cfg


def binary_vectors_in_suffix(number_rows, number_last_entries):
    """Verbatim copy of reproduce_fig3_RL_zc11_asced_5G_LDPC.py's helper:
    yields all 2**number_last_entries - 1 NONZERO patterns for the last
    `number_last_entries` entries of a length-`number_rows` vector (prefix
    zero)."""
    v = [0] * number_rows
    for mask in range(1, 1 << number_last_entries):
        for i in range(number_last_entries):
            v[number_rows - number_last_entries + i] = (mask >> i) & 1
        yield v.copy()


def get_row_segments(candidate_rows, block_offset, Zc, split_pattern, num_used_blocks):
    """Mirrors create_asced_config's block/segment loop, returning the
    actual appended-row segments (not yet turned into BP_configs) so hard
    and soft construction can share the IDENTICAL row selection."""
    segments = []  # list of (block_idx, rows) with rows shape (Delta_seg, n)
    for block in range(num_used_blocks):
        assert block + block_offset <= number_additional_cyclic_blocks
        row_block = candidate_rows[
            (block + block_offset) * Zc : (block + block_offset + 1) * Zc
        ]
        if len(split_pattern) == 0:
            row_segments = [row_block]
        else:
            row_segments = np.split(row_block, split_pattern)
        for rows in row_segments:
            segments.append((block, rows))
    return segments


def create_hard_asced_config(H, candidate_rows, block_offset, Zc, split_pattern, num_used_blocks):
    """Hard aSCED baseline -- a near-verbatim copy of
    reproduce_fig3_RL_zc11_asced_5G_LDPC.py's create_asced_config with
    use_all_zero=False hardcoded (the project's standard asced48/384
    behavior: every path, zero affine offset AND all nonzero suffix
    patterns, i.e. all 2**Delta_seg patterns per segment)."""
    path_configs = []
    for block, rows in get_row_segments(candidate_rows, block_offset, Zc, split_pattern, num_used_blocks):
        H_aux = gf2(np.vstack((H, rows)))
        m_aux = H_aux.shape[0]
        rank_difference = rows.shape[0]

        cfg = ccl.BP_config(H_aux)
        _apply_bp_defaults(cfg)
        path_configs.append(cfg)  # zero affine offset

        for offset in binary_vectors_in_suffix(m_aux, rank_difference):
            cfg = ccl.BP_config(H_aux)
            _apply_bp_defaults(cfg)
            cfg.affine_offset = offset
            path_configs.append(cfg)

    print("Simulated num. hard aSCED paths:", len(path_configs))
    return ccl.Ensemble_config(H, path_configs)


def make_soft_bp_config(H_aux_int, Delta_seg, pattern, bias, H0):
    """H_aux_int: (m_aux, n) int ndarray (NOT yet extended with syndrome-VN
    columns). Builds H_soft = [H_aux | Delta_seg-column identity block],
    identity placed in the LAST Delta_seg rows (the appended segment's own
    rows) -- matching the already-validated Delta=1 degree-1 construction,
    just repeated Delta_seg times. One BP_config per full Delta_seg-bit
    pattern."""
    m_aux = H_aux_int.shape[0]
    ident_block = np.zeros((m_aux, Delta_seg), dtype=int)
    ident_block[m_aux - Delta_seg:, :] = np.eye(Delta_seg, dtype=int)
    H_soft = np.hstack((H_aux_int, ident_block)).astype(int)

    cfg = ccl.BP_config(gf2(H_soft))
    _apply_bp_defaults(cfg)
    cfg.use_avns = True
    cfg.num_syndrome_vns = Delta_seg
    cfg.syndrome_vn_pattern = list(pattern)
    cfg.syndrome_vn_bias = float(bias)
    cfg.set_H0(gf2(H0))
    return cfg


def create_soft_asced_config(H, candidate_rows, block_offset, Zc, split_pattern, num_used_blocks, bias):
    """Soft-aSCED generalization of create_asced_config: SAME block/segment
    loop and row selection as create_hard_asced_config (via
    get_row_segments), SAME per-segment path count (2**Delta_seg, matching
    hard exactly) -- only the per-pattern path construction differs (soft
    syndrome VNs instead of a hard-baked affine_offset). `bias` is one
    shared scalar swept across ALL segments/patterns in this first pass."""
    path_configs = []
    for block, rows in get_row_segments(candidate_rows, block_offset, Zc, split_pattern, num_used_blocks):
        H_aux_int = np.array(gf2(np.vstack((H, rows)))).astype(int)
        Delta_seg = rows.shape[0]
        for pattern in itertools.product((0, 1), repeat=Delta_seg):
            path_configs.append(make_soft_bp_config(H_aux_int, Delta_seg, pattern, bias, H))

    print("Simulated num. soft aSCED paths:", len(path_configs))
    return ccl.Ensemble_config(H, path_configs)


VARIANTS = {
    "asced48": dict(split_pattern=SPLIT_PATTERN_48, num_used_blocks=2),
    "asced384": dict(split_pattern=SPLIT_PATTERN_384, num_used_blocks=4),
}


def expected_num_paths(split_pattern, num_used_blocks):
    segs = get_row_segments(candidate_rows, block_offset, Zc, split_pattern, num_used_blocks)
    return sum(2 ** rows.shape[0] for _, rows in segs), [rows.shape[0] for _, rows in segs]


# --------------------------------------------------------------------------
# Shared helpers: random-message/channel simulation (verbatim convention
# from prototype_soft_asced_single_row.py), Wilson CI, ensemble-sim runner.
# --------------------------------------------------------------------------

NEAR_ZERO = 1e-10  # matches channel_coding_lib's safe::NEAR_ZERO erasure-LLR magnitude


def awgn_llrs(codeword_bits, snr_db, rng):
    """BPSK + AWGN channel LLRs for the FULL (un-punctured) n-length
    codeword, punctured positions overwritten by a near-zero erasure LLR of
    random sign -- matches Simulation_Env::simulate_transmission. Sign
    convention: LLR < 0 => hard bit 1."""
    r = k / n_transmitted
    snr_lin = 10 ** (snr_db / 10.0)
    L_c = 4.0 * r * snr_lin
    sigma = 1.0 / np.sqrt(2.0 * r * snr_lin)
    bpsk = 1.0 - 2.0 * np.asarray(codeword_bits, dtype=float)
    noisy = bpsk + rng.normal(0.0, sigma, size=len(codeword_bits))
    llrs = noisy * L_c
    for i in message_bit_pucturing:
        llrs[i] = NEAR_ZERO if rng.integers(0, 2) == 0 else -NEAR_ZERO
    return llrs.tolist()


def encode_random_message(rng):
    msg = rng.integers(0, 2, size=k)
    codeword = np.array(gf2(msg) @ gf2(G)).astype(np.int64)
    return msg, codeword


def wilson_interval(num_events, num_trials, z=1.96):
    """Wilson score interval for a binomial proportion (95% by default).
    Falls back to the "rule of three" upper bound when num_events == 0."""
    if num_trials == 0:
        return (float("nan"), float("nan"))
    if num_events == 0:
        return (0.0, 3.0 / num_trials)
    phat = num_events / num_trials
    denom = 1 + z**2 / num_trials
    center = (phat + z**2 / (2 * num_trials)) / denom
    half = z * np.sqrt(phat * (1 - phat) / num_trials + z**2 / (4 * num_trials**2)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def _lookup_by_snr(d, snr):
    for k_str, v in d.items():
        if abs(float(k_str) - snr) < 1e-9:
            return v
    return None


def run_sim(ensemble_cfg, snr_points, save_dir, target_errors=200, max_transmissions=int(2e6)):
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
    effort, latency, counts = {}, {}, {}

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


# --------------------------------------------------------------------------
# Generalized sanity checks (Delta>1): same four checks as
# prototype_soft_asced_single_row.py, generalized for a multi-bit pattern
# per segment. All operate on ONE representative segment at a time (passed
# in as `rows`, an ndarray of shape (Delta_seg, n)).
# --------------------------------------------------------------------------

def sanity_check_1_bias_applied(rows, label):
    """bias=0 must collapse EVERY pattern (incl. all-erasure) to the same
    decision; a large bias must make different full patterns diverge, and
    changing a SINGLE bit of the pattern (holding the rest fixed) must also
    be able to change the outcome (checks per-VN independence, not just
    whole-pattern on/off)."""
    print(f"\n--- Sanity check 1: bias is really applied ({label}, Delta_seg={rows.shape[0]}) ---")
    Delta_seg = rows.shape[0]
    H_aux_int = np.array(gf2(np.vstack((H, rows)))).astype(int)
    rng = np.random.default_rng(101)
    _, codeword = encode_random_message(rng)
    llrs = awgn_llrs(codeword, snr_db=0.5, rng=rng)

    all_patterns = list(itertools.product((0, 1), repeat=Delta_seg))
    erasure_pattern = tuple([-1] * Delta_seg)

    def decode_with(pattern, bias):
        cfg = make_soft_bp_config(H_aux_int, Delta_seg, pattern, bias, H)
        dec = ccl.build_decoder(cfg)
        dec.decode(llrs)
        return tuple(dec.get_codeword())

    bias0_results = {p: decode_with(p, 0.0) for p in all_patterns}
    bias0_results[erasure_pattern] = decode_with(erasure_pattern, 0.0)
    bias0_collapses = len(set(bias0_results.values())) == 1

    zero_pat = tuple([0] * Delta_seg)
    one_pat = tuple([1] * Delta_seg)
    full_flip_differs = decode_with(zero_pat, 30.0) != decode_with(one_pat, 30.0)

    # Single-bit flip (holding the rest at 0) must ALSO be able to change
    # the decision relative to the all-zero pattern -- confirms each
    # syndrome-VN column has its own independent effect, not just "any
    # guess vs. its total opposite".
    single_bit_pattern = tuple([1] + [0] * (Delta_seg - 1))
    single_bit_differs_from_zero = decode_with(zero_pat, 30.0) != decode_with(single_bit_pattern, 30.0)

    print(f"  bias=0 collapses ALL {len(all_patterns)} patterns + erasure to the same decision: {bias0_collapses}")
    print(f"  bias=30, all-zero vs all-one pattern differ: {full_flip_differs}")
    print(f"  bias=30, single-bit-flip vs all-zero pattern differ: {single_bit_differs_from_zero}")
    return bias0_collapses, full_flip_differs, single_bit_differs_from_zero


def sanity_check_2_hard_limit_equivalence(rows, label, num_samples=1200, snr_db=2.0):
    """On samples with a KNOWN true per-segment pattern, soft with
    pattern=true bits, bias=30 (saturating) vs. hard aSCED's path using
    that SAME correct full guess should (nearly) agree. Reports a Wilson
    95% CI on the DISAGREEMENT rate (not just a raw agreement fraction --
    the Delta=1 review's explicit correction), both unconditionally and
    restricted to the subset where hard did NOT converge (the only place
    saturation and a hard-fixed guess could plausibly differ)."""
    print(f"\n--- Sanity check 2: hard-limit equivalence ({label}, Delta_seg={rows.shape[0]}, N={num_samples}) ---")
    Delta_seg = rows.shape[0]
    H_aux_int = np.array(gf2(np.vstack((H, rows)))).astype(int)
    m_aux = H_aux_int.shape[0]
    rng = np.random.default_rng(102)

    agreements = 0
    hard_not_converged_count = 0
    hard_not_converged_agreements = 0
    both_converged_count = 0
    both_converged_agreements = 0

    for _ in range(num_samples):
        _, codeword = encode_random_message(rng)
        true_pattern = tuple(int(b) for b in (rows @ codeword) % 2)
        llrs = awgn_llrs(codeword, snr_db=snr_db, rng=rng)

        hard_offset = np.zeros(m_aux, dtype=int)
        hard_offset[m_aux - Delta_seg:] = true_pattern
        hard_cfg = ccl.BP_config(H_aux_int)
        _apply_bp_defaults(hard_cfg)
        hard_cfg.affine_offset = list(hard_offset)
        hard_dec = ccl.build_decoder(hard_cfg)
        hard_dec.decode(llrs)
        hard_cw = hard_dec.get_codeword()

        soft_cfg = make_soft_bp_config(H_aux_int, Delta_seg, true_pattern, 30.0, H)
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
            hard_not_converged_count += 1
            if agree:
                hard_not_converged_agreements += 1

    disagreements = num_samples - agreements
    agreement_rate = agreements / num_samples
    disagree_ci = wilson_interval(disagreements, num_samples)
    cond_rate = both_converged_agreements / both_converged_count if both_converged_count else float("nan")
    hnc_rate = (
        hard_not_converged_agreements / hard_not_converged_count if hard_not_converged_count else float("nan")
    )
    hnc_disagree_ci = wilson_interval(hard_not_converged_count - hard_not_converged_agreements, hard_not_converged_count)

    print(f"  overall decision agreement: {agreement_rate:.4f} ({agreements}/{num_samples})")
    print(f"  => disagreement rate 95% CI: [{disagree_ci[0]:.5f}, {disagree_ci[1]:.5f}] ({disagreements} observed disagreements)")
    print(f"  agreement given BOTH converged: {cond_rate:.4f} ({both_converged_agreements}/{both_converged_count})")
    print(f"  agreement given HARD DID NOT CONVERGE: {hnc_rate:.4f} ({hard_not_converged_agreements}/{hard_not_converged_count}); "
          f"disagreement 95% CI [{hnc_disagree_ci[0]:.5f}, {hnc_disagree_ci[1]:.5f}]")

    return dict(
        agreement_rate=agreement_rate, disagreements=disagreements, num_samples=num_samples,
        disagreement_ci=list(disagree_ci), conditional_agreement_rate=cond_rate,
        both_converged_count=both_converged_count,
        hard_not_converged_count=hard_not_converged_count,
        hard_not_converged_agreement_rate=hnc_rate,
        hard_not_converged_disagreement_ci=list(hnc_disagree_ci),
    )


def sanity_check_3_reused_decoder(rows, label, num_samples=50, snr_db=2.0, bias=10.0):
    """Decoding the same samples once on a reused decoder object (with an
    unrelated sample decoded in between) and once with a fresh decoder per
    sample must give byte-for-byte identical results."""
    print(f"\n--- Sanity check 3: reused decoder vs. fresh decoder ({label}) ---")
    Delta_seg = rows.shape[0]
    H_aux_int = np.array(gf2(np.vstack((H, rows)))).astype(int)
    pattern = tuple([0] * Delta_seg)
    rng = np.random.default_rng(103)
    samples = []
    for _ in range(num_samples):
        _, codeword = encode_random_message(rng)
        samples.append(awgn_llrs(codeword, snr_db=snr_db, rng=rng))

    cfg = make_soft_bp_config(H_aux_int, Delta_seg, pattern, bias, H)
    reused_decoder = ccl.build_decoder(cfg)
    reused_results = []
    for llrs in samples:
        reused_decoder.decode(llrs)
        reused_results.append((tuple(reused_decoder.get_codeword()), reused_decoder.converged))

    fresh_results = []
    for llrs in samples:
        fresh_cfg = make_soft_bp_config(H_aux_int, Delta_seg, pattern, bias, H)
        fresh_decoder = ccl.build_decoder(fresh_cfg)
        fresh_decoder.decode(llrs)
        fresh_results.append((tuple(fresh_decoder.get_codeword()), fresh_decoder.converged))

    identical = reused_results == fresh_results
    print(f"  reused-decoder results identical to fresh-decoder-per-sample: {identical}")
    return identical


def sanity_check_4_effort_by_hamming(rows, label, bias=10.0, num_samples=300, snr_db=2.0):
    """Per-path effort/convergence mechanism check, generalized to Delta>1:
    a multi-bit guess can be PARTIALLY wrong, so bucket by Hamming distance
    between the guess pattern and the TRUE pattern (0=fully correct,
    Delta_seg=fully wrong), separately for hard and soft, and report mean
    BP iterations / convergence fraction per bucket."""
    Delta_seg = rows.shape[0]
    print(f"\n--- Sanity check 4: effort/convergence by Hamming distance ({label}, Delta_seg={Delta_seg}, N={num_samples}) ---")
    H_aux_int = np.array(gf2(np.vstack((H, rows)))).astype(int)
    m_aux = H_aux_int.shape[0]
    all_patterns = list(itertools.product((0, 1), repeat=Delta_seg))
    rng = np.random.default_rng(104)

    hard_buckets = {d: ([], []) for d in range(Delta_seg + 1)}
    soft_buckets = {d: ([], []) for d in range(Delta_seg + 1)}

    for _ in range(num_samples):
        _, codeword = encode_random_message(rng)
        true_pattern = tuple(int(b) for b in (rows @ codeword) % 2)
        llrs = awgn_llrs(codeword, snr_db=snr_db, rng=rng)

        for pattern in all_patterns:
            dist = sum(a != b for a, b in zip(pattern, true_pattern))

            offset = np.zeros(m_aux, dtype=int)
            offset[m_aux - Delta_seg:] = pattern
            hard_cfg = ccl.BP_config(H_aux_int)
            _apply_bp_defaults(hard_cfg)
            hard_cfg.affine_offset = list(offset)
            hdec = ccl.build_decoder(hard_cfg)
            hdec.decode(llrs)
            hard_buckets[dist][0].append(hdec.decoding_effort)
            hard_buckets[dist][1].append(bool(hdec.converged))

            soft_cfg = make_soft_bp_config(H_aux_int, Delta_seg, pattern, bias, H)
            sdec = ccl.build_decoder(soft_cfg)
            sdec.decode(llrs)
            soft_buckets[dist][0].append(sdec.decoding_effort)
            soft_buckets[dist][1].append(bool(sdec.converged))

    summary = {}
    for dist in range(Delta_seg + 1):
        he, hc = hard_buckets[dist]
        se, sc = soft_buckets[dist]
        row = dict(
            n=len(he),
            hard_mean_effort=sum(he) / len(he) if he else float("nan"),
            hard_converged_frac=sum(hc) / len(hc) if hc else float("nan"),
            soft_mean_effort=sum(se) / len(se) if se else float("nan"),
            soft_converged_frac=sum(sc) / len(sc) if sc else float("nan"),
        )
        summary[dist] = row
        print(f"  Hamming dist {dist}: n={row['n']:5d}  "
              f"hard mean_effort={row['hard_mean_effort']:6.2f} conv={row['hard_converged_frac']:.3f}  |  "
              f"soft mean_effort={row['soft_mean_effort']:6.2f} conv={row['soft_converged_frac']:.3f}")
    return summary


def run_sanity_suite(variant_name, num_samples_hard_limit=1200, num_samples_hamming=300):
    """Runs all 4 generalized sanity checks on the FIRST segment of the
    given variant (asced48: Delta_seg=2; asced384: Delta_seg=5)."""
    params = VARIANTS[variant_name]
    segs = get_row_segments(candidate_rows, block_offset, Zc, params["split_pattern"], params["num_used_blocks"])
    block, rows = segs[0]
    label = f"{variant_name} segment 0 (block {block}, Delta_seg={rows.shape[0]})"

    results = {}
    results["bias_applied"] = sanity_check_1_bias_applied(rows, label)
    results["hard_limit_equivalence"] = sanity_check_2_hard_limit_equivalence(rows, label, num_samples=num_samples_hard_limit)
    results["reused_decoder_matches_fresh"] = sanity_check_3_reused_decoder(rows, label)
    results["effort_by_hamming_distance"] = sanity_check_4_effort_by_hamming(rows, label, num_samples=num_samples_hamming)
    return results


# --------------------------------------------------------------------------
# Headline hard-vs-soft comparison + bias sweep
# --------------------------------------------------------------------------

BIAS_SWEEP = {
    "asced48": [1, 2, 3, 5, 10, 20, 30, 35],
    # Informed by the asced48 result (bias 1-3 measurably underperforms
    # hard's FER there -- unlike Delta=1 where bias=1 already matched --
    # so bias=5 is included as the empirically-found "sweet spot", plus a
    # low-bias point to document the same underperformance generalizes.
    "asced384": [1, 5, 10, 20, 30],
}

SNR_POINTS = {
    "asced48": [1.0, 2.0, 3.0],
    "asced384": [1.0, 1.5, 2.0],
}

TARGET_ERRORS = {"asced48": 200, "asced384": 150}


def run_variant_sweep(variant_name, results_root=None):
    params = VARIANTS[variant_name]
    snr_points = SNR_POINTS[variant_name]
    bias_values = BIAS_SWEEP[variant_name]
    target_errors = TARGET_ERRORS[variant_name]
    if results_root is None:
        results_root = f"RESULTS/soft_asced_delta_gt1_{variant_name}"
    # Only clean the sim-output subdirs this function (re)writes -- NOT the
    # whole results_root, which may also hold sanity_checks.json from a
    # separate `sanity` run on the same variant.
    shutil.rmtree(os.path.join(results_root, "hard"), ignore_errors=True)
    for bias in BIAS_SWEEP[variant_name]:
        shutil.rmtree(os.path.join(results_root, f"soft_bias{bias}"), ignore_errors=True)

    print("\n" + "=" * 70)
    print(f"{variant_name}: hard aSCED vs. soft aSCED (bias swept)")
    print("=" * 70)

    hard_cfg = create_hard_asced_config(H, candidate_rows, block_offset, Zc, **params)
    hard_fer, hard_effort, hard_latency, hard_counts = run_sim(
        hard_cfg, snr_points, os.path.join(results_root, "hard"), target_errors=target_errors,
    )
    for snr in snr_points:
        print(f"  hard FER @ {snr}dB: {_fer_ci_str(hard_fer, hard_counts, snr)}")
    print(f"  hard avg ensemble effort: {hard_effort}")
    print(f"  hard avg ensemble latency: {hard_latency}")

    soft_results = {}
    for bias in bias_values:
        soft_cfg = create_soft_asced_config(H, candidate_rows, block_offset, Zc, bias=bias, **params)
        fer, effort, latency, counts = run_sim(
            soft_cfg, snr_points, os.path.join(results_root, f"soft_bias{bias}"), target_errors=target_errors,
        )
        soft_results[bias] = dict(fer=fer, effort=effort, latency=latency, counts=counts)
        for snr in snr_points:
            print(f"  soft(bias={bias}) FER @ {snr}dB: {_fer_ci_str(fer, counts, snr)}")
        print(f"  soft(bias={bias}) avg ensemble effort: {effort}")
        print(f"  soft(bias={bias}) avg ensemble latency: {latency}")

    return dict(
        variant=variant_name,
        snr_points=snr_points,
        target_errors=target_errors,
        hard=dict(fer=hard_fer, effort=hard_effort, latency=hard_latency, counts=hard_counts),
        soft=soft_results,
    )


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    mode = sys.argv[1]
    variant = sys.argv[2] if len(sys.argv) > 2 else "asced48"

    if mode == "sanity":
        sanity = run_sanity_suite(variant)
        out_path = f"RESULTS/soft_asced_delta_gt1_{variant}/sanity_checks.json"
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(sanity, f, indent=2, default=str)
        print(f"\nWrote sanity checks to {out_path}")
    elif mode == "sweep":
        sweep = run_variant_sweep(variant)
        out_path = f"RESULTS/soft_asced_delta_gt1_{variant}/summary.json"
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as f:
            json.dump(sweep, f, indent=2, default=str)
        print(f"\nWrote sweep summary to {out_path}")
    else:
        for name, params in VARIANTS.items():
            total, deltas = expected_num_paths(**params)
            print(f"{name}: segments Delta_seg={deltas}, total paths={total}")


if __name__ == "__main__":
    main()
