# pylint: disable=invalid-name
"""Quick scouting step for the Zc=6 greedy block-selection search: find the
SNR operating points achieving (roughly) FER=1e-1 and FER=1e-3 for a
"reasonable baseline config" (subsequent blocks, i.e. block_offset=0, split3
[Delta=3 sub-splitting], num_used_blocks=4 - matches the reference script's
"64_split3" variant). Uses the reference script's own search_target_fer()
binary-search helper (copied verbatim) with a CHEAP evaluation budget
(target_errors=50, max_transmissions=2e5) since this is only for picking two
ranking operating points for the greedy chain, not a production number.

Unlike the Zc=11 greedy search (which already had a completed full_parallel
production sweep to read 1.8dB/3.3dB operating points off), there is no
existing cheap-budget full_parallel baseline for Zc=6, so this step is
required (see sweeps/STATUS.md and the task description).

Usage:
    python sweeps/scout_zc6_fer_targets.py <n_simul: 78|180>

Writes sweeps/zc6_scout_<n_simul>.json:
    {"n_simul": .., "fer_0.1": {"snr": .., "fer": ..}, "fer_0.001": {"snr": .., "fer": ..}}
"""

import json
import os
import sys

import numpy as np
import galois

gf2 = galois.GF2

import channel_code_lib2

if len(sys.argv) != 2:
    raise ValueError("Usage: python sweeps/scout_zc6_fer_targets.py <n_simul: 78|180>")

n_simul = int(sys.argv[1])
if n_simul not in (78, 180):
    raise ValueError(f"Unsupported n_simul={n_simul}, expected 78 or 180")

SPLIT_PATTERN = [3]  # Delta=3, same as the greedy search itself
NUM_USED_BLOCKS = 4  # "mid-size L" baseline, matches the "64_split3" variant
BLOCK_OFFSET = 0  # "subsequent blocks" i.e. the reference script's baseline

bg_vn = 12
bg_cn = 4
Zc = 6

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
use_all_zero = False

# Rough starting guesses per code, informed by the existing natural-order
# baseline data at /home/pj9034/aSCED/RESULTS/fig_scatter_zc6_n={78,180}/
# (spot-checked before writing this script - split=3,numblocks=4 required_snr
# is ~5.07dB for n=78, ~2.99dB for n=180 at target_fer=1e-3).
START_SNR = {
    78: {"fer_0.1": 3.0, "fer_0.001": 5.0},
    180: {"fer_0.1": 1.0, "fer_0.001": 2.9},
}

# ---- Same PCM construction as reproduce_fig5_scatter_plot_zc6_5G_LDPC.py ----
code = np.load("Codes/TCOM_aSCED/5G_zc=6/max_rank_5G_zc=6.npz")
H_full = code["h"].astype(int)

number_vn_simul = n_simul + 2 * Zc
number_vns_start = bg_vn * Zc + 2 * Zc
harq_part = (number_vn_simul - number_vns_start) // Zc
m = Zc * bg_cn + harq_part * Zc
number_vns = number_vns_start + harq_part * Zc
assert number_vn_simul == number_vns

H = gf2(H_full[:m, :number_vns])
G = np.array(gf2(H[:m, :number_vns]).null_space()).astype(int)
k, n = G.shape
assert n == number_vns
print(f"Simulating number vns={n}, n={number_vns - len(message_bit_pucturing)}, k={k}")

candidate_rows = H_full[m:, : number_vns_start + harq_part * Zc]
number_additional_cyclic_blocks = candidate_rows.shape[0] // Zc
assert number_additional_cyclic_blocks * Zc == candidate_rows.shape[0]

enc_cfg = channel_code_lib2.PCM_Encoder_config(H, k, n)


def binary_vectors_in_suffix(number_rows, number_last_entries):
    v = [0] * number_rows
    for mask in range(1, 1 << number_last_entries):
        for i in range(number_last_entries):
            v[number_rows - number_last_entries + i] = (mask >> i) & 1
        yield v.copy()


def create_bp_config(H_aux):
    cfg = channel_code_lib2.BP_config(H_aux)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    return cfg


def create_asced_config(H, candidate_rows, block_offset, Zc, split_pattern, num_used_blocks, use_all_zero):
    """Verbatim contiguous-block construction (block_offset=0, blocks
    0..num_used_blocks-1) - the "subsequent blocks" baseline convention.
    """
    path_configs = []
    for block in range(num_used_blocks):
        assert block + block_offset <= number_additional_cyclic_blocks
        row_block = candidate_rows[(block + block_offset) * Zc : (block + block_offset + 1) * Zc]
        row_segments = [row_block] if len(split_pattern) == 0 else np.split(row_block, split_pattern)
        for rows in row_segments:
            H_aux = gf2(np.vstack((H, rows)))
            m_aux = H_aux.shape[0]
            path_configs.append(create_bp_config(H_aux))
            if not use_all_zero:
                rank_difference = rows.shape[0]
                for offset in binary_vectors_in_suffix(m_aux, rank_difference):
                    cfg = create_bp_config(H_aux)
                    cfg.affine_offset = offset
                    path_configs.append(cfg)
    print("Simulated num. aSCED paths:", len(path_configs))
    return channel_code_lib2.Ensemble_config(H, path_configs), len(path_configs)


def run_at_snr(sim, snr):
    sim.get_error_rates([snr])
    return sim.error_rates["FER-SNR"][snr]


def search_target_fer(
    sim,
    start_snr,
    target_fer,
    tolerance=0.02,
    initial_step=0.2,
    min_step=0.0001,
    max_iter=40,
):
    """Copied verbatim from reproduce_fig5_scatter_plot_zc6_5G_LDPC.py."""
    lower_fer = target_fer * (1 - tolerance)
    upper_fer = target_fer * (1 + tolerance)

    def evaluate(snr):
        fer = run_at_snr(sim, snr)
        print(f"SNR = {snr:.3f} dB, FER = {fer:.3e}")
        return fer

    snr = start_snr
    fer = evaluate(snr)
    if lower_fer <= fer <= upper_fer:
        return snr, fer

    step = initial_step
    if fer > target_fer:
        low_snr = snr
        low_fer = fer
        while True:
            snr += step
            fer = evaluate(snr)
            if lower_fer <= fer <= upper_fer:
                return snr, fer
            if fer < target_fer:
                high_snr = snr
                high_fer = fer
                break
    else:
        high_snr = snr
        high_fer = fer
        while True:
            snr -= step
            fer = evaluate(snr)
            if lower_fer <= fer <= upper_fer:
                return snr, fer
            if fer > target_fer:
                low_snr = snr
                low_fer = fer
                break

    for _ in range(max_iter):
        if high_snr - low_snr < min_step:
            print("stop due to too small step")
            break
        mid_snr = 0.5 * (low_snr + high_snr)
        mid_fer = evaluate(mid_snr)
        if lower_fer <= mid_fer <= upper_fer:
            return mid_snr, mid_fer
        if mid_fer > target_fer:
            low_snr = mid_snr
            low_fer = mid_fer
        else:
            high_snr = mid_snr
            high_fer = mid_fer

    if abs(low_fer - target_fer) < abs(high_fer - target_fer):
        return low_snr, low_fer
    return high_snr, high_fer


ensemble_cfg, ensemble_size = create_asced_config(
    H=H,
    candidate_rows=candidate_rows,
    block_offset=BLOCK_OFFSET,
    Zc=Zc,
    split_pattern=SPLIT_PATTERN,
    num_used_blocks=NUM_USED_BLOCKS,
    use_all_zero=use_all_zero,
)

results = {"n_simul": n_simul, "ensemble_size": ensemble_size, "num_used_blocks": NUM_USED_BLOCKS}

for key, target_fer in (("fer_0.1", 0.1), ("fer_0.001", 0.001)):
    sim = channel_code_lib2.Simulation_Env(k, n, "all")
    sim.target_errors = 50
    sim.max_transmissions = int(2e5)
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus is not None:
        sim.num_threads = int(slurm_cpus)
    sim.auto_save = auto_save
    sim.save_dir = f"RESULTS/greedy_search_zc6/n{n_simul}/scout/{key}"
    sim.puncturing(message_bit_pucturing)
    sim.init(enc_cfg, ensemble_cfg, use_all_zero)

    required_snr, accepted_fer = search_target_fer(sim, START_SNR[n_simul][key], target_fer)
    print(f"n_simul={n_simul} target_fer={target_fer}: snr={required_snr:.4f}, fer={accepted_fer:.4e}")
    results[key] = {"snr": required_snr, "fer": accepted_fer}

out_path = f"sweeps/zc6_scout_{n_simul}.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)
print("wrote", out_path)
print(json.dumps(results, indent=2))
