# pylint: disable=invalid-name
"""Production validation script for the GREEDY block-selection search
(sweeps/greedy_block_search.py + sweeps/run_greedy_search.py).

Same RL/QC-block aSCED ensemble construction and CLI shape as
reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py, but instead of using
`num_used_blocks` *contiguous* blocks starting at block_offset=0 ("the
subsequent blocks" baseline), it uses the specific, greedy-discovered list of
block indices recorded in `sweeps/greedy_search_state_<variant>.json`
("fixed_blocks", in discovery order: index 0 was picked first / is the most
valuable).

Path order matters here (unlike in the search driver, which only cared about
ensemble membership): `create_asced_config` below iterates `block_indices` in
the exact order given (NOT sorted), so for `fixed_sequential` the resulting
path index order naturally starts with the most-valuable (first-picked)
block's paths first - this is the "intrinsic ordering" insight from the
greedy search (see STATUS.md).

Only supports the two ensembles the greedy search was run for: asced48
(split pattern [2,4,6,8], 2 blocks) and asced384 (split pattern [5], 4
blocks).

Usage: identical to the sequential script:
    python reproduce_fig3_RL_zc11_asced_5G_LDPC_greedy.py \
        <variant: asced48|asced384> <n_simul> <snr_start> <snr_end> \
        <selector_type: full_parallel|fixed_sequential|syndrome_sequential> \
        <members_per_group> <target_num_converged (ignored for full_parallel)>
"""

import json
import os
import sys

import numpy as np

import galois

gf2 = galois.GF2

import channel_code_lib2

if len(sys.argv) != 8:
    raise ValueError(
        "Usage: python reproduce_fig3_RL_zc11_asced_5G_LDPC_greedy.py "
        "<variant: asced48|asced384> <n_simul> <snr_start> <snr_end> "
        "<selector_type: full_parallel|fixed_sequential|syndrome_sequential> "
        "<members_per_group> <target_num_converged (ignored for full_parallel)>"
    )

decoder_variant = sys.argv[1].lower()
n_simul = int(sys.argv[2])
snr_start = float(sys.argv[3])
snr_end = float(sys.argv[4])
selector_type = sys.argv[5].lower()
members_per_group = int(sys.argv[6])
target_num_converged = int(sys.argv[7])

# Include the end point
sim_regime = np.arange(snr_start, snr_end + 0.25, 0.5)

remove = 4  # k=66 code rather than k=110
bg_vn = 12 - remove  # bg2
bg_cn = 4  # bg2
Zc = 11

# "greedy"/"opt" in the path so this never collides with the baseline sweep's
# RESULTS/fig_x_zc11_r4_seq_<selector>_mpg<mpg>_n<n_simul>/ directories.
results_dir = f"RESULTS/fig_x_zc11_r{remove}_seq_greedy_{selector_type}_mpg{members_per_group}_n{n_simul}"

SPLIT_PATTERNS = {
    "asced48": [2, 4, 6, 8],  # sub-blocks of size 2,2,2,2,3 (Delta in {2,3})
    "asced384": [5],  # sub-blocks of size 5,6 (Delta in {5,6})
}
STATE_FILES = {
    "asced48": "sweeps/greedy_search_state_asced48.json",
    "asced384": "sweeps/greedy_search_state_asced384.json",
}
EXPECTED_NUM_BLOCKS = {"asced48": 2, "asced384": 4}

if decoder_variant not in SPLIT_PATTERNS:
    raise ValueError(f"Unknown variant '{decoder_variant}', expected one of {list(SPLIT_PATTERNS)}")

split_pattern = SPLIT_PATTERNS[decoder_variant]

with open(STATE_FILES[decoder_variant]) as f:
    _state = json.load(f)
block_indices = _state["fixed_blocks"]  # greedy discovery order, 1st pick first
if len(block_indices) != EXPECTED_NUM_BLOCKS[decoder_variant]:
    raise ValueError(
        f"{STATE_FILES[decoder_variant]} has {len(block_indices)} fixed_blocks, "
        f"expected {EXPECTED_NUM_BLOCKS[decoder_variant]} - greedy search for "
        f"{decoder_variant} may not have finished."
    )
print(f"Using greedy-discovered block order for {decoder_variant}: {block_indices}")

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
use_all_zero = False

# Same PCM construction as reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py.
code = np.load("Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz")
H_full = code["h"].astype(int)
H_full = np.delete(H_full, np.arange(6 * Zc, 10 * Zc), axis=1)

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
for b in block_indices:
    assert 0 <= b < number_additional_cyclic_blocks, f"block index {b} out of range"


def binary_vectors_in_suffix(number_rows, number_last_entries):
    v = [0] * number_rows
    for mask in range(1, 1 << number_last_entries):
        for i in range(number_last_entries):
            v[number_rows - number_last_entries + i] = (mask >> i) & 1
        yield v.copy()


if not use_all_zero:
    enc_cfg = channel_code_lib2.PCM_Encoder_config(H, k, n)


def create_bp_config(H_aux):
    """Create a BP configuration with the default decoder settings."""
    cfg = channel_code_lib2.BP_config(H_aux)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    return cfg


def create_asced_config(H, candidate_rows, block_indices, Zc, split_pattern, use_all_zero):
    """Construct all decoder paths for an aSCED configuration built from the
    greedy-discovered LIST of block indices, in discovery order (unlike the
    search driver, order here is kept exactly as given, not sorted, because
    it determines fixed_sequential's path index order).
    """
    batches = []

    for block in block_indices:
        row_block = candidate_rows[block * Zc : (block + 1) * Zc]

        row_segments = [row_block] if len(split_pattern) == 0 else np.split(row_block, split_pattern)

        for rows in row_segments:
            H_aux = gf2(np.vstack((H, rows)))
            m_aux = H_aux.shape[0]

            batch = [create_bp_config(H_aux)]  # zero affine offset

            if not use_all_zero:
                rank_difference = rows.shape[0]
                for offset in binary_vectors_in_suffix(m_aux, rank_difference):
                    cfg = create_bp_config(H_aux)
                    cfg.affine_offset = offset
                    batch.append(cfg)

            batches.append(batch)

    path_configs = [cfg for batch in batches for cfg in batch]
    print("Simulated num. aSCED paths:", len(path_configs), "in", len(batches), "batches")

    return channel_code_lib2.Ensemble_config(H, path_configs), len(path_configs)


def run_asced(save_dir_name):
    ensemble_cfg, ensemble_size = create_asced_config(
        H=H,
        candidate_rows=candidate_rows,
        block_indices=block_indices,
        Zc=Zc,
        split_pattern=split_pattern,
        use_all_zero=use_all_zero,
    )

    if selector_type != "full_parallel":
        ensemble_cfg.set_mConvergedConfig(target_num_converged)

        if selector_type == "fixed_sequential":
            ensemble_cfg.set_fixed_sequential(members_per_group)
        elif selector_type == "syndrome_sequential":
            ensemble_cfg.set_syndrome_sequential(members_per_group)
        else:
            raise ValueError(f"Unknown selector_type '{selector_type}'")

    print("starting", save_dir_name)

    sim = channel_code_lib2.Simulation_Env(k, n, "all")
    # Production budget, matching the baseline sweep (reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py).
    sim.target_errors = 200
    sim.max_transmissions = int(2e6)
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus is not None:
        sim.num_threads = int(slurm_cpus)
    sim.auto_save = auto_save
    # One directory per (config, SNR) - see STATUS.md's race-condition note.
    snr_tag = f"{sim_regime[0]:g}"
    sim.save_dir = results_dir + "/" + save_dir_name + "/" + f"snr_{snr_tag}"

    sim.puncturing(message_bit_pucturing)

    if use_all_zero:
        sim.all_zero_init(ensemble_cfg)
    else:
        sim.init(enc_cfg, ensemble_cfg, use_all_zero)

    sim.get_error_rates(sim_regime)

    print(save_dir_name, "finished")
    print(sim.error_rates["FER-SNR"])

    return sim.error_rates["FER-SNR"]


if selector_type == "full_parallel":
    save_name = f"{decoder_variant}_{selector_type}"
else:
    save_name = f"{decoder_variant}_{selector_type}_mpg{members_per_group}_target{target_num_converged}"
FER_results = {save_name: run_asced(save_name)}
