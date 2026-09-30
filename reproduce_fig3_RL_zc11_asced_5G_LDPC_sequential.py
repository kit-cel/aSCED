# pylint: disable=invalid-name
"""Sequential/mconverged variant of reproduce_fig3_RL_zc11_asced_5G_LDPC.py.

Same RL/QC-block aSCED ensemble construction (aSCED-48, aSCED-384, ...) for
the 5G LDPC code, extended with a DecoderSelector (grouping the ensemble's
paths into sequentially-executed groups) and an mconverged stopping policy,
to measure how much decoding effort can be saved relative to full-parallel
decoding at matched FER.
"""

import numpy as np
import os
import sys

import galois

gf2 = galois.GF2

import channel_code_lib2

if len(sys.argv) != 8:
    raise ValueError(
        "Usage: python reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py "
        "<decoder_variant> <n_simul> <snr_start> <snr_end> "
        "<selector_type: full_parallel|fixed_sequential|random_sequential|syndrome_sequential> "
        "<members_per_group> <target_fraction_converged>"
    )

decoder_variant = sys.argv[1].lower()
n_simul = int(sys.argv[2])
snr_start = float(sys.argv[3])
snr_end = float(sys.argv[4])
selector_type = sys.argv[5].lower()
members_per_group = int(sys.argv[6])
target_fraction_converged = float(sys.argv[7])

# Include the end point
sim_regime = np.arange(snr_start, snr_end + 0.25, 0.5)

remove = 4  # k=66 code rather than k=110
bg_vn = 12 - remove  # bg2
bg_cn = 4  # bg2
Zc = 11

results_dir = f"RESULTS/fig_x_zc11_r{remove}_seq_{selector_type}_mpg{members_per_group}_n{n_simul}"

splitting_pattern = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 4, 6, 8], [5]]


def create_bp_config(H_aux):
    """Create a BP configuration with the default decoder settings."""
    cfg = channel_code_lib2.BP_config(H_aux)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    return cfg


VARIANTS = {
    "asced22": (splitting_pattern[0], 1),
    "asced44": (splitting_pattern[0], 2),
    "asced88": (splitting_pattern[0], 4),
    "asced24": (splitting_pattern[1], 1),
    "asced48": (splitting_pattern[1], 2),
    "asced96_4batch": (splitting_pattern[1], 4),
    "asced192_8batch": (splitting_pattern[1], 8),
    "asced96": (splitting_pattern[2], 1),
    "asced192": (splitting_pattern[2], 2),
    "asced288": (splitting_pattern[2], 3),
    "asced384": (splitting_pattern[2], 4),
    "asced2048": ([], 1),
}

if decoder_variant not in VARIANTS:
    raise ValueError(f"Unknown decoder variant '{decoder_variant}'")

split_pattern, num_used_blocks = VARIANTS[decoder_variant]

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
use_all_zero = False

if use_all_zero:
    results_dir += "_AZ"

# Same PCM construction as reproduce_fig3_RL_zc11_asced_5G_LDPC.py.
code = np.load("Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz")
H_full = code["h"].astype(int)
H_full = np.delete(H_full, np.arange(6 * Zc, 10 * Zc), axis=1)

block_offset = 0
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


def binary_vectors_in_suffix(number_rows, number_last_entries):
    v = [0] * number_rows
    for mask in range(1, 1 << number_last_entries):
        for i in range(number_last_entries):
            v[number_rows - number_last_entries + i] = (mask >> i) & 1
        yield v.copy()


if not use_all_zero:
    enc_cfg = channel_code_lib2.PCM_Encoder_config(H, k, n)


def create_asced_config(H, candidate_rows, block_offset, Zc, split_pattern, num_used_blocks, use_all_zero):
    """Construct all decoder paths for one aSCED configuration.

    Returns the path configs grouped by "batch" (one row segment's zero-offset
    path plus all its affine-offset siblings), since that's the natural,
    physically-meaningful grouping: every path in a batch shares the same
    Tanner graph (same H_aux), differing only in the affine check-node offset.
    """
    batches = []

    for block in range(num_used_blocks):
        assert block + block_offset <= number_additional_cyclic_blocks
        row_block = candidate_rows[(block + block_offset) * Zc : (block + block_offset + 1) * Zc]

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


def run_asced(save_dir_name, split_pattern, num_used_blocks):
    ensemble_cfg, ensemble_size = create_asced_config(
        H=H,
        candidate_rows=candidate_rows,
        block_offset=block_offset,
        Zc=Zc,
        split_pattern=split_pattern,
        num_used_blocks=num_used_blocks,
        use_all_zero=use_all_zero,
    )

    if selector_type != "full_parallel":
        target_num_converged = max(1, round(target_fraction_converged * ensemble_size))
        ensemble_cfg.set_mConvergedConfig(target_num_converged)

        if selector_type == "fixed_sequential":
            ensemble_cfg.set_fixed_sequential(members_per_group)
        elif selector_type == "random_sequential":
            ensemble_cfg.set_random_sequential(members_per_group)
        elif selector_type == "syndrome_sequential":
            ensemble_cfg.set_syndrome_sequential(members_per_group)
        else:
            raise ValueError(f"Unknown selector_type '{selector_type}'")

    print("starting", save_dir_name)

    sim = channel_code_lib2.Simulation_Env(k, n, "all")
    # Bounded for the SLURM sweep: at high SNR, target_errors=1000 with the
    # original max_transmissions=3e8 can take weeks (the paper's own aSCED-384
    # data hits the same wall at 4.5 dB, collecting only 13 FEs). Capping
    # max_transmissions bounds per-task worst-case runtime; low-FER points will
    # simply collect fewer than target_errors FEs.
    sim.target_errors = 200
    sim.max_transmissions = int(2e6)
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus is not None:
        sim.num_threads = int(slurm_cpus)
    sim.auto_save = auto_save
    sim.save_dir = results_dir + "/" + save_dir_name

    sim.puncturing(message_bit_pucturing)

    if use_all_zero:
        sim.all_zero_init(ensemble_cfg)
    else:
        sim.init(enc_cfg, ensemble_cfg, use_all_zero)

    sim.get_error_rates(sim_regime)

    print(save_dir_name, "finished")
    print(sim.error_rates["FER-SNR"])

    return sim.error_rates["FER-SNR"]


save_name = f"{decoder_variant}_{selector_type}_mpg{members_per_group}_target{target_fraction_converged}"
FER_results = {save_name: run_asced(save_name, split_pattern, num_used_blocks)}
