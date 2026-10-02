# pylint: disable=invalid-name
"""Sweep BP_config.norm_factor (the MSA normalization constant) over a grid
of dyadic values, holding everything else at the established defaults
(scheduling_type="flooding", max_iterations=32, early_stopping=True,
cn_update_type="msa").

Reuses the exact H/PCM construction and create_bp_config()/create_asced_config()
patterns from reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py (RL/QC-block
5G LDPC C_5G(132,66), Zc=11), but only varies norm_factor - no selector,
stopping-policy, or scheduling changes (those are separate studies). Three
decoder variants:

  nmsa      - standalone NMSA: plain BP on the base PCM H (no ensemble).
  asced48   - aSCED-48 ensemble, full_parallel (splitting_pattern [2,4,6,8],
              2 blocks).
  asced384  - aSCED-384 ensemble, full_parallel (splitting_pattern [5],
              4 blocks).

Usage:
    python alpha_sweep.py <variant: nmsa|asced48|asced384> <n_simul> \
        <snr_start> <snr_end> <norm_factor>
"""

import numpy as np
import os
import sys

import galois

gf2 = galois.GF2

import channel_code_lib2

if len(sys.argv) != 6:
    raise ValueError(
        "Usage: python alpha_sweep.py <variant: nmsa|asced48|asced384> "
        "<n_simul> <snr_start> <snr_end> <norm_factor>"
    )

decoder_variant = sys.argv[1].lower()
n_simul = int(sys.argv[2])
snr_start = float(sys.argv[3])
snr_end = float(sys.argv[4])
norm_factor = float(sys.argv[5])

# Include the end point
sim_regime = np.arange(snr_start, snr_end + 0.25, 0.5)

remove = 4  # k=66 code rather than k=110
bg_vn = 12 - remove  # bg2
bg_cn = 4  # bg2
Zc = 11

# "alpha" in the path keeps this clearly distinguishable from the sequential/
# mconverged sweep's RESULTS directories.
results_dir = f"RESULTS/fig_x_zc11_r{remove}_alpha_sweep_n{n_simul}/{decoder_variant}/alpha_{norm_factor:g}"

splitting_pattern = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 4, 6, 8], [5]]


def create_bp_config(H_aux):
    """Create a BP configuration with the established defaults, varying only
    norm_factor."""
    cfg = channel_code_lib2.BP_config(H_aux)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = norm_factor
    return cfg


VARIANTS = {
    "asced48": (splitting_pattern[1], 2),
    "asced384": (splitting_pattern[2], 4),
}

if decoder_variant != "nmsa" and decoder_variant not in VARIANTS:
    raise ValueError(f"Unknown decoder variant '{decoder_variant}' (expected nmsa|asced48|asced384)")

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
use_all_zero = False

# Same PCM construction as reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py.
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

print(f"Simulating number vns={n}, n={number_vns - len(message_bit_pucturing)}, k={k}, norm_factor={norm_factor}")

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
    """Construct all decoder paths for one aSCED configuration (full_parallel
    only - no DecoderSelector/stopping policy applied here, see module
    docstring)."""
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


def make_sim(save_dir_name):
    sim = channel_code_lib2.Simulation_Env(k, n, "all")
    # Overridable via env vars only for quick local sanity checks; the real
    # sweep (sbatch array) never sets these, so it always uses the
    # established defaults below.
    sim.target_errors = int(os.environ.get("ALPHA_SWEEP_TARGET_ERRORS", 200))
    sim.max_transmissions = int(os.environ.get("ALPHA_SWEEP_MAX_TRANSMISSIONS", int(2e6)))
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus is not None:
        sim.num_threads = int(slurm_cpus)
    sim.auto_save = auto_save
    # One directory per (config, SNR) - avoids the concurrent-write race on
    # the C++ side's non-atomic JSON read-modify-write documented in
    # sweeps/STATUS.md (multiple SNR tasks of the same config launched in the
    # same SLURM array wave).
    snr_tag = f"{sim_regime[0]:g}"
    sim.save_dir = results_dir + "/" + save_dir_name + "/" + f"snr_{snr_tag}"
    sim.puncturing(message_bit_pucturing)
    return sim


def run_nmsa():
    nmsa_config = create_bp_config(H)
    sim = make_sim("nmsa")
    if use_all_zero:
        sim.all_zero_init(nmsa_config)
    else:
        sim.init(enc_cfg, nmsa_config, use_all_zero)
    sim.get_error_rates(sim_regime)
    print("nmsa finished")
    print(sim.error_rates["FER-SNR"])
    return sim.error_rates["FER-SNR"]


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

    print("starting", save_dir_name)

    sim = make_sim(save_dir_name)

    if use_all_zero:
        sim.all_zero_init(ensemble_cfg)
    else:
        sim.init(enc_cfg, ensemble_cfg, use_all_zero)

    sim.get_error_rates(sim_regime)

    print(save_dir_name, "finished")
    print(sim.error_rates["FER-SNR"])

    return sim.error_rates["FER-SNR"]


if decoder_variant == "nmsa":
    FER_results = {"nmsa": run_nmsa()}
else:
    split_pattern, num_used_blocks = VARIANTS[decoder_variant]
    save_name = f"{decoder_variant}_full_parallel"
    FER_results = {save_name: run_asced(save_name, split_pattern, num_used_blocks)}
