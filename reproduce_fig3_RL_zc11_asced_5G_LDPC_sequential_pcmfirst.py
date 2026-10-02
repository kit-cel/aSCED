# pylint: disable=invalid-name
""""PCM-first" variant of reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py.

Adds two new ensemble variants, aSCED-49 and aSCED-385, which are exactly
the existing aSCED-48 / aSCED-384 ensembles plus one extra decoding path on
the plain original PCM `H` (no appended auxiliary rows, hence a single path
rather than a batch of affine-offset siblings -- there is no rank deficiency
to offset against). That extra path is prepended at path index 0, i.e. before
all the aSCED batches.

Why index 0: `FixedSequentialSelector` (see channel-code-lib2's
`FixedSequentialSelector.cpp`) partitions the ensemble into *contiguous*
index ranges, in order: group i covers indices
[i*members_per_group, (i+1)*members_per_group). Putting the PCM-only path at
index 0 guarantees it lands in the very first group under `fixed_sequential`,
regardless of members_per_group. The hypothesis: the plain-PCM path often
converges fast and cleanly on its own, so scheduling it first should let the
m-converged stopping policy trigger earlier, reducing both
average_ensemble_effort and average_ensemble_latency relative to plain
aSCED-48/384 -- and since Ensemble_Decoder::ML_criterion() just considers
one more candidate codeword, FER should be <= the plain ensemble's at every
SNR (never worse).

This is deliberately scoped to `fixed_sequential` only (verified by reading
the C++ selector source, not assumed) -- `syndrome_sequential` reorders
dynamically per received word, so "pin path 0 to the first group" would
require a real C++ "pin" feature, out of scope here. `full_parallel` is
included too since it ignores ordering/selectors entirely and is used as
the correctness/FER baseline.

members_per_group=7 divides both new ensemble sizes evenly: 49 = 7x7,
385 = 7x55 (384+1).
"""

import numpy as np
import os
import sys

import galois

gf2 = galois.GF2

import channel_code_lib2

if len(sys.argv) != 8:
    raise ValueError(
        "Usage: python reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential_pcmfirst.py "
        "<decoder_variant> <n_simul> <snr_start> <snr_end> "
        "<selector_type: full_parallel|fixed_sequential|random_sequential|syndrome_sequential> "
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

# "pcmfirst" in the path keeps this clearly distinguishable from the plain
# sequential script's RESULTS/fig_x_zc11_r4_seq_<selector>_mpg<N>_n<N>/ tree.
results_dir = f"RESULTS/fig_x_zc11_r{remove}_seq_pcmfirst_{selector_type}_mpg{members_per_group}_n{n_simul}"

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


# (split_pattern, num_used_blocks, prepend_pcm_path)
# prepend_pcm_path=True adds one extra create_bp_config(H) path (plain PCM,
# no auxiliary rows) at index 0, before all the aSCED batches.
VARIANTS = {
    "asced22": (splitting_pattern[0], 1, False),
    "asced44": (splitting_pattern[0], 2, False),
    "asced88": (splitting_pattern[0], 4, False),
    "asced24": (splitting_pattern[1], 1, False),
    "asced48": (splitting_pattern[1], 2, False),
    "asced96_4batch": (splitting_pattern[1], 4, False),
    "asced192_8batch": (splitting_pattern[1], 8, False),
    "asced96": (splitting_pattern[2], 1, False),
    "asced192": (splitting_pattern[2], 2, False),
    "asced288": (splitting_pattern[2], 3, False),
    "asced384": (splitting_pattern[2], 4, False),
    "asced2048": ([], 1, False),
    # New "PCM-first" variants: aSCED-48/384 + one plain-PCM path at index 0.
    "asced49": (splitting_pattern[1], 2, True),
    "asced385": (splitting_pattern[2], 4, True),
}

if decoder_variant not in VARIANTS:
    raise ValueError(f"Unknown decoder variant '{decoder_variant}'")

split_pattern, num_used_blocks, prepend_pcm_path = VARIANTS[decoder_variant]

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
use_all_zero = False

if use_all_zero:
    results_dir += "_AZ"

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


def create_asced_config(
    H, candidate_rows, block_offset, Zc, split_pattern, num_used_blocks, use_all_zero, prepend_pcm_path=False
):
    """Construct all decoder paths for one aSCED configuration.

    Returns the path configs grouped by "batch" (one row segment's zero-offset
    path plus all its affine-offset siblings), since that's the natural,
    physically-meaningful grouping: every path in a batch shares the same
    Tanner graph (same H_aux), differing only in the affine check-node offset.

    If prepend_pcm_path, a single extra path decoding the plain original PCM
    `H` (no appended auxiliary rows, hence no affine-offset siblings -- there
    is no rank deficiency to offset against) is prepended as its own batch,
    at index 0 of the resulting path list -- so it is the very first path of
    the very first group under `fixed_sequential`, whatever
    `members_per_group` is.
    """
    batches = []

    if prepend_pcm_path:
        batches.append([create_bp_config(H)])

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
        prepend_pcm_path=prepend_pcm_path,
    )

    if selector_type != "full_parallel":
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
    # One directory per (config, SNR): sim_regime here is always a single SNR
    # point (snr_start==snr_end), and these JSON outputs are merged by the
    # C++ side via a non-atomic read-modify-write - concurrent SLURM array
    # tasks sharing a save_dir (e.g. all SNR points of one config launched in
    # the same wave) can race and silently drop entries. A unique directory
    # per run sidesteps that entirely and is better provenance anyway.
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
FER_results = {save_name: run_asced(save_name, split_pattern, num_used_blocks)}
