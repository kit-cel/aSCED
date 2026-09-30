# pylint: disable=invalid-name
"""Small smoke test for the sequential DecoderSelector + mconverged stopping
policy path (AED with random/fixed sequential grouping), sized to run in
seconds on 4 threads. Not a numerical accuracy check - just exercises the
code paths under valgrind / repeated runs to catch crashes and leaks."""

import numpy as np
import galois

gf2 = galois.GF2

import channel_code_lib2

from Codes.generate_5G_LDPC import (
    generate_5G_LDPC,
    get_final_matrices_and_message_bit_pucturing,
)

n_ = 132
k_ = 66
H, p, s, Z, BG = generate_5G_LDPC(2, k_, n_)
H, G, k, n, message_bit_pucturing = get_final_matrices_and_message_bit_pucturing(H, s, p)
G = gf2(G)

use_all_zero = False
enc_cfg = channel_code_lib2.PCM_Encoder_config(H, k, n)

# tiny run: 1 SNR point, few target errors, hard cap on transmissions
sim_regime = np.array([1.0])
target_errors = 5
max_transmissions = 200
num_threads = 4
target_fraction_converged_path = 0.25
members_per_group = 2


def quasi_cyclic_permutation_vector(length, block_size=11):
    permuted_indices = np.arange(length)
    for start in range(0, length, block_size):
        end = min(start + block_size, length)
        block_indices = permuted_indices[start:end]
        if len(block_indices) == block_size:
            permuted_indices[start:end] = np.roll(block_indices, 1)
    return permuted_indices


def run_aed(random_sequential: bool):
    permutation = quasi_cyclic_permutation_vector(n, Z)
    H_del = H[1:, :]
    shifted_permutations = [permutation]

    undercomplete_bp_config = [channel_code_lib2.BP_config(H_del)]
    undercomplete_bp_config[0].early_stopping = True
    undercomplete_bp_config[0].max_iterations = 32
    undercomplete_bp_config[0].cn_update_type = "msa"
    undercomplete_bp_config[0].scheduling_type = "flooding"
    undercomplete_bp_config[0].norm_factor = 0.75
    for i in range(1, 11):
        shifted_permutations.append(shifted_permutations[i - 1][permutation])

    for per in shifted_permutations:
        assert np.all(gf2(H) @ gf2(G[:, per]).T == 0)

    processing_config = channel_code_lib2.Automorphism_config(shifted_permutations)
    ensemble_decoder_config = channel_code_lib2.Ensemble_config(
        H, undercomplete_bp_config, processing_config
    )
    ensemble_decoder_config.set_mConvergedConfig(
        int(np.ceil(target_fraction_converged_path * 11))
    )

    if random_sequential:
        ensemble_decoder_config.set_random_sequential(members_per_group)
    else:
        ensemble_decoder_config.set_fixed_sequential(members_per_group)

    sim = channel_code_lib2.Simulation_Env(k, n, "all")
    sim.target_errors = target_errors
    sim.max_transmissions = max_transmissions
    sim.num_threads = num_threads
    sim.auto_save = False

    sim.puncturing(message_bit_pucturing)
    sim.init(enc_cfg, ensemble_decoder_config, use_all_zero)
    sim.get_error_rates(sim_regime)
    label = "random_sequential" if random_sequential else "fixed_sequential"
    print(f"[{label}] FER-SNR:", sim.error_rates["FER-SNR"])


if __name__ == "__main__":
    run_aed(random_sequential=True)
    run_aed(random_sequential=False)
    print("smoke test finished OK")
