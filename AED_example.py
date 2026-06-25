import numpy as np
from time import time

import galois

gf2 = galois.GF2

import channel_code_lib2

from Codes.generate_5G_LDPC import (
    generate_5G_LDPC,
    get_final_matrices_and_message_bit_pucturing,
)
from Codes.generate_RM import generate_RM
from Codes.overcomplete import overcomplete
from Codes.read_AList import read_AList
import matplotlib.pyplot as plt

import show_results

# import data


# - Generate 5G LDPC code
n_ = 132
k_ = 66
H, p, s, Z, BG = generate_5G_LDPC(2, k_, n_, return_lifting_size=True)

H, G, k, n, message_bit_pucturing = get_final_matrices_and_message_bit_pucturing(
    H, s, p
)

use_all_zero = True

if not use_all_zero:
    enc_cfg = channel_code_lib2.PCM_Encoder_config(H, k, n)


## First setup interprets AED as MBBP instanciated with shifted parity-check matrices obtained by cyclically permuting the columns of the original parity-check matrix.
## Should yield the same performance as AED using same permutations
## Thereby nice check if both work


def quasi_cyclic_permutation_vector(length, block_size=11):
    permuted_indices = np.arange(length)
    for start in range(0, length, block_size):
        end = min(start + block_size, length)
        block_indices = permuted_indices[start:end]
        if len(block_indices) == block_size:
            permuted_indices[start:end] = np.roll(block_indices, 1)
    return permuted_indices


print(quasi_cyclic_permutation_vector(n, Z))

permutation = quasi_cyclic_permutation_vector(n, Z)

H_del = H[1:, :]
print(H_del.shape)
shifted_permutations = [permutation]

undercomplete_bp_config = [channel_code_lib2.BP_config(H_del)]

for i in range(1, 11):
    shifted_permutations.append(shifted_permutations[i - 1][permutation])
    print(shifted_permutations[i])


for per in shifted_permutations:
    assert np.all(gf2(H) @ G[:, per].T == 0)

processing_config = channel_code_lib2.Automorphism_config(shifted_permutations)

##PCM provided to ensemble config used for ML in the list
# ensemble_decoder_config = channel_code_lib2.Ensemble_config(
#     H, configs, processing_config
# )

ensemble_decoder_config = channel_code_lib2.Ensemble_config(
    H, undercomplete_bp_config, processing_config
)

# identical to

# ensemble_decoder_config = channel_code_lib2.Ensemble_config(
#     H, configs
# ) since Identity_config is default
sim = channel_code_lib2.Simulation_Env( k, n, "all")


# cfg.H = H

sim.puncturing(message_bit_pucturing)
if not use_all_zero:
    sim.init(enc_cfg, ensemble_decoder_config, use_all_zero)
else:
    sim.all_zero_init(ensemble_decoder_config)

sim.get_error_rates(np.linspace(1, 3, 7))
FER_AED = sim.error_rates["FER-SNR"]

bp_config = channel_code_lib2.BP_config(H)

sim_bp = channel_code_lib2.Simulation_Env( k, n, "all")

sim_bp.puncturing(message_bit_pucturing)
if not use_all_zero:
    sim_bp.init(enc_cfg, bp_config, use_all_zero)
else:
    sim_bp.all_zero_init(bp_config)

sim_bp.get_error_rates(np.linspace(1, 3, 7))
FER_bp = sim_bp.error_rates["FER-SNR"]


print(FER_AED)
show_results.plot_error_rates((FER_AED, "AED"), (FER_bp, "BP"))
