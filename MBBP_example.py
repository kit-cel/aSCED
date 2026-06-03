import numpy as np
from time import time

import galois

gf2 = galois.GF2

import channel_code_lib2

from Codes.generate_5G_LDPC import generate_5G_LDPC
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


G = gf2(H).null_space()
k, n = G.shape

print(quasi_cyclic_permutation_vector(n, Z))

permutation = quasi_cyclic_permutation_vector(n, Z)

H_del = H[1:, :]
shifted_PCMs = [H_del]

configs = [channel_code_lib2.BP_config(H_del)]

for i in range(1, 11):
    shifted_PCMs.append(shifted_PCMs[i - 1][:, permutation])
    configs.append(channel_code_lib2.BP_config(shifted_PCMs[i - 1]))


processing_config = channel_code_lib2.Identity_config()

##PCM provided to ensemble config used for ML in the list
# ensemble_decoder_config = channel_code_lib2.Ensemble_config(
#     H, configs, processing_config
# )

ensemble_decoder_config = channel_code_lib2.Ensemble_config(H, configs)
# identical to

# ensemble_decoder_config = channel_code_lib2.Ensemble_config(
#     H, configs
# ) since Identity_config is default
sim = channel_code_lib2.Simulation_Env(H, k, n, "all")

# cfg.H = H

sim.use_all_zero_codeword = True
sim.puncturing(p)
sim.shortening(s)
# sim.Z = Z
# sim.set_ensemble_decoding('SED', 8)
sim.init(ensemble_decoder_config)

sim.get_error_rates(np.linspace(1, 3, 7))

FER = sim.error_rates["FER-SNR"]

print(FER)


bp_config = channel_code_lib2.BP_config(H)

sim_bp = channel_code_lib2.Simulation_Env(H, k, n, "all")

sim_bp.use_all_zero_codeword = True
sim_bp.puncturing(p)
sim_bp.shortening(s)
sim_bp.init(bp_config)
sim_bp.get_error_rates(np.linspace(1, 3, 7))


FER_bp = sim_bp.error_rates["FER-SNR"]

show_results.plot_error_rates((FER, "MBBP"), (FER_bp, "BP"))
