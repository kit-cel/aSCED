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

flag_aed=True # if true simulates AED-11

flag_nmsa=True # if true simulates NMSA-32

flag_nmsa352=True #if true simulates NMSA-352


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

if(flag_nmsa):

    nmsa_config=channel_code_lib2.BP_config(H)
    nmsa_config.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    nmsa_config.max_iterations = 32  # set maximum number of BP iterations; default is 32
    nmsa_config.cn_update_type = "msa"
    nmsa_config.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding
    nmsa_config.norm_factor=0.75

    sim_nmsa = channel_code_lib2.Simulation_Env(H, k, n, "all")
    sim_nmsa.use_all_zero_codeword = True
    sim_nmsa.puncturing(p)
    sim_nmsa.shortening(s)
    sim_nmsa.init(nmsa_config)
    sim_nmsa.get_error_rates(np.linspace(1, 4, 7))

    FER_nmsa = sim_nmsa.error_rates["FER-SNR"]
    print(FER_nmsa)
    print("NMSA 32 finished")  

if(flag_nmsa352):

    nmsa_config_352=channel_code_lib2.BP_config(H)
    nmsa_config_352.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    nmsa_config_352.max_iterations = 352  # set maximum number of BP iterations; default is 32
    nmsa_config_352.cn_update_type = "msa" 
    nmsa_config_352.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding
    nmsa_config_352.norm_factor=0.75

    sim_nmsa_352 = channel_code_lib2.Simulation_Env(H, k, n, "all")
    sim_nmsa_352.use_all_zero_codeword = True
    sim_nmsa_352.puncturing(p)
    sim_nmsa_352.shortening(s)
    sim_nmsa_352.init(nmsa_config)
    sim_nmsa_352.get_error_rates(np.linspace(1, 4, 7))

    FER_nmsa_352 = sim_nmsa_352.error_rates["FER-SNR"]
    print(FER_nmsa_352)
    print("NMSA 352 finished")



if(flag_aed):
    permutation = quasi_cyclic_permutation_vector(n, Z)

    H_del = H[1:, :]
    shifted_permutations = [permutation]

    undercomplete_bp_config = [channel_code_lib2.BP_config(H_del)]

    undercomplete_bp_config[0]=channel_code_lib2.BP_config(H)
    undercomplete_bp_config[0].early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    undercomplete_bp_config[0].max_iterations = 32  # set maximum number of BP iterations; default is 32
    undercomplete_bp_config[0].cn_update_type = "msa"  # Check node update rule (msa, spa, spa_phi); default is spa)
    undercomplete_bp_config[0].scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding
    undercomplete_bp_config[0].norm_factor=0.75
    for i in range(1, 11):
        shifted_permutations.append(shifted_permutations[i - 1][permutation])


    for per in shifted_permutations:
        assert np.all(gf2(H) @ G[:, per].T == 0)

    processing_config = channel_code_lib2.Automorphism_config(shifted_permutations)



    ensemble_decoder_config = channel_code_lib2.Ensemble_config(
        H, undercomplete_bp_config, processing_config
    )
    sim_aed = channel_code_lib2.Simulation_Env(H, k, n, "all")


    # cfg.H = H

    sim_aed.use_all_zero_codeword = True
    sim_aed.puncturing(p)
    sim_aed.shortening(s)
    sim_aed.init(ensemble_decoder_config)

    sim_aed.get_error_rates(np.linspace(1, 3.6, 7))
    FER_AED = sim_aed.error_rates["FER-SNR"]


    print(FER_AED)
show_results.plot_error_rates((FER_AED, "AED"), (FER_nmsa, "NMSA 32"), (FER_nmsa_352, "NMSA 352"))
