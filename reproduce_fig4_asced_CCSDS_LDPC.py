import numpy as np
from time import time
import os

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
n_ = 256
k_ = 128
Z = 32
H = np.array(read_AList("Codes/CCSDS/256_128.txt"), dtype=int)

use_all_zero = False  # SCED requires false or will result in to good performance since az covered by all paths

flag_aed = True  # if true simulates AED-11

flag_spa = True  # if true simulates spa-32

flag_asced_17 = True

flag_asced_31 = True  # if true simulate aSCED-11

plot_using_tex = True


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

if flag_spa:

    spa_config = channel_code_lib2.BP_config(H)
    spa_config.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    spa_config.max_iterations = 32  # set maximum number of BP iterations; default is 32
    spa_config.cn_update_type = "spa"
    spa_config.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding

    sim_spa = channel_code_lib2.Simulation_Env(H, k, n, "all")
    sim_spa.use_all_zero_codeword = use_all_zero
    sim_spa.init(spa_config)
    sim_spa.get_error_rates(np.linspace(1, 4, 7))

    FER_spa = sim_spa.error_rates["FER-SNR"]
    print(FER_spa)
    print("SPA 32 finished")


if flag_aed:
    permutation = quasi_cyclic_permutation_vector(n, Z)

    H_del = H[1:, :]
    shifted_permutations = [permutation]

    undercomplete_bp_config = [channel_code_lib2.BP_config(H_del)]

    undercomplete_bp_config[0].early_stopping = True
    undercomplete_bp_config[0].max_iterations = 32
    undercomplete_bp_config[0].cn_update_type = "spa"
    undercomplete_bp_config[0].scheduling_type = "flooding"
    for i in range(1, Z):
        shifted_permutations.append(shifted_permutations[i - 1][permutation])

    for per in shifted_permutations:
        assert np.all(gf2(H) @ G[:, per].T == 0)

    processing_config = channel_code_lib2.Automorphism_config(shifted_permutations)

    ensemble_decoder_config = channel_code_lib2.Ensemble_config(
        H, undercomplete_bp_config, processing_config
    )
    sim_aed = channel_code_lib2.Simulation_Env(H, k, n, "all")

    # cfg.H = H

    sim_aed.use_all_zero_codeword = use_all_zero
    sim_aed.init(ensemble_decoder_config)

    sim_aed.get_error_rates(np.linspace(1, 4, 7))
    FER_AED = sim_aed.error_rates["FER-SNR"]

    print(FER_AED)
    print("AED finished")


if flag_asced_17:
    # load pcms

    asced_17_path_configs = []
    asced_17_pcms = []
    for i in range(8):
        asced_17_pcms.append(
            np.load(f"Codes/TCOM_aSCED/CCSDS_affine_ensemble/CCSDS_path_{i}" + ".npy")
        )
        # path 0 is in fact the original pcm and has shape 88,154
        # all others generate two paths, one with affine_offset and one without
        asced_17_path_configs.append(channel_code_lib2.BP_config(asced_17_pcms[-1]))
        asced_17_path_configs[-1].early_stopping = True
        asced_17_path_configs[-1].max_iterations = 32
        asced_17_path_configs[-1].cn_update_type = "spa"
        asced_17_path_configs[-1].scheduling_type = "flooding"
        if i > 0 and not use_all_zero:
            # 2nd (i.e. affine) path of the batch
            affine_offset = np.zeros(asced_17_pcms[-1].shape[0], dtype=int)
            affine_offset[-1] = 1

            asced_17_path_configs.append(channel_code_lib2.BP_config(asced_17_pcms[-1]))
            asced_17_path_configs[-1].early_stopping = True
            asced_17_path_configs[-1].max_iterations = 32
            asced_17_path_configs[-1].cn_update_type = "spa"
            asced_17_path_configs[-1].scheduling_type = "flooding"
            asced_17_path_configs[-1].affine_offset = affine_offset
    print("Simulated num. aSCED paths:", len(asced_17_path_configs))
    asced_17_config = channel_code_lib2.Ensemble_config(H, asced_17_path_configs)
    sim_asced_17 = channel_code_lib2.Simulation_Env(H, k, n, "all")

    # cfg.H = H

    sim_asced_17.use_all_zero_codeword = use_all_zero

    sim_asced_17.init(asced_17_config)

    sim_asced_17.get_error_rates(np.linspace(1, 4, 7))
    FER_aSCED_17 = sim_asced_17.error_rates["FER-SNR"]
    print("aSCED finished")

    print(FER_aSCED_17)

if flag_asced_17:
    # load pcms

    asced_31_path_configs = []
    asced_31_pcms = []
    for i in range(15):
        asced_31_pcms.append(
            np.load(f"Codes/TCOM_aSCED/CCSDS_affine_ensemble/CCSDS_path_{i}" + ".npy")
        )
        # path 0 is in fact the original pcm and has shape 88,154
        # all others generate two paths, one with affine_offset and one without
        asced_31_path_configs.append(channel_code_lib2.BP_config(asced_31_pcms[-1]))
        asced_31_path_configs[-1].early_stopping = True
        asced_31_path_configs[-1].max_iterations = 32
        asced_31_path_configs[-1].cn_update_type = "spa"
        asced_31_path_configs[-1].scheduling_type = "flooding"
        if i > 0 and not use_all_zero:
            # 2nd (i.e. affine) path of the batch
            affine_offset = np.zeros(asced_31_pcms[-1].shape[0], dtype=int)
            affine_offset[-1] = 1

            asced_31_path_configs.append(channel_code_lib2.BP_config(asced_31_pcms[-1]))
            asced_31_path_configs[-1].early_stopping = True
            asced_31_path_configs[-1].max_iterations = 32
            asced_31_path_configs[-1].cn_update_type = "spa"
            asced_31_path_configs[-1].scheduling_type = "flooding"
            asced_31_path_configs[-1].affine_offset = affine_offset
    print("Simulated num. aSCED paths:", len(asced_31_path_configs))
    asced_31_config = channel_code_lib2.Ensemble_config(H, asced_31_path_configs)
    sim_asced_31 = channel_code_lib2.Simulation_Env(H, k, n, "all")

    # cfg.H = H

    sim_asced_31.use_all_zero_codeword = use_all_zero

    sim_asced_31.init(asced_31_config)

    sim_asced_31.get_error_rates(np.linspace(1, 4, 7))
    FER_aSCED_31 = sim_asced_31.error_rates["FER-SNR"]
    print("aSCED finished")

    print(FER_aSCED_31)


if plot_using_tex:
    show_results.plot_error_rates(
        (FER_aSCED_31, "aSCED-31"),
        (FER_aSCED_17, "aSCED-17"),
        (FER_AED, "AED"),
        (FER_spa, "spa 32"),
    )
