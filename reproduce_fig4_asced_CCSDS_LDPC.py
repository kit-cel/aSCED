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
n, k, H = read_AList("Codes/CCSDS/256_128.txt")

G = gf2(H).null_space()

k, n = G.shape

print(n, k)

use_all_zero = False  # Currently only all-zero since bug in encode of ccsds 256,128


if not use_all_zero:
    g_enc_cfg = channel_code_lib2.G_Encoder_config(G, k, n)

flag_aed = True  # if true simulates AED-11

flag_spa = True  # if true simulates spa-32

flag_asced_17 = True

flag_asced_31 = True  # if true simulate aSCED-11

plot_using_tex = True
# np.linspace(1, 4,7 )

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

print(k, n)

if flag_spa:
    spa_config = channel_code_lib2.BP_config(H)
    spa_config.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    spa_config.max_iterations = 32  # set maximum number of BP iterations; default is 32
    spa_config.cn_update_type = "spa"
    spa_config.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding

    sim_spa = channel_code_lib2.Simulation_Env(k, n, "all")

    if not use_all_zero:
        sim_spa.init(g_enc_cfg, spa_config, use_all_zero)
    else:
        sim_spa.all_zero_init(spa_config)

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
    sim_aed = channel_code_lib2.Simulation_Env(k, n, "all")

    # cfg.H = H

    if not use_all_zero:
        sim_aed.init(g_enc_cfg, ensemble_decoder_config, use_all_zero)
    else:
        sim_aed.all_zero_init(ensemble_decoder_config)

    sim_aed.get_error_rates(np.linspace(1, 4, 7))
    FER_AED = sim_aed.error_rates["FER-SNR"]

    print(FER_AED)
    print("AED finished")


if flag_asced_17:
    # load pcms

    asced_17_path_configs = []
    asced_17_pcms = []
    for i in range(9):
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
    sim_asced_17 = channel_code_lib2.Simulation_Env(k, n, "all")

    # cfg.H = H

    if not use_all_zero:
        sim_asced_17.init(g_enc_cfg, asced_17_config, use_all_zero)
    else:
        sim_asced_17.all_zero_init(asced_17_config)

    sim_asced_17.get_error_rates(np.linspace(1, 4, 7))
    FER_aSCED_17 = sim_asced_17.error_rates["FER-SNR"]
    print("aSCED finished")

    print(FER_aSCED_17)

if flag_asced_31:
    # load pcms

    asced_31_path_configs = []
    asced_31_pcms = []
    for i in range(16):
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
    sim_asced_31 = channel_code_lib2.Simulation_Env(k, n, "all")

    # cfg.H = H

    if not use_all_zero:
        sim_asced_31.init(g_enc_cfg, asced_31_config, use_all_zero)
    else:
        sim_asced_31.all_zero_init(asced_31_config)
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


# 256 128
# 128 256
# Using all-zero codeword assumption
# {1.0: 0.8366606170598911, 1.5: 0.584354628422425, 2.0: 0.281539080305671, 2.5: 0.08950442811381194, 3.0: 0.02025435704192181, 3.5: 0.003055046381158696, 4.0: 0.0002910863849150182}
# SPA 32 finished
# Using all-zero codeword assumption
# {1.0: 0.8027976799727056, 1.5: 0.500794551645857, 2.0: 0.22726751519321178, 2.5: 0.06373572820941242, 3.0: 0.009498272795584411, 3.5: 0.00083249924723978, 4.0: 3.823989919435125e-05}
# AED finished
# Simulated num. aSCED paths: 8
# Using all-zero codeword assumption
# aSCED finished
# {1.0: 0.7831164497831165, 1.5: 0.4828653548529737, 2.0: 0.20079443892750745, 2.5: 0.054939400868968674, 3.0: 0.008216598113695468, 3.5: 0.000827142822315791, 4.0: 4.632158102043766e-05}
# Simulated num. aSCED paths: 15
# Using all-zero codeword assumption
# aSCED finished
# {1.0: 0.7730040595399188, 1.5: 0.4
