import numpy as np
from time import time
import os

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

flag_aed = True  # if true simulates AED-11

flag_nmsa = True  # if true simulates NMSA-32

flag_nmsa352 = True  # if true simulates NMSA-352

flag_sced = True

flag_asced = True  # if true simulate aSCED-11

plot_using_tex = False


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



if flag_nmsa:

    nmsa_config = channel_code_lib2.BP_config(H)
    nmsa_config.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    nmsa_config.max_iterations = (
        32  # set maximum number of BP iterations; default is 32
    )
    nmsa_config.cn_update_type = "msa"
    nmsa_config.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding
    nmsa_config.norm_factor = 0.75

    sim_nmsa = channel_code_lib2.Simulation_Env(H, k, n, "all")
    sim_nmsa.puncturing(message_bit_pucturing)
    if not use_all_zero:
        sim_nmsa.use_all_zero_codeword = use_all_zero
        sim_nmsa.init(enc_cfg, nmsa_config)
    else:
        sim_nmsa.all_zero_init(nmsa_config)
    sim_nmsa.get_error_rates(np.linspace(1, 4, 7))

    FER_nmsa = sim_nmsa.error_rates["FER-SNR"]
    print(FER_nmsa)
    print("NMSA 32 finished")

if flag_nmsa352:

    nmsa_config_352 = channel_code_lib2.BP_config(H)
    nmsa_config_352.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    nmsa_config_352.max_iterations = (
        352  # set maximum number of BP iterations; default is 32
    )
    nmsa_config_352.cn_update_type = "msa"
    nmsa_config_352.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding
    nmsa_config_352.norm_factor = 0.75

    sim_nmsa_352 = channel_code_lib2.Simulation_Env(H, k, n, "all")
    sim_nmsa_352.puncturing(message_bit_pucturing)
    if not use_all_zero:
        sim_nmsa_352.use_all_zero_codeword = use_all_zero
        sim_nmsa_352.init(enc_cfg, nmsa_config_352)
    else:
        sim_nmsa_352.all_zero_init(nmsa_config_352) 
    
    sim_nmsa_352.get_error_rates(np.linspace(1, 4, 7))

    FER_nmsa_352 = sim_nmsa_352.error_rates["FER-SNR"]
    print(FER_nmsa_352)
    print("NMSA 352 finished")


if flag_aed:
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
        assert np.all(gf2(H) @ G[:, per].T == 0)

    processing_config = channel_code_lib2.Automorphism_config(shifted_permutations)

    ensemble_decoder_config = channel_code_lib2.Ensemble_config(
        H, undercomplete_bp_config, processing_config
    )
    sim_aed = channel_code_lib2.Simulation_Env(H, k, n, "all")

    # cfg.H = H


    sim_aed.puncturing(message_bit_pucturing)
    if not use_all_zero:
        sim_aed.use_all_zero_codeword = use_all_zero
        sim_aed.init(enc_cfg, ensemble_decoder_config)
    else:
        sim_aed.all_zero_init(ensemble_decoder_config) 

    sim_aed.get_error_rates(np.linspace(1, 4, 7))
    FER_AED = sim_aed.error_rates["FER-SNR"]

    print(FER_AED)
    print("AED finished")


if flag_sced:
    # load pcms

    sced_path_configs = [channel_code_lib2.BP_config(H)]
    sced_pcms = []
    for i in range(10):
        sced_pcms.append(
            np.load(f"Codes/TCOM_aSCED/5G_sced_ensemble/5G_aux_sced_path_{i}" + ".npy")
        )
        sced_path_configs.append(channel_code_lib2.BP_config(sced_pcms[-1]))

        # path 0 is in fact the original pcm and has shape 88,154
        # all others generate two paths, one with affine_offset and one without

    for cfg in sced_path_configs:
        cfg.early_stopping = True
        cfg.max_iterations = 32
        cfg.cn_update_type = "msa"
        cfg.scheduling_type = "flooding"
        cfg.norm_factor = 0.75

    print("Simulated num. sced paths:", len(sced_path_configs))
    sced_config = channel_code_lib2.Ensemble_config(H, sced_path_configs)
    sim_sced = channel_code_lib2.Simulation_Env(H, k, n, "all")

    # cfg.H = H

    sim_sced.puncturing(message_bit_pucturing)
    if not use_all_zero:
        sim_sced.use_all_zero_codeword = use_all_zero
        sim_sced.init(enc_cfg, sced_config)
    else:
        sim_sced.all_zero_init(sced_config) 
    sim_sced.get_error_rates(np.linspace(1, 4, 7))
    FER_sced = sim_sced.error_rates["FER-SNR"]
    print("SCED finished")

    print(FER_sced)

if flag_asced:
    # load pcms

    asced_path_configs = []
    asced_pcms = []
    for i in range(6):
        asced_pcms.append(
            np.load(f"Codes/TCOM_aSCED/5G_affine_ensemble/5G_path_{i}" + ".npy")
        )
        # path 0 is in fact the original pcm and has shape 88,154
        # all others generate two paths, one with affine_offset and one without
        asced_path_configs.append(channel_code_lib2.BP_config(asced_pcms[-1]))
        asced_path_configs[-1].early_stopping = True
        asced_path_configs[-1].max_iterations = 32
        asced_path_configs[-1].cn_update_type = "msa"
        asced_path_configs[-1].scheduling_type = "flooding"
        asced_path_configs[-1].norm_factor = 0.75
        if i > 0 and not use_all_zero:
            # 2nd (i.e. affine) path of the batch
            affine_offset = np.zeros(asced_pcms[-1].shape[0], dtype=int)
            affine_offset[-1] = 1

            asced_path_configs.append(channel_code_lib2.BP_config(asced_pcms[-1]))
            asced_path_configs[-1].early_stopping = True
            asced_path_configs[-1].max_iterations = 32
            asced_path_configs[-1].cn_update_type = "msa"
            asced_path_configs[-1].scheduling_type = "flooding"
            asced_path_configs[-1].norm_factor = 0.75
            asced_path_configs[-1].affine_offset = affine_offset
    print("Simulated num. aSCED paths:", len(asced_path_configs))
    asced_config = channel_code_lib2.Ensemble_config(H, asced_path_configs)
    sim_asced = channel_code_lib2.Simulation_Env(H, k, n, "all")

    # cfg.H = H

    sim_asced.puncturing(message_bit_pucturing)
    if not use_all_zero:
        sim_asced.use_all_zero_codeword = use_all_zero
        sim_asced.init(enc_cfg, asced_config)
    else:
        sim_asced.all_zero_init(asced_config) 
    sim_asced.get_error_rates(np.linspace(1, 4, 7))
    FER_aSCED = sim_asced.error_rates["FER-SNR"]
    print("aSCED finished")

    print(FER_aSCED)


if plot_using_tex:
    show_results.plot_error_rates(
        (FER_aSCED, "aSCED"),
        (FER_sced, "SCED"),
        (FER_AED, "AED"),
        (FER_nmsa, "NMSA 32"),
        (FER_nmsa_352, "NMSA 352"),
    )
