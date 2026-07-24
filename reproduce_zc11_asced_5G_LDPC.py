# pylint: disable=invalid-name
"""Reproduce Figure x results for ASCED on 5G LDPC codes."""

import numpy as np
from time import time
import os

import galois

gf2 = galois.GF2

import channel_code_lib2


import matplotlib.pyplot as plt

results_dir = "RESULTS/fig_x_zc11"

bg_vn = 12  ##bg2
bg_cn = 4  ##bg2
Zc = 11

sim_regime = np.linspace(2, 3.5, 4)

n_simul = 143  # or increase in stepzsizes of 11 e.g. 165

splitting_pattern = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [5], [2, 4, 6, 8]]


flag_nmsa = True  # if true simulates NMSA-32

flag_aed = True  # if true simulates AED-11

##splitting_pattern [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
flag_aed_asced_22 = True  # 1 batch of 11*2**1
flag_aed_asced_44 = True  # 2 batches of 11*2**1
flag_aed_asced_88 = True  # 4 batches of 11*2**1

##splitting pattern [2, 4, 6, 8]
flag_asced_24 = True  # 2**2*4+2**3
flag_asced_48 = True  # 2 batches of 2**2*4+2**3
flag_asced_96_4batch = True  # 4 batches of 2**2*4+2**3

##splitting_pattern [5]
flag_asced_96 = True  # 2**5+2**6
flag_asced_192 = True  # 2batches of 2**5+2**6

flag_asced_2048 = True  #  2**11

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
bool_emulate_stopping = False
target_fraction_coverged_path = 0.5


use_all_zero = False


# we start first with code max_rank_5G_zc=11
code = np.load("Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz")
# this is the full BG2 PCM lifted with Zc=11
H_full = code["h"].astype(int)


block_offset = 0  # first batch, uses block 0, 2nd batch, uses block and so on


number_vn_simul = n_simul + 2 * Zc

number_vns_start = (
    bg_vn * 11
    + 2 * Zc  # is +11 not standardized? #+22 to have highest rate one from standard
)  # this here is only starting n; the code we transmit will be extendend by harq_part*Zc bits,


harq_part = (number_vn_simul - number_vns_start) // Zc

print(harq_part)

m = Zc * bg_cn + harq_part * Zc

number_vns = number_vns_start + harq_part * Zc

assert number_vn_simul == number_vns

H = gf2(H_full[:m, :number_vns])  # take until current m

G = np.array(gf2(H[:m, :number_vns]).null_space()).astype(int)

k, n = G.shape

assert n == number_vns

print(f"Simulating number vns={n}, n={number_vns-len(message_bit_pucturing)},k={k} ")

candidate_rows = H_full[m:, : number_vns_start + harq_part * Zc]  # future possible rows

number_additional_cyclic_blocks = candidate_rows.shape[0] // Zc

assert number_additional_cyclic_blocks * Zc == candidate_rows.shape[0]


def binary_vectors_in_suffix(number_rows, number_last_entries):
    v = [0] * number_rows

    for mask in range(1, 1 << number_last_entries):
        # overwrite the last m entries
        for i in range(number_last_entries):
            v[number_rows - number_last_entries + i] = (mask >> i) & 1
        yield v.copy()


if not use_all_zero:
    enc_cfg = channel_code_lib2.PCM_Encoder_config(H, k, n)


if flag_nmsa:

    nmsa_config = channel_code_lib2.BP_config(H)
    nmsa_config.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    nmsa_config.max_iterations = (
        32  # set maximum number of BP iterations; default is 32
    )
    nmsa_config.cn_update_type = "msa"
    nmsa_config.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding
    nmsa_config.norm_factor = 0.75

    sim_nmsa = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_nmsa.auto_save = auto_save
    sim_nmsa.save_dir = results_dir + "/nmsa"
    sim_nmsa.puncturing(message_bit_pucturing)
    if not use_all_zero:
        sim_nmsa.init(enc_cfg, nmsa_config, use_all_zero)
    else:
        sim_nmsa.all_zero_init(nmsa_config)
    sim_nmsa.get_error_rates(sim_regime)

    FER_nmsa = sim_nmsa.error_rates["FER-SNR"]
    print(FER_nmsa)
    print("NMSA 32 finished")


if flag_aed_asced_22:
    pattern = splitting_pattern[0]

    num_used_blocks = 1

    ##TODO: append all used blocks once to H and assert that increased rank == number appended rows

    # construct pcms
    asced_22_path_configs = []
    asced_22_pcms = []
    for i in range(num_used_blocks):
        assert i + block_offset <= number_additional_cyclic_blocks
        row_block = candidate_rows[
            (i + block_offset) * Zc : (i + block_offset + 1) * Zc
        ]
        # load pcms
        row_segments = np.split(row_block, pattern)
        print(len(row_segments))

        for rows in row_segments:
            ##due to previous assert,
            rank_difference = rows.shape[0]
            H_aux = gf2(np.vstack((H, rows)))
            m_aux, _ = H_aux.shape
            asced_22_path_configs.append(channel_code_lib2.BP_config(H_aux))
            asced_22_path_configs[-1].early_stopping = True
            asced_22_path_configs[-1].max_iterations = 32
            asced_22_path_configs[-1].cn_update_type = "msa"
            asced_22_path_configs[-1].scheduling_type = "flooding"
            asced_22_path_configs[-1].norm_factor = 0.75
            if not use_all_zero:
                ##TODO offset takes on all binary vectors of length rank_difference
                all_affine_offsets = binary_vectors_in_suffix(m_aux, rank_difference)

                for offset in all_affine_offsets:
                    asced_22_path_configs.append(channel_code_lib2.BP_config(H_aux))
                    asced_22_path_configs[-1].early_stopping = True
                    asced_22_path_configs[-1].max_iterations = 32
                    asced_22_path_configs[-1].cn_update_type = "msa"
                    asced_22_path_configs[-1].scheduling_type = "flooding"
                    asced_22_path_configs[-1].norm_factor = 0.75
                    asced_22_path_configs[-1].affine_offset = offset

    print("Simulated num. aSCED paths:", len(asced_22_path_configs))
    asced_22_config = channel_code_lib2.Ensemble_config(H, asced_22_path_configs)

    sim_22_asced = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_22_asced.auto_save = auto_save
    sim_22_asced.save_dir = results_dir + "/aSCED_22_split0_batch1"
    # cfg.H = H

    sim_22_asced.puncturing(message_bit_pucturing)
    if not use_all_zero:
        sim_22_asced.init(enc_cfg, asced_22_config, use_all_zero)
    else:
        sim_22_asced.all_zero_init(asced_22_config)
    sim_22_asced.get_error_rates(sim_regime)
    FER_aSCED_22 = sim_22_asced.error_rates["FER-SNR"]
    print("aSCED finished")

    print(FER_aSCED_22)


if flag_aed_asced_44:
    pattern = splitting_pattern[0]
    num_used_blocks = 2

    ##TODO: append all used blocks once to H and assert that increased rank == number appended rows

    # construct pcms
    asced_44_path_configs = []
    asced_44_pcms = []
    for i in range(num_used_blocks):
        assert i + block_offset <= number_additional_cyclic_blocks
        row_block = candidate_rows[
            (i + block_offset) * Zc : (i + block_offset + 1) * Zc
        ]
        # load pcms
        row_segments = np.split(row_block, pattern)
        for rows in row_segments:
            ##due to previous assert,
            rank_difference = rows.shape[0]
            H_aux = gf2(np.vstack((H, rows)))
            m_aux, _ = H_aux.shape
            asced_44_path_configs.append(channel_code_lib2.BP_config(H_aux))
            asced_44_path_configs[-1].early_stopping = True
            asced_44_path_configs[-1].max_iterations = 32
            asced_44_path_configs[-1].cn_update_type = "msa"
            asced_44_path_configs[-1].scheduling_type = "flooding"
            asced_44_path_configs[-1].norm_factor = 0.75
            if not use_all_zero:
                ##TODO offset takes on all binary vectors of length rank_difference
                all_affine_offsets = binary_vectors_in_suffix(m_aux, rank_difference)
                for offset in all_affine_offsets:
                    asced_44_path_configs.append(channel_code_lib2.BP_config(H_aux))
                    asced_44_path_configs[-1].early_stopping = True
                    asced_44_path_configs[-1].max_iterations = 32
                    asced_44_path_configs[-1].cn_update_type = "msa"
                    asced_44_path_configs[-1].scheduling_type = "flooding"
                    asced_44_path_configs[-1].norm_factor = 0.75
                    asced_44_path_configs[-1].affine_offset = offset

    print("Simulated num. aSCED paths:", len(asced_44_path_configs))
    asced_44_config = channel_code_lib2.Ensemble_config(H, asced_44_path_configs)

    sim_44_asced = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_44_asced.auto_save = auto_save
    sim_44_asced.save_dir = results_dir + "/aSCED_44_split0_batch1"
    # cfg.H = H

    sim_44_asced.puncturing(message_bit_pucturing)
    if not use_all_zero:
        sim_44_asced.init(enc_cfg, asced_44_config, use_all_zero)
    else:
        sim_44_asced.all_zero_init(asced_44_config)
    sim_44_asced.get_error_rates(sim_regime)
    FER_aSCED_44 = sim_44_asced.error_rates["FER-SNR"]
    print("aSCED finished")

    print(FER_aSCED_44)
