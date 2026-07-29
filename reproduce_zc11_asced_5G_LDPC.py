# pylint: disable=invalid-name
"""Reproduce Figure x results for ASCED on 5G LDPC codes."""

import numpy as np
from time import time
import os
import sys

import galois

gf2 = galois.GF2

import channel_code_lib2


import matplotlib.pyplot as plt

if len(sys.argv) != 5:
    raise ValueError("Usage: python simulate.py <decoder_variant> <n_simul>")

decoder_variant = sys.argv[1].lower()
n_simul = int(sys.argv[2])

snr_start = float(sys.argv[3])
snr_end = float(sys.argv[4])

# Include the end point
sim_regime = np.arange(snr_start, snr_end + 0.25, 0.5)



remove = 4

results_dir = f"RESULTS/fig_x_zc11_r{remove}"
# für k=66
# wirf bg_vns 6-10 raus --> checke ob nmsa korrekt mit in paper
bg_vn = 12  ##bg2
bg_cn = 4  ##bg2
Zc = 11

bg_vn=bg_vn-remove

# sim_regime = np.linspace(2, 4.5, 6)

# n_simul = 143  # or increase in stepzsizes of 11 e.g. 165

results_dir += f"n={n_simul}"

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


flag_nmsa = False  # if true simulates NMSA-32

flag_aed = False  # if true simulates AED-11

##splitting_pattern [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
flag_aed_asced_22 = False  # 1 batch of 11*2**1
flag_aed_asced_44 = False  # 2 batches of 11*2**1
flag_aed_asced_88 = False  # 4 batches of 11*2**1

##splitting pattern [2, 4, 6, 8]
flag_asced_24 = False  # 2**2*4+2**3
flag_asced_48 = False  # 2 batches of 2**2*4+2**3
flag_asced_96_4batch = False  # 4 batches of 2**2*4+2**3

##splitting_pattern [5]
flag_asced_96 = False  # 2**5+2**6
flag_asced_192 = False  # 2batches of 2**5+2**6

flag_asced_2048 = False  #  2**11


if decoder_variant == "nmsa":
    flag_nmsa = True

elif decoder_variant == "aed":
    flag_aed = True

elif decoder_variant == "asced22":
    flag_aed_asced_22 = True

elif decoder_variant == "asced44":
    flag_aed_asced_44 = True

elif decoder_variant == "asced88":
    flag_aed_asced_88 = True

elif decoder_variant == "asced24":
    flag_asced_24 = True

elif decoder_variant == "asced48":
    flag_asced_48 = True

elif decoder_variant == "asced96_4batch":
    flag_asced_96_4batch = True

elif decoder_variant == "asced96":
    flag_asced_96 = True

elif decoder_variant == "asced192":
    flag_asced_192 = True

elif decoder_variant == "asced2048":
    flag_asced_2048 = True

else:
    raise ValueError(f"Unknown decoder variant '{decoder_variant}'")

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
bool_emulate_stopping = False

target_fraction_coverged_path = 0.5


use_all_zero = True

if use_all_zero:
    results_dir += "_AZ"


# we start first with code max_rank_5G_zc=11
code = np.load("Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz")
# this is the full BG2 PCM lifted with Zc=11
H_full = code["h"].astype(int)



H_full = np.delete(H_full, np.arange(6*Zc, 10*Zc), axis=1)

block_offset = 0  # first batch, uses block 0, 2nd batch, uses block and so on


number_vn_simul = n_simul + 2 * Zc

number_vns_start = (
    bg_vn * Zc
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


def create_asced_config(
    H,
    candidate_rows,
    block_offset,
    Zc,
    split_pattern,
    num_used_blocks,
    use_all_zero,
):
    """Construct all decoder paths for one ASCED configuration."""

    path_configs = []

    for block in range(num_used_blocks):

        assert block + block_offset <= number_additional_cyclic_blocks

        row_block = candidate_rows[
            (block + block_offset) * Zc : (block + block_offset + 1) * Zc
        ]

        # no split -> use the whole block
        if len(split_pattern) == 0:
            row_segments = [row_block]
        else:
            row_segments = np.split(row_block, split_pattern)

        for rows in row_segments:

            H_aux = gf2(np.vstack((H, rows)))
            m_aux = H_aux.shape[0]

            # zero affine offset
            path_configs.append(create_bp_config(H_aux))

            if not use_all_zero:

                rank_difference = rows.shape[0]

                for offset in binary_vectors_in_suffix(
                    m_aux,
                    rank_difference,
                ):
                    cfg = create_bp_config(H_aux)
                    cfg.affine_offset = offset
                    path_configs.append(cfg)

    print("Simulated num. aSCED paths:", len(path_configs))

    return channel_code_lib2.Ensemble_config(H, path_configs)


def run_asced(
    save_dir_name,
    split_pattern,
    num_used_blocks,
):
    ensemble_cfg = create_asced_config(
        H=H,
        candidate_rows=candidate_rows,
        block_offset=block_offset,
        Zc=Zc,
        split_pattern=split_pattern,
        num_used_blocks=num_used_blocks,
        use_all_zero=use_all_zero,
    )
    print("starting", save_dir_name)

    sim = channel_code_lib2.Simulation_Env(k, n, "all")

    sim.auto_save = auto_save
    sim.save_dir = results_dir + "/" + save_dir_name

    sim.puncturing(message_bit_pucturing)

    if use_all_zero:
        sim.all_zero_init(ensemble_cfg)
    else:
        sim.init(enc_cfg, ensemble_cfg, use_all_zero)

    sim.get_error_rates(sim_regime)

    print(save_dir_name, "finished")
    print(sim.error_rates["FER-SNR"])

    return sim.error_rates["FER-SNR"]


if flag_nmsa:
    nmsa_config = channel_code_lib2.BP_config(H)
    nmsa_config.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    nmsa_config.max_iterations = (
        32  # set maximum number of BP iterations; default is 32  )
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
    print(FER_nmsa)  #
    print("NMSA 32 finished")

experiments = [
    # splitting pattern [1,2,3,...,10]
    (
        flag_aed_asced_22,
        "aSCED_22_split0_batch1",
        splitting_pattern[0],
        1,
    ),
    (
        flag_aed_asced_44,
        "aSCED_44_split0_batch2",
        splitting_pattern[0],
        2,
    ),
    (
        flag_aed_asced_88,
        "aSCED_88_split0_batch4",
        splitting_pattern[0],
        4,
    ),
    # splitting pattern [2,4,6,8]
    (
        flag_asced_24,
        "aSCED_24_split2_batch1",
        splitting_pattern[2],
        1,
    ),
    (
        flag_asced_48,
        "aSCED_48_split2_batch2",
        splitting_pattern[2],
        2,
    ),
    (
        flag_asced_96_4batch,
        "aSCED_96_split2_batch4",
        splitting_pattern[2],
        4,
    ),
    # splitting pattern [5]
    (
        flag_asced_96,
        "aSCED_96_split1_batch1",
        splitting_pattern[1],
        1,
    ),
    (
        flag_asced_192,
        "aSCED_192_split1_batch2",
        splitting_pattern[1],
        2,
    ),
    # no split
    (
        flag_asced_2048,
        "aSCED_2048",
        [],
        1,
    ),
]

FER_results = {}

for enabled, save_name, split_pattern, num_blocks in experiments:

    if not enabled:
        continue

    FER_results[save_name] = run_asced(
        save_dir_name=save_name,
        split_pattern=split_pattern,
        num_used_blocks=num_blocks,
    )
