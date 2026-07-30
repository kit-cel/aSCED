## Goal is to create a plot for a certain target FER: Plot number of paths vs required SNR


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
import pandas as pd

if len(sys.argv) != 5:
    raise ValueError("Usage: python simulate.py <decoder_variant> <n_simul>")

decoder_variant = sys.argv[1].lower()
n_simul = int(sys.argv[2])

snr_start = float(sys.argv[3])
target_fer = float(sys.argv[4])


# n_simul = 180  # min= 78 or increase in stepzsizes of 6 upto  276 (Max supported is 300, however, some asced require 4 block à 6 rows)


results_dir = "RESULTS/fig_scatter_zc6"


use_all_zero = False

bg_vn = 12  ##bg2
bg_cn = 4  ##bg2
Zc = 6

# sim_regime = np.linspace(4, 4.5, 2)


# for maj rev: n=78 sim_regime = np.linspace(1, 5, 9)
# for maj rev: n=180 sim_regime = np.linspace(1, 3.5, 6)
# for maj rev: n=276

results_dir += f"_n={n_simul}"

splitting_pattern = [[1, 2, 3, 4, 5], [2, 4], [3]]


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

##splitting_pattern [1, 2, 3, 4, 5]
flag_aed_asced_12_split0 = False  # 1 batch of 6*2**1
flag_aed_asced_24_split0 = False  # 2 batches of  6*2**1
flag_aed_asced_48_split0 = False  # 4 batches of  6*2**1

##splitting pattern [2, 4]
flag_asced_12_split1 = False  # 3*(2**2)
flag_asced_24_split1 = False  # 2 batches of 3*(2**2)
flag_asced_48_split1 = False  # 4 batches of 3*(2**2)


##splitting pattern [3]

flag_asced_16 = False  # 2*(2**3)
flag_asced_48 = False  # 3*(2**3)
flag_asced_64 = False  # 4*(2**3)

##no spliitng
flag_asced_64_nosplit = False  # 2**6
flag_asced_128 = False  # 2batches of 2**6

if decoder_variant == "nmsa":
    flag_nmsa = True

elif decoder_variant == "aed":
    flag_aed = True

# splitting pattern [1,2,3,4,5]
elif decoder_variant == "12_split0":
    flag_aed_asced_12_split0 = True

elif decoder_variant == "24_split0":
    flag_aed_asced_24_split0 = True

elif decoder_variant == "48_split0":
    flag_aed_asced_48_split0 = True

# splitting pattern [2,4]
elif decoder_variant == "12_split1":
    flag_asced_12_split1 = True

elif decoder_variant == "24_split1":
    flag_asced_24_split1 = True

elif decoder_variant == "48_split1":
    flag_asced_48_split1 = True

# splitting pattern [3]
elif decoder_variant == "16":
    flag_asced_16 = True

elif decoder_variant == "48":
    flag_asced_48 = True

elif decoder_variant == "64":
    flag_asced_64 = True

# no splitting
elif decoder_variant == "64_nosplit":
    flag_asced_64_nosplit = True

elif decoder_variant == "128":
    flag_asced_128 = True

else:
    raise ValueError(f"Unknown decoder variant '{decoder_variant}'")

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
bool_emulate_stopping = False

target_fraction_coverged_path = 0.5


if use_all_zero:
    results_dir += "_AZ"


# we start first with code max_rank_5G_zc=11
code = np.load("Codes/TCOM_aSCED/5G_zc=6/max_rank_5G_zc=6.npz")
# this is the full BG2 PCM lifted with Zc=11
H_full = code["h"].astype(int)

print(H_full.shape)


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

    return channel_code_lib2.Ensemble_config(H, path_configs), len(path_configs)


def create_asced_sim(
    save_dir_name,
    split_pattern,
    num_used_blocks,
):
    ensemble_cfg, num_paths = create_asced_config(
        H=H,
        candidate_rows=candidate_rows,
        block_offset=block_offset,
        Zc=Zc,
        split_pattern=split_pattern,
        num_used_blocks=num_used_blocks,
        use_all_zero=use_all_zero,
    )

    sim = channel_code_lib2.Simulation_Env(k, n, "all")

    sim.auto_save = auto_save
    sim.save_dir = os.path.join(results_dir, save_dir_name)

    sim.puncturing(message_bit_pucturing)

    if use_all_zero:
        sim.all_zero_init(ensemble_cfg)
    else:
        sim.init(enc_cfg, ensemble_cfg, use_all_zero)

    return sim, num_paths


def run_at_snr(sim, snr):

    sim.get_error_rates([snr])
    return sim.error_rates["FER-SNR"][snr]


def search_target_fer(
    sim,
    start_snr,
    target_fer,
    tolerance=0.03,      # ±5%
    initial_step=0.3,
    min_step=0.001,
    max_iter=30,
):
    """
    Find the SNR required to achieve the target FER.

    Parameters
    ----------
    sim : Simulation_Env
    start_snr : float
        Initial guess.
    target_fer : float
    tolerance : float
        Relative tolerance. 0.05 means ±5%.
    initial_step : float
        Initial search step while bracketing.
    min_step : float
        Stop once the SNR interval becomes smaller than this.
    """

    lower_fer = target_fer * (1 - tolerance)
    upper_fer = target_fer * (1 + tolerance)

    def evaluate(snr):
        fer = run_at_snr(sim, snr)
        print(f"SNR = {snr:.3f} dB, FER = {fer:.3e}")
        return fer

    # -------------------------------------------------------
    # Evaluate starting point
    # -------------------------------------------------------

    snr = start_snr
    fer = evaluate(snr)

    if lower_fer <= fer <= upper_fer:
        return snr, fer

    # -------------------------------------------------------
    # Bracketing phase
    # -------------------------------------------------------

    step = initial_step

    if fer > target_fer:
        # Need higher SNR
        low_snr = snr
        low_fer = fer

        while True:
            snr += step
            fer = evaluate(snr)

            if lower_fer <= fer <= upper_fer:
                return snr, fer

            if fer < target_fer:
                high_snr = snr
                high_fer = fer
                break

    else:
        # Need lower SNR
        high_snr = snr
        high_fer = fer

        while True:
            snr -= step
            fer = evaluate(snr)

            if lower_fer <= fer <= upper_fer:
                return snr, fer

            if fer > target_fer:
                low_snr = snr
                low_fer = fer
                break

    # -------------------------------------------------------
    # Binary search
    # -------------------------------------------------------

    for _ in range(max_iter):

        if high_snr - low_snr < min_step:
            print("stop due to too small step")
            break

        mid_snr = 0.5 * (low_snr + high_snr)
        mid_fer = evaluate(mid_snr)

        if lower_fer <= mid_fer <= upper_fer:
            return mid_snr, mid_fer

        if mid_fer > target_fer:
            # Too many errors -> increase SNR
            low_snr = mid_snr
            low_fer = mid_fer
        else:
            # Too few errors -> decrease SNR
            high_snr = mid_snr
            high_fer = mid_fer

    # -------------------------------------------------------
    # Return whichever endpoint is closer
    # -------------------------------------------------------

    if abs(low_fer - target_fer) < abs(high_fer - target_fer):
        return low_snr, low_fer

    return high_snr, high_fer



experiments = [
    # splitting pattern [1,2,3,4,5]
    (
        flag_aed_asced_12_split0,
        "aSCED_12_split0_batch1",
        splitting_pattern[0],
        1,
    ),
    (
        flag_aed_asced_24_split0,
        "aSCED_24_split0_batch2",
        splitting_pattern[0],
        2,
    ),
    (
        flag_aed_asced_48_split0,
        "aSCED_48_split0_batch4",
        splitting_pattern[0],
        4,
    ),
    # splitting pattern [2,4]
    (
        flag_asced_12_split1,
        "aSCED_12_split1_batch1",
        splitting_pattern[1],
        1,
    ),
    (
        flag_asced_24_split1,
        "aSCED_24_split1_batch2",
        splitting_pattern[1],
        2,
    ),
    (
        flag_asced_48_split1,
        "aSCED_48_split1_batch4",
        splitting_pattern[1],
        4,
    ),
    # splitting pattern [3]
    (
        flag_asced_16,
        "aSCED_16_split2_batch1",
        splitting_pattern[2],
        1,
    ),
    (
        flag_asced_48,
        "aSCED_48_split2_batch3",
        splitting_pattern[2],
        3,
    ),
    (
        flag_asced_64,
        "aSCED_64_split2_batch4",
        splitting_pattern[2],
        4,
    ),
    # no splitting
    (
        flag_asced_64_nosplit,
        "aSCED_64_nosplit_batch1",
        [],
        1,
    ),
    (
        flag_asced_128,
        "aSCED_128_nosplit_batch2",
        [],
        2,
    ),
]

summary = []


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

    
    required_snr, accepted_fer = search_target_fer(
            sim_nmsa,
            snr_start,
            target_fer,
        )

    summary.append(
        {
            "decoder": "nmsa",
            "split_pattern": str(""),
            "num_blocks": 0,
            "num_paths": 1,
            "required_snr": required_snr,
            "accepted_fer": accepted_fer,
        }
    )

for enabled, save_name, split_pattern, num_blocks in experiments:

    if not enabled:
        continue

    sim, num_paths = create_asced_sim(
        save_name,
        split_pattern,
        num_blocks,
    )

    required_snr, accepted_fer = search_target_fer(
        sim,
        snr_start,
        target_fer,
    )

    summary.append(
        {
            "decoder": save_name,
            "split_pattern": str(split_pattern),
            "num_blocks": num_blocks,
            "num_paths": num_paths,
            "required_snr": required_snr,
            "accepted_fer": accepted_fer,
        }
    )


df = pd.DataFrame(summary)

summary_dir = os.path.join(results_dir, "summary")

os.makedirs(summary_dir, exist_ok=True)

df.to_csv(
    os.path.join(summary_dir, "paths_vs_required_snr.csv"),
    index=False,
)