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
flag_aed_asced_2_split1_subsplit1 = False
flag_aed_asced_4_split1_subsplit2 = False
flag_aed_asced_6_split1_subsplit3 = False
flag_aed_asced_8_split1_subsplit4 = False
flag_aed_asced_10_split1_subsplit5 = False
flag_aed_asced_12_split1 = False  # 1 batch of 6*2**1
flag_aed_asced_24_split1 = False  # 2 batches of  6*2**1
flag_aed_asced_48_split1 = False  # 4 batches of  6*2**1
flag_aed_asced_60_split1 = False  # 5 batches of  6*2**1
flag_aed_asced_96_split1 = False


##splitting pattern [2, 4]
# split2 subsplits
flag_asced_4_split2_subsplit1 = False
flag_asced_8_split2_subsplit2 = False
flag_asced_12_split2 = False  # 3*(2**2)
flag_asced_24_split2 = False  # 2 batches of 3*(2**2)
flag_asced_48_split2 = False  # 4 batches of 3*(2**2)
flag_asced_96_split2 = False


##splitting pattern [3]
# split3 subsplits
flag_asced_8_split3_subsplit1 = False
flag_asced_16_split3 = False  # 2*(2**3)
flag_asced_48_split3 = False  # 3*(2**3)
flag_asced_64_split3 = False  # 4*(2**3)
flag_asced_96_split3 = False
flag_asced_128_split3 = False

##no spliitng
flag_asced_64_nosplit = False  # 2**6
flag_asced_128_nosplit = False  # 2batches of 2**6
flag_asced_192_nosplit = False

if decoder_variant == "nmsa":
    flag_nmsa = True

elif decoder_variant == "aed":
    flag_aed = True

# ----------------------------------------------------------
# splitting pattern [1,2,3,4,5]
# ----------------------------------------------------------
elif decoder_variant == "2_split1_subsplit1":
    flag_aed_asced_2_split1_subsplit1 = True

elif decoder_variant == "4_split1_subsplit2":
    flag_aed_asced_4_split1_subsplit2 = True

elif decoder_variant == "6_split1_subsplit3":
    flag_aed_asced_6_split1_subsplit3 = True

elif decoder_variant == "8_split1_subsplit4":
    flag_aed_asced_8_split1_subsplit4 = True

elif decoder_variant == "10_split1_subsplit5":
    flag_aed_asced_10_split1_subsplit5 = True

elif decoder_variant == "12_split1":
    flag_aed_asced_12_split1 = True

elif decoder_variant == "24_split1":
    flag_aed_asced_24_split1 = True

elif decoder_variant == "48_split1":
    flag_aed_asced_48_split1 = True

elif decoder_variant == "60_split1":
    flag_aed_asced_60_split1 = True

elif decoder_variant == "96_split1":
    flag_aed_asced_96_split1 = True

# ----------------------------------------------------------
# splitting pattern [2,4]
# ----------------------------------------------------------
elif decoder_variant == "4_split2_subsplit1":
    flag_asced_4_split2_subsplit1 = True

elif decoder_variant == "8_split2_subsplit2":
    flag_asced_8_split2_subsplit2 = True

elif decoder_variant == "12_split2":
    flag_asced_12_split2 = True

elif decoder_variant == "24_split2":
    flag_asced_24_split2 = True

elif decoder_variant == "48_split2":
    flag_asced_48_split2 = True

elif decoder_variant == "96_split2":
    flag_asced_96_split2 = True

# ----------------------------------------------------------
# splitting pattern [3]
# ----------------------------------------------------------
elif decoder_variant == "8_split3_subsplit1":
    flag_asced_8_split3_subsplit1 = True

elif decoder_variant == "16_split3":
    flag_asced_16_split3 = True

elif decoder_variant == "48_split3":
    flag_asced_48_split3 = True

elif decoder_variant == "64_split3":
    flag_asced_64_split3 = True

elif decoder_variant == "96_split3":
    flag_asced_96_split3 = True

elif decoder_variant == "128_split3":
    flag_asced_128_split3 = True
# ----------------------------------------------------------
# no splitting
# ----------------------------------------------------------
elif decoder_variant == "64_nosplit":
    flag_asced_64_nosplit = True

elif decoder_variant == "128_nosplit":
    flag_asced_128_nosplit = True

elif decoder_variant == "192_nosplit":
    flag_asced_192_nosplit = True

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
    subsplit=None,
):
    """Construct all decoder paths for one ASCED configuration."""

    path_configs = []

    for block in range(num_used_blocks):
        print(block)
        print(block_offset)
        print(number_additional_cyclic_blocks)
        assert block + block_offset <= number_additional_cyclic_blocks

        row_block = candidate_rows[
            (block + block_offset) * Zc : (block + block_offset + 1) * Zc
        ]

        # no split -> use the whole block
        if len(split_pattern) == 0:
            row_segments = [row_block]
        else:
            row_segments = np.split(row_block, split_pattern)

        if subsplit is not None:
            assert 1 <= subsplit <= len(row_segments)
            row_segments = row_segments[:subsplit]

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
    tolerance=0.03,  # ±5%
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
        flag_aed_asced_2_split1_subsplit1,
        "aSCED_12_split1_batch1_subsplit1",
        splitting_pattern[0],
        1,
        1,
    ),
    (
        flag_aed_asced_4_split1_subsplit2,
        "aSCED_12_split1_batch1_subsplit2",
        splitting_pattern[0],
        1,
        2,
    ),
    (
        flag_aed_asced_6_split1_subsplit3,
        "aSCED_12_split1_batch1_subsplit3",
        splitting_pattern[0],
        1,
        3,
    ),
    (
        flag_aed_asced_8_split1_subsplit4,
        "aSCED_12_split1_batch1_subsplit4",
        splitting_pattern[0],
        1,
        4,
    ),
    (
        flag_aed_asced_10_split1_subsplit5,
        "aSCED_12_split1_batch1_subsplit5",
        splitting_pattern[0],
        1,
        5,
    ),
    (
        flag_aed_asced_12_split1,
        "aSCED_12_split1_batch1",
        splitting_pattern[0],
        1,
        None,
    ),
    (
        flag_aed_asced_24_split1,
        "aSCED_24_split1_batch2",
        splitting_pattern[0],
        2,
        None,
    ),
    (
        flag_aed_asced_48_split1,
        "aSCED_48_split1_batch4",
        splitting_pattern[0],
        4,
        None,
    ),
    (
        flag_aed_asced_60_split1,
        "aSCED_60_split1_batch5",
        splitting_pattern[0],
        5,
        None,
    ),
    (
        flag_aed_asced_96_split1,
        "aSCED_96_split1_batch8",
        splitting_pattern[0],
        8,
        None,
    ),
    # splitting pattern [2,4]
    (
        flag_asced_4_split2_subsplit1,
        "aSCED_12_split2_batch1_subsplit1",
        splitting_pattern[1],
        1,
        1,
    ),
    (
        flag_asced_8_split2_subsplit2,
        "aSCED_12_split2_batch1_subsplit2",
        splitting_pattern[1],
        1,
        2,
    ),
    (
        flag_asced_12_split2,
        "aSCED_12_split2_batch1",
        splitting_pattern[1],
        1,
        None,
    ),
    (
        flag_asced_24_split2,
        "aSCED_24_split2_batch2",
        splitting_pattern[1],
        2,
        None,
    ),
    (
        flag_asced_48_split2,
        "aSCED_48_split2_batch4",
        splitting_pattern[1],
        4,
        None,
    ),
    (
        flag_asced_96_split2,
        "aSCED_96_split2_batch8",
        splitting_pattern[1],
        8,
        None,
    ),
    # splitting pattern [3]
    (
        flag_asced_8_split3_subsplit1,
        "aSCED_16_split3_batch1_subsplit1",
        splitting_pattern[2],
        1,
        1,
    ),
    (
        flag_asced_16_split3,
        "aSCED_16_split3_batch1",
        splitting_pattern[2],
        1,
        None,
    ),
    (
        flag_asced_48_split3,
        "aSCED_48_split3_batch3",
        splitting_pattern[2],
        3,
        None,
    ),
    (
        flag_asced_64_split3,
        "aSCED_64_split3_batch4",
        splitting_pattern[2],
        4,
        None,
    ),
    (
        flag_asced_96_split3,
        "aSCED_96_split3_batch6",
        splitting_pattern[2],
        6,
        None,
    ),
    (
        flag_asced_128_split3,
        "aSCED_128_split3_batch8",
        splitting_pattern[2],
        8,
        None,
    ),
    # no splitting
    (
        flag_asced_64_nosplit,
        "aSCED_64_nosplit_batch1",
        [],
        1,
        None,
    ),
    (
        flag_asced_128_nosplit,
        "aSCED_128_nosplit_batch2",
        [],
        2,
        None,
    ),
    (
        flag_asced_192_nosplit,
        "aSCED_192_nosplit_batch3",
        [],
        3,
        None,
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
            "split": 0,
            "numblocks": 0,
            "numpaths": 1,
            "requiredsnr": required_snr,
            "acceptedfer": accepted_fer,
        }
    )


def quasi_cyclic_permutation_vector(length, block_size=11):
    permuted_indices = np.arange(length)
    for start in range(0, length, block_size):
        end = min(start + block_size, length)
        block_indices = permuted_indices[start:end]
        if len(block_indices) == block_size:
            permuted_indices[start:end] = np.roll(block_indices, 1)
    return permuted_indices


if flag_aed:
    permutation = quasi_cyclic_permutation_vector(n, Zc)
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
    sim_aed = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_aed.auto_save = auto_save
    sim_aed.save_dir = results_dir + f"/aed{Zc}"
    sim_aed.puncturing(message_bit_pucturing)
    if not use_all_zero:
        sim_aed.init(enc_cfg, ensemble_decoder_config, use_all_zero)
    else:
        sim_aed.all_zero_init(ensemble_decoder_config)

    required_snr, accepted_fer = search_target_fer(
        sim_aed,
        snr_start,
        target_fer,
    )
    summary.append(
        {
            "split": 5,
            "numblocks": 1,
            "numpaths": Zc,
            "requiredsnr": required_snr,
            "acceptedfer": accepted_fer,
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

    if "split1" in save_name:
        split = 1
    elif "split2" in save_name:
        split = 2
    elif "split3" in save_name:
        split = 3
    elif "nosplit" in save_name:
        split = 4
    else:
        split = -1

    summary.append(
        {
            "split": split,
            "numblocks": num_blocks,
            "numpaths": num_paths,
            "requiredsnr": required_snr,
            "acceptedfer": accepted_fer,
        }
    )


df = pd.DataFrame(summary)

summary_dir = os.path.join(results_dir, "summary")

os.makedirs(summary_dir, exist_ok=True)

summary_path = os.path.join(summary_dir, "paths_vs_required_snr.dat")

new_df = pd.DataFrame(summary)

if os.path.exists(summary_path):
    old_df = pd.read_csv(summary_path, sep=" ")

    # Remove already existing identical split/path combinations
    new_df_keys = new_df[["split", "num_paths"]]
    old_df = old_df.merge(
        new_df_keys,
        on=["split", "num_paths"],
        how="left",
        indicator=True,
    )

    old_df = old_df[old_df["_merge"] == "left_only"]
    old_df = old_df.drop(columns=["_merge"])

    new_df = pd.concat([old_df, new_df], ignore_index=True)


new_df = new_df.sort_values(["split", "num_paths"])

new_df.to_csv(
    summary_path,
    sep=" ",
    index=False,
)

print(new_df)
print()
print(f"Saved summary to {summary_path}")
