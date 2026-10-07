## Goal: Plot number of paths vs required SNR, GREEDY block-selection variant.
#
# Production validation script for the Zc=6 greedy block-selection search
# (sweeps/greedy_block_search_zc6.py + sweeps/run_greedy_search_zc6.py).
#
# Identical structure/CLI to reproduce_fig5_scatter_plot_zc6_5G_LDPC.py (the
# reference script this was forked from), but instead of using
# `num_used_blocks` *contiguous* blocks starting at `block_offset=0` ("the
# subsequent blocks" baseline), it uses the specific, greedy-discovered LIST
# of block indices recorded in sweeps/greedy_search_state_zc6_n<n_simul>.json
# ("fixed_blocks", in discovery order: index 0 was picked first / is the most
# valuable), for whichever split pattern / num_blocks / subsplit the variant
# calls for. The greedy chain itself was only ever evaluated under the split3
# (Delta=3) sub-splitting pattern (a pure block-SELECTION search); the same
# discovered order is then reused as-is across split1/split2/split3/nosplit
# here, per the task's explicit scope ("use the first L blocks of the
# discovered chain in place of blocks 0..L-1").
#
# Deliberate deviation from the reference script (flagged per the task
# description): the reference's create_asced_sim() never sets
# sim.max_transmissions (unbounded), which risks an unbounded run at
# low-FER operating points when combined with search_target_fer()'s binary
# search. This script adds max_transmissions=int(2e6), same convention used
# elsewhere in this effort (e.g. the Zc=11 production sweeps).
#
# pylint: disable=invalid-name
"""Reproduce Figure 5-style results for ASCED on the Zc=6 5G LDPC codes,
using greedy-discovered block order instead of the natural "subsequent
blocks" order."""

import json
import numpy as np
import os
import sys

import galois

gf2 = galois.GF2

import channel_code_lib2

import pandas as pd

if len(sys.argv) != 5:
    raise ValueError("Usage: python reproduce_fig5_scatter_plot_zc6_5G_LDPC_greedy.py <decoder_variant> <n_simul> <snr_start> <target_fer>")

decoder_variant = sys.argv[1].lower()
n_simul = int(sys.argv[2])

snr_start = float(sys.argv[3])
target_fer = float(sys.argv[4])

if n_simul not in (78, 180):
    raise ValueError(f"Unsupported n_simul={n_simul}, expected 78 or 180")

# Same results_dir ROOT as the baseline (reference) script, so the final
# summary .dat file lands in the SAME summary/ folder (parallel file, see
# below) - but per-decoder detailed run output goes to a "greedy" subdir to
# avoid colliding with the baseline sweep's identically-named leaf
# directories (the same collision class fixed for the Zc=11 greedy search,
# see STATUS.md's "Fixed a latent collision bug" note).
results_dir_base = "RESULTS/fig_scatter_zc6"
results_dir_base += f"_n={n_simul}"
results_dir = os.path.join(results_dir_base, "greedy")

bg_vn = 12  ##bg2
bg_cn = 4  ##bg2
Zc = 6

splitting_pattern = [[1, 2, 3, 4, 5], [2, 4], [3]]

STATE_FILE = f"sweeps/greedy_search_state_zc6_n{n_simul}.json"
EXPECTED_CHAIN_LEN = 8  # see task scope: chain depth up to 8 blocks

with open(STATE_FILE) as f:
    _state = json.load(f)
GREEDY_BLOCK_ORDER = _state["fixed_blocks"]  # greedy discovery order, 1st pick first
if len(GREEDY_BLOCK_ORDER) < EXPECTED_CHAIN_LEN:
    raise ValueError(
        f"{STATE_FILE} has only {len(GREEDY_BLOCK_ORDER)} fixed_blocks, expected "
        f"{EXPECTED_CHAIN_LEN} - greedy search for n_simul={n_simul} may not have finished."
    )
print(f"Using greedy-discovered block order for n_simul={n_simul}: {GREEDY_BLOCK_ORDER}")


def create_bp_config(H_aux):
    """Create a BP configuration with the default decoder settings."""
    cfg = channel_code_lib2.BP_config(H_aux)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    return cfg


flag_nmsa = False
flag_aed = False

##splitting_pattern [1, 2, 3, 4, 5]
flag_aed_asced_2_split1_subsplit1 = False
flag_aed_asced_4_split1_subsplit2 = False
flag_aed_asced_6_split1_subsplit3 = False
flag_aed_asced_8_split1_subsplit4 = False
flag_aed_asced_10_split1_subsplit5 = False
flag_aed_asced_12_split1 = False
flag_aed_asced_24_split1 = False
flag_aed_asced_36_split1 = False
flag_aed_asced_48_split1 = False
flag_aed_asced_60_split1 = False
flag_aed_asced_72_split1 = False
flag_aed_asced_84_split1 = False
flag_aed_asced_96_split1 = False

##splitting pattern [2, 4]
flag_asced_4_split2_subsplit1 = False
flag_asced_8_split2_subsplit2 = False
flag_asced_12_split2 = False
flag_asced_24_split2 = False
flag_asced_36_split2 = False
flag_asced_48_split2 = False
flag_asced_60_split2 = False
flag_asced_72_split2 = False
flag_asced_84_split2 = False
flag_asced_96_split2 = False

##splitting pattern [3]
flag_asced_8_split3_subsplit1 = False
flag_asced_16_split3 = False
flag_asced_32_split3 = False
flag_asced_48_split3 = False
flag_asced_64_split3 = False
flag_asced_80_split3 = False
flag_asced_96_split3 = False
flag_asced_112_split3 = False
flag_asced_128_split3 = False

##no spliitng
flag_asced_64_nosplit = False
flag_asced_128_nosplit = False
flag_asced_192_nosplit = False

if decoder_variant == "nmsa":
    flag_nmsa = True
elif decoder_variant == "aed":
    flag_aed = True
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
elif decoder_variant == "36_split1":
    flag_aed_asced_36_split1 = True
elif decoder_variant == "48_split1":
    flag_aed_asced_48_split1 = True
elif decoder_variant == "60_split1":
    flag_aed_asced_60_split1 = True
elif decoder_variant == "72_split1":
    flag_aed_asced_72_split1 = True
elif decoder_variant == "84_split1":
    flag_aed_asced_84_split1 = True
elif decoder_variant == "96_split1":
    flag_aed_asced_96_split1 = True
elif decoder_variant == "4_split2_subsplit1":
    flag_asced_4_split2_subsplit1 = True
elif decoder_variant == "8_split2_subsplit2":
    flag_asced_8_split2_subsplit2 = True
elif decoder_variant == "12_split2":
    flag_asced_12_split2 = True
elif decoder_variant == "24_split2":
    flag_asced_24_split2 = True
elif decoder_variant == "36_split2":
    flag_asced_36_split2 = True
elif decoder_variant == "48_split2":
    flag_asced_48_split2 = True
elif decoder_variant == "60_split2":
    flag_asced_60_split2 = True
elif decoder_variant == "72_split2":
    flag_asced_72_split2 = True
elif decoder_variant == "84_split2":
    flag_asced_84_split2 = True
elif decoder_variant == "96_split2":
    flag_asced_96_split2 = True
elif decoder_variant == "8_split3_subsplit1":
    flag_asced_8_split3_subsplit1 = True
elif decoder_variant == "16_split3":
    flag_asced_16_split3 = True
elif decoder_variant == "32_split3":
    flag_asced_32_split3 = True
elif decoder_variant == "48_split3":
    flag_asced_48_split3 = True
elif decoder_variant == "64_split3":
    flag_asced_64_split3 = True
elif decoder_variant == "80_split3":
    flag_asced_80_split3 = True
elif decoder_variant == "96_split3":
    flag_asced_96_split3 = True
elif decoder_variant == "112_split3":
    flag_asced_112_split3 = True
elif decoder_variant == "128_split3":
    flag_asced_128_split3 = True
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

use_all_zero = False

# we start first with code max_rank_5G_zc=6
code = np.load("Codes/TCOM_aSCED/5G_zc=6/max_rank_5G_zc=6.npz")
H_full = code["h"].astype(int)

print(H_full.shape)

number_vn_simul = n_simul + 2 * Zc
number_vns_start = bg_vn * Zc + 2 * Zc
harq_part = (number_vn_simul - number_vns_start) // Zc
print(harq_part)
m = Zc * bg_cn + harq_part * Zc
number_vns = number_vns_start + harq_part * Zc
assert number_vn_simul == number_vns

H = gf2(H_full[:m, :number_vns])
G = np.array(gf2(H[:m, :number_vns]).null_space()).astype(int)
k, n = G.shape
assert n == number_vns
print(f"Simulating number vns={n}, n={number_vns-len(message_bit_pucturing)},k={k} ")

candidate_rows = H_full[m:, : number_vns_start + harq_part * Zc]
number_additional_cyclic_blocks = candidate_rows.shape[0] // Zc
assert number_additional_cyclic_blocks * Zc == candidate_rows.shape[0]

for b in GREEDY_BLOCK_ORDER:
    assert 0 <= b < number_additional_cyclic_blocks, f"block index {b} out of range"


def binary_vectors_in_suffix(number_rows, number_last_entries):
    v = [0] * number_rows
    for mask in range(1, 1 << number_last_entries):
        for i in range(number_last_entries):
            v[number_rows - number_last_entries + i] = (mask >> i) & 1
        yield v.copy()


if not use_all_zero:
    enc_cfg = channel_code_lib2.PCM_Encoder_config(H, k, n)


def create_asced_config(
    H,
    candidate_rows,
    block_order,
    Zc,
    split_pattern,
    num_used_blocks,
    use_all_zero,
    subsplit=None,
):
    """Construct all decoder paths for one ASCED configuration, using the
    first `num_used_blocks` entries of `block_order` (the greedy-discovered
    chain, IN ORDER - order matters for fixed_sequential path indexing, see
    the module docstring) in place of the reference script's contiguous
    `block_offset=0, blocks 0..num_used_blocks-1` selection.
    """

    path_configs = []

    assert num_used_blocks <= len(block_order), (
        f"requested {num_used_blocks} blocks but greedy chain only has {len(block_order)}"
    )

    for i in range(num_used_blocks):
        block = block_order[i]
        assert 0 <= block <= number_additional_cyclic_blocks

        row_block = candidate_rows[block * Zc : (block + 1) * Zc]

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


def create_asced_sim(save_dir_name, split_pattern, num_used_blocks, subsplit=None):
    ensemble_cfg, num_paths = create_asced_config(
        H=H,
        candidate_rows=candidate_rows,
        block_order=GREEDY_BLOCK_ORDER,
        Zc=Zc,
        split_pattern=split_pattern,
        num_used_blocks=num_used_blocks,
        use_all_zero=use_all_zero,
        subsplit=subsplit,
    )

    sim = channel_code_lib2.Simulation_Env(k, n, "all")
    sim.target_errors = 1000
    # Deliberate deviation from the reference script (see module docstring):
    # cap max_transmissions so search_target_fer() can't run unboundedly at
    # low-FER operating points.
    sim.max_transmissions = int(2e6)
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus is not None:
        sim.num_threads = int(slurm_cpus)

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
    tolerance=0.02,
    initial_step=0.2,
    min_step=0.0001,
    max_iter=40,
):
    """Copied verbatim from reproduce_fig5_scatter_plot_zc6_5G_LDPC.py - same
    binary-search methodology, matching Fig. 5's actual methodology exactly
    (per the task's explicit choice), modulo the max_transmissions cap set
    above on the sim itself."""

    lower_fer = target_fer * (1 - tolerance)
    upper_fer = target_fer * (1 + tolerance)

    def evaluate(snr):
        fer = run_at_snr(sim, snr)
        print(f"SNR = {snr:.3f} dB, FER = {fer:.3e}")
        return fer

    snr = start_snr
    fer = evaluate(snr)

    if lower_fer <= fer <= upper_fer:
        return snr, fer

    step = initial_step

    if fer > target_fer:
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

    for _ in range(max_iter):

        if high_snr - low_snr < min_step:
            print("stop due to too small step")
            break

        mid_snr = 0.5 * (low_snr + high_snr)
        mid_fer = evaluate(mid_snr)

        if lower_fer <= mid_fer <= upper_fer:
            return mid_snr, mid_fer

        if mid_fer > target_fer:
            low_snr = mid_snr
            low_fer = mid_fer
        else:
            high_snr = mid_snr
            high_fer = mid_fer

    if abs(low_fer - target_fer) < abs(high_fer - target_fer):
        return low_snr, low_fer

    return high_snr, high_fer


# 6th element of each tuple is the SPLIT LABEL (1/2/3/4 for
# split1/split2/split3/nosplit respectively), tracked explicitly from the
# splitting_pattern index used at construction time rather than re-derived
# from `save_name` after the fact. This fixes a real pre-existing bug
# (inherited faithfully from the original reference script,
# reproduce_fig5_scatter_plot_zc6_5G_LDPC.py, lines ~887-900 there) where the
# classifier did `if "split1" in save_name: split = 1 elif "split2" in
# save_name: ... elif "split3" in save_name: ...` -- since "subsplit1"
# trivially CONTAINS the substring "split1", every "..._subsplitN" variant
# whose save_name happens to end in "subsplit1" (e.g.
# "aSCED_12_split2_batch1_subsplit1", "aSCED_16_split3_batch1_subsplit1")
# was silently misclassified as split=1 and hit the FIRST branch regardless
# of its true pattern. Because the summary-file merge dedupes/overwrites by
# (split, numpaths), this caused a genuine key collision with the real
# split=1 variants of the same numpaths (4_split2_subsplit1 numpaths=4 collides with
# the real 4_split1_subsplit2; 8_split3_subsplit1 numpaths=8 collides with
# the real 8_split1_subsplit4) -- not a race, a deterministic
# misclassification that silently drops one of the two colliding rows.
experiments = [
    (flag_aed_asced_2_split1_subsplit1, "aSCED_2_split1_batch1_subsplit1", splitting_pattern[0], 1, 1, 1),
    (flag_aed_asced_4_split1_subsplit2, "aSCED_4_split1_batch1_subsplit2", splitting_pattern[0], 1, 2, 1),
    (flag_aed_asced_6_split1_subsplit3, "aSCED_6_split1_batch1_subsplit3", splitting_pattern[0], 1, 3, 1),
    (flag_aed_asced_8_split1_subsplit4, "aSCED_8_split1_batch1_subsplit4", splitting_pattern[0], 1, 4, 1),
    (flag_aed_asced_10_split1_subsplit5, "aSCED_10_split1_batch1_subsplit5", splitting_pattern[0], 1, 5, 1),
    (flag_aed_asced_12_split1, "aSCED_12_split1_batch1", splitting_pattern[0], 1, None, 1),
    (flag_aed_asced_24_split1, "aSCED_24_split1_batch2", splitting_pattern[0], 2, None, 1),
    (flag_aed_asced_36_split1, "aSCED_36_split1_batch3", splitting_pattern[0], 3, None, 1),
    (flag_aed_asced_48_split1, "aSCED_48_split1_batch4", splitting_pattern[0], 4, None, 1),
    (flag_aed_asced_60_split1, "aSCED_60_split1_batch5", splitting_pattern[0], 5, None, 1),
    (flag_aed_asced_72_split1, "aSCED_72_split1_batch6", splitting_pattern[0], 6, None, 1),
    (flag_aed_asced_84_split1, "aSCED_84_split1_batch7", splitting_pattern[0], 7, None, 1),
    (flag_aed_asced_96_split1, "aSCED_96_split1_batch8", splitting_pattern[0], 8, None, 1),
    (flag_asced_4_split2_subsplit1, "aSCED_12_split2_batch1_subsplit1", splitting_pattern[1], 1, 1, 2),
    (flag_asced_8_split2_subsplit2, "aSCED_12_split2_batch1_subsplit2", splitting_pattern[1], 1, 2, 2),
    (flag_asced_12_split2, "aSCED_12_split2_batch1", splitting_pattern[1], 1, None, 2),
    (flag_asced_24_split2, "aSCED_24_split2_batch2", splitting_pattern[1], 2, None, 2),
    (flag_asced_36_split2, "aSCED_36_split2_batch3", splitting_pattern[1], 3, None, 2),
    (flag_asced_48_split2, "aSCED_48_split2_batch4", splitting_pattern[1], 4, None, 2),
    (flag_asced_60_split2, "aSCED_60_split2_batch5", splitting_pattern[1], 5, None, 2),
    (flag_asced_72_split2, "aSCED_72_split2_batch6", splitting_pattern[1], 6, None, 2),
    (flag_asced_84_split2, "aSCED_84_split2_batch7", splitting_pattern[1], 7, None, 2),
    (flag_asced_96_split2, "aSCED_96_split2_batch8", splitting_pattern[1], 8, None, 2),
    (flag_asced_8_split3_subsplit1, "aSCED_16_split3_batch1_subsplit1", splitting_pattern[2], 1, 1, 3),
    (flag_asced_16_split3, "aSCED_16_split3_batch1", splitting_pattern[2], 1, None, 3),
    (flag_asced_32_split3, "aSCED_16_split3_batch2", splitting_pattern[2], 2, None, 3),
    (flag_asced_48_split3, "aSCED_48_split3_batch3", splitting_pattern[2], 3, None, 3),
    (flag_asced_64_split3, "aSCED_64_split3_batch4", splitting_pattern[2], 4, None, 3),
    (flag_asced_80_split3, "aSCED_80_split3_batch5", splitting_pattern[2], 5, None, 3),
    (flag_asced_96_split3, "aSCED_96_split3_batch6", splitting_pattern[2], 6, None, 3),
    (flag_asced_112_split3, "aSCED_112_split3_batch6", splitting_pattern[2], 7, None, 3),
    (flag_asced_128_split3, "aSCED_128_split3_batch8", splitting_pattern[2], 8, None, 3),
    (flag_asced_64_nosplit, "aSCED_64_nosplit_batch1", [], 1, None, 4),
    (flag_asced_128_nosplit, "aSCED_128_nosplit_batch2", [], 2, None, 4),
    (flag_asced_192_nosplit, "aSCED_192_nosplit_batch3", [], 3, None, 4),
]

summary = []

if flag_nmsa or flag_aed:
    raise ValueError("nmsa/aed are not part of the block-selection search; use the baseline (non-greedy) script for those.")

for enabled, save_name, split_pattern, num_blocks, subsplit, split in experiments:

    if not enabled:
        continue

    sim, num_paths = create_asced_sim(save_name, split_pattern, num_blocks, subsplit)

    required_snr, accepted_fer = search_target_fer(
        sim,
        snr_start,
        target_fer,
    )

    # split label comes directly from the experiments tuple now (see the
    # comment above the experiments list) -- NOT re-derived from save_name
    # substring matching, which was the source of the collision bug.
    subsplit_out = 1 if subsplit is not None else 0

    summary.append(
        {
            "split": split,
            "numblocks": num_blocks,
            "numpaths": num_paths,
            "requiredsnr": required_snr,
            "acceptedfer": accepted_fer,
            "subsplit": subsplit_out,
        }
    )

df = pd.DataFrame(summary)

summary_dir = os.path.join(results_dir_base, "summary")

os.makedirs(summary_dir, exist_ok=True)

# Parallel (not overwriting) filename - the "_greedy" suffix distinguishes
# this from the baseline (natural-order) n{n_simul}_paths_vs_required_snr.dat
# in the same summary/ folder, per the task's explicit layout.
summary_path = os.path.join(summary_dir, f"n{n_simul}_paths_vs_required_snr_greedy.dat")

new_df = pd.DataFrame(summary)

if os.path.exists(summary_path):
    old_df = pd.read_csv(summary_path, sep=" ")

    new_df_keys = new_df[["split", "numpaths"]]
    old_df = old_df.merge(
        new_df_keys,
        on=["split", "numpaths"],
        how="left",
        indicator=True,
    )

    old_df = old_df[old_df["_merge"] == "left_only"]
    old_df = old_df.drop(columns=["_merge"])

    new_df = pd.concat([old_df, new_df], ignore_index=True)

new_df = new_df.sort_values(["split", "numpaths"])

new_df.to_csv(
    summary_path,
    sep=" ",
    index=False,
)

print(new_df)
print()
print(f"Saved summary to {summary_path}")
