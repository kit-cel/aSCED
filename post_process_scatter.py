#!/usr/bin/env python3

import json
import os
import pandas as pd

##############################################################################
# Configuration
##############################################################################

results_dir = "RESULTS/fig_scatter_zc6_n=78"

target_fer = 1e-3

# Change this if your FER file has another name
fer_filename = "FER.json"

##############################################################################

summary = []

# Number of paths for every experiment
num_paths_lookup = {
    "aSCED_12_split0_batch1": 12,
    "aSCED_24_split0_batch2": 24,
    "aSCED_48_split0_batch4": 48,
    "aSCED_12_split1_batch1": 12,
    "aSCED_24_split1_batch2": 24,
    "aSCED_48_split1_batch4": 48,
    "aSCED_16_split2_batch1": 16,
    "aSCED_48_split2_batch3": 48,
    "aSCED_64_split2_batch4": 64,
    "aSCED_64_nosplit_batch1": 64,
    "aSCED_128_nosplit_batch2": 128,
    "nmsa": 1,
}

for decoder in sorted(os.listdir(results_dir)):

    decoder_dir = os.path.join(results_dir, decoder)

    if not os.path.isdir(decoder_dir):
        continue

    fer_path = os.path.join(decoder_dir, fer_filename)

    if not os.path.exists(fer_path):
        print(f"Skipping {decoder}: no FER file found.")
        continue

    with open(fer_path, "r") as f:
        fer_dict = json.load(f)

    # Convert SNR keys to float
    fer_data = [(float(snr), fer) for snr, fer in fer_dict.items()]

    # Find FER closest to target
    required_snr, accepted_fer = min(
        fer_data,
        key=lambda x: abs(x[1] - target_fer),
    )

    summary.append(
        {
            "decoder": decoder,
            "num_paths": num_paths_lookup.get(decoder, -1),
            "required_snr": required_snr,
            "accepted_fer": accepted_fer,
            "target_fer": target_fer,
        }
    )

summary_df = pd.DataFrame(summary)

summary_df = summary_df.sort_values("num_paths")

summary_dir = os.path.join(results_dir, "summary")
os.makedirs(summary_dir, exist_ok=True)

summary_path = os.path.join(summary_dir, "paths_vs_required_snr.csv")

summary_df.to_csv(summary_path, index=False)

print(summary_df)
print()
print(f"Saved summary to {summary_path}")