"""Manifest for the two follow-up sweeps on asced48:

1. Parallelism (group-size) trade-off: members_per_group swept from 1 (full
   sequential, best complexity) up to the ensemble size (full parallel, best
   latency), ordering fixed to syndrome_sequential, target_num_converged
   fixed to 6 so group size is the only varying axis.
2. Ordering comparison: members_per_group fixed to 4, selector swept over
   fixed_sequential/random_sequential (syndrome_sequential@mpg=4 is already
   covered by sweep 1 and is not duplicated here).

Both reuse the same SNR grid as sequential_rl_manifest.csv (1.0-4.0dB), which
already brackets both FER=1e-1 (~2.15dB) and FER=1e-3 (~3.7dB) for asced48,
per the completed full_parallel baseline.

mpg=48 (full parallel) is NOT included here - it's equivalent to the already
-completed asced48_full_parallel run, which is reused directly instead of
rerunning it.
"""

import csv

VARIANT = "asced48"
N_SIMUL = 132  # n_simul IS the final transmitted length n directly; 110 was wrong
# (simulated the higher-rate C_5G(110,66) code instead of C_5G(132,66)), see STATUS.md
SNR_POINTS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
TARGET_NUM_CONVERGED = 6

rows = []

# Sweep 1: parallelism trade-off (group size), syndrome ordering fixed.
for members_per_group in [1, 2, 4, 8, 16, 24]:
    for snr in SNR_POINTS:
        rows.append(
            {
                "variant": VARIANT,
                "n_simul": N_SIMUL,
                "snr": snr,
                "selector": "syndrome_sequential",
                "members_per_group": members_per_group,
                "target_num_converged": TARGET_NUM_CONVERGED,
            }
        )

# Sweep 2: ordering comparison, group size fixed to 4.
for selector in ["fixed_sequential", "random_sequential"]:
    for snr in SNR_POINTS:
        rows.append(
            {
                "variant": VARIANT,
                "n_simul": N_SIMUL,
                "snr": snr,
                "selector": selector,
                "members_per_group": 4,
                "target_num_converged": TARGET_NUM_CONVERGED,
            }
        )

with open("sweeps/parallelism_ordering_manifest.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} tasks to sweeps/parallelism_ordering_manifest.csv")
