"""Generate the parameter grid manifest for the "PCM-first" aSCED-49/385
sweep (reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential_pcmfirst.py).

aSCED-49/385 = aSCED-48/384 + one extra plain-PCM decoding path prepended at
path index 0. Scoped to `fixed_sequential` only (plus `full_parallel` as the
FER/correctness baseline) since the "pin path 0 to the first group" guarantee
only holds for the static, index-range-based FixedSequentialSelector -
`syndrome_sequential` reorders dynamically per received word and is out of
scope here (would need a real C++ "pin" feature).

members_per_group=7 divides both new ensemble sizes evenly: 49 = 7x7,
385 = 7x55 (384+1) - so the PCM-first path is guaranteed in the very first
7-member group regardless of variant.

target_num_converged values {2, 6} reused as-is from the plain
aSCED-48/384 sweep's calibration (sweeps/generate_sequential_rl_manifest.py) -
same ensemble family/sizes (49≈48, 385≈384), so the same absolute targets
should be in a comparable, reachable range.
"""

import csv

VARIANTS = ["asced49", "asced385"]
SELECTOR = "fixed_sequential"
SNR_POINTS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
MEMBERS_PER_GROUP = 7
TARGET_NUM_CONVERGED_VALUES = [2, 6]
N_SIMUL = 132  # see STATUS.md's n_simul=110 CRITICAL BUG section - 132 is required

rows = []
for variant in VARIANTS:
    # full_parallel ignores target_num_converged entirely - one run per
    # (variant, snr) is enough, not one per target value.
    for snr in SNR_POINTS:
        rows.append(
            {
                "variant": variant,
                "n_simul": N_SIMUL,
                "snr": snr,
                "selector": "full_parallel",
                "members_per_group": MEMBERS_PER_GROUP,
                "target_num_converged": TARGET_NUM_CONVERGED_VALUES[0],  # unused
            }
        )

    for target in TARGET_NUM_CONVERGED_VALUES:
        for snr in SNR_POINTS:
            rows.append(
                {
                    "variant": variant,
                    "n_simul": N_SIMUL,
                    "snr": snr,
                    "selector": SELECTOR,
                    "members_per_group": MEMBERS_PER_GROUP,
                    "target_num_converged": target,
                }
            )

with open("sweeps/pcmfirst_manifest.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} tasks to sweeps/pcmfirst_manifest.csv")
