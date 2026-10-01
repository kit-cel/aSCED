"""Generate the parameter grid manifest for the sequential/mconverged RL aSCED
sweep (reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py), so the SLURM array
job is reproducible from the manifest alone.

Run this once before submitting sequential_rl_array.sbatch; re-run it if the
grid changes (the sbatch script always reads whatever is on disk).

target_num_converged values were picked from the first (target=0.5-fraction)
sweep's measured average_number_converged_path (3.5-10 for aSCED-48, 4.8-8.7
for aSCED-384 over 1.0-4.0 dB) - 0.5*ensemble_size (24 / 192) was never
reachable, so mconverged never triggered. 2 and 6 bracket the observed range:
2 should trigger stopping even at low SNR, 6 only once SNR is more favorable.
"""

import csv

VARIANTS = ["asced48", "asced384"]
SEQUENTIAL_SELECTORS = ["fixed_sequential", "syndrome_sequential"]
SNR_POINTS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
MEMBERS_PER_GROUP = 8
TARGET_NUM_CONVERGED_VALUES = [2, 6]
N_SIMUL = 132  # reproduces the exact C_5G(132,66) code used in Fig. 3 (n_simul IS
# the final transmitted length n directly, not an offset - n_simul=110 was wrong,
# it simulated the higher-rate C_5G(110,66) code instead; see STATUS.md)

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

    for selector in SEQUENTIAL_SELECTORS:
        for target in TARGET_NUM_CONVERGED_VALUES:
            for snr in SNR_POINTS:
                rows.append(
                    {
                        "variant": variant,
                        "n_simul": N_SIMUL,
                        "snr": snr,
                        "selector": selector,
                        "members_per_group": MEMBERS_PER_GROUP,
                        "target_num_converged": target,
                    }
                )

with open("sweeps/sequential_rl_manifest.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} tasks to sweeps/sequential_rl_manifest.csv")
