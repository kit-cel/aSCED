"""Generate the parameter grid manifest for the sequential/mconverged RL aSCED
sweep (reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py), so the SLURM array
job is reproducible from the manifest alone.

Run this once before submitting sequential_rl_array.sbatch; re-run it if the
grid changes (the sbatch script always reads whatever is on disk).
"""

import csv

VARIANTS = ["asced48", "asced384"]
SELECTORS = ["full_parallel", "fixed_sequential", "syndrome_sequential"]
SNR_POINTS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
MEMBERS_PER_GROUP = 8
TARGET_FRACTION_CONVERGED = 0.5
N_SIMUL = 110  # reproduces the exact C_5G(132,66) code used in Fig. 3

rows = []
for variant in VARIANTS:
    for selector in SELECTORS:
        for snr in SNR_POINTS:
            rows.append(
                {
                    "variant": variant,
                    "n_simul": N_SIMUL,
                    "snr": snr,
                    "selector": selector,
                    "members_per_group": MEMBERS_PER_GROUP,
                    "target_fraction_converged": TARGET_FRACTION_CONVERGED,
                }
            )

with open("sweeps/sequential_rl_manifest.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} tasks to sweeps/sequential_rl_manifest.csv")
