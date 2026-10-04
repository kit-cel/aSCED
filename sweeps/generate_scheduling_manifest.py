"""Generate the parameter grid manifest for the pure scheduling/iteration-count
study (sweeps/scheduling_study.py), so the SLURM array job is reproducible
from the manifest alone.

Grid: 3 variants {nmsa, asced48, asced384} x 3 scheduling configs
{flooding, row_layered_natural, row_layered_appended_first} x 2 max_iterations
{32, 16} x 7 SNR points (1.0-4.0dB step 0.5) = 126 tasks.

Run this once before submitting scheduling_study_array.sbatch; re-run it if
the grid changes (the sbatch script always reads whatever is on disk).
"""

import csv

VARIANTS = ["nmsa", "asced48", "asced384"]
SCHEDULING_CONFIGS = ["flooding", "row_layered_natural", "row_layered_appended_first"]
MAX_ITERATIONS_VALUES = [32, 16]
SNR_POINTS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
N_SIMUL = 132  # n_simul IS the final transmitted length n directly, not an
# offset - n_simul=110 was wrong (simulates the higher-rate C_5G(110,66) code
# instead of Fig. 3's C_5G(132,66)); see sweeps/STATUS.md's "CRITICAL BUG" section.

rows = []
for variant in VARIANTS:
    for scheduling_config in SCHEDULING_CONFIGS:
        for max_iterations in MAX_ITERATIONS_VALUES:
            for snr in SNR_POINTS:
                rows.append(
                    {
                        "variant": variant,
                        "n_simul": N_SIMUL,
                        "snr": snr,
                        "scheduling_config": scheduling_config,
                        "max_iterations": max_iterations,
                    }
                )

with open("sweeps/scheduling_manifest.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} tasks to sweeps/scheduling_manifest.csv")
