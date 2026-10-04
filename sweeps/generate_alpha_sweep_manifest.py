"""Generate the parameter grid manifest for the NMSA normalization-constant
(BP_config.norm_factor) sweep (sweeps/alpha_sweep.py), so the SLURM array job
is reproducible from the manifest alone.

Sweeps norm_factor over 9 dyadic values (each a sum of negative powers of two
- hardware-friendly fixed-point constants) across three decoders: standalone
NMSA, and aSCED-48/aSCED-384 full_parallel ensembles. Everything else is held
at the established defaults (scheduling_type="flooding", max_iterations=32,
early_stopping=True) - this sweep varies only norm_factor, no selector/
stopping-policy/scheduling changes (those are separate studies happening
elsewhere).

Run this once before submitting the sbatch array job; re-run it if the grid
changes (the sbatch script always reads whatever is on disk).
"""

import csv

VARIANTS = ["nmsa", "asced48", "asced384"]
NORM_FACTOR_VALUES = [0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.0]
SNR_POINTS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
N_SIMUL = 132  # reproduces the exact C_5G(132,66) code used in Fig. 3 - see
# sweeps/STATUS.md's "CRITICAL BUG" section (n_simul=110 silently simulates a
# different, wrong code).

rows = []
for variant in VARIANTS:
    for norm_factor in NORM_FACTOR_VALUES:
        for snr in SNR_POINTS:
            rows.append(
                {
                    "variant": variant,
                    "n_simul": N_SIMUL,
                    "snr": snr,
                    "norm_factor": norm_factor,
                }
            )

with open("sweeps/alpha_sweep_manifest.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} tasks to sweeps/alpha_sweep_manifest.csv")
