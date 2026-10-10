"""Generate the manifest for recalibrating `target_num_converged` at a fixed
`members_per_group=8`, for the aSCED-48/384 ensembles
(reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py).

Sweep A (sweeps/generate_sequential_rl_manifest.py) reused target_num_converged
values {2, 6} as-is from an earlier, wrong-code (n_simul=110) calibration,
explicitly flagged as "not re-derived yet" (see sweeps/STATUS.md). This script
re-derives a proper sweep from Sweep A's own measured
`average_number_converged_path` on the corrected (n_simul=132) code, read
directly from RESULTS/fig_x_zc11_r4_seq_full_parallel_mpg8_n132/
<variant>_full_parallel/snr_*/decoder_stats.json:

  asced48:  4.66 (1.0dB) -> 10.02 (4.0dB), saturating ~10 by 3.0dB
  asced384: 5.76 (1.0dB) ->  8.73 (3.5dB), slightly down to 8.68 at 4.0dB

TARGET_NUM_CONVERGED_VALUES brackets this range: well below (should trigger
stopping even at low SNR), within, and above (should rarely/never trigger,
converging toward full_parallel behavior) -- includes the original {2, 6}
for direct comparability with Sweep A's existing data.

Does NOT regenerate full_parallel rows (target-independent, already have
that data from Sweep A) -- only the two sequential selectors that actually
honor members_per_group/target_num_converged.
"""

import csv

VARIANTS = ["asced48", "asced384"]
SEQUENTIAL_SELECTORS = ["fixed_sequential", "syndrome_sequential"]
SNR_POINTS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
MEMBERS_PER_GROUP = 8
TARGET_NUM_CONVERGED_VALUES = [1, 2, 3, 4, 6, 8, 10, 12, 16, 24]
N_SIMUL = 132  # C_5G(132,66), the corrected code (n_simul IS the final
# transmitted length n directly, not an offset; see sweeps/STATUS.md)

rows = []
for variant in VARIANTS:
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

with open("sweeps/target_converged_recalibration_manifest.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} tasks to sweeps/target_converged_recalibration_manifest.csv")
