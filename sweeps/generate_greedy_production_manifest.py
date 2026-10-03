"""Generate the manifest for the production-quality validation sweep of the
GREEDY block-selection search's winning ensembles
(reproduce_fig3_RL_zc11_asced_5G_LDPC_greedy.py), comparing them against the
existing baseline "subsequent blocks" sweep
(RESULTS/fig_x_zc11_r4_seq_<selector>_mpg8_n132/{asced48,asced384}_*).

Structure mirrors sweeps/generate_sequential_rl_manifest.py exactly (same
grid), just a different results_dir prefix (results_dir is computed inside
the greedy script itself, not passed on the CLI) so it doesn't collide with
the baseline sweep's output directories.
"""

import csv

VARIANTS = ["asced48", "asced384"]
SEQUENTIAL_SELECTORS = ["fixed_sequential", "syndrome_sequential"]
SNR_POINTS = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
MEMBERS_PER_GROUP = 8
TARGET_NUM_CONVERGED_VALUES = [2, 6]
N_SIMUL = 132

rows = []
for variant in VARIANTS:
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

with open("sweeps/greedy_production_manifest.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} tasks to sweeps/greedy_production_manifest.csv")
