"""Generate the manifest for the production validation sweep of the Zc=6
GREEDY block-selection search (reproduce_fig5_scatter_plot_zc6_5G_LDPC_greedy.py),
mirroring reproduce_fig5_scatter_plot_zc6_5G_LDPC.py's `experiments` list
structure (split1 L=1..8 + 5 subsplit variants, split2 L=1..8, split3 L=1..8,
nosplit L=1..3) for both n_simul in {78, 180}, to compare against the existing
natural-order ("subsequent blocks") baseline at
RESULTS/fig_scatter_zc6_n={78,180}/summary/n{78,180}_paths_vs_required_snr.dat.

snr_start per code is a single reasonable starting guess (search_target_fer
brackets/binary-searches from there regardless of variant) informed by the
existing baseline data's mid-range required_snr for target_fer=1e-3.
"""

import csv

VARIANTS = [
    "2_split1_subsplit1",
    "4_split1_subsplit2",
    "6_split1_subsplit3",
    "8_split1_subsplit4",
    "10_split1_subsplit5",
    "12_split1",
    "24_split1",
    "36_split1",
    "48_split1",
    "60_split1",
    "72_split1",
    "84_split1",
    "96_split1",
    "4_split2_subsplit1",
    "8_split2_subsplit2",
    "12_split2",
    "24_split2",
    "36_split2",
    "48_split2",
    "60_split2",
    "72_split2",
    "84_split2",
    "96_split2",
    "8_split3_subsplit1",
    "16_split3",
    "32_split3",
    "48_split3",
    "64_split3",
    "80_split3",
    "96_split3",
    "112_split3",
    "128_split3",
    "64_nosplit",
    "128_nosplit",
    "192_nosplit",
]

assert len(VARIANTS) == 35

N_SIMUL_SNR_START = {78: 5.1, 180: 2.95}
TARGET_FER = 1e-3

rows = []
for n_simul, snr_start in N_SIMUL_SNR_START.items():
    for variant in VARIANTS:
        rows.append(
            {
                "decoder_variant": variant,
                "n_simul": n_simul,
                "snr_start": snr_start,
                "target_fer": TARGET_FER,
            }
        )

with open("sweeps/scatter_zc6_greedy_manifest.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} tasks to sweeps/scatter_zc6_greedy_manifest.csv")
