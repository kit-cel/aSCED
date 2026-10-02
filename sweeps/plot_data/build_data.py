import json, csv, math
from pathlib import Path

ROOT = Path("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/RESULTS")
OUT = Path("/tmp/claude-246141/-home-pj9034-aSCED--claude-worktrees-sim-orchestrator-6e0478/442a2e43-0413-4065-902f-1e00d80ea3b2/scratchpad")
SNR_GRID = ["1", "1.5", "2", "2.5", "3", "3.5", "4"]

import re
PATTERN = re.compile(
    r"^(?P<variant>asced\d+)_(?P<selector>full_parallel|fixed_sequential|random_sequential|syndrome_sequential)"
    r"(?:_mpg(?P<mpg>[\d.]+))?(?:_target(?P<target>[\d.]+))?$"
)

def load(p):
    try:
        with open(p) as f:
            return json.load(f)
    except FileNotFoundError:
        return {}

def get_config_data(cdir: Path):
    flat_fer = load(cdir / "FER.json")
    flat_stats = load(cdir / "stats.json")
    flat_dstats = load(cdir / "decoder_stats.json")
    out = {}
    for snr in SNR_GRID:
        src_fer, src_stats, src_dstats = flat_fer, flat_stats, flat_dstats
        snr_dir = cdir / f"snr_{snr}"
        if (snr_dir / "FER.json").exists():
            src_fer = load(snr_dir / "FER.json")
            src_stats = load(snr_dir / "stats.json")
            src_dstats = load(snr_dir / "decoder_stats.json")
        if snr not in src_fer:
            continue
        out[snr] = dict(
            fer=src_fer.get(snr),
            bit_errors=src_stats.get("bit_errors", {}).get(snr),
            frame_errors=src_stats.get("frame_errors", {}).get(snr),
            trials=src_stats.get("trials", {}).get(snr),
            effort=src_dstats.get("average_ensemble_effort", {}).get(snr),
            latency=src_dstats.get("average_ensemble_latency", {}).get(snr),
            num_converged=src_dstats.get("average_number_converged_path", {}).get(snr),
            num_active=src_dstats.get("average_number_active_path", {}).get(snr),
            threshold=src_dstats.get("average_threshold", {}).get(snr),
        )
    return out

configs = {}
for cdir in sorted(ROOT.glob("*_n132/*/")):
    name = cdir.name.rstrip("/")
    m = PATTERN.match(name)
    if not m:
        continue
    g = m.groupdict()
    if g["target"] == "0.5":
        continue
    data = get_config_data(cdir)
    if not data:
        continue
    key = (g["variant"], g["selector"], g["mpg"], g["target"])
    configs[key] = data

# ---------------------------------------------------------------- raw numbers
raw_rows = []
for (variant, selector, mpg, target), data in sorted(configs.items()):
    for snr, d in sorted(data.items(), key=lambda kv: float(kv[0])):
        raw_rows.append(dict(
            variant=variant, selector=selector,
            members_per_group=mpg if mpg else ("-" if selector != "full_parallel" else "all"),
            target_num_converged=target if target else "-",
            snr_db=snr, fer=d["fer"], bit_errors=d["bit_errors"], frame_errors=d["frame_errors"],
            trials=d["trials"], average_ensemble_effort=d["effort"], average_ensemble_latency=d["latency"],
            average_number_converged_path=d["num_converged"], average_number_active_path=d["num_active"],
            average_threshold=d["threshold"],
        ))

with open(OUT / "raw_numbers.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(raw_rows[0].keys()))
    w.writeheader()
    w.writerows(raw_rows)
json.dump(raw_rows, open(OUT / "raw_numbers.json", "w"), indent=1)
print(f"raw rows: {len(raw_rows)}")

# ---------------------------------------------------------------- literature (Fig. 3, fig3_literature_data.tex)
def pts(pairs):
    return [[x, y] for x, y in pairs]

literature = {
    "NMSA": pts([(1.0,5.467980e-01),(1.5,3.553223e-01),(2.0,1.798365e-01),(2.5,7.278189e-02),
                 (3.0,2.111038e-02),(3.5,4.332920e-03),(4.0,7.442976e-04),(4.5,1.148381e-04),(5.0,1.290430e-05)]),
    "NMSA (Imax=352)": pts([(1.0,5.141388e-01),(1.5,2.551020e-01),(2.0,1.264223e-01),(2.5,4.380201e-02),
                 (3.0,1.090275e-02),(3.5,2.144956e-03),(4.0,3.055329e-04),(4.5,4.053300e-05),(5.0,4.830237e-06)]),
    "AED-11": pts([(1.0,4.926108e-01),(1.5,2.998501e-01),(2.0,1.362398e-01),(2.5,4.621072e-02),
                 (3.0,1.192677e-02),(3.5,2.103359e-03),(4.0,2.896100e-04),(4.5,3.245119e-05),(5.0,3.464770e-06)]),
    "SCED-11": pts([(1.0,4.464286e-01),(1.5,2.743484e-01),(2.0,1.147447e-01),(2.5,4.145078e-02),
                 (3.0,8.008329e-03),(3.5,1.445713e-03),(4.0,1.989583e-04),(4.5,1.514168e-05),(5.0,1.494075e-06)]),
    "SCED-43": pts([(1.0,3.960396e-01),(1.5,2.242152e-01),(2.0,7.855460e-02),(2.5,2.734108e-02),
                 (3.0,5.499945e-03),(3.5,6.574903e-04),(4.0,6.791694e-05),(4.5,5.026574e-06)]),
    "aSCED-11": pts([(1.0,5.122951e-01),(1.5,2.467917e-01),(2.0,1.062473e-01),(2.5,3.351206e-02),
                 (3.0,7.658470e-03),(3.5,1.160389e-03),(4.0,1.392763e-04),(4.5,1.281535e-05),(5.0,1.007620e-06)]),
    "aSCED-48": pts([(1.0,0.3571035747021082),(1.5,0.18488151184881513),(2.0,0.06884949348769899),
                 (2.5,0.019854137385318563),(3.0,0.0036196964605847448),(3.5,0.000491052420269672),
                 (4.0,4.4328988197752285e-05),(4.5,2.965920104771095e-06),(5.0,1.3379296624288492e-07)]),
    "aSCED-384": pts([(1.0,0.23617896631524812),(1.5,0.10047980361526444),(2.0,0.030786870657493925),
                 (2.5,0.006419873173834153),(3.0,0.0009913115828252789),(3.5,9.581785784534925e-05),
                 (4.0,7.385968565865376e-06)]),
}
reference_osd4 = pts([(1.00,1.120e-01),(1.50,3.609e-02),(2.00,9.891e-03),(2.50,1.623e-03),(3.00,2.514e-04)])

json.dump({"literature": literature, "reference_osd4": reference_osd4}, open(OUT / "literature.json", "w"), indent=1)
print("literature curves:", list(literature.keys()))
