"""NMSA norm_factor (alpha) sweep: 9 dyadic alpha values x 3 variants (nmsa,
asced48, asced384 full_parallel) x 7 SNR. See sweeps/STATUS.md's "NMSA
norm_factor sweep" section. Headline: alpha=0.75 (the long-standing default)
is at/near the optimum for all three variants.

Root: RESULTS/fig_x_zc11_r4_alpha_sweep_n132/<variant>/alpha_<X>/<save_name>/snr_<snr>/
  save_name = "nmsa" for nmsa, "<variant>_full_parallel" for asced48/384.
  NOTE: alpha=1.0 is stored as "alpha_1" (Python %g formatting), not
  "alpha_1.0" -- confirmed via `ls`, not assumed.
"""
import json, csv, re
from pathlib import Path

ROOT = Path("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/RESULTS/fig_x_zc11_r4_alpha_sweep_n132")
OUT = Path(__file__).parent
SNR_GRID = ["1", "1.5", "2", "2.5", "3", "3.5", "4"]
ALPHAS = ["0.5", "0.5625", "0.625", "0.6875", "0.75", "0.8125", "0.875", "0.9375", "1"]
VARIANTS = ["nmsa", "asced48", "asced384"]

def load(p):
    try:
        with open(p) as f:
            return json.load(f)
    except FileNotFoundError:
        return {}

def save_name_for(variant):
    return "nmsa" if variant == "nmsa" else f"{variant}_full_parallel"

def get_config_data(cdir: Path):
    out = {}
    for snr in SNR_GRID:
        snr_dir = cdir / f"snr_{snr}"
        src_fer = load(snr_dir / "FER.json") if (snr_dir / "FER.json").exists() else load(cdir / "FER.json")
        src_stats = load(snr_dir / "stats.json") if (snr_dir / "stats.json").exists() else load(cdir / "stats.json")
        src_dstats = load(snr_dir / "decoder_stats.json") if (snr_dir / "decoder_stats.json").exists() else load(cdir / "decoder_stats.json")
        if snr not in src_fer:
            continue
        out[snr] = dict(fer=src_fer.get(snr), frame_errors=src_stats.get("frame_errors", {}).get(snr),
                         trials=src_stats.get("trials", {}).get(snr),
                         effort=src_dstats.get("average_ensemble_effort", {}).get(snr))
    return out

configs = {}  # (variant, alpha) -> {snr: {...}}
for variant in VARIANTS:
    for alpha in ALPHAS:
        cdir = ROOT / variant / f"alpha_{alpha}" / save_name_for(variant)
        data = get_config_data(cdir)
        if data:
            configs[(variant, alpha)] = data

print(f"configs found: {len(configs)} / {len(VARIANTS) * len(ALPHAS)} expected")
manifest_rows = list(csv.DictReader(open("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/sweeps/alpha_sweep_manifest.csv")))
print(f"alpha_sweep_manifest.csv: {len(manifest_rows)} rows (expect 189)")

raw_rows = []
for (variant, alpha), data in sorted(configs.items(), key=lambda kv: (kv[0][0], float(kv[0][1]))):
    for snr, d in sorted(data.items(), key=lambda kv: float(kv[0])):
        raw_rows.append(dict(variant=variant, norm_factor=alpha, snr_db=snr, fer=d["fer"],
                              frame_errors=d["frame_errors"], trials=d["trials"], average_ensemble_effort=d["effort"]))
with open(OUT / "alpha_raw.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(raw_rows[0].keys()))
    w.writeheader()
    w.writerows(raw_rows)
print(f"alpha raw rows: {len(raw_rows)} (expect 3 variants x 9 alphas x 7 SNR = 189)")

# ---- FER(alpha) curves, one line per SNR, per variant ----
def curve_for_variant(variant):
    lines = []
    for snr in SNR_GRID:
        pts = []
        for alpha in ALPHAS:
            d = configs.get((variant, alpha), {}).get(snr)
            if d is not None:
                pts.append([float(alpha), d["fer"], d["frame_errors"]])
        lines.append({"snr": float(snr), "pts": pts})
    return lines

alpha_data = {v: curve_for_variant(v) for v in VARIANTS}
json.dump({"alphas": [float(a) for a in ALPHAS], "variants": alpha_data}, open(OUT / "alpha_data.json", "w"), indent=1)

# ---- headline check: argmin FER over alpha, per variant per SNR ----
print("\nargmin-FER alpha per (variant, SNR) -- headline: should cluster near 0.75")
best_counts = {}
for variant in VARIANTS:
    for snr in SNR_GRID:
        best = None
        for alpha in ALPHAS:
            d = configs.get((variant, alpha), {}).get(snr)
            if d is None:
                continue
            if best is None or d["fer"] < best[1]:
                best = (alpha, d["fer"])
        if best:
            print(f"  {variant} @ {snr}dB: best alpha={best[0]} (FER={best[1]:.4g})")
            best_counts[best[0]] = best_counts.get(best[0], 0) + 1
print("best-alpha histogram:", dict(sorted(best_counts.items(), key=lambda kv: -kv[1])))
