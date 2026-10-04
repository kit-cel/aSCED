"""Greedy RL/QC-block selection: 'baseline (subsequent blocks)' vs 'greedy blocks'
production comparison for aSCED-48/384, across full_parallel / fixed_sequential /
syndrome_sequential. See sweeps/greedy_block_search_results.md for the search
methodology and sweeps/STATUS.md's "Greedy RL/QC-block selection search" section.

Baseline root: RESULTS/fig_x_zc11_r4_seq_<selector>_mpg8_n132/
Greedy root:   RESULTS/fig_x_zc11_r4_seq_greedy_<selector>_mpg8_n132/
Both share the exact same leaf-directory naming convention
(<variant>_<selector>[_mpg8_target<N>]), which is why build_data.py's baseline
extraction has to explicitly guard against the "greedy" root tag -- see the
comment there. This script reads both root trees independently and explicitly,
so there is no risk of the two colliding here.
"""
import json, csv, re
from pathlib import Path

ROOT = Path("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/RESULTS")
OUT = Path(__file__).parent
SNR_GRID = ["1", "1.5", "2", "2.5", "3", "3.5", "4"]

PATTERN = re.compile(
    r"^(?P<variant>asced\d+)_(?P<selector>full_parallel|fixed_sequential|syndrome_sequential)"
    r"(?:_mpg(?P<mpg>[\d.]+))?(?:_target(?P<target>[\d.]+))?$"
)

def load(p):
    try:
        with open(p) as f:
            return json.load(f)
    except FileNotFoundError:
        return {}

def get_config_data(cdir: Path):
    out = {}
    for snr in SNR_GRID:
        snr_dir = cdir / f"snr_{snr}"
        src_fer = load(snr_dir / "FER.json") if (snr_dir / "FER.json").exists() else load(cdir / "FER.json")
        src_stats = load(snr_dir / "stats.json") if (snr_dir / "stats.json").exists() else load(cdir / "stats.json")
        src_dstats = load(snr_dir / "decoder_stats.json") if (snr_dir / "decoder_stats.json").exists() else load(cdir / "decoder_stats.json")
        if snr not in src_fer:
            continue
        out[snr] = dict(
            fer=src_fer.get(snr),
            frame_errors=src_stats.get("frame_errors", {}).get(snr),
            trials=src_stats.get("trials", {}).get(snr),
            effort=src_dstats.get("average_ensemble_effort", {}).get(snr),
            latency=src_dstats.get("average_ensemble_latency", {}).get(snr),
        )
    return out

def load_tree(roots):
    configs = {}
    for root in roots:
        for cdir in sorted(root.glob("*/")):
            m = PATTERN.match(cdir.name.rstrip("/"))
            if not m:
                continue
            g = m.groupdict()
            data = get_config_data(cdir)
            if not data:
                continue
            key = (g["variant"], g["selector"], g["mpg"], g["target"])
            configs[key] = data
    return configs

SELECTORS = "full_parallel|fixed_sequential|random_sequential|syndrome_sequential"
BASELINE_ROOT_RE = re.compile(rf"^fig_x_zc11_r4_seq_(?:{SELECTORS})_mpg8_n132$")
GREEDY_ROOT_RE = re.compile(rf"^fig_x_zc11_r4_seq_greedy_(?:{SELECTORS})_mpg8_n132$")
all_roots = sorted(ROOT.glob("fig_x_zc11_r4_seq_*_mpg8_n132"))
baseline = load_tree([r for r in all_roots if BASELINE_ROOT_RE.match(r.name)])
greedy = load_tree([r for r in all_roots if GREEDY_ROOT_RE.match(r.name)])

print("baseline configs:", sorted(baseline.keys()))
print("greedy configs:  ", sorted(greedy.keys()))

manifest_rows = list(csv.DictReader(open("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/sweeps/greedy_production_manifest.csv")))
print(f"greedy_production_manifest.csv: {len(manifest_rows)} rows")
expected = set()
for r in manifest_rows:
    target = r["target_num_converged"] if r["selector"] != "full_parallel" else None
    expected.add((r["variant"], r["selector"], r["members_per_group"], target, r["snr"].rstrip("0").rstrip(".") if "." in r["snr"] else r["snr"]))

def pts4(data_dict, snr_list):
    return [[float(s), data_dict[s]["fer"], data_dict[s]["frame_errors"], data_dict[s]["effort"]] for s in snr_list if s in data_dict]

def rows_for(cfgs, variant, selector, target):
    key = (variant, selector, "8", target)
    return cfgs.get(key, {})

def setting(variant):
    curves = []
    for grp_name, cfgs in (("Baseline blocks", baseline), ("Greedy blocks", greedy)):
        fp = rows_for(cfgs, variant, "full_parallel", None)
        curves.append({"name": f"{grp_name.split()[0].lower()} · full parallel", "group": grp_name,
                        "visible": True, "pts": pts4(fp, SNR_GRID)})
        for target, vis in (("2", True), ("6", False)):
            fx = rows_for(cfgs, variant, "fixed_sequential", target)
            curves.append({"name": f"{grp_name.split()[0].lower()} · fixed order · target={target}",
                            "group": grp_name, "visible": vis, "pts": pts4(fx, SNR_GRID)})
            sy = rows_for(cfgs, variant, "syndrome_sequential", target)
            curves.append({"name": f"{grp_name.split()[0].lower()} · syndrome order · target={target}",
                            "group": grp_name, "visible": vis, "pts": pts4(sy, SNR_GRID)})
    return {"xlabel": "Eb/N0 [dB]", "xlog": False, "ylabel": "FER", "curves": curves}

greedy_data = {"asced48": setting("asced48"), "asced384": setting("asced384")}
json.dump(greedy_data, open(OUT / "greedy_data.json", "w"), indent=1)

raw_rows = []
for grp_name, cfgs in (("baseline", baseline), ("greedy", greedy)):
    for (variant, selector, mpg, target), data in sorted(cfgs.items()):
        for snr, d in sorted(data.items(), key=lambda kv: float(kv[0])):
            raw_rows.append(dict(source=grp_name, variant=variant, selector=selector, members_per_group=mpg,
                                  target_num_converged=target or "-", snr_db=snr, fer=d["fer"],
                                  frame_errors=d["frame_errors"], trials=d["trials"],
                                  average_ensemble_effort=d["effort"], average_ensemble_latency=d["latency"]))
with open(OUT / "greedy_raw.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(raw_rows[0].keys()))
    w.writeheader()
    w.writerows(raw_rows)
print(f"greedy raw rows: {len(raw_rows)} (expect 2 variants x 2 sources x (1 full_parallel + 2 fixed targets + 2 syndrome targets) x 7 SNR = 140)")
