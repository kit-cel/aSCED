"""Scheduling / iteration-count study: BP scheduling_type {flooding,
row_layered_natural, row_layered_appended_first} x max_iterations {32, 16},
for 3 variants (nmsa, asced48, asced384), full_parallel only.
See sweeps/STATUS.md's "Scheduling/iteration-count study" section.

Root: RESULTS/scheduling_study_zc11_r4_n132/<variant>_<scheduling_config>_maxiter<N>/snr_<X>/
"""
import json, csv, re
from pathlib import Path

ROOT = Path("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/RESULTS/scheduling_study_zc11_r4_n132")
OUT = Path(__file__).parent
SNR_GRID = ["1", "1.5", "2", "2.5", "3", "3.5", "4"]

PATTERN = re.compile(
    r"^(?P<variant>nmsa|asced48|asced384)_(?P<config>flooding|row_layered_natural|row_layered_appended_first)"
    r"_maxiter(?P<maxiter>\d+)$"
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
        out[snr] = dict(fer=src_fer.get(snr), frame_errors=src_stats.get("frame_errors", {}).get(snr),
                         trials=src_stats.get("trials", {}).get(snr),
                         effort=src_dstats.get("average_ensemble_effort", {}).get(snr),
                         latency=src_dstats.get("average_ensemble_latency", {}).get(snr))
    return out

configs = {}
for cdir in sorted(ROOT.glob("*/")):
    m = PATTERN.match(cdir.name.rstrip("/"))
    if not m:
        continue
    g = m.groupdict()
    data = get_config_data(cdir)
    if not data:
        continue
    configs[(g["variant"], g["config"], g["maxiter"])] = data

print("configs found:", len(configs), "/ 18 expected")
manifest_rows = list(csv.DictReader(open("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/sweeps/scheduling_manifest.csv")))
print(f"scheduling_manifest.csv: {len(manifest_rows)} rows (expect 126)")

def pts4(data_dict, snr_list):
    return [[float(s), data_dict[s]["fer"], data_dict[s]["frame_errors"], data_dict[s]["effort"]] for s in snr_list if s in data_dict]

SCHED_CONFIGS = ["flooding", "row_layered_natural", "row_layered_appended_first"]
GROUP_LABEL = {"flooding": "flooding", "row_layered_natural": "row-layered (natural)",
               "row_layered_appended_first": "row-layered (appended-first)"}

def setting(variant):
    curves = []
    for sc in SCHED_CONFIGS:
        for maxiter in ("32", "16"):  # order matters: 32 first => solid, 16 second => dash (see curveTraces dash cycling)
            data = configs.get((variant, sc, maxiter), {})
            curves.append({"name": f"{GROUP_LABEL[sc]} · maxiter={maxiter}", "group": GROUP_LABEL[sc],
                            "visible": True, "pts": pts4(data, SNR_GRID)})
    return {"xlabel": "Eb/N0 [dB]", "xlog": False, "ylabel": "FER", "curves": curves}

scheduling_data = {v: setting(v) for v in ("nmsa", "asced48", "asced384")}
json.dump(scheduling_data, open(OUT / "scheduling_data.json", "w"), indent=1)

raw_rows = []
for (variant, sc, maxiter), data in sorted(configs.items()):
    for snr, d in sorted(data.items(), key=lambda kv: float(kv[0])):
        raw_rows.append(dict(variant=variant, scheduling_config=sc, max_iterations=maxiter, snr_db=snr,
                              fer=d["fer"], frame_errors=d["frame_errors"], trials=d["trials"],
                              average_ensemble_effort=d["effort"], average_ensemble_latency=d["latency"]))
with open(OUT / "scheduling_raw.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(raw_rows[0].keys()))
    w.writeheader()
    w.writerows(raw_rows)
print(f"scheduling raw rows: {len(raw_rows)} (expect 18 configs x 7 SNR = 126)")

# headline check: row_layered_natural@16 vs flooding@32, nmsa, a few SNRs
for snr in ("1.5", "2", "2.5"):
    fl32 = configs.get(("nmsa", "flooding", "32"), {}).get(snr, {}).get("fer")
    rl16 = configs.get(("nmsa", "row_layered_natural", "16"), {}).get(snr, {}).get("fer")
    print(f"  nmsa @ {snr}dB: flooding/32it={fl32}  row_layered_natural/16it={rl16}")
