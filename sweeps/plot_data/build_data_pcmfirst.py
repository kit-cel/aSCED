"""PCM-first variants (aSCED-49/385) vs. plain aSCED-48/384 baseline.

**Confounded as of this build** (see sweeps/STATUS.md): the production pcmfirst
sweep used members_per_group=7 (required so the prepended plain-PCM path at
index 0 divides evenly: 49=7x7, 385=7x55), while the existing baseline uses
members_per_group=8. Job 535646 (the mpg=7, no-PCM-first de-confounding
control) was launched to isolate the two effects but had not finished at
build time (`sacct -j 535646` showed 16/28 COMPLETED, 10 RUNNING, 1 PENDING)
-- so this script emits the confounded comparison only, and build_html.py
must show a clearly visible note about it. Re-run this script (and
build_html.py) once job 535646 completes to add the clean, de-confounded
comparison.

Baseline root:  RESULTS/fig_x_zc11_r4_seq_<selector>_mpg8_n132/ (asced48/384)
PCM-first root: RESULTS/fig_x_zc11_r4_seq_pcmfirst_<selector>_mpg7_n132/ (asced49/385)
"""
import json, csv, re
from pathlib import Path

ROOT = Path("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/RESULTS")
OUT = Path(__file__).parent
SNR_GRID = ["1", "1.5", "2", "2.5", "3", "3.5", "4"]

PATTERN = re.compile(
    r"^(?P<variant>asced\d+)_(?P<selector>full_parallel|fixed_sequential)"
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
        out[snr] = dict(fer=src_fer.get(snr), frame_errors=src_stats.get("frame_errors", {}).get(snr),
                         trials=src_stats.get("trials", {}).get(snr),
                         effort=src_dstats.get("average_ensemble_effort", {}).get(snr),
                         latency=src_dstats.get("average_ensemble_latency", {}).get(snr))
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
            configs[(g["variant"], g["selector"], g["mpg"], g["target"])] = data
    return configs

SELECTORS = "full_parallel|fixed_sequential"
BASE_RE = re.compile(rf"^fig_x_zc11_r4_seq_(?:{SELECTORS})_mpg8_n132$")
PCMFIRST_RE = re.compile(rf"^fig_x_zc11_r4_seq_pcmfirst_(?:{SELECTORS})_mpg7_n132$")
all_roots = sorted(ROOT.glob("fig_x_zc11_r4_seq_*"))
baseline = load_tree([r for r in all_roots if BASE_RE.match(r.name)])
pcmfirst = load_tree([r for r in all_roots if PCMFIRST_RE.match(r.name)])

print("baseline configs:", sorted(baseline.keys()))
print("pcmfirst configs:", sorted(pcmfirst.keys()))

manifest_rows = list(csv.DictReader(open("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/sweeps/pcmfirst_manifest.csv")))
print(f"pcmfirst_manifest.csv: {len(manifest_rows)} rows (expect 42)")

def pts4(data_dict, snr_list):
    return [[float(s), data_dict[s]["fer"], data_dict[s]["frame_errors"], data_dict[s]["effort"]] for s in snr_list if s in data_dict]

# family: ("48", baseline variant, pcmfirst variant), ("384", ...)
FAMILIES = [("48", "asced48", "asced49"), ("384", "asced384", "asced385")]

def setting(base_variant, pcm_variant):
    curves = []
    for grp_name, cfgs, variant in (("No PCM-first (mpg=8, baseline)", baseline, base_variant),
                                     ("PCM-first (mpg=7)", pcmfirst, pcm_variant)):
        tag = "baseline" if "baseline" in grp_name else "pcmfirst"
        fp = cfgs.get((variant, "full_parallel", None, None), {})
        curves.append({"name": f"{tag} · full parallel", "group": grp_name, "visible": True, "pts": pts4(fp, SNR_GRID)})
        for target, vis in (("2", True), ("6", False)):
            fx = cfgs.get((variant, "fixed_sequential", "8" if tag == "baseline" else "7", target), {})
            curves.append({"name": f"{tag} · fixed order · target={target}", "group": grp_name,
                            "visible": vis, "pts": pts4(fx, SNR_GRID)})
    return {"xlabel": "Eb/N0 [dB]", "xlog": False, "ylabel": "FER", "curves": curves}

pcmfirst_data = {f"K{k}": setting(base_v, pcm_v) for k, base_v, pcm_v in FAMILIES}
json.dump(pcmfirst_data, open(OUT / "pcmfirst_data.json", "w"), indent=1)

raw_rows = []
for grp_name, cfgs in (("baseline_mpg8", baseline), ("pcmfirst_mpg7", pcmfirst)):
    for (variant, selector, mpg, target), data in sorted(cfgs.items()):
        for snr, d in sorted(data.items(), key=lambda kv: float(kv[0])):
            raw_rows.append(dict(source=grp_name, variant=variant, selector=selector, members_per_group=mpg,
                                  target_num_converged=target or "-", snr_db=snr, fer=d["fer"],
                                  frame_errors=d["frame_errors"], trials=d["trials"],
                                  average_ensemble_effort=d["effort"], average_ensemble_latency=d["latency"]))
with open(OUT / "pcmfirst_raw.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(raw_rows[0].keys()))
    w.writeheader()
    w.writerows(raw_rows)
print(f"pcmfirst raw rows: {len(raw_rows)} (expect 2 families x (1 full_parallel + 2 fixed targets) x 7 SNR x 2 sources = 84)")
