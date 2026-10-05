"""PCM-first variants (aSCED-49/385) vs. plain aSCED-48/384.

De-confounded (job 535646 complete, 28/28): three groups are shown so the
PCM-first effect can be read off cleanly, separate from the group-size
effect (members_per_group=7 was required for the PCM-first variants so the
prepended plain-PCM path at index 0 divides evenly: 49=7x7, 385=7x55):
  - "baseline (mpg=8)": the original asced48/384, for reference.
  - "no PCM-first (mpg=7 control)": plain asced48/384 at mpg=7, job 535646 -
    isolates the pure group-size effect (mpg=8 -> mpg=7) alone.
  - "PCM-first (mpg=7)": asced49/385, job 534062 - the PCM-first effect
    alone is this curve vs. the mpg=7 control (both at mpg=7), NOT vs. the
    mpg=8 baseline.

Baseline root: RESULTS/fig_x_zc11_r4_seq_<selector>_mpg8_n132/ (asced48/384)
Control root:  RESULTS/fig_x_zc11_r4_seq_<selector>_mpg7_n132/ (asced48/384, job 535646)
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
CONTROL_RE = re.compile(rf"^fig_x_zc11_r4_seq_(?:{SELECTORS})_mpg7_n132$")
PCMFIRST_RE = re.compile(rf"^fig_x_zc11_r4_seq_pcmfirst_(?:{SELECTORS})_mpg7_n132$")
all_roots = sorted(ROOT.glob("fig_x_zc11_r4_seq_*"))
baseline = load_tree([r for r in all_roots if BASE_RE.match(r.name)])
control = load_tree([r for r in all_roots if CONTROL_RE.match(r.name)])
pcmfirst = load_tree([r for r in all_roots if PCMFIRST_RE.match(r.name)])

print("baseline configs:", sorted(baseline.keys()))
print("control (mpg=7, no PCM) configs:", sorted(control.keys()))
print("pcmfirst configs:", sorted(pcmfirst.keys()))

manifest_rows = list(csv.DictReader(open("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/sweeps/pcmfirst_manifest.csv")))
print(f"pcmfirst_manifest.csv: {len(manifest_rows)} rows (expect 42)")
control_manifest_rows = list(csv.DictReader(open("/home/pj9034/aSCED/.claude/worktrees/sim-orchestrator-6e0478/sweeps/mpg7_control_manifest.csv")))
print(f"mpg7_control_manifest.csv: {len(control_manifest_rows)} rows (expect 28)")

def pts4(data_dict, snr_list):
    return [[float(s), data_dict[s]["fer"], data_dict[s]["frame_errors"], data_dict[s]["effort"]] for s in snr_list if s in data_dict]

# family: ("48", baseline variant, pcmfirst variant), ("384", ...)
FAMILIES = [("48", "asced48", "asced49"), ("384", "asced384", "asced385")]

GROUPS = (
    ("baseline (mpg=8)", baseline, "8", "baseline", False),
    ("no PCM-first (mpg=7 control)", control, "7", "mpg7 control", False),
    ("PCM-first (mpg=7)", pcmfirst, "7", "PCM-first", True),
)

def setting(base_variant, pcm_variant):
    curves = []
    for grp_name, cfgs, mpg, tag, is_pcmfirst in GROUPS:
        variant = pcm_variant if is_pcmfirst else base_variant
        fp = cfgs.get((variant, "full_parallel", None, None), {})
        if fp:
            curves.append({"name": f"{tag} · full parallel", "group": grp_name, "visible": True, "pts": pts4(fp, SNR_GRID)})
        for target, vis in (("2", True), ("6", False)):
            fx = cfgs.get((variant, "fixed_sequential", mpg, target), {})
            if fx:
                curves.append({"name": f"{tag} · fixed order · target={target}", "group": grp_name,
                                "visible": vis, "pts": pts4(fx, SNR_GRID)})
    return {"xlabel": "Eb/N0 [dB]", "xlog": False, "ylabel": "FER", "curves": curves}

pcmfirst_data = {f"K{k}": setting(base_v, pcm_v) for k, base_v, pcm_v in FAMILIES}
json.dump(pcmfirst_data, open(OUT / "pcmfirst_data.json", "w"), indent=1)

raw_rows = []
for source, cfgs in (("baseline_mpg8", baseline), ("control_mpg7_no_pcm", control), ("pcmfirst_mpg7", pcmfirst)):
    for (variant, selector, mpg, target), data in sorted(cfgs.items()):
        for snr, d in sorted(data.items(), key=lambda kv: float(kv[0])):
            raw_rows.append(dict(source=source, variant=variant, selector=selector, members_per_group=mpg,
                                  target_num_converged=target or "-", snr_db=snr, fer=d["fer"],
                                  frame_errors=d["frame_errors"], trials=d["trials"],
                                  average_ensemble_effort=d["effort"], average_ensemble_latency=d["latency"]))
with open(OUT / "pcmfirst_raw.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(raw_rows[0].keys()))
    w.writeheader()
    w.writerows(raw_rows)
print(f"pcmfirst raw rows: {len(raw_rows)}")
