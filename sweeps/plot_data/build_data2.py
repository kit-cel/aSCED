import json, math
from pathlib import Path

OUT = Path("/tmp/claude-246141/-home-pj9034-aSCED--claude-worktrees-sim-orchestrator-6e0478/442a2e43-0413-4065-902f-1e00d80ea3b2/scratchpad")
raw = json.load(open(OUT / "raw_numbers.json"))
lit = json.load(open(OUT / "literature.json"))

def rows_for(variant, selector, mpg, target):
    out = [r for r in raw if r["variant"] == variant and r["selector"] == selector
           and str(r["members_per_group"]) == str(mpg) and str(r["target_num_converged"]) == str(target)]
    return sorted(out, key=lambda r: float(r["snr_db"]))

def pts4(rows):  # [snr, fer, frame_errors, effort]  (effort used in place of cpu_seconds_total)
    return [[float(r["snr_db"]), r["fer"], r["frame_errors"], r["average_ensemble_effort"]] for r in rows]

# ===================================================================== Sweep A
def setting(variant):
    fp = rows_for(variant, "full_parallel", "all", "-")
    fx2 = rows_for(variant, "fixed_sequential", "8", "2")
    fx6 = rows_for(variant, "fixed_sequential", "8", "6")
    sy2 = rows_for(variant, "syndrome_sequential", "8", "2")
    sy6 = rows_for(variant, "syndrome_sequential", "8", "6")
    curves = [
        {"name": "full parallel (K=%s)" % variant[5:], "group": "Baseline", "visible": True, "pts": pts4(fp)},
        {"name": "fixed order · mpg=8 · target=2", "group": "Fixed order", "visible": True, "pts": pts4(fx2)},
        {"name": "fixed order · mpg=8 · target=6", "group": "Fixed order", "visible": False, "pts": pts4(fx6)},
        {"name": "syndrome order · mpg=8 · target=2", "group": "Syndrome order", "visible": True, "pts": pts4(sy2)},
        {"name": "syndrome order · mpg=8 · target=6", "group": "Syndrome order", "visible": False, "pts": pts4(sy6)},
    ]
    return {
        "xlabel": "Eb/N0 [dB]", "xlog": False, "ylabel": "FER",
        "curves": curves,
        "reference": {"name": "OSD-4, non-ML", "pts": lit["reference_osd4"]},
        "paper": lit["literature"],
    }

sweepA = {"asced48": setting("asced48"), "asced384": setting("asced384")}
json.dump(sweepA, open(OUT / "sweepA_data.json", "w"), indent=1)
print("Sweep A settings:", list(sweepA.keys()))
for k, v in sweepA.items():
    for c in v["curves"]:
        print(" ", k, c["name"], len(c["pts"]), "pts")

# ===================================================================== interpolation helpers
def interp_at_fer(rows, target_fer):
    """Find SNR where FER == target_fer (log-linear interp, bracket from grid),
    then return (snr*, effort*, latency*, min_frame_errors_in_bracket) or None if out of range."""
    xs = [float(r["snr_db"]) for r in rows]
    ys = [math.log10(r["fer"]) for r in rows]
    tlog = math.log10(target_fer)
    # ys is decreasing in xs (FER falls with SNR); find bracket
    for i in range(len(xs) - 1):
        y0, y1 = ys[i], ys[i + 1]
        if (y0 >= tlog >= y1) or (y1 >= tlog >= y0):
            if y1 == y0:
                t = 0.0
            else:
                t = (tlog - y0) / (y1 - y0)
            snr_star = xs[i] + t * (xs[i + 1] - xs[i])
            e0, e1 = rows[i]["average_ensemble_effort"], rows[i + 1]["average_ensemble_effort"]
            l0, l1 = rows[i]["average_ensemble_latency"], rows[i + 1]["average_ensemble_latency"]
            effort_star = 10 ** (math.log10(e0) + t * (math.log10(e1) - math.log10(e0)))
            latency_star = 10 ** (math.log10(l0) + t * (math.log10(l1) - math.log10(l0)))
            fe_min = min(rows[i]["frame_errors"], rows[i + 1]["frame_errors"])
            return dict(snr=snr_star, effort=effort_star, latency=latency_star, frame_errors=fe_min)
    return None

TARGETS = [0.1, 0.001]

# ===================================================================== Sweep B part 1: group-size trade-off
mpg_values = ["1", "2", "4", "8", "16", "24"]
groupsize_curves = []
for mpg in mpg_values:
    rows = rows_for("asced48", "syndrome_sequential", mpg, "6")
    groupsize_curves.append({"mpg": int(mpg), "selector": "syndrome_sequential", "rows": rows})
# mpg = 48 endpoint: reuse full_parallel (one group of the whole ensemble)
fp48 = rows_for("asced48", "full_parallel", "all", "-")
groupsize_curves.append({"mpg": 48, "selector": "full_parallel", "rows": fp48})
groupsize_curves.sort(key=lambda c: c["mpg"])

groupsize_result = {"mpg_values": [c["mpg"] for c in groupsize_curves], "operating_points": {}}
for target in TARGETS:
    pts = []
    for c in groupsize_curves:
        r = interp_at_fer(c["rows"], target)
        if r is None:
            print(f"  [groupsize] mpg={c['mpg']} target_fer={target}: OUT OF RANGE, skipped")
            continue
        pts.append(dict(mpg=c["mpg"], **r))
    groupsize_result["operating_points"][str(target)] = pts
json.dump(groupsize_result, open(OUT / "sweepB_groupsize.json", "w"), indent=1)
print("\nSweep B groupsize operating points:")
print(json.dumps(groupsize_result, indent=1))

# ===================================================================== Sweep B part 2: ordering @ mpg=4
ordering_curves = [
    {"selector": "fixed_sequential", "rows": rows_for("asced48", "fixed_sequential", "4", "6")},
    {"selector": "random_sequential", "rows": rows_for("asced48", "random_sequential", "4", "6")},
    {"selector": "syndrome_sequential", "rows": rows_for("asced48", "syndrome_sequential", "4", "6")},
]
ordering_result = {"operating_points": {}}
for target in TARGETS:
    pts = []
    for c in ordering_curves:
        r = interp_at_fer(c["rows"], target)
        if r is None:
            print(f"  [ordering] selector={c['selector']} target_fer={target}: OUT OF RANGE, skipped")
            continue
        pts.append(dict(selector=c["selector"], **r))
    ordering_result["operating_points"][str(target)] = pts
json.dump(ordering_result, open(OUT / "sweepB_ordering.json", "w"), indent=1)
print("\nSweep B ordering operating points:")
print(json.dumps(ordering_result, indent=1))
