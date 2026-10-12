"""Summarize the target_num_converged recalibration sweep (job 545344) at
members_per_group=8, for asced48/384 x {fixed_sequential, syndrome_sequential},
against the full_parallel (no stopping) baseline.

Reads sweeps/plot_data/raw_numbers.json (built by sweeps/plot_data/build_data.py,
which already scans the shared RESULTS/fig_x_zc11_r4_seq_*_mpg8_n132/ tree) and
writes sweeps/target_converged_recalibration_results.md.
"""

import json
from pathlib import Path

HERE = Path(__file__).parent
rows = json.load(open(HERE / "plot_data" / "raw_numbers.json"))

VARIANTS = ["asced48", "asced384"]
SELECTORS = ["fixed_sequential", "syndrome_sequential"]
TARGETS = ["1", "2", "3", "4", "6", "8", "10", "12", "16", "24"]
SNR_POINTS = ["1", "1.5", "2", "2.5", "3", "3.5", "4"]

by_key = {}
for r in rows:
    if r["members_per_group"] not in ("8", "all"):
        continue
    key = (r["variant"], r["selector"], r["members_per_group"], r["target_num_converged"], r["snr_db"])
    by_key[key] = r


def get(variant, selector, mpg, target, snr):
    return by_key.get((variant, selector, mpg, target, snr))


def fmt_fer(d):
    if d is None or d.get("fer") is None:
        return "-"
    fe, tr = d.get("frame_errors"), d.get("trials")
    tag = "" if (fe or 0) >= 50 else "~"  # flag low-sample points
    return f"{tag}{d['fer']:.4g}"


def fmt_num(x, nd=2):
    return "-" if x is None else f"{x:.{nd}f}"


lines = []
lines.append("# target_num_converged recalibration: results (job 545344)\n")
lines.append(
    "Sweep A reused `target_num_converged={2,6}` from an old, pre-`n_simul`-fix "
    "calibration without re-deriving it on the corrected code (flagged at the "
    "time, see `sweeps/STATUS.md`). This sweep reruns both sequential "
    "selectors at `members_per_group=8` over "
    f"`target_num_converged in {{{','.join(TARGETS)}}}` (280 tasks, all "
    "completed, 0 failures), for asced48 and asced384, against the "
    "`full_parallel` (no early stopping) baseline. `target_errors=200`, "
    "`max_transmissions=2e6` per point (same as Sweep A).\n"
)
lines.append(
    "Low-sample FER points (<50 frame errors observed) are prefixed `~` -- "
    "read those as directional only, not precise.\n"
)

lines.append("## Key findings\n")
lines.append(
    "- **asced48**: FER is essentially insensitive to `target_num_converged` "
    "across the whole swept range {1..24} for both selectors -- even "
    "`target_num_converged=1` (stop after the very first converged path) "
    "matches `full_parallel` within noise at every SNR. Effort/latency savings "
    "are largest at low target values and shrink monotonically back toward the "
    "`full_parallel` baseline as target increases past ~10 (the natural "
    "average-converged-path count at high SNR). **`target_num_converged=1` or "
    "`2` looks close to optimal for asced48** -- old Sweep A's `6` was "
    "unnecessarily conservative (far higher effort/latency than `1`/`2` for no "
    "measurable FER benefit).\n"
)
lines.append(
    "- **asced384 is different and this matters**: `target_num_converged=1` "
    "shows a clear, well-sampled (300+ frame errors, not a low-sample "
    "artifact) FER FLOOR at moderate-to-high SNR -- e.g. fixed_sequential @ "
    "4.0dB: FER=0.0163 at target=1 vs. full_parallel's ~1.1e-05 (~1500x "
    "worse); @ 3.0dB: 0.026 vs. 0.00097 (~27x worse). `target=2` is much "
    "better but still measurably worse than full_parallel at high SNR (@ "
    "4.0dB: 2.3e-04 vs ~1.1e-05, ~21x worse; 332 frame errors, well-sampled). "
    "**`target=3` is roughly where asced384 recovers to matching full_parallel "
    "within noise**, and `target>=4` is clearly safe. **This means Sweep A's "
    "asced384 results at `target_num_converged=2` likely understate FER "
    "performance relative to full_parallel/large-target configurations at "
    "high SNR -- flag this if that comparison is used anywhere (e.g. the "
    "published plot).**\n"
)
lines.append(
    "- Intuition: with only 48 total paths, even 1-2 converged candidates "
    "are enough to reliably include the true codeword; with 384 paths, "
    "stopping after just 1-2 converged candidates risks locking in a "
    "spurious-but-valid low-weight codeword before the genuinely correct one "
    "(among the larger pool) has had a chance to converge too -- a real, "
    "ensemble-size-dependent effect, not a code-family artifact.\n"
)
lines.append(
    "- Recommendation: recalibrate to **asced48: target_num_converged=2** "
    "(small safety margin over the apparently-fine target=1, negligible extra "
    "cost), **asced384: target_num_converged=4** (comfortably past the "
    "target=3 recovery point, still far cheaper than full_parallel or the old "
    "target=6). Both selectors (fixed_sequential, syndrome_sequential) behave "
    "similarly, no strong reason to prefer one from this sweep alone.\n"
)

for variant in VARIANTS:
    lines.append(f"\n## {variant}\n")
    fp_row = {snr: get(variant, "full_parallel", "all", "-", snr) for snr in SNR_POINTS}
    lines.append(
        "**full_parallel baseline** (no stopping) -- FER / effort / latency:\n"
    )
    lines.append("| SNR | " + " | ".join(SNR_POINTS) + " |")
    lines.append("|---|" + "---|" * len(SNR_POINTS))
    lines.append("| FER | " + " | ".join(fmt_fer(fp_row[s]) for s in SNR_POINTS) + " |")
    lines.append(
        "| effort | "
        + " | ".join(fmt_num(fp_row[s]["average_ensemble_effort"] if fp_row[s] else None) for s in SNR_POINTS)
        + " |"
    )
    lines.append(
        "| latency | "
        + " | ".join(fmt_num(fp_row[s]["average_ensemble_latency"] if fp_row[s] else None) for s in SNR_POINTS)
        + " |"
    )
    lines.append(
        "| avg #converged (of "
        + ("48" if variant == "asced48" else "384")
        + ") | "
        + " | ".join(fmt_num(fp_row[s]["average_number_converged_path"] if fp_row[s] else None) for s in SNR_POINTS)
        + " |"
    )

    for selector in SELECTORS:
        lines.append(f"\n### {variant}, `{selector}`, mpg=8\n")
        lines.append("FER by target_num_converged x SNR:\n")
        lines.append("| target | " + " | ".join(SNR_POINTS) + " |")
        lines.append("|---|" + "---|" * len(SNR_POINTS))
        for target in TARGETS:
            cells = [fmt_fer(get(variant, selector, "8", target, s)) for s in SNR_POINTS]
            lines.append(f"| {target} | " + " | ".join(cells) + " |")
        lines.append(f"| full_parallel (ref) | " + " | ".join(fmt_fer(fp_row[s]) for s in SNR_POINTS) + " |")

        lines.append("\neffort (avg ensemble BP-iteration cost) by target_num_converged x SNR:\n")
        lines.append("| target | " + " | ".join(SNR_POINTS) + " |")
        lines.append("|---|" + "---|" * len(SNR_POINTS))
        for target in TARGETS:
            cells = [
                fmt_num(d["average_ensemble_effort"] if (d := get(variant, selector, "8", target, s)) else None)
                for s in SNR_POINTS
            ]
            lines.append(f"| {target} | " + " | ".join(cells) + " |")
        lines.append(
            f"| full_parallel (ref) | "
            + " | ".join(fmt_num(fp_row[s]["average_ensemble_effort"] if fp_row[s] else None) for s in SNR_POINTS)
            + " |"
        )

        lines.append("\nlatency (avg per-group max iteration cost) by target_num_converged x SNR:\n")
        lines.append("| target | " + " | ".join(SNR_POINTS) + " |")
        lines.append("|---|" + "---|" * len(SNR_POINTS))
        for target in TARGETS:
            cells = [
                fmt_num(d["average_ensemble_latency"] if (d := get(variant, selector, "8", target, s)) else None)
                for s in SNR_POINTS
            ]
            lines.append(f"| {target} | " + " | ".join(cells) + " |")
        lines.append(
            f"| full_parallel (ref) | "
            + " | ".join(fmt_num(fp_row[s]["average_ensemble_latency"] if fp_row[s] else None) for s in SNR_POINTS)
            + " |"
        )

out_path = HERE / "target_converged_recalibration_results.md"
out_path.write_text("\n".join(lines) + "\n")
print(f"Wrote {out_path}")
