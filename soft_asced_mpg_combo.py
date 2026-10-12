# pylint: disable=invalid-name
"""Soft-aSCED combined with members-per-group (mpg) / MConvergedPolicy
sequential early-stopping, following up on the explicit open question left
at the end of both ``sweeps/soft_asced_single_row_results.md`` ("not yet
combined with members_per_group/MConvergedPolicy") and
``sweeps/soft_asced_delta_gt1_results.md`` (asced384's full-parallel
"latency" barely improves at FER-matching bias, because the ensemble's
single worst path dominates the max-over-paths metric).

Imports construction helpers directly from ``soft_asced_delta_gt1.py``
(importable: all its CLI logic is gated behind ``if __name__ ==
"__main__":``) rather than duplicating them. Does NOT modify
``soft_asced_delta_gt1.py``, ``reproduce_fig3_RL_zc11_asced_5G_LDPC.py``, or
``reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py``.

Usage:
    uv run python soft_asced_mpg_combo.py natural asced48
    uv run python soft_asced_mpg_combo.py natural asced384
    uv run python soft_asced_mpg_combo.py sweep  asced48
    uv run python soft_asced_mpg_combo.py sweep  asced384
    uv run python soft_asced_mpg_combo.py hard_extra asced384   # targets 48/96 only
"""

import json
import os
import sys

import soft_asced_delta_gt1 as base
from soft_asced_delta_gt1 import (
    H,
    candidate_rows,
    block_offset,
    Zc,
    VARIANTS,
    SNR_POINTS,
    create_soft_asced_config,
    create_hard_asced_config,
    run_sim,
    wilson_interval,
    _fer_ci_str,
    _lookup_by_snr,
)

BIAS = 10  # fixed throughout -- the already-validated FER-matching bias for
# both variants per sweeps/soft_asced_delta_gt1_results.md.
MPG = 8  # matches the existing hard recalibration's mpg, for direct
# comparability (sweeps/STATUS.md's "target_num_converged recalibration
# sweep, mpg=8" section).
SELECTOR = "fixed_sequential"  # see "Design decision: selector" below --
# syndrome_sequential is confirmed broken for soft paths (segfault), not
# used here. fixed_sequential is the primary, required selector per scope.

TARGET_ERRORS = {"asced48": 200, "asced384": 150}
MAX_TRANSMISSIONS = int(2e6)

# Soft target_num_converged sweep, chosen AFTER measuring each variant's
# natural average_number_converged_path at bias=10 full-parallel (see
# `measure_natural` below / the results doc's "Step 1" section): soft
# converges far more readily than hard (asced48 natural: 13.3 -> 44.3 over
# 1-3dB, out of 48; asced384 natural: 45.5 -> 140.8 over 1-2dB, out of 384 --
# vs. hard's 4.7-10.0 (asced48) / 5.8-8.7 (asced384) from the existing
# recalibration, sweeps/target_converged_recalibration_results.md). Both
# sweeps below reuse hard's existing recalibration grid {1,2,3,4,6,8,10,12,
# 16,24} wherever possible, for a FREE, already-run matched-target hard
# comparison (no rerun needed) -- see the results doc's hardcoded
# HARD_RECAL table. asced384 additionally tests two higher values (48, 96)
# specifically because its measured natural average (45-141) is well above
# hard's whole existing grid, so the FER-recovery point (if any, within a
# sane budget) could lie beyond target=24; hard is re-run fresh ONLY at
# these two new values (mode "hard_extra"), since they are not in the
# existing recalibration data.
TARGET_SWEEP = {
    "asced48": [1, 2, 4, 8, 16, 24],
    "asced384": [4, 8, 16, 24, 48, 96],
}
HARD_NEW_TARGETS = {"asced48": [], "asced384": [48, 96]}

RESULTS_ROOT = "RESULTS/soft_asced_mpg_combo"


def _decoder_stats(save_dir):
    path = os.path.join(save_dir, "decoder_stats.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def measure_natural(variant_name, out_path=None):
    """Full-parallel (no mpg/MConverged) soft ensemble at bias=10, measuring
    `average_number_converged_path` from decoder_stats.json -- the natural
    convergence rate that informs the target_num_converged sweep below."""
    params = VARIANTS[variant_name]
    snr_points = SNR_POINTS[variant_name]
    save_dir = os.path.join(RESULTS_ROOT, variant_name, "natural_fullparallel_bias10")

    soft_cfg = create_soft_asced_config(
        H, candidate_rows, block_offset, Zc, bias=BIAS, **params
    )
    fer, effort, latency, counts = run_sim(
        soft_cfg, snr_points, save_dir,
        target_errors=TARGET_ERRORS[variant_name], max_transmissions=MAX_TRANSMISSIONS,
    )
    stats = _decoder_stats(save_dir)
    avg_converged = stats.get("average_number_converged_path", {})

    result = dict(variant=variant_name, bias=BIAS, fer=fer, effort=effort,
                  latency=latency, counts=counts, avg_converged=avg_converged)
    print(f"\n=== {variant_name}: natural avg_number_converged_path (soft, full_parallel, bias={BIAS}) ===")
    for snr in snr_points:
        ac = _lookup_by_snr(avg_converged, snr)
        total_paths = sum(2 ** rows.shape[0] for _, rows in base.get_row_segments(
            candidate_rows, block_offset, Zc, params["split_pattern"], params["num_used_blocks"]))
        print(f"  SNR={snr}dB: avg_converged={ac:.2f} / {total_paths} paths "
              f"({100*ac/total_paths:.1f}%)   FER={_fer_ci_str(fer, counts, snr)}")

    if out_path is None:
        out_path = os.path.join(RESULTS_ROOT, variant_name, "natural_summary.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2, default=str)
    return result


def run_soft_mpg(variant_name, target, mpg=MPG, selector=SELECTOR):
    params = VARIANTS[variant_name]
    snr_points = SNR_POINTS[variant_name]
    save_dir = os.path.join(RESULTS_ROOT, variant_name, f"soft_bias{BIAS}_{selector}_mpg{mpg}_target{target}")

    soft_cfg = create_soft_asced_config(
        H, candidate_rows, block_offset, Zc, bias=BIAS, **params
    )
    soft_cfg.set_mConvergedConfig(target)
    if selector == "fixed_sequential":
        soft_cfg.set_fixed_sequential(mpg)
    else:
        raise ValueError(f"selector '{selector}' not supported here (see Design decisions)")

    fer, effort, latency, counts = run_sim(
        soft_cfg, snr_points, save_dir,
        target_errors=TARGET_ERRORS[variant_name], max_transmissions=MAX_TRANSMISSIONS,
    )
    stats = _decoder_stats(save_dir)
    avg_converged = stats.get("average_number_converged_path", {})
    return dict(variant=variant_name, target=target, mpg=mpg, selector=selector,
                fer=fer, effort=effort, latency=latency, counts=counts,
                avg_converged=avg_converged)


def run_hard_mpg(variant_name, target, mpg=MPG, selector=SELECTOR):
    params = VARIANTS[variant_name]
    snr_points = SNR_POINTS[variant_name]
    save_dir = os.path.join(RESULTS_ROOT, variant_name, f"hard_{selector}_mpg{mpg}_target{target}")

    hard_cfg = create_hard_asced_config(
        H, candidate_rows, block_offset, Zc, **params
    )
    hard_cfg.set_mConvergedConfig(target)
    if selector == "fixed_sequential":
        hard_cfg.set_fixed_sequential(mpg)
    else:
        raise ValueError(f"selector '{selector}' not supported here (see Design decisions)")

    fer, effort, latency, counts = run_sim(
        hard_cfg, snr_points, save_dir,
        target_errors=TARGET_ERRORS[variant_name], max_transmissions=MAX_TRANSMISSIONS,
    )
    return dict(variant=variant_name, target=target, mpg=mpg, selector=selector,
                fer=fer, effort=effort, latency=latency, counts=counts)


def sweep_soft(variant_name):
    results = {}
    for target in TARGET_SWEEP[variant_name]:
        print(f"\n--- soft {variant_name}, target_num_converged={target}, mpg={MPG}, {SELECTOR} ---")
        r = run_soft_mpg(variant_name, target)
        results[target] = r
        for snr in SNR_POINTS[variant_name]:
            print(f"  FER @ {snr}dB: {_fer_ci_str(r['fer'], r['counts'], snr)}")
        print(f"  effort: {r['effort']}")
        print(f"  latency: {r['latency']}")
        print(f"  avg_converged: {r['avg_converged']}")

    out_path = os.path.join(RESULTS_ROOT, variant_name, "soft_sweep_summary.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    return results


def sweep_hard_new_targets(variant_name):
    """Fresh hard+mpg runs ONLY at the target values added for soft beyond
    the existing recalibration grid (see HARD_NEW_TARGETS)."""
    results = {}
    for target in HARD_NEW_TARGETS[variant_name]:
        print(f"\n--- hard {variant_name}, target_num_converged={target} (NEW, not in existing recalibration grid), mpg={MPG}, {SELECTOR} ---")
        r = run_hard_mpg(variant_name, target)
        results[target] = r
        for snr in SNR_POINTS[variant_name]:
            print(f"  FER @ {snr}dB: {_fer_ci_str(r['fer'], r['counts'], snr)}")
        print(f"  effort: {r['effort']}")
        print(f"  latency: {r['latency']}")

    out_path = os.path.join(RESULTS_ROOT, variant_name, "hard_new_targets_summary.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    return results


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    mode, variant = sys.argv[1], sys.argv[2]

    if mode == "natural":
        measure_natural(variant)
    elif mode == "sweep":
        sweep_soft(variant)
    elif mode == "hard_extra":
        sweep_hard_new_targets(variant)
    else:
        raise ValueError(f"Unknown mode '{mode}'")


if __name__ == "__main__":
    main()
