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
    uv run python soft_asced_mpg_combo.py hard_extra asced384   # targets 48/96 only, fixed_sequential

    # syndrome_sequential (added once channel-code-lib2's pybind11 dangling-
    # reference segfault was fixed upstream, commit 64884c1 -- see
    # sweeps/soft_asced_syndrome_sequential_results.md):
    uv run python soft_asced_mpg_combo.py sweep_syn asced48   # soft, full TARGET_SWEEP grid
    uv run python soft_asced_mpg_combo.py sweep_syn asced384
    uv run python soft_asced_mpg_combo.py hard_syn   asced48   # hard, full TARGET_SWEEP grid (fresh)
    uv run python soft_asced_mpg_combo.py hard_syn   asced384
"""

import json
import os
import sys

import numpy as np

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
SELECTOR = "fixed_sequential"  # default/primary selector, still supported
# as before. "syndrome_sequential" is now ALSO supported (the segfault
# documented below is fixed upstream -- channel-code-lib2 commit 64884c1,
# see sweeps/STATUS.md's "syndrome_sequential segfault: fixed" section --
# this worktree has been rebased onto that fix).

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


def build_soft_syndrome_mappings(variant_name):
    """Build the per-path mapping-matrix list for
    ``Ensemble_config.set_syndrome_sequential(mpg, mappings)`` on a SOFT
    ensemble.

    Design decision (flagged): one all-zero ``(Delta_seg, n)`` matrix per
    path, where ``Delta_seg`` is THAT path's own segment's row count (soft
    paths within a single variant are not all the same Delta_seg -- e.g.
    asced48 mixes Delta_seg=2 and Delta_seg=3 segments) and ``n`` is the
    physical codeword length (``H.shape[1]``, the ensemble's reference PCM
    column count -- NOT the per-path soft PCM's wider n+Delta_seg column
    count, which is what SyndromeSequentialSelector's shape validation
    actually expects of the mapping's column dimension).

    Chosen because it makes soft's syndrome-weight ranking heuristic for
    ordering groups exactly equivalent to hard's: per
    SyndromeSequentialSelector.cpp's ``compute_syndrome_weight``, the
    mapping produces the auxiliary-node bits appended to the physical
    codeword's hard-decision before computing ``H_path @ [hd | aux]``; an
    all-zero mapping forces ``aux=0`` unconditionally, i.e. the guessed
    syndrome-VN pattern is always treated as zero when ranking -- exactly
    mirroring hard's own ``compute_syndrome_weight``, which (per
    sweeps/soft_asced_mpg_combo_results.md's Step 0) already ignores
    ``affine_offset`` entirely when computing its own syndrome weight.
    Verified non-degenerate (real segfault-free, real per-codeword
    reordering, not a silent no-op) via a direct smoke test before trusting
    any sweep built on this -- see this task's final report.

    Mirrors create_soft_asced_config's own segment/pattern loop EXACTLY
    (same get_row_segments call, same per-segment `2**Delta_seg` multiplier)
    so the returned mapping list lines up index-for-index with the
    path_configs list that Ensemble_config was built from.
    """
    params = VARIANTS[variant_name]
    n = H.shape[1]
    mappings = []
    for _, rows in base.get_row_segments(
        candidate_rows, block_offset, Zc, params["split_pattern"], params["num_used_blocks"]
    ):
        Delta_seg = rows.shape[0]
        zero_mapping = np.zeros((Delta_seg, n), dtype=int)
        for _ in range(2 ** Delta_seg):
            mappings.append(zero_mapping)
    return mappings


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
    elif selector == "syndrome_sequential":
        mappings = build_soft_syndrome_mappings(variant_name)
        soft_cfg.set_syndrome_sequential(mpg, mappings)
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
    elif selector == "syndrome_sequential":
        # No mapping needed: hard aSCED's per-path PCM (H_aux = vstack(H,
        # rows)) has the SAME column count as the reference H (the guessed
        # suffix is baked into affine_offset, not extra columns) -- it is
        # not a "generalized PCM" in SyndromeSequentialSelector's sense, so
        # the existing no-mapping convenience overload (used throughout the
        # rest of this project) applies unchanged.
        hard_cfg.set_syndrome_sequential(mpg)
    else:
        raise ValueError(f"selector '{selector}' not supported here (see Design decisions)")

    fer, effort, latency, counts = run_sim(
        hard_cfg, snr_points, save_dir,
        target_errors=TARGET_ERRORS[variant_name], max_transmissions=MAX_TRANSMISSIONS,
    )
    return dict(variant=variant_name, target=target, mpg=mpg, selector=selector,
                fer=fer, effort=effort, latency=latency, counts=counts)


def sweep_soft(variant_name, selector=SELECTOR):
    results = {}
    for target in TARGET_SWEEP[variant_name]:
        print(f"\n--- soft {variant_name}, target_num_converged={target}, mpg={MPG}, {selector} ---")
        r = run_soft_mpg(variant_name, target, selector=selector)
        results[target] = r
        for snr in SNR_POINTS[variant_name]:
            print(f"  FER @ {snr}dB: {_fer_ci_str(r['fer'], r['counts'], snr)}")
        print(f"  effort: {r['effort']}")
        print(f"  latency: {r['latency']}")
        print(f"  avg_converged: {r['avg_converged']}")

    suffix = "" if selector == "fixed_sequential" else f"_{selector}"
    out_path = os.path.join(RESULTS_ROOT, variant_name, f"soft_sweep_summary{suffix}.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    return results


def sweep_hard_all_targets(variant_name, selector):
    """Fresh hard+mpg runs at EVERY target in TARGET_SWEEP (not just the
    targets beyond the existing fixed_sequential recalibration grid, unlike
    `sweep_hard_new_targets` below).

    Required for `syndrome_sequential` specifically, for two independent
    reasons (either alone would justify it): (1) the task's explicit
    instruction not to assume `fixed_sequential`'s hard numbers transfer,
    since `syndrome_sequential` reorders paths dynamically per-codeword
    rather than using `fixed_sequential`'s static index order; (2) this
    worktree's (gitignored) RESULTS/ tree does not contain job 545344's raw
    per-SNR frame-error/trial counts for ANY selector -- only the published
    point estimates in sweeps/target_converged_recalibration_results.md --
    so real Wilson CIs for hard+syndrome_sequential need a fresh run
    regardless of selector-equivalence questions.
    """
    results = {}
    for target in TARGET_SWEEP[variant_name]:
        print(f"\n--- hard {variant_name}, target_num_converged={target}, mpg={MPG}, {selector} ---")
        r = run_hard_mpg(variant_name, target, selector=selector)
        results[target] = r
        for snr in SNR_POINTS[variant_name]:
            print(f"  FER @ {snr}dB: {_fer_ci_str(r['fer'], r['counts'], snr)}")
        print(f"  effort: {r['effort']}")
        print(f"  latency: {r['latency']}")

    out_path = os.path.join(RESULTS_ROOT, variant_name, f"hard_sweep_summary_{selector}.json")
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
    elif mode == "sweep_syn":
        # Soft, syndrome_sequential, full TARGET_SWEEP grid.
        sweep_soft(variant, selector="syndrome_sequential")
    elif mode == "hard_syn":
        # Hard, syndrome_sequential, full TARGET_SWEEP grid (fresh -- see
        # sweep_hard_all_targets's docstring for why this can't reuse the
        # fixed_sequential-era matched-target shortcut).
        sweep_hard_all_targets(variant, selector="syndrome_sequential")
    else:
        raise ValueError(f"Unknown mode '{mode}'")


if __name__ == "__main__":
    main()
