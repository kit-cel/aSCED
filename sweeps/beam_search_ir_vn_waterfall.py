# pylint: disable=invalid-name
"""First (preliminary/exploratory) waterfall comparison for the beam-search
LDPC decoder's IR-VN splitting extension.

Compares, on the aSCED-48/NMSA reference's IDENTICAL base code (BG2,
Zc=11, n_simul=132, k=66 -- see Codes/build_ir_vn_extended_pcm.py):

  1. Beam_Search_Decoder on the IR-VN-EXTENDED PCM (H_ext = base H + 2
     candidate HARQ blocks appended as new check rows, each with its own
     new degree-1 IR-VN column -- see Codes/build_ir_vn_extended_pcm.py),
     with ALL columns (state + IR) eligible for splitting, across 4
     configs: beam_width in {16, 32} x max_iterations in {32 (capped,
     same total budget as the NMSA/asced-48 references), 1000 (generous --
     empirically confirmed below NOT to hang/take unreasonable wall-clock
     time per trial; see calibration note at the bottom of this docstring)}.
  2. A FRESH NMSA-32 baseline, run at the SAME (reduced) SNR grid/frame
     budget as the beam-search configs above (NOT the existing
     RESULTS/.../nmsa/FER.json, which was run at a different
     target_errors/SNR range) -- on the BASE (non-extended) H, matching
     reproduce_fig3_RL_zc11_asced_5G_LDPC.py's `flag_nmsa` block
     (cn_update_type=msa, norm_factor=0.75, max_iterations=32,
     early_stopping=True, scheduling_type=flooding).
  3. The existing aSCED-48 TCOM-paper reference curve, REUSED (not rerun --
     it's an expensive, already-validated reference at target_errors=1000):
     /home/pj9034/aSCED/RESULTS/fig_x_zc11_r4n=132/aSCED_48_split1_batch2/FER.json

All three use the SAME real-codeword (PCM/G encoder) simulation mode as
reproduce_fig3_RL_zc11_asced_5G_LDPC.py (`use_all_zero = False` there --
confirmed by reading that script), i.e. `sim.init(enc_cfg, cfg, False)`,
NOT the all-zero-codeword mode.

PRELIMINARY STATISTICAL BUDGET, NOT matched to the asced-48 reference's
much larger one (target_errors=1000, essentially uncapped transmissions):
this script uses `target_errors=200` with `max_transmissions` CAPPED at a
few thousand per SNR point (see TARGET_ERRORS/MAX_TRANSMISSIONS below) so
the whole sweep completes in minutes rather than hours -- high-SNR points
may show noisier (or zero) error counts as a result; this is expected and
reported as such (hollow/flagged points, per PLOTTING_STYLE.md).

SNR grid: 2.0-4.0 dB in 0.5 dB steps (5 points), overlapping the asced-48
reference's own 1.0-5.0 dB / 0.5 dB grid.

For every beam-search config, also reports the split-location histogram
(Beam_Search_Decoder::vn_split_counts, read from decoder_stats.json via
the new Beam_Search_Decoder_Stats plumbing) bucketed into:
  (1) the first 2*Zc punctured state VNs (standard 5G systematic puncturing)
  (2) the remaining, normally-transmitted state VNs
  (3) the appended IR-VN columns

Calibration note (max_iterations=1000 "generous" budget, done interactively
before committing to this script -- see session report): beam_width in
{16, 32}, 500 trials at SNR in {2.0..4.0}dB took ~4-7 seconds total (8
threads) at every point -- no blowup, no hang. 1000 was chosen as
iters_per_round(5) * (n_ext(176) + settle margin), rounded to a clean
number, comfortably above the max possible split depth (176 eligible VNs).

IMPORTANT (discovered while building this script -- see session report):
``Simulation_Env::simulate_transmission`` declares its per-thread scratch
buffer as ``thread_local TrialContext local(n);`` -- C++ ``thread_local``
function-statics are constructed ONCE per OS thread for the lifetime of
the process, not re-constructed on a later call with a different `n`. Since
this codebase's OpenMP runtime reuses its thread pool (same OS threads)
across separate ``#pragma omp parallel`` regions, running a SECOND
``Simulation_Env`` with a DIFFERENT `n` (e.g. the base H, n=154) in the
SAME process after an earlier one (e.g. H_ext, n=176) silently reuses the
FIRST call's stale, wrong-sized thread-local buffer and crashes every
worker thread with an uncaught ``BP_Decoder::decode: transmitted_bits > n``
(or ``!= n``) ``std::logic_error`` (OpenMP can't propagate it across the
parallel region, so the whole process aborts via ``std::terminate``) --
reproduced directly in this session, confirmed to disappear when the same
config is run alone in a fresh process. This is a genuine latent bug in
the C++ library (out of this task's explicit scope: it only authorizes
adding the per-VN split counter there, not broader Simulation_Env
refactors), triggered only by driving differently-sized codes from one
long-lived process -- which is exactly what a naive "run all 5 configs in
one script" loop does here (4 beam-search configs on H_ext, n=176, then
NMSA on the base H, n=154).
WORKAROUND (this script): each config (4 beam-search + NMSA) runs in its
OWN subprocess (`--run <name>`); the no-argument invocation is an
orchestrator that launches them one at a time and merges the resulting
per-config JSON files. This sidesteps the bug with no changes to the C++
library, at the cost of re-paying each subprocess's Python/library import
and PCM-construction overhead (a few seconds, negligible next to the
simulation time itself).

Usage:
    uv run python sweeps/beam_search_ir_vn_waterfall.py             # orchestrator: runs everything
    uv run python sweeps/beam_search_ir_vn_waterfall.py --run NAME  # one config, in-process (internal use)
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
os.chdir(REPO_ROOT)  # Codes/build_ir_vn_extended_pcm.py's default npz path is CWD-relative,
                     # matching every other reproduce_*/sweeps script's convention of being
                     # run from the repo root.

import channel_code_lib2 as ccl
from Codes.build_ir_vn_extended_pcm import build_ir_vn_code
RESULTS_DIR = REPO_ROOT / "RESULTS" / "beam_search_ir_vn_waterfall"
PLOT_DATA_DIR = REPO_ROOT / "sweeps" / "plot_data"
ASCED48_FER_PATH = Path(
    "/home/pj9034/aSCED/RESULTS/fig_x_zc11_r4n=132/aSCED_48_split1_batch2/FER.json"
)  # NOTE: absolute path to the main aSCED checkout, NOT this worktree -- this
# expensive (target_errors=1000) reference run lives there, outside any
# worktree-local RESULTS/ dir (read-only; never rerun, never modified).

SNR_GRID = [2.0, 2.5, 3.0, 3.5, 4.0]
TARGET_ERRORS = 200          # modest/preliminary -- NOT the asced-48 reference's target_errors=1000
MAX_TRANSMISSIONS = 3000     # capped per SNR point so the whole sweep finishes in minutes, not hours
NUM_THREADS = 8
NUM_AD_BLOCKS = 2            # mirrors asced-48's own num_used_blocks=2 candidate-row footprint 1:1
BLOCK_OFFSET = 0


def run_one_point(sim, snr):
    """Runs get_error_rates for a SINGLE SNR point (so we can time it and
    read its own decoder_stats.json snapshot cleanly), returns
    (fer, errors, trials, wall_seconds)."""
    t0 = time.time()
    sim.get_error_rates([snr])
    wall = time.time() - t0
    fer = sim.error_rates["FER-SNR"][snr]
    # Reconstruct trial/error counts from stats.json (written by auto_save).
    with open(os.path.join(sim.save_dir, "stats.json")) as f:
        stats = json.load(f)
    key = str(snr) if str(snr) in stats["trials"] else next(
        k for k in stats["trials"] if abs(float(k) - snr) < 1e-9
    )
    trials = stats["trials"][key]
    errors = stats["frame_errors"][key]
    return fer, errors, trials, wall


def run_beam_config(name, H_ext, k, n_ext, combined_punct, beam_width, max_iterations,
                     number_vns, Zc):
    """Runs one beam-search config across SNR_GRID, returns (pts, split_hist, wall_total)."""
    save_dir = str(RESULTS_DIR / name)
    os.makedirs(save_dir, exist_ok=True)

    enc_cfg = ccl.PCM_Encoder_config(H_ext, k, n_ext)

    cfg = ccl.Beam_Search_config(H_ext)
    cfg.early_stopping = True
    cfg.max_iterations = max_iterations
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    cfg.iters_per_round = 5
    cfg.beam_width = beam_width
    cfg.target_converged = 1
    # eligible_vns left empty (default) -- ALL columns (state + IR) eligible,
    # per the owner's explicit request to allow splitting everywhere and
    # just track where it happens.
    cfg.affine_offset = np.zeros(H_ext.shape[0], dtype=int)

    sim = ccl.Simulation_Env(k, n_ext, "all")
    sim.target_errors = TARGET_ERRORS
    sim.max_transmissions = MAX_TRANSMISSIONS
    sim.num_threads = NUM_THREADS
    sim.auto_save = True
    sim.save_dir = save_dir
    sim.puncturing(combined_punct)
    sim.init(enc_cfg, cfg, False)  # real-codeword mode, matching the asced-48/NMSA references

    pts = []
    wall_total = 0.0
    split_hist_per_snr = {}
    for snr in SNR_GRID:
        fer, errors, trials, wall = run_one_point(sim, snr)
        wall_total += wall
        pts.append([snr, fer, errors, wall * NUM_THREADS])
        print(f"  [{name}] SNR={snr}: FER={fer:.4g} ({errors}/{trials} trials), wall={wall:.1f}s")

        with open(os.path.join(save_dir, "decoder_stats.json")) as f:
            dstats = json.load(f)
        key = str(snr) if str(snr) in dstats.get("vn_split_counts", {}) else next(
            k for k in dstats["vn_split_counts"] if abs(float(k) - snr) < 1e-9
        )
        counts = np.array(dstats["vn_split_counts"][key], dtype=np.int64)
        split_hist_per_snr[snr] = {
            "punctured_state_vns": int(counts[:2 * Zc].sum()),
            "other_state_vns": int(counts[2 * Zc:number_vns].sum()),
            "ir_vns": int(counts[number_vns:].sum()),
            "total": int(counts.sum()),
        }

    agg = {
        "punctured_state_vns": sum(v["punctured_state_vns"] for v in split_hist_per_snr.values()),
        "other_state_vns": sum(v["other_state_vns"] for v in split_hist_per_snr.values()),
        "ir_vns": sum(v["ir_vns"] for v in split_hist_per_snr.values()),
    }
    agg["total"] = agg["punctured_state_vns"] + agg["other_state_vns"] + agg["ir_vns"]

    return pts, {"per_snr": split_hist_per_snr, "aggregate": agg}, wall_total


def run_nmsa_baseline(H, k, n):
    """Fresh NMSA-32 baseline on the BASE (non-extended) H, matching
    reproduce_fig3_RL_zc11_asced_5G_LDPC.py's flag_nmsa block exactly, at
    the SAME SNR grid / frame budget as the beam-search runs."""
    Zc = 11
    save_dir = str(RESULTS_DIR / "nmsa_fresh")
    os.makedirs(save_dir, exist_ok=True)

    enc_cfg = ccl.PCM_Encoder_config(H, k, n)

    cfg = ccl.BP_config(H)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    cfg.affine_offset = np.zeros(H.shape[0], dtype=int)

    sim = ccl.Simulation_Env(k, n, "all")
    sim.target_errors = TARGET_ERRORS
    sim.max_transmissions = MAX_TRANSMISSIONS
    sim.num_threads = NUM_THREADS
    sim.auto_save = True
    sim.save_dir = save_dir
    sim.puncturing(np.arange(2 * Zc, dtype=int))
    sim.init(enc_cfg, cfg, False)

    pts = []
    wall_total = 0.0
    for snr in SNR_GRID:
        fer, errors, trials, wall = run_one_point(sim, snr)
        wall_total += wall
        pts.append([snr, fer, errors, wall * NUM_THREADS])
        print(f"  [nmsa_fresh] SNR={snr}: FER={fer:.4g} ({errors}/{trials} trials), wall={wall:.1f}s")

    return pts, wall_total


def load_asced48_reference():
    with open(ASCED48_FER_PATH) as f:
        d = json.load(f)
    return sorted([[float(k), v] for k, v in d.items()], key=lambda p: p[0])


BEAM_CONFIGS = [
    ("beam16_capped32", 16, 32),
    ("beam32_capped32", 32, 32),
    ("beam16_generous", 16, 1000),
    ("beam32_generous", 32, 1000),
]
ALL_RUN_NAMES = [name for name, _, _ in BEAM_CONFIGS] + ["nmsa_fresh"]


def run_single_config(name):
    """Runs exactly ONE config (one beam-search config or nmsa_fresh) and
    writes its own small JSON (pts + split_hist, if applicable) to
    RESULTS_DIR/<name>/run_result.json. Always invoked as a FRESH
    subprocess by the orchestrator (see module docstring: thread_local
    TrialContext sizing bug) -- never call this twice in one process with
    a different code size."""
    print(f"Building IR-VN-extended PCM (num_ad_blocks={NUM_AD_BLOCKS})...")
    result = build_ir_vn_code(num_ad_blocks=NUM_AD_BLOCKS, block_offset=BLOCK_OFFSET)
    H, H_ext = result["H"], result["H_ext"]
    k, n_ext, Delta = result["k"], result["n_ext"], result["Delta"]
    number_vns = result["number_vns"]
    Zc = result["Zc"]
    combined_punct = result["combined_puncturing"]

    beam_lookup = {n: (bw, mi) for n, bw, mi in BEAM_CONFIGS}
    os.makedirs(RESULTS_DIR / name, exist_ok=True)
    out_path = RESULTS_DIR / name / "run_result.json"

    if name == "nmsa_fresh":
        pts, wall_total = run_nmsa_baseline(H, k, number_vns)
        payload = {"name": name, "pts": pts, "wall_total": wall_total}
    elif name in beam_lookup:
        beam_width, max_iterations = beam_lookup[name]
        print(f"Running {name} (beam_width={beam_width}, max_iterations={max_iterations})...")
        pts, split_hist, wall_total = run_beam_config(
            name, H_ext, k, n_ext, combined_punct, beam_width, max_iterations,
            number_vns, Zc,
        )
        payload = {
            "name": name, "beam_width": beam_width, "max_iterations": max_iterations,
            "pts": pts, "split_hist": split_hist, "wall_total": wall_total,
        }
    else:
        raise ValueError(f"Unknown config name {name!r} (expected one of {ALL_RUN_NAMES})")

    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"  {name} wall_total={payload['wall_total']:.1f}s -> {out_path}")

    # Also stash code/combined-puncturing metadata once (idempotent, tiny) so
    # the orchestrator can read it back without rebuilding the PCM itself.
    meta_path = RESULTS_DIR / "code_meta.json"
    if not meta_path.exists():
        with open(meta_path, "w") as f:
            json.dump({
                "Zc": Zc, "number_vns": int(number_vns), "n_ext": int(n_ext),
                "Delta": int(Delta), "k": int(k), "num_ad_blocks": NUM_AD_BLOCKS,
                "combined_puncturing": combined_punct.tolist(),
            }, f, indent=2)


def main_orchestrate():
    """No-arg entry point: launches each config as its OWN subprocess (see
    module docstring for why), then merges the per-config JSON files into
    the final sweeps/plot_data/beam_search_ir_vn_data.json."""
    t_start = time.time()
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(PLOT_DATA_DIR, exist_ok=True)

    for name in ALL_RUN_NAMES:
        print(f"=== subprocess: --run {name} ===", flush=True)
        t0 = time.time()
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--run", name],
            cwd=str(REPO_ROOT), check=True,
        )
        print(f"=== {name} subprocess finished in {time.time() - t0:.1f}s ===", flush=True)

    with open(RESULTS_DIR / "code_meta.json") as f:
        code_meta = json.load(f)

    curves = []
    split_histograms = {}
    for name, beam_width, max_iterations in BEAM_CONFIGS:
        with open(RESULTS_DIR / name / "run_result.json") as f:
            r = json.load(f)
        curves.append({
            "name": f"beam_search bw={beam_width} max_it={max_iterations}",
            "group": "Beam search",
            "visible": True,
            "pts": r["pts"],
        })
        split_histograms[name] = r["split_hist"]

    with open(RESULTS_DIR / "nmsa_fresh" / "run_result.json") as f:
        r = json.load(f)
    curves.append({
        "name": "NMSA-32 (fresh)",
        "group": "Baseline",
        "visible": True,
        "pts": r["pts"],
    })

    asced48_pts = load_asced48_reference()

    data = {
        "beam_search_ir_vn_zc11_r4n132": {
            "xlabel": "Eb/N0 [dB]",
            "xlog": False,
            "ylabel": "FER",
            "curves": curves,
            "reference": {
                "name": "aSCED-48 (TCOM ref., target_errors=1000, NOT rerun)",
                "pts": asced48_pts,
            },
            "paper": {},
        },
        "split_histograms": split_histograms,
        "meta": {
            "note": "Preliminary/exploratory sweep: target_errors="
                    f"{TARGET_ERRORS}, max_transmissions={MAX_TRANSMISSIONS} per SNR point "
                    "(NOT matched to the asced-48 reference's target_errors=1000 budget).",
            **code_meta,
            "snr_grid": SNR_GRID,
            "wall_clock_total_seconds": time.time() - t_start,
        },
    }

    out_path = PLOT_DATA_DIR / "beam_search_ir_vn_data.json"
    with open(out_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Wrote {out_path}")
    print(f"TOTAL wall-clock: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", default=None, help=f"Run one config in-process: one of {ALL_RUN_NAMES}")
    args = parser.parse_args()
    if args.run:
        run_single_config(args.run)
    else:
        main_orchestrate()
