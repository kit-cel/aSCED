"""Orchestrate the greedy RL/QC-block selection search for one aSCED ensemble
(asced48 or asced384), end to end: generate each round's manifest, submit a
SLURM array job (sweeps/run_greedy_search_array.sbatch), poll sacct until the
whole array is terminal, read back the per-candidate FER results, rank
candidates at each SNR, combine ranks (sum, lower is better), pick the
round's winner, fix it, and move to the next round. Resumable: round state is
persisted to sweeps/greedy_search_state_<ensemble>.json after every round, so
a restart continues from the first incomplete round.

Usage:
    python sweeps/run_greedy_search.py <ensemble_name: asced48|asced384> <num_rounds> [array_concurrency]

This is meant to be run as a single long-lived background process (polls
sacct every ~60s internally) - submit it and wait, don't poll it externally.
"""

import csv
import json
import os
import subprocess
import sys
import time
from collections import Counter

if len(sys.argv) not in (3, 4):
    raise ValueError(
        "Usage: python sweeps/run_greedy_search.py <ensemble_name: asced48|asced384> "
        "<num_rounds> [array_concurrency]"
    )

ENSEMBLE = sys.argv[1].lower()
NUM_ROUNDS = int(sys.argv[2])
ARRAY_CONCURRENCY = int(sys.argv[3]) if len(sys.argv) > 3 else 10

if ENSEMBLE not in ("asced48", "asced384"):
    raise ValueError(f"Unknown ensemble_name '{ENSEMBLE}'")

SNR_POINTS = [1.8, 3.3]
ALL_BLOCKS = list(range(34))
STATE_PATH = f"sweeps/greedy_search_state_{ENSEMBLE}.json"

TERMINAL_STATES_PREFIXES = (
    "COMPLETED",
    "FAILED",
    "CANCELLED",
    "TIMEOUT",
    "OUT_OF_MEMORY",
    "NODE_FAIL",
    "DEADLINE",
    "BOOT_FAIL",
)


def snr_tag(s):
    return f"{s:g}"


def fixed_tag(blocks):
    return "none" if not blocks else "-".join(str(b) for b in sorted(blocks))


def result_path(fixed_blocks, candidate, snr):
    return (
        f"RESULTS/greedy_search/{ENSEMBLE}/round_results/"
        f"fixed_{fixed_tag(fixed_blocks)}_cand_{candidate}_snr_{snr_tag(snr)}.json"
    )


def load_state():
    if os.path.exists(STATE_PATH):
        with open(STATE_PATH) as f:
            return json.load(f)
    return {"ensemble": ENSEMBLE, "fixed_blocks": [], "rounds": []}


def save_state(state):
    with open(STATE_PATH, "w") as f:
        json.dump(state, f, indent=2)
    print(f"[state saved] {STATE_PATH}")


def write_manifest(path, fixed_blocks, candidates):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ensemble_name", "fixed_blocks", "candidate_block", "snr"])
        for c in candidates:
            for snr in SNR_POINTS:
                w.writerow([ENSEMBLE, fixed_tag(fixed_blocks), c, snr])


def submit_array(manifest_path, n_rows):
    cmd = [
        "sbatch",
        "--parsable",
        f"--array=1-{n_rows}%{min(ARRAY_CONCURRENCY, n_rows)}",
        "sweeps/run_greedy_search_array.sbatch",
        manifest_path,
    ]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True)
    jobid = out.stdout.strip()
    print(f"Submitted {manifest_path} ({n_rows} tasks) -> job {jobid}")
    return jobid


def poll_job(jobid):
    while True:
        res = subprocess.run(
            ["sacct", "-j", jobid, "--format=State", "-X", "--noheader"],
            capture_output=True,
            text=True,
        )
        lines = [l.strip() for l in res.stdout.splitlines() if l.strip()]
        counts = Counter(lines)
        print(f"[poll] job {jobid}: {dict(counts)}")
        non_terminal = [
            l for l in lines if not any(l.startswith(p) for p in TERMINAL_STATES_PREFIXES)
        ]
        if lines and not non_terminal:
            print(f"Job {jobid} fully terminal: {dict(counts)}")
            return counts
        time.sleep(60)


def read_results(fixed_blocks, candidates):
    data = {}
    for c in candidates:
        per_snr = {}
        for snr in SNR_POINTS:
            p = result_path(fixed_blocks, c, snr)
            if os.path.exists(p):
                with open(p) as f:
                    per_snr[snr] = json.load(f)["fer"]
            else:
                per_snr[snr] = None
        data[c] = per_snr
    return data


def rank_candidates(results):
    """Per-SNR rank (1=best/lowest FER), then combined = sum of per-SNR ranks."""
    ranking = {c: {} for c in results}
    for snr in SNR_POINTS:
        vals = [(c, results[c][snr]) for c in results if results[c][snr] is not None]
        vals.sort(key=lambda x: x[1])
        for rank, (c, _fer) in enumerate(vals, start=1):
            ranking[c][snr] = rank
    combined = {}
    for c in results:
        rs = [ranking[c].get(snr) for snr in SNR_POINTS]
        combined[c] = None if any(r is None for r in rs) else sum(rs)
    return ranking, combined


def run_round(round_idx, fixed_blocks):
    candidates = [b for b in ALL_BLOCKS if b not in fixed_blocks]
    manifest_path = f"sweeps/greedy_round_manifest_{ENSEMBLE}_r{round_idx}.csv"
    write_manifest(manifest_path, fixed_blocks, candidates)
    jobid = submit_array(manifest_path, len(candidates) * len(SNR_POINTS))
    poll_job(jobid)

    results = read_results(fixed_blocks, candidates)
    missing = [(c, snr) for c in candidates for snr in SNR_POINTS if results[c][snr] is None]
    retry_jobid = None
    if missing:
        print(f"Round {round_idx}: {len(missing)} missing result(s), retrying once: {missing}")
        retry_manifest = f"sweeps/greedy_round_manifest_{ENSEMBLE}_r{round_idx}_retry.csv"
        with open(retry_manifest, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["ensemble_name", "fixed_blocks", "candidate_block", "snr"])
            for c, snr in missing:
                w.writerow([ENSEMBLE, fixed_tag(fixed_blocks), c, snr])
        retry_jobid = submit_array(retry_manifest, len(missing))
        poll_job(retry_jobid)
        results = read_results(fixed_blocks, candidates)

    ranking, combined = rank_candidates(results)
    still_missing = [c for c, v in combined.items() if v is None]
    usable = {c: v for c, v in combined.items() if v is not None}
    if not usable:
        raise RuntimeError(f"Round {round_idx}: no usable results at all! results={results}")
    if still_missing:
        print(f"WARNING round {round_idx}: candidates with missing data, excluded from pick: {still_missing}")

    winner = min(usable, key=lambda c: (usable[c], c))
    print(f"Round {round_idx} winner: block {winner} (combined_rank={usable[winner]})")

    return {
        "round": round_idx,
        "fixed_blocks_before": list(fixed_blocks),
        "candidates": candidates,
        "fer": {str(c): results[c] for c in candidates},
        "rank": {str(c): ranking.get(c) for c in candidates},
        "combined_rank": {str(c): combined.get(c) for c in candidates},
        "still_missing": still_missing,
        "winner": winner,
        "jobid": jobid,
        "retry_jobid": retry_jobid,
    }


def main():
    state = load_state()
    fixed_blocks = state["fixed_blocks"]
    start_round = len(state["rounds"]) + 1
    print(f"Starting {ENSEMBLE} greedy search from round {start_round} (fixed so far: {fixed_blocks})")

    for round_idx in range(start_round, NUM_ROUNDS + 1):
        record = run_round(round_idx, fixed_blocks)
        state["rounds"].append(record)
        fixed_blocks = fixed_blocks + [record["winner"]]
        state["fixed_blocks"] = fixed_blocks
        save_state(state)

    print(f"DONE. Final greedy block order for {ENSEMBLE}: {fixed_blocks}")


if __name__ == "__main__":
    main()
