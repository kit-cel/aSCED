# pylint: disable=invalid-name
"""Per-candidate evaluation driver for the greedy RL/QC-block selection search
on the Zc=6 codes C_5G(78,60) and C_5G(180,60) (Fig. 5 of the paper).

Generalizes `create_asced_config()` from
reproduce_fig5_scatter_plot_zc6_5G_LDPC.py (which always uses
`num_used_blocks` *contiguous* blocks starting at `block_offset=0`) to an
arbitrary LIST of block indices into the candidate RL/QC-block pool (37
blocks for n_simul=78, 20 for n_simul=180 - verified via sweeps/STATUS.md's
greedy-search task description and a direct shape check). Reuses the exact
same H load / candidate_rows / affine-offset path construction logic as the
reference script so results are directly comparable.

Evaluated using ONLY the split3 (Delta=3) sub-splitting pattern
(splitting_pattern[2] = [3] in the reference script - splits each 6-row block
into two Delta=3 sub-batches) - this is a pure block-SELECTION search, not a
re-optimization of within-block sub-batch selection.

This script is invoked once per (n_simul, fixed_blocks, candidate_block, snr)
row of a round's manifest, with a CHEAP evaluation budget (target_errors=50,
max_transmissions=2e5) - for ranking candidates only, not final production
numbers.

Usage:
    python sweeps/greedy_block_search_zc6.py <n_simul: 78|180> \
        <fixed_blocks: comma-separated ints, or "none"> <candidate_block: int> \
        <snr>

Writes a small JSON result file to a deterministic path (independent of fixed
block order - the ensemble is a SET of blocks, order doesn't affect the
Tanner graph/FER, only the greedy *discovery order* we record separately):
    RESULTS/greedy_search_zc6/n<n_simul>/round_results/fixed_<a-b-c>_cand_<c>_snr_<snr>.json
("fixed_none" when the fixed set is empty, i.e. round 1).

Detailed per-run decoder stats (for provenance, not read by the orchestrator)
go to a separate unique save_dir per (n_simul, block set, snr), same
one-dir-per-(config,SNR) convention used throughout this effort, to avoid the
known concurrent-JSON-write race (see STATUS.md).
"""

import json
import os
import sys

import numpy as np
import galois

gf2 = galois.GF2

import channel_code_lib2

if len(sys.argv) != 5:
    raise ValueError(
        "Usage: python sweeps/greedy_block_search_zc6.py "
        "<n_simul: 78|180> <fixed_blocks: comma-separated ints or 'none'> "
        "<candidate_block: int> <snr>"
    )

n_simul = int(sys.argv[1])
fixed_blocks_arg = sys.argv[2].strip().lower()
candidate_block = int(sys.argv[3])
snr = float(sys.argv[4])

if fixed_blocks_arg in ("none", ""):
    fixed_blocks = []
else:
    fixed_blocks = [int(x) for x in fixed_blocks_arg.split("-") if x != ""]

if n_simul not in (78, 180):
    raise ValueError(f"Unsupported n_simul={n_simul}, expected 78 or 180")

SPLIT_PATTERN = [3]  # Delta=3 sub-splitting, per the task's explicit scope

bg_vn = 12  # bg2
bg_cn = 4  # bg2
Zc = 6

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
use_all_zero = False

EXPECTED_NUM_BLOCKS = {78: 37, 180: 20}

# ---- Same PCM construction as reproduce_fig5_scatter_plot_zc6_5G_LDPC.py ----
code = np.load("Codes/TCOM_aSCED/5G_zc=6/max_rank_5G_zc=6.npz")
H_full = code["h"].astype(int)

number_vn_simul = n_simul + 2 * Zc
number_vns_start = bg_vn * Zc + 2 * Zc
harq_part = (number_vn_simul - number_vns_start) // Zc
m = Zc * bg_cn + harq_part * Zc
number_vns = number_vns_start + harq_part * Zc
assert number_vn_simul == number_vns

H = gf2(H_full[:m, :number_vns])
G = np.array(gf2(H[:m, :number_vns]).null_space()).astype(int)
k, n = G.shape
assert n == number_vns

candidate_rows = H_full[m:, : number_vns_start + harq_part * Zc]
number_additional_cyclic_blocks = candidate_rows.shape[0] // Zc
assert number_additional_cyclic_blocks * Zc == candidate_rows.shape[0]
assert number_additional_cyclic_blocks == EXPECTED_NUM_BLOCKS[n_simul], (
    f"Expected {EXPECTED_NUM_BLOCKS[n_simul]} candidate blocks for n_simul={n_simul}, "
    f"got {number_additional_cyclic_blocks} - did the code/config change upstream?"
)

for b in fixed_blocks + [candidate_block]:
    assert 0 <= b < number_additional_cyclic_blocks, f"block index {b} out of range [0,{number_additional_cyclic_blocks})"

enc_cfg = channel_code_lib2.PCM_Encoder_config(H, k, n)


def binary_vectors_in_suffix(number_rows, number_last_entries):
    v = [0] * number_rows
    for mask in range(1, 1 << number_last_entries):
        for i in range(number_last_entries):
            v[number_rows - number_last_entries + i] = (mask >> i) & 1
        yield v.copy()


def create_bp_config(H_aux):
    """Create a BP configuration with the default decoder settings."""
    cfg = channel_code_lib2.BP_config(H_aux)
    cfg.early_stopping = True
    cfg.max_iterations = 32
    cfg.cn_update_type = "msa"
    cfg.scheduling_type = "flooding"
    cfg.norm_factor = 0.75
    return cfg


def create_asced_config(H, candidate_rows, block_indices, Zc, split_pattern, use_all_zero):
    """Construct all decoder paths for an aSCED configuration built from an
    ARBITRARY list of block indices (generalization of the reference script's
    `create_asced_config`, which only supports num_used_blocks contiguous
    blocks starting at block_offset).
    """
    batches = []

    for block in block_indices:
        row_block = candidate_rows[block * Zc : (block + 1) * Zc]

        row_segments = [row_block] if len(split_pattern) == 0 else np.split(row_block, split_pattern)

        for rows in row_segments:
            H_aux = gf2(np.vstack((H, rows)))
            m_aux = H_aux.shape[0]

            batch = [create_bp_config(H_aux)]  # zero affine offset

            if not use_all_zero:
                rank_difference = rows.shape[0]
                for offset in binary_vectors_in_suffix(m_aux, rank_difference):
                    cfg = create_bp_config(H_aux)
                    cfg.affine_offset = offset
                    batch.append(cfg)

            batches.append(batch)

    path_configs = [cfg for batch in batches for cfg in batch]

    return channel_code_lib2.Ensemble_config(H, path_configs), len(path_configs)


block_set = sorted(set(fixed_blocks + [candidate_block]))
block_set_tag = "-".join(str(b) for b in block_set)
fixed_tag = "none" if not fixed_blocks else "-".join(str(b) for b in sorted(fixed_blocks))
snr_tag = f"{snr:g}"

results_root = f"RESULTS/greedy_search_zc6/n{n_simul}"
result_json_path = (
    f"{results_root}/round_results/fixed_{fixed_tag}_cand_{candidate_block}_snr_{snr_tag}.json"
)
os.makedirs(os.path.dirname(result_json_path), exist_ok=True)

ensemble_cfg, ensemble_size = create_asced_config(
    H=H,
    candidate_rows=candidate_rows,
    block_indices=block_set,
    Zc=Zc,
    split_pattern=SPLIT_PATTERN,
    use_all_zero=use_all_zero,
)

print(
    f"[n_simul={n_simul}] fixed={fixed_blocks} candidate={candidate_block} "
    f"-> block_set={block_set}, ensemble_size={ensemble_size}, snr={snr}"
)

sim = channel_code_lib2.Simulation_Env(k, n, "all")
sim.target_errors = 50
sim.max_transmissions = int(2e5)
slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
if slurm_cpus is not None:
    sim.num_threads = int(slurm_cpus)
sim.auto_save = auto_save
sim.save_dir = f"{results_root}/detailed/blocks_{block_set_tag}/snr_{snr_tag}"

sim.puncturing(message_bit_pucturing)
sim.init(enc_cfg, ensemble_cfg, use_all_zero)

sim.get_error_rates(np.array([snr]))

fer_snr = sim.error_rates["FER-SNR"]
print("FER-SNR raw:", fer_snr)
if isinstance(fer_snr, dict):
    assert len(fer_snr) == 1, f"expected exactly 1 SNR entry, got {fer_snr}"
    fer = float(next(iter(fer_snr.values())))
else:
    fer = float(fer_snr[0])

print(f"finished block_set={block_set} snr={snr} -> FER={fer}")

with open(result_json_path, "w") as f:
    json.dump(
        {
            "n_simul": n_simul,
            "fixed_blocks": sorted(fixed_blocks),
            "candidate_block": candidate_block,
            "block_set": block_set,
            "ensemble_size": ensemble_size,
            "snr": snr,
            "fer": fer,
        },
        f,
        indent=2,
    )

print("wrote", result_json_path)
