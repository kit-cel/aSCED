# pylint: disable=invalid-name
"""Pure scheduling/iteration-count study: flooding vs. row-layered BP
scheduling (two visiting-order variants), decoupled from the sequential-
selector/stopping work in `reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py`.

Same RL/QC-block aSCED ensemble construction (aSCED-48, aSCED-384) for the 5G
LDPC code C_5G(132,66)/Zc=11, reused from that script, but:
  - only the `full_parallel` selector is exercised (no fixed_sequential/
    syndrome_sequential/mconverged stopping here - one clean variable at a
    time);
  - a third variant, standalone "nmsa" (BP directly on the base PCM `H`, no
    aSCED ensemble at all), isolates the pure scheduling effect from any
    ensemble complexity;
  - every BP_config's `scheduling_type`/`max_iterations` and (for row-layered)
    `schedule`/`nodes_in_layer` are driven by this script's CLI args instead
    of being hardcoded to flooding/32 as in the original scripts.

Row-layered schedule construction (see sweeps/STATUS.md and the task prompt
for the full channel-code-lib2 source derivation this is based on):
  - `BP_Row_Layered::decode()` (src/Decoder/BP/BP_Scheduling/BP_Row_Layered.cpp)
    iterates `for layer in schedule: for j in nodes_in_layer[layer]: ...`
    once per BP iteration - `schedule` is a *permutation of layer indices*
    (visiting order), `nodes_in_layer[layer]` is that layer's list of check-
    node (H_aux row) indices. `BP_Row_Layered::init()` copies both verbatim
    from `BP_config` - there is no auto-generation from `Z` (confirmed: `Z`
    is bound to Python but never read anywhere in BP_Decoder.cpp/
    BP_Row_Layered.cpp/Decoder_Factory.h - it's a dead field).
  - For a given H_aux = vstack(H, appended_rows) (H's `m` rows are always a
    multiple of Zc by construction - m = Zc*(bg_cn + harq_part)), one layer
    per Zc=11-row block of H's original rows (row order 0, 1, 2, ...), plus
    (if there are appended rows, i.e. every aSCED path, but NOT the standalone
    "nmsa" variant which has no appended rows) one final layer for the
    appended-rows block (already exactly one contiguous row-group per batch
    by construction - Delta in {2,3} for aSCED-48, {5,6} for aSCED-384).
  - "natural": layers visited in matrix row order (H's blocks first, appended
    last). "appended_first": appended-rows layer visited first, then H's
    blocks in original order. For "nmsa" (no appended layer to reorder) the
    two variants coincide by construction - run anyway for a uniform 3x3x2x7
    grid, not because a real difference is expected there.
"""

import numpy as np
import os
import sys

import galois

gf2 = galois.GF2

import channel_code_lib2

if len(sys.argv) != 7:
    raise ValueError(
        "Usage: python scheduling_study.py "
        "<decoder_variant: nmsa|asced48|asced384> <n_simul> <snr_start> <snr_end> "
        "<scheduling_config: flooding|row_layered_natural|row_layered_appended_first> "
        "<max_iterations>"
    )

decoder_variant = sys.argv[1].lower()
n_simul = int(sys.argv[2])
snr_start = float(sys.argv[3])
snr_end = float(sys.argv[4])
scheduling_config = sys.argv[5].lower()
max_iterations = int(sys.argv[6])

SCHEDULING_CONFIGS = ("flooding", "row_layered_natural", "row_layered_appended_first")
if scheduling_config not in SCHEDULING_CONFIGS:
    raise ValueError(f"Unknown scheduling_config '{scheduling_config}', must be one of {SCHEDULING_CONFIGS}")

# Include the end point
sim_regime = np.arange(snr_start, snr_end + 0.25, 0.5)

remove = 4  # k=66 code rather than k=110
bg_vn = 12 - remove  # bg2
bg_cn = 4  # bg2
Zc = 11

results_dir = f"RESULTS/scheduling_study_zc11_r{remove}_n{n_simul}"

splitting_pattern = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [2, 4, 6, 8], [5]]

ASCED_VARIANTS = {
    "asced48": (splitting_pattern[1], 2),
    "asced384": (splitting_pattern[2], 4),
}


def build_row_layered_layout(m_H, Zc, m_appended, order):
    """Build (nodes_in_layer, schedule) for row-layered scheduling on an
    H_aux = vstack(H, appended_rows) with H.shape[0] == m_H and
    appended_rows.shape[0] == m_appended (0 for the standalone nmsa variant,
    which has no appended rows at all).

    One layer per Zc-row block of H's rows, in original row order, plus (if
    m_appended > 0) one final layer for the appended-rows block. `order` is
    "natural" (H's blocks first, appended last) or "appended_first" (appended
    layer first, then H's blocks in original order).
    """
    assert m_H % Zc == 0, f"H's row count {m_H} is not a multiple of Zc={Zc}"
    num_h_layers = m_H // Zc

    nodes_in_layer = [list(range(b * Zc, (b + 1) * Zc)) for b in range(num_h_layers)]
    has_appended = m_appended > 0
    if has_appended:
        nodes_in_layer.append(list(range(m_H, m_H + m_appended)))

    if order == "natural" or not has_appended:
        schedule = list(range(len(nodes_in_layer)))
    elif order == "appended_first":
        appended_layer_idx = len(nodes_in_layer) - 1
        schedule = [appended_layer_idx] + list(range(num_h_layers))
    else:
        raise ValueError(f"Unknown order '{order}'")

    return nodes_in_layer, schedule


def create_bp_config(H_aux, m_H, m_appended):
    """Create a BP configuration, applying this run's scheduling_config/
    max_iterations. m_H/m_appended are used only to build the row-layered
    layout (ignored for flooding)."""
    cfg = channel_code_lib2.BP_config(H_aux)
    cfg.early_stopping = True
    cfg.max_iterations = max_iterations
    cfg.cn_update_type = "msa"
    cfg.norm_factor = 0.75

    if scheduling_config == "flooding":
        cfg.scheduling_type = "flooding"
    else:
        order = "natural" if scheduling_config == "row_layered_natural" else "appended_first"
        nodes_in_layer, schedule = build_row_layered_layout(m_H, Zc, m_appended, order)
        cfg.scheduling_type = "row_layered"
        cfg.nodes_in_layer = nodes_in_layer
        cfg.schedule = schedule
        cfg.Z = Zc  # unused by BP_Row_Layered (dead field, confirmed from source) - set for provenance only

    return cfg


if decoder_variant not in ("nmsa",) and decoder_variant not in ASCED_VARIANTS:
    raise ValueError(f"Unknown decoder variant '{decoder_variant}'")

message_bit_pucturing = np.arange(2 * Zc, dtype=int)
auto_save = True
use_all_zero = False

# Same PCM construction as reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py.
code = np.load("Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz")
H_full = code["h"].astype(int)
H_full = np.delete(H_full, np.arange(6 * Zc, 10 * Zc), axis=1)

block_offset = 0
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

print(f"Simulating number vns={n}, n={number_vns - len(message_bit_pucturing)}, k={k}")

candidate_rows = H_full[m:, : number_vns_start + harq_part * Zc]
number_additional_cyclic_blocks = candidate_rows.shape[0] // Zc
assert number_additional_cyclic_blocks * Zc == candidate_rows.shape[0]

m_H = H.shape[0]
assert m_H == m


def binary_vectors_in_suffix(number_rows, number_last_entries):
    v = [0] * number_rows
    for mask in range(1, 1 << number_last_entries):
        for i in range(number_last_entries):
            v[number_rows - number_last_entries + i] = (mask >> i) & 1
        yield v.copy()


if not use_all_zero:
    enc_cfg = channel_code_lib2.PCM_Encoder_config(H, k, n)


def create_asced_config(H, candidate_rows, block_offset, Zc, split_pattern, num_used_blocks, use_all_zero):
    """Construct all decoder paths for one aSCED configuration (full_parallel
    only - no selector/stopping policy is set here, see module docstring).

    Mirrors reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py's
    create_asced_config() exactly for the H_aux/batch construction, but each
    path's BP_config gets this run's scheduling_config/max_iterations (every
    path in a batch shares the same H_aux, hence the same row-layered
    nodes_in_layer/schedule).
    """
    batches = []

    for block in range(num_used_blocks):
        assert block + block_offset <= number_additional_cyclic_blocks
        row_block = candidate_rows[(block + block_offset) * Zc : (block + block_offset + 1) * Zc]

        row_segments = [row_block] if len(split_pattern) == 0 else np.split(row_block, split_pattern)

        for rows in row_segments:
            H_aux = gf2(np.vstack((H, rows)))
            m_aux = H_aux.shape[0]
            m_appended = rows.shape[0]

            batch = [create_bp_config(H_aux, m_H, m_appended)]  # zero affine offset

            if not use_all_zero:
                rank_difference = rows.shape[0]
                for offset in binary_vectors_in_suffix(m_aux, rank_difference):
                    cfg = create_bp_config(H_aux, m_H, m_appended)
                    cfg.affine_offset = offset
                    batch.append(cfg)

            batches.append(batch)

    path_configs = [cfg for batch in batches for cfg in batch]
    print("Simulated num. aSCED paths:", len(path_configs), "in", len(batches), "batches")

    return channel_code_lib2.Ensemble_config(H, path_configs), len(path_configs)


def make_sim(save_dir_name):
    sim = channel_code_lib2.Simulation_Env(k, n, "all")
    sim.target_errors = 200
    sim.max_transmissions = int(2e6)
    slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm_cpus is not None:
        sim.num_threads = int(slurm_cpus)
    sim.auto_save = auto_save
    # Per-(config, SNR) directory - sidesteps the known concurrent-write race
    # on the C++ side's non-atomic save-dir JSON read-modify-write (see
    # sweeps/STATUS.md's "Race condition found and fixed" section).
    snr_tag = f"{sim_regime[0]:g}"
    sim.save_dir = results_dir + "/" + save_dir_name + "/" + f"snr_{snr_tag}"
    sim.puncturing(message_bit_pucturing)
    return sim


def run_nmsa(save_dir_name):
    nmsa_config = create_bp_config(H, m_H, 0)
    sim = make_sim(save_dir_name)
    if use_all_zero:
        sim.all_zero_init(nmsa_config)
    else:
        sim.init(enc_cfg, nmsa_config, use_all_zero)
    sim.get_error_rates(sim_regime)
    print(save_dir_name, "finished")
    print(sim.error_rates["FER-SNR"])
    return sim.error_rates["FER-SNR"]


def run_asced(save_dir_name, split_pattern, num_used_blocks):
    ensemble_cfg, ensemble_size = create_asced_config(
        H=H,
        candidate_rows=candidate_rows,
        block_offset=block_offset,
        Zc=Zc,
        split_pattern=split_pattern,
        num_used_blocks=num_used_blocks,
        use_all_zero=use_all_zero,
    )

    print("starting", save_dir_name)
    sim = make_sim(save_dir_name)

    if use_all_zero:
        sim.all_zero_init(ensemble_cfg)
    else:
        sim.init(enc_cfg, ensemble_cfg, use_all_zero)

    sim.get_error_rates(sim_regime)

    print(save_dir_name, "finished")
    print(sim.error_rates["FER-SNR"])

    return sim.error_rates["FER-SNR"]


save_name = f"{decoder_variant}_{scheduling_config}_maxiter{max_iterations}"

if decoder_variant == "nmsa":
    FER_results = {save_name: run_nmsa(save_name)}
else:
    split_pattern, num_used_blocks = ASCED_VARIANTS[decoder_variant]
    FER_results = {save_name: run_asced(save_name, split_pattern, num_used_blocks)}
