"""
Minimal usage example for the new Beam_Search_config / Beam_Search_Decoder
(channel-code-lib2, branch beam_search). Mirrors bp_example.py's structure.

See channel-code-lib2's BEAM_SEARCH_PHASE0_PLAN.md for the full design
writeup of what this is doing underneath (chunk -> split -> score -> prune
adaptive ensemble decoding).

NOTE: this needs `uv sync` to have picked up channel-code-lib2's
`beam_search` branch (see this repo's pyproject.toml) -- i.e. it only works
once that branch has been pushed and synced, same as any other dependency
bump.
"""
import numpy as np

import channel_code_lib2

from Codes.generate_5G_LDPC import (
    generate_5G_LDPC,
    get_final_matrices_and_message_bit_pucturing,
)

import show_results

# - Generate 5G LDPC code (same toy setup as bp_example.py)
n_ = 132
k_ = 66
H, p, s, Z, BG = generate_5G_LDPC(2, k_, n_)
H, G, k, n, message_puncturing = get_final_matrices_and_message_bit_pucturing(H, s, p)

## Beam-search config: same BP_config fields (H, early_stopping, max_iterations,
## cn_update_type, scheduling_type, norm_factor, affine_offset) PLUS the 5
## beam-search hyperparameters.
beam_cfg = channel_code_lib2.Beam_Search_config(H)
beam_cfg.early_stopping = True
beam_cfg.max_iterations = 50          # total iteration budget across all rounds/splits
beam_cfg.cn_update_type = "spa"
beam_cfg.scheduling_type = "flooding"  # only supported value in this phase -- see plan doc

beam_cfg.iters_per_round = 10          # BP iterations per round before a split is attempted
beam_cfg.beam_width = 16               # max active (non-frozen) paths
beam_cfg.target_converged = 1          # stop as soon as 1 path converges
beam_cfg.warm_start = True             # child/parent continue from live BP state at split time
# eligible_vns left empty => all VNs eligible (fine here: this code has no
# punctured/shortened bits left after get_final_matrices_and_message_bit_pucturing
# strips them out of H -- for a code that still has punctured columns in H,
# set beam_cfg.eligible_vns explicitly, excluding them).

beam_cfg.affine_offset = np.zeros(H.shape[0], dtype=int)

sim_beam = channel_code_lib2.Simulation_Env(k, n, "all")
sim_beam.puncturing(message_puncturing)
sim_beam.all_zero_init(beam_cfg)
sim_beam.get_error_rates(np.linspace(1, 4, 7))
FER_beam = sim_beam.error_rates["FER-SNR"]
print("Beam search FER:", FER_beam)

## Plain BP baseline, identical BP-level settings, for comparison.
bp_cfg = channel_code_lib2.BP_config(H)
bp_cfg.early_stopping = True
bp_cfg.max_iterations = 50
bp_cfg.cn_update_type = "spa"
bp_cfg.scheduling_type = "flooding"
bp_cfg.affine_offset = np.zeros(H.shape[0], dtype=int)

sim_bp = channel_code_lib2.Simulation_Env(k, n, "all")
sim_bp.puncturing(message_puncturing)
sim_bp.all_zero_init(bp_cfg)
sim_bp.get_error_rates(np.linspace(1, 4, 7))
FER_bp = sim_bp.error_rates["FER-SNR"]
print("Plain BP FER:", FER_bp)

show_results.plot_error_rates((FER_beam, "Beam search"), (FER_bp, "Plain BP"))
