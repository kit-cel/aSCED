# pylint: disable=invalid-name
"""Build an IR-VN-extended PCM for the aSCED-48 TCOM-paper code family
(BG2, Zc=11) -- beam-search-decoder "IR VN splitting" work.

Design (see the beam_search task write-up this module was built for):
Incremental-Redundancy (IR) VNs are treated as ORDINARY extended-PCM
columns, exactly like puncturing -- NOT via this codebase's existing
``num_syndrome_vns``/``syndrome_vn_pattern``/``H0`` "soft-aSCED syndrome VN"
mechanism (see ``soft_asced_delta_gt1.make_soft_bp_config`` in
``/home/pj9034/aSCED/.claude/worktrees/soft_asced_wireless/soft_asced_delta_gt1.py``,
lines ~172-191). That mechanism checks convergence only against the
ORIGINAL un-extended PCM ``H0`` -- it deliberately EXCLUDES its guess
columns from the real validity check. IR VNs here are REAL codeword
columns: the only convergence check is the real, full one,
``H_ext @ codeword == 0`` (via ``Simulation_Env``/``PCM_Encoder_config``
built from ``H_ext`` directly) -- there is no separate "H0".

This module only REUSES that other function's row/column-append PATTERN
(append ``Delta`` new check rows, then one new degree-1 VN column per new
row via an identity block) -- it does not import or touch any of that
function's BP_config / num_syndrome_vns / H0 fields.

Base-PCM construction is copied verbatim (same column slicing convention)
from ``reproduce_fig3_RL_zc11_asced_5G_LDPC.py`` (this project's asced-48/
NMSA reference script) so the IR-VN-extended code and the asced-48/NMSA
baselines being compared against are built from the IDENTICAL base H:

    H = H_full[:m, :number_vns]                      # base PCM
    candidate_rows = H_full[m:, :number_vns]          # HARQ "future" rows

``H_full`` is loaded from
``Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz`` (key ``"h"``): the full
BG2 PCM lifted with Zc=11, 6 "parity-puncturing" base-graph VN columns
already deleted (``np.delete(H_full, np.arange(6*Zc, 10*Zc), axis=1)``).
Rows beyond ``m`` are candidate HARQ/IR check rows, sliced in whole
Zc-sized blocks (see ``apply_harq_pcm_shortening.py``'s docstring for the
same ``H_full[:m,:]``/``H_full[m:,:]`` convention, and
``generate_5G_LDPC.py`` / that docstring confirming
``galois.GF2(H).null_space()`` doesn't care how H was built -- it
generalizes unmodified to this IR-VN-extended H too).

Default ``num_ad_blocks=2`` deliberately mirrors the asced-48 reference's
own candidate-row footprint 1:1 (``num_used_blocks=2`` in
``create_asced_config``/``get_row_segments``, same ``Zc=11``,
``block_offset=0`` default) -- i.e. this module's default IR-VN extension
consumes EXACTLY the same 22 candidate rows (blocks 0-1) that asced-48
already uses for its ensemble splitting, just turned into 22 new REAL
degree-1 VN columns instead of ensemble/affine-offset guesses. With the
combined puncturing set (2*Zc systematic + Delta=22 IR-VN columns), the
transmitted length works out to n_ext - len(combined_puncturing) = 176 - 44
= 132 -- the SAME transmitted length as the asced-48/NMSA references'
n_simul=132, so FER-vs-SNR curves are directly comparable (same channel
uses per codeword); the IR-VN columns add decoder-side Tanner-graph
structure (extra, currently-erased check constraints) without changing
what's actually sent over the channel.

Run this file directly for a self-contained sanity check (shapes,
H_ext @ G.T == 0, k unchanged, combined-puncturing non-contiguity).
"""

import numpy as np
import galois

gf2 = galois.GF2

# ---- aSCED-48/NMSA reference code-family parameters (reproduce_fig3_RL_zc11_asced_5G_LDPC.py) ----
DEFAULT_ZC = 11
DEFAULT_REMOVE = 4
DEFAULT_BG_VN = 12  # bg2
DEFAULT_BG_CN = 4  # bg2
DEFAULT_N_SIMUL = 132  # "k=66" code (reproduce_fig3's remove=4 deployment)
DEFAULT_CODE_PATH = "Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz"


def load_base_pcm_and_candidates(
    n_simul=DEFAULT_N_SIMUL,
    remove=DEFAULT_REMOVE,
    bg_vn=DEFAULT_BG_VN,
    bg_cn=DEFAULT_BG_CN,
    Zc=DEFAULT_ZC,
    code_path=DEFAULT_CODE_PATH,
):
    """Reproduce reproduce_fig3_RL_zc11_asced_5G_LDPC.py's exact base-PCM
    slice (``H_full[:m, :number_vns]``) and its trailing HARQ
    candidate-row block (``H_full[m:, :number_vns]``), so IR-VN extension
    starts from the IDENTICAL base code as the asced-48/NMSA references
    (n=154, k=66, m=88 for the default n_simul=132).

    Returns
    -------
    H : np.ndarray[int], shape (m, number_vns)
    candidate_rows : np.ndarray[int], shape (num_candidate_blocks*Zc, number_vns)
    number_vns : int
    num_candidate_blocks : int
    """
    bg_vn_eff = bg_vn - remove

    code = np.load(code_path)
    H_full = code["h"].astype(int)
    # Delete the 4 "parity-puncturing" base-graph VN columns not used by
    # this project's remove=4 deployment (verbatim from the reference script).
    H_full = np.delete(H_full, np.arange(6 * Zc, 10 * Zc), axis=1)

    number_vn_simul = n_simul + 2 * Zc
    number_vns_start = bg_vn_eff * Zc + 2 * Zc
    harq_part = (number_vn_simul - number_vns_start) // Zc
    m = Zc * bg_cn + harq_part * Zc
    number_vns = number_vns_start + harq_part * Zc
    assert number_vn_simul == number_vns, (
        f"n_simul={n_simul} does not land on a whole number of Zc={Zc} "
        "HARQ blocks beyond number_vns_start -- pick a different n_simul."
    )

    H = np.array(gf2(H_full[:m, :number_vns])).astype(int)
    candidate_rows = H_full[m:, :number_vns]
    num_candidate_blocks = candidate_rows.shape[0] // Zc
    assert num_candidate_blocks * Zc == candidate_rows.shape[0]

    return H, candidate_rows, number_vns, num_candidate_blocks


def build_ir_vn_extended_pcm(H, candidate_rows, Zc, num_ad_blocks, block_offset=0):
    """Extend H with ``num_ad_blocks*Zc`` new check rows (taken from
    ``candidate_rows[block_offset*Zc : (block_offset+num_ad_blocks)*Zc]``)
    and one new degree-1 IR-VN column per new row (identity block),
    following the row-append + identity-column-append PATTERN of
    ``soft_asced_delta_gt1.make_soft_bp_config`` -- but building a PLAIN
    PCM only (no ``num_syndrome_vns``/``H0``/``BP_config`` fields; see
    module docstring for why).

    The new IR-VN columns are REAL codeword columns: there is no held-out
    "original H0" here, unlike the syndrome-VN mechanism -- the real
    convergence/validity check is simply ``H_ext @ codeword == 0``.

    Parameters
    ----------
    H : array-like, shape (m, n)
    candidate_rows : array-like, shape (num_candidate_blocks*Zc, n)
    Zc : int
    num_ad_blocks : int
        Number of Zc-sized candidate-row blocks to append as new check rows.
    block_offset : int, default 0
        Which candidate block to start from (0 = the first/lowest-index
        HARQ block, matching ``create_asced_config``'s ``block_offset``
        convention).

    Returns
    -------
    H_ext : np.ndarray[int], shape (m + Delta, n + Delta)
    ir_vn_indices : np.ndarray[int], shape (Delta,)
        Column indices of the newly appended IR-VN columns in H_ext
        (``np.arange(n, n + Delta)``).
    Delta : int
        Number of new rows/columns appended (``num_ad_blocks * Zc``).
    """
    H = np.asarray(H, dtype=int)
    n = H.shape[1]

    num_blocks_avail = candidate_rows.shape[0] // Zc
    if block_offset < 0 or num_ad_blocks <= 0:
        raise ValueError("block_offset must be >= 0 and num_ad_blocks must be > 0.")
    if block_offset + num_ad_blocks > num_blocks_avail:
        raise ValueError(
            f"Requested candidate blocks [{block_offset}, {block_offset + num_ad_blocks}) "
            f"exceed available blocks (0, {num_blocks_avail})."
        )

    new_rows = np.asarray(
        candidate_rows[block_offset * Zc : (block_offset + num_ad_blocks) * Zc],
        dtype=int,
    )
    Delta = new_rows.shape[0]
    assert Delta == num_ad_blocks * Zc

    H_aux = np.vstack((H, new_rows)).astype(int)
    m_aux = H_aux.shape[0]

    # Identity block: new column i connects ONLY to new row i (degree-1),
    # placed in the LAST Delta rows (the just-appended segment's own rows) --
    # same placement convention as soft_asced_delta_gt1.make_soft_bp_config.
    ident_block = np.zeros((m_aux, Delta), dtype=int)
    ident_block[m_aux - Delta :, :] = np.eye(Delta, dtype=int)
    H_ext = np.hstack((H_aux, ident_block)).astype(int)

    ir_vn_indices = np.arange(n, n + Delta)
    return H_ext, ir_vn_indices, Delta


def compute_generator_matrix(H_ext):
    """``G`` via ``galois.GF2(H_ext).null_space()`` -- the same approach
    already verified (this session, HARQ-shortening work) to generalize
    regardless of how H was built (QC base-graph lifting, a precomputed
    max-rank HARQ slice, or -- here -- an IR-VN extension); see
    ``apply_harq_pcm_shortening.py``'s docstring and ``generate_5G_LDPC.py``.

    Returns (G, k, n).
    """
    G = np.array(gf2(H_ext).null_space()).astype(int)
    k, n = G.shape
    return G, k, n


def build_combined_puncturing(Zc, ir_vn_indices):
    """Combined puncturing index set: the standard leading ``2*Zc``
    systematic puncturing block PLUS the trailing IR-VN columns, as ONE
    non-contiguous array passed together to ``sim.puncturing(...)`` -- this
    is exactly the case that needed the ``fix/noncontiguous-puncturing``
    merge (see task write-up)."""
    message_bit_puncturing = np.arange(2 * Zc, dtype=int)
    return np.concatenate(
        [message_bit_puncturing, np.asarray(ir_vn_indices, dtype=int)]
    )


def build_ir_vn_code(
    num_ad_blocks,
    block_offset=0,
    n_simul=DEFAULT_N_SIMUL,
    remove=DEFAULT_REMOVE,
    bg_vn=DEFAULT_BG_VN,
    bg_cn=DEFAULT_BG_CN,
    Zc=DEFAULT_ZC,
    code_path=DEFAULT_CODE_PATH,
):
    """Convenience one-shot: base PCM -> IR-VN extension -> G -> combined
    puncturing. Returns a dict with every quantity a caller (e.g. the
    waterfall sweep script) needs."""
    H, candidate_rows, number_vns, num_candidate_blocks = load_base_pcm_and_candidates(
        n_simul=n_simul, remove=remove, bg_vn=bg_vn, bg_cn=bg_cn, Zc=Zc, code_path=code_path
    )
    H_ext, ir_vn_indices, Delta = build_ir_vn_extended_pcm(
        H, candidate_rows, Zc, num_ad_blocks, block_offset=block_offset
    )
    G, k, n_ext = compute_generator_matrix(H_ext)
    combined_puncturing = build_combined_puncturing(Zc, ir_vn_indices)
    n_transmitted = n_ext - len(combined_puncturing)

    return dict(
        H=H,
        candidate_rows=candidate_rows,
        number_vns=number_vns,
        num_candidate_blocks=num_candidate_blocks,
        H_ext=H_ext,
        G=G,
        k=k,
        n_ext=n_ext,
        Delta=Delta,
        ir_vn_indices=ir_vn_indices,
        state_vn_indices=np.arange(number_vns),  # original (pre-extension) VN columns
        combined_puncturing=combined_puncturing,
        n_transmitted=n_transmitted,
        Zc=Zc,
    )


if __name__ == "__main__":
    result = build_ir_vn_code(num_ad_blocks=2, block_offset=0)
    H_ext, G = result["H_ext"], result["G"]
    k, n_ext, Delta = result["k"], result["n_ext"], result["Delta"]

    print(f"H (base) shape:        {result['H'].shape}")
    print(f"candidate_rows shape:  {result['candidate_rows'].shape} "
          f"({result['num_candidate_blocks']} Zc-blocks available)")
    print(f"H_ext shape:           {H_ext.shape}")
    print(f"G shape:               {G.shape}  (k={k}, n_ext={n_ext})")
    print(f"Delta (new rows/cols): {Delta}")
    print(f"ir_vn_indices:         {result['ir_vn_indices']}")
    print(f"combined_puncturing:   {result['combined_puncturing']}")
    print(f"n_transmitted:         {result['n_transmitted']}")

    # ---- Sanity checks ----
    check = np.array(gf2(H_ext) @ gf2(G).T)
    assert np.all(check == 0), "H_ext @ G.T != 0 -- invalid generator matrix!"
    assert k == result["number_vns"] - result["H"].shape[0], "k changed vs. the base code (should be invariant)!"
    assert n_ext == result["number_vns"] + Delta
    # Combined puncturing set must be non-contiguous for a non-trivial Delta>0
    # (this is the exact scenario that needed fix/noncontiguous-puncturing).
    gaps = np.diff(result["combined_puncturing"])
    assert np.any(gaps > 1), "expected a non-contiguous combined puncturing set"
    print("\nAll sanity checks passed.")
