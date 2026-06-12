import pandas as pd
import numpy as np
import galois

import pathlib

codes_base_dir = pathlib.Path(__file__).parent.resolve()


def generate_5G_LDPC(BG_no, K_target, N_target, return_lifting_size=False):

    # Check input sizes
    K_max = 8448 if BG_no == 1 else 3840
    if K_target > K_max:
        raise ValueError("Max. message size is %d" % K_max)

    # Get base graph table and sizes
    if BG_no == 1:
        file_name = codes_base_dir / "5G_LDPC_BaseGraphs/5G_LDPC_BG1.xls"
        M_b = 46
        N_b = 68
        K_b = 22  # ??? why only 22 out of [19,22]
    elif BG_no == 2:
        file_name = codes_base_dir / "5G_LDPC_BaseGraphs/5G_LDPC_BG2.xls"
        M_b = 42
        N_b = 52
        if K_target > 640:
            K_b = 10
        elif K_target > 560:
            K_b = 9
        elif K_target > 192:
            K_b = 8
        else:  # ??? why no 7 out of [6,10]
            K_b = 6
    else:
        raise ValueError("Invalid input for BG_no (Base graph number)")

    # Get lifting size Z
    # (min. value in set to fullfill K_b*Z >= K_target)
    all_Z = [
        [2, 4, 8, 16, 32, 64, 128, 256],
        [3, 6, 12, 24, 48, 96, 192, 384],
        [5, 10, 20, 40, 80, 160, 320],
        [7, 14, 28, 56, 112, 224],
        [9, 18, 36, 72, 144, 288],
        [11, 22, 44, 88, 176, 352],
        [13, 26, 52, 104, 208],
        [15, 30, 60, 120, 240],
    ]
    i_best = -1  # Index of lifting size set [0,7]
    j_best = -1  # Index of element in set
    Z_best = 385  # (Initial value > max. Z)
    for i in range(len(all_Z)):
        for j in range(len(all_Z[i])):
            Z = all_Z[i][j]
            if K_b * Z >= K_target and Z < Z_best:
                i_best = i
                j_best = j
                Z_best = Z
    i = i_best
    j = j_best
    Z = Z_best

    # Get base graph from table
    a = [2, 3, 5, 7, 9, 11, 13, 15][i]  # Needed to adress sheet
    BG = pd.read_excel(
        file_name, "Set %d (a = %d)" % (i + 1, a), header=None
    ).to_numpy()
    non_zero = np.where(BG != -1)
    BG[non_zero] = BG[non_zero] % Z

    # Cut out unused message part in the middle
    BG = np.hstack((BG[:, :K_b], BG[:, -M_b:]))

    # Lifting
    M_b, N_b = BG.shape
    H = np.zeros([M_b * Z, N_b * Z], dtype=np.uint8)
    for i in range(M_b):
        for j in range(N_b):
            if BG[i, j] != -1:
                H[i * Z : (i + 1) * Z, j * Z : (j + 1) * Z] = np.roll(
                    np.eye(Z, dtype=np.uint8), BG[i, j], axis=1
                )

    # Rate matching
    num_shortened_bits = K_b * Z - K_target
    initial_N_b = int(np.ceil(N_target / Z)) + 2

    # Adjust N_b in case additional bits are needed
    num_excess_bits = -((initial_N_b - 2) * Z - (N_target + num_shortened_bits))
    if num_excess_bits > 0:
        N_b = initial_N_b + int(np.ceil(num_excess_bits / Z))
    else:
        N_b = initial_N_b

    num_punctured_bits = (N_b - 2) * Z - (N_target + num_shortened_bits)

    M_b = N_b - K_b

    if num_punctured_bits > 0:
        # Message and parity bit puncturing
        P = np.concatenate(
            [
                np.arange(2 * Z),
                N_b * Z - num_punctured_bits + np.arange(num_punctured_bits),
            ]
        )
    else:
        # Message bit puncturing only
        P = np.arange(2 * Z)

    if num_shortened_bits > 0:
        # Message bit shortening
        S = K_target + np.arange(num_shortened_bits)
    else:
        # No shortening
        S = np.array([], dtype=int)

    # Get parity check matrix
    H = H[: M_b * Z, : N_b * Z]

    if return_lifting_size:
        return H, P, S, Z, BG[:3, K_b]
    else:
        return H, P, S


def get_final_matrices_and_message_bit_pucturing(H, s, p):
    gf2 = galois.GF2

    # to ensure ascending order
    # s might be empty
    # p will contain at least message bit puncturing
    s = np.sort(s)
    p = np.sort(p)

    # Initial generator matrix
    G = gf2(H).null_space()
    k, n = G.shape

    # Indices of punctured parity bits and message bits
    parity_bit_puncturing = p[p >= k]
    message_bit_pucturing = p[p < k]

    # Sanity checks
    if len(s) > 0:
        assert s[-1] < k
    if len(parity_bit_puncturing) > 0:
        assert parity_bit_puncturing[0] > s[-1]

        # Check that punctured parity bits are the last columns of H
        num_parity_punctured = len(parity_bit_puncturing)
        expected = np.arange(n - num_parity_punctured, n)
        assert np.array_equal(
            parity_bit_puncturing, expected
        ), "Parity-bit puncturing must correspond to the last columns of H"

        # Remove last columns and corresponding rows
        H = H[:, : n - num_parity_punctured]

    # Remove message-bit punctured columns
    if len(s) > 0:
        keep_cols = np.setdiff1d(np.arange(H.shape[1]), s)
        H = H[:, keep_cols]

    # Recompute generator matrix after puncturing
    G = gf2(H).null_space()
    k, n = G.shape

    return H, np.array(G), k, n, message_bit_pucturing
