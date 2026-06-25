import galois
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

gf2 = galois.GF2


def get_last_one_row_indices(matrix, last_row=True):
    """
    For each column in the matrix, find the row index of the last '1' entry.
    Return the set of row indices that are always the last '1' in their respective columns.
    Assert that there are no duplicate entries in the last one indices.
    """

    matrix = np.array(matrix)

    _, delta = matrix.shape

    last_one_indices = []

    # We are iterating from 0th to last column-> hence last one indices are already sorted such that rows of the will be in the correct order
    for col in range(delta):
        rows_with_one = np.where(matrix[:, col] == 1)[0]
        # print(matrix[rows_with_one])
        # Remove indices where the corresponding row has more than 1 entry in the matrix
        # Otherwise,
        rows_with_one = [row for row in rows_with_one if np.sum(matrix[row, :]) == 1]

        assert len(rows_with_one) > 0
        if last_row:
            last_one_indices.append(rows_with_one[-1])
        else:
            # check if randomly choosing works aswell; last row should always work
            last_one_indices.append(np.random.choice(rows_with_one))
    # Assert no duplicate entries in last_one_indices
    if len(last_one_indices) != len(set(last_one_indices)):
        raise ValueError("Duplicate entries found in last one indices!")
    # Return unique indices that are always the last one in their columns
    return np.array(last_one_indices)


def find_mapping_from_codeword_onto_avn(matrix, n_original):
    """
    given matrix, find T such that x in Csubcode is mappend onto Tx=a such that matrix @ (x|a)=0

    """

    if n_original == matrix.shape[1]:
        print("no mapping needed, n_original equals number of columns in matrix")
        return np.array([[-1]])
    elif n_original > matrix.shape[1]:
        raise ValueError(
            "n_original should not be larger than the number of columns in the matrix."
        )

    row_reduced_matrix = gf2(matrix[:, ::-1])
    row_reduced_matrix = row_reduced_matrix.row_reduce()
    row_reduced_matrix = gf2(row_reduced_matrix[:, ::-1])

    last_indices = get_last_one_row_indices(row_reduced_matrix[:, n_original:], True)
    # fixxed bug here, replaced matrix by row_reduced_matrix
    T = row_reduced_matrix[last_indices, :n_original]
    return T


def get_affine_offset_structured_PCMs(G_original, extended_H_subcode, expect_rank=1):
    """
    Given: Generator matrix G_original of a code and an extended parity check matrix of a subcode, expect_rank = rank defficiency of the subcode
    Returns: affine offset for aSCED of the extended subcode PCM
    Meanwhile: performs some test to check for consistency
    Currently: only implemented for subcodes of dim k-1
    """
    k, n = G_original.shape

    T = find_mapping_from_codeword_onto_avn(extended_H_subcode, n_original=n)

    # Careful, extended G might still be corrupted by supoerposition with codeword that is not in the subcode
    # hence, yields many columns with equal syndrome
    # However, by eliminating the influece (e.g. g0+g_ns from g1+g_ns) to g1+g_ns -(g0+g_ns) = g1+g0 without g_ns we are happy
    extended_G = np.hstack((gf2(G_original), gf2(G_original) @ gf2(T).T))

    syndrome_matrix = gf2(extended_H_subcode) @ gf2(extended_G.T)

    assert not np.all(
        syndrome_matrix == 0
    ), "syndrome_matrix should not be all zero; otherweise, wouldn't be a subcode pcm"

    if np.linalg.matrix_rank(syndrome_matrix) != expect_rank:
        # This check currently assumes the subcode is of dimension k-1.

        raise ValueError("Syndrome matrix rank is not = expect rank")
        # if this occurs, the rank defficiency is not as expected; either the subcode is not of the expected dimension or the PCM is not correct

    # perform gaussian elimination to find the basis of the column space of the syndrome matrix
    # those will be the affine offsets
    # Transpose, row reduce, transpose back to get column reduced form
    affine_offsets = (gf2(syndrome_matrix.T).row_reduce()).T

    nonzero_col_indices = np.where(~np.all(affine_offsets == 0, axis=0))[0]
    # Extract nonzero columns; no duplicate columns should be present after col_reduce
    affine_offsets = affine_offsets[:, nonzero_col_indices]

    assert (
        affine_offsets.shape[1] == expect_rank
    ), "There should be exactly expect_rank unique nonzero column type in the syndrome matrix."
    # now, given basis vectors of the column space, we can form all affine offsets by taking all linear combinations
    # using itertools.product to generate all combinations of 0 and 1 for expect_rank length

    all_combinations = list(product([0, 1], repeat=expect_rank))
    list_of_affine_offsets = []
    for combination in all_combinations:
        offset = gf2.Zeros(syndrome_matrix.shape[0])
        for i in range(expect_rank):
            if combination[i] == 1:
                offset += affine_offsets[:, i]
        # do not append all-zero offset
        if not np.all(offset == 0):
            list_of_affine_offsets.append(offset)

    return list_of_affine_offsets
