import numpy as np
from itertools import combinations



def overcomplete(H):
    '''
    Makes given parity check matrix overcomplete by adding all possible parity
    check equations of minimum weight.
    '''

    m, n = H.shape

    min_weight = n + 1                  # Minimum weight (Init out of bound)
    min_weight_parity_checks = []       # Parity check equations of minimum weight

    # Generate all linear combinations of rows in parity check matrix
    for j in range(1, m+1):
        for comb in combinations(range(m), j):
            parity_check = np.sum(H[list(comb), :], axis=0) % 2

            # Check for minimum weight
            weight = np.count_nonzero(parity_check)

            if weight == 0:
                # Don't use all-zero row
                continue

            if weight < min_weight:
                # Update minimum weight and add parity check equation
                min_weight = weight
                min_weight_parity_checks = [parity_check]

            if weight == min_weight:
                # Add parity check equation if not a duplicate
                if not any(np.array_equal(parity_check, existing) for existing in min_weight_parity_checks):
                    min_weight_parity_checks.append(parity_check)
                
    return np.stack(min_weight_parity_checks, axis=0)   # Overcomplete parity check matrix
