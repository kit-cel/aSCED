import numpy as np
import scipy.special as sp
from itertools import combinations, product



def get_generator_matrix(m, r):
    '''
    Returns generator matrix of (m,r)-Reed-Muller code.
    '''

    # Get all possible binary vectors of length m
    binary_inputs = list(product([1, 0], repeat=m))

    # Get monomials
    monomials = []
    for deg in reversed(range(r + 1)):
        monomials.extend(combinations(range(m), deg))

    G = []  # Generator matrix
    for monomial in monomials:
        row = []
        for input in binary_inputs:
            val = 1
            for i in monomial:
                val &= input[i]
            row.append(val)
        G.append(row)

    return np.array(G, dtype=int)



def generate_RM(m, r):
    '''
    Generates parity check matrix of (m,r)-Reed-Muller code and returns it
    alongside its parameters.
    '''

    # Get code parameters
    n = 2**m   
    k = sum(sp.comb(m, i, exact=True) for i in range(r + 1))

    # Get parity check matrix (via dual code)
    H = get_generator_matrix(m, m - r - 1)

    return n, k, H
