import numpy as np



def read_AList(path):
    '''
    Get parity check matrix & parameters from AList format
    '''

    with open(path, 'r') as file:

        lines = file.readlines()

        # get parameters
        n = int(lines[0].strip().split()[0])
        m = int(lines[0].strip().split()[1])
        k = n - m

        # init all-zero matrix of size m x n
        H = np.zeros((m, n))

        # indices of non-zero entires in rows
        row_indices = []
        for row in range(m):
            row_indices.append([int(ind) for ind in lines[4 + n + row].strip().split()])

        # add non-zero entries to matrix
        for row in range(m):
            for ind in row_indices[row]:
                H[row][ind-1] = 1

        return n, k, H.astype(int)
