# pylint: disable=invalid-name
"""Reproduce Figure x results for ASCED on 5G LDPC codes."""

import numpy as np
from time import time
import os

import galois

gf2 = galois.GF2

import channel_code_lib2

from Codes.generate_5G_LDPC import (
    generate_5G_LDPC,
    get_final_matrices_and_message_bit_pucturing,
)
from Codes.overcomplete import overcomplete
from Codes.read_AList import read_AList
import matplotlib.pyplot as plt

import show_results

results_dir = "RESULTS/fig_x"

bg_vn = 12
n_simul = 132  # or increase in stepzsizes of 11 e.g. 165

# we start first with code max_rank_5G_zc=11
code = np.load("Codes/TCOM_aSCED/5G_zc=11/max_rank_5G_zc=11.npz")
# this is the full BG2 PCM lifted with Zc=11
H_full = code["h"].astype(int)


plt.spy(H_full)
plt.show()
block = 0

Zc = 11

number_vn_simul = n_simul + 2 * Zc

number_vns_start = (
    bg_vn * 11
    + 2 * Zc  # is +11 not standardized? #+22 to have highest rate one from standard
)  # this here is only starting n; the code we transmit will be extendend by harq_part*Zc bits,


harq_part = (number_vn_simul - number_vns_start) // Zc

print(harq_part)

m = 44 + harq_part * Zc

number_vns = number_vns_start + harq_part * Zc

assert number_vn_simul == number_vns

H = gf2(H_full[:m, :number_vns])  # take until current m

G_sub = np.array(gf2(H[:m, :number_vns]).null_space()).astype(int)


plt.spy(H[:m, :number_vns])
plt.show()


candidate_rows = H_full[m:, : number_vns_start + harq_part * Zc]  # future possible rows

next_quasi_cyclic_block = candidate_rows[block * Zc : (block + 1) * Zc]
number_additional_cyclic_blocks = candidate_rows.shape[0] // Zc

assert number_additional_cyclic_blocks * Zc == candidate_rows.shape[0]


print(number_additional_cyclic_blocks)
