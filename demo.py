import numpy as np 
from time import time

import galois
gf2 = galois.GF2

import channel_code_lib2

from Codes.generate_5G_LDPC import generate_5G_LDPC
from Codes.generate_RM import generate_RM
from Codes.overcomplete import overcomplete
from Codes.read_AList import read_AList

import show_results
#import data



# - Generate 5G LDPC code
n_ = 132
k_ = 66
H, p, s, Z, BG = generate_5G_LDPC(2, k_, n_, return_lifting_size=True)
G = gf2(H).null_space()
k,n = G.shape

# - Read code from npy
#H = np.load('Codes/RM/PCM_II_eta_max_20250507_142613.npy')
#G = gf2(H).null_space()
#k,n = G.shape

# - Generate Reed-Muller code
#n, k, H = generate_RM(5, 2)

# - Read code from AList
#n, k, H = read_AList('Codes/CCSDS/CCSDS_ldpc_n128_k64.alist.txt')

# Make parity check matrix overcomplete
#H = overcomplete(H)



sim = channel_code_lib2.Simulation_Env('bp-spa', H, k, n, 'all')
sim.use_all_zero_codeword = False
sim.puncturing(p)
sim.shortening(s)
sim.Z = Z
#sim.set_ensemble_decoding('SED', 8)
sim.init()

sim.get_error_rates(8, np.linspace(1, 4, 7))

FER = sim.error_rates['FER-SNR']
show_results.plot_error_rates((FER, 'BP'))
