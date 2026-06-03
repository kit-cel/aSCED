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

# import data


# - Generate 5G LDPC code
n_ = 132
k_ = 66
H, p, s, Z, BG = generate_5G_LDPC(2, k_, n_, return_lifting_size=True)
G = gf2(H).null_space()
k, n = G.shape



cfg = channel_code_lib2.BP_config(H)


sim = channel_code_lib2.Simulation_Env(H, k, n, "all")

# cfg.H = H

sim.use_all_zero_codeword = False
sim.puncturing(p)
sim.shortening(s)
# sim.Z = Z
# sim.set_ensemble_decoding('SED', 8)
sim.init(cfg)

sim.get_error_rates(np.linspace(1, 4, 7))

FER = sim.error_rates["FER-SNR"]

print(FER)
show_results.plot_error_rates((FER, "MBBP"))
