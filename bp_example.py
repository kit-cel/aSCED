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


## This explains the possibilities of configuring BP decoding;

# Each BP decoder requires a PCM H, all other parameters have default values
# The PCM can be overcomplete (i.e. rank(H)>n-k)
# TODO: The PCM may use auxiliary variable nodes (see cfg.use_avns()
# and for literature for instance "Iterative Decoding of Linear Block Codes: A Parity-Check Orthogonalization Approach")
# In this case, num_columns(H)>H. Those auxilary variable nodes behave similarly to punctured nodes, i.e., the respective LLRs are set to 0
#
# However, they are a decoder property and hence not set as punctured bits but handled by the BP decoder
#
cfg = channel_code_lib2.BP_config(H)


## # Typically, we recommend setting each paramter. However, for convenience we set some default parameters ##

# To avoid silently using, the cfg requires the parameter use_avns which if true allows to use a PCM with avns
# cfg.use_avns(False)


cfg.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
cfg.max_iterations = 32  # set maximum number of BP iterations; default is 32
cfg.cn_update_type = (
    "spa"  # Check node update rule (nmsa, spa, spa_phi); default is spa
)
cfg.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding


cfg.norm_factor(0.75) #set normalization constant used for normalized min sum or normalized sum product; **default is 1**;


## Affine offset used for aSCED; For an affine offset z_a, BP aims at solving H@x_hat=z_a
# in this case, early stopping is activated as soon as H@x_hat=z_a
# Input is a binary vector with length = num rows of PCMs
# Hence, it shares similarities with syndrome-based BP decoding; however, z_a here is fixed

cfg.affine_offset = np.zeros(H.shape[0], dtype=int)




# cfg.schedule(np.array)
# cfg.nodes_in_layer
# cfg.Z(int) #Lifting factor of QC LDPC Code; enables on-the-fly generation of QC automorphism or usable for scheduling
sim = channel_code_lib2.Simulation_Env(H, k, n, "all")

# cfg.H = H

sim.use_all_zero_codeword = True
sim.puncturing(p)
sim.shortening(s)
# sim.Z = Z
# sim.set_ensemble_decoding('SED', 8)
sim.init(cfg)

sim.get_error_rates(np.linspace(1, 4, 7))

FER = sim.error_rates["FER-SNR"]

print(FER)
show_results.plot_error_rates((FER, "MBBP"))
