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


n = 7
k = 4
H = np.array(
    [
        [1, 0, 0, 1, 1, 0, 1],
        [0, 1, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 1],
    ]
)

G = gf2(H).null_space()

G[0]=G[0]+G[1]
G[1]=G[1]+G[2]

cfg = channel_code_lib2.BP_config(H)

enc_cfg = channel_code_lib2.G_Encoder_config(G, k, n)

cfg.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
cfg.max_iterations = 20  # set maximum number of BP iterations; default is 32
cfg.cn_update_type = "msa"  # Check node update rule (msa, spa, spa_phi); default is spa
cfg.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding


cfg.norm_factor = 0.75  # set normalization constant used for normalized min sum or normalized sum product; **default is 1**;



cfg.affine_offset = np.zeros(H.shape[0], dtype=int)


sim = channel_code_lib2.Simulation_Env(H, k, n, "all")

sim.use_all_zero_codeword = False

# sim.all_zero_init(cfg)

sim.init(enc_cfg, cfg)

sim.get_error_rates(np.linspace(1, 4, 7))

FER = sim.error_rates["FER-SNR"]

print(FER)
show_results.plot_error_rates((FER, "MBBP"))
