import numpy as np
from time import time

import galois

gf2 = galois.GF2

import channel_code_lib2


from Codes.overcomplete import overcomplete
from Codes.read_AList import read_AList

import show_results

n_ = 7
k_ = 4

# PCM of the 7,4 hamming code
H = np.array([[1, 0, 0, 1, 1, 0, 1], [0, 1, 0, 0, 1, 1, 1], [0, 0, 1, 1, 0, 1, 1]])
G = gf2(H).null_space()
k, n = G.shape


# Here a singl avn and acn pair can remove all 4-cycles
H_avn = np.array(
    [
        [1, 0, 0, 1, 1, 0, 1, 0],
        [0, 1, 0, 0, 1, 0, 0, 1],
        [0, 0, 1, 1, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 1, 1, 1],  # acn row couples vn 6 and 7 into avn 8
    ]
)

## This explains the possibilities of configuring BP decoding;

# Each BP decoder requires a PCM H, all other parameters have default values
# The PCM can be overcomplete (i.e. rank(H)>n-k)
# The PCM may use auxiliary variable nodes (see cfg.use_avns()
# and for literature for instance "Iterative Decoding of Linear Block Codes: A Parity-Check Orthogonalization Approach")
# In this case, num_columns(H)>H. Those auxilary variable nodes behave similarly to punctured nodes, i.e., the respective LLRs are set to 0
#
# However, they are a decoder property and hence not set as punctured bits but handled by the BP decoder
#
cfg = channel_code_lib2.BP_config(H)

cfg.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
cfg.max_iterations = 32  # set maximum number of BP iterations; default is 32
cfg.cn_update_type = (
    "spa"  # Check node update rule (nmsa, spa, spa_phi); default is spa
)
cfg.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding
# cfg.norm_factor(0.75) #set normalization constant used for nmsa; default is 0.75;


cfg_avns = channel_code_lib2.BP_config(H_avn)
cfg_avns.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
cfg_avns.max_iterations = 32  # set maximum number of BP iterations; default is 32
cfg_avns.cn_update_type = (
    "spa"  # Check node update rule (nmsa, spa, spa_phi); default is spa
)
cfg_avns.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding

cfg_avns.use_avns = True


sim = channel_code_lib2.Simulation_Env(H, k, n, "all")
sim_avn = channel_code_lib2.Simulation_Env(H, k, n, "all")


sim.use_all_zero_codeword = True
sim_avn.use_all_zero_codeword = True

# sim.Z = Z
# sim.set_ensemble_decoding('SED', 8)
sim.init(cfg)
sim_avn.init(cfg)

sim.get_error_rates(np.linspace(1, 4, 7))
sim_avn.get_error_rates(np.linspace(1, 4, 7))


FER = sim.error_rates["FER-SNR"]
FER_avn = sim_avn.error_rates["FER-SNR"]

show_results.plot_error_rates((FER, "BP"), (FER_avn, "BP with AVNs"))

# the following causes would cause an error



# cfg_wrong = channel_code_lib2.BP_config(H_avn)
# ##use avns not set!

# sim_wrong = channel_code_lib2.Simulation_Env(H, k, n, "all")w
# sim_wrong.init(cfg_wrong)
# sim_wrong.get_error_rates(np.linspace(1, 4, 7))
