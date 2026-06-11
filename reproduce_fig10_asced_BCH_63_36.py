import numpy as np
from time import time
from pathlib import Path
import galois

gf2 = galois.GF2

import channel_code_lib2


from Codes.read_AList import read_AList
import matplotlib.pyplot as plt

import show_results

from affine_helpers import get_affine_offset_structured_PCMs

use_all_zero = False  # Currently only all-zero since bug in encode of ccsds 256,128

sim_regime = np.linspace(2, 2.5, 2)

norm_const = 0.5
max_iter = 20

flag_1min = True  # if true simulates AED-11

flag_ssPCM2 = True  # if true simulates spa-32

flag_asced_6 = True

flag_asced_30 = True  # if true simulate aSCED-11

plot_using_tex = True

simulate_affine = False

n, k, H = read_AList("Codes/BCH63_36/BCH_63_36.alist")

G = gf2(H).null_space()


use_all_zero = True  # Currently only all-zero since bug in encode


k, n = G.shape


print(k, n)


if flag_1min:
    _, _, H_1min = read_AList("Codes/BCH63_36/BCH_63_36_1min.alist")

    ##CAREFUL; PCM PROVIDED BY RPTU ONLY GIVES THE PCM OF AN EQUIVALENT CODE; ONLY VALID PLOT FOR AZ SIM
    print("Only equivalent!", np.all(gf2(H_1min) @ G.T == 0))
    msa_1min_config = channel_code_lib2.BP_config(H_1min)

    msa_1min_config.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    msa_1min_config.max_iterations = max_iter
    msa_1min_config.cn_update_type = "msa"
    msa_1min_config.norm_factor = norm_const
    msa_1min_config.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding

    sim_msa_1min = channel_code_lib2.Simulation_Env(H, k, n, "all")
    sim_msa_1min.use_all_zero_codeword = use_all_zero

    sim_msa_1min.init(msa_1min_config)

    sim_msa_1min.get_error_rates(sim_regime)

    FER_1min = sim_msa_1min.error_rates["FER-SNR"]
    print(FER_1min)
    print("H1min finished")


if flag_ssPCM2:

    H_ssPCM2 = np.load("Codes/BCH63_36/ssPCM2_20250515_162503_228x139.npy").astype(int)
    msa_ssPCM2_config = channel_code_lib2.BP_config(H_ssPCM2)
    msa_ssPCM2_config.use_avns = True

    msa_ssPCM2_config.early_stopping = True
    msa_ssPCM2_config.max_iterations = max_iter
    msa_ssPCM2_config.cn_update_type = "msa"
    msa_ssPCM2_config.norm_factor = norm_const
    msa_ssPCM2_config.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding

    sim_msa_ssPCM2 = channel_code_lib2.Simulation_Env(H, k, n, "all")
    sim_msa_ssPCM2.use_all_zero_codeword = use_all_zero

    sim_msa_ssPCM2.init(msa_ssPCM2_config)

    sim_msa_ssPCM2.get_error_rates(sim_regime)

    FER_ssPCM2 = sim_msa_ssPCM2.error_rates["FER-SNR"]
    print(FER_ssPCM2)
    print("HssPCM2 finished")

if flag_asced_6:
    parent_folder = Path("Codes/BCH63_36/aSCED-6")

    asced_path_configs = []
    for subdir in parent_folder.iterdir():
        if subdir.is_dir():
            # get the single file inside the subdirectory to setup batch
            file_path = next(subdir.iterdir())
            subcode_ssPCM = np.load(file_path)
            asced_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))
            if not use_all_zero or simulate_affine:
                offsets = get_affine_offset_structured_PCMs(
                    G_original=G,
                    extended_H_subcode=gf2(subcode_ssPCM),
                    expect_rank=1,
                )
                asced_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))
                asced_path_configs[-1].affine_offset = offsets[0]
        for cfg in asced_path_configs:
            cfg.use_avns = True
            cfg.early_stopping = True
            cfg.max_iterations = max_iter
            cfg.cn_update_type = "msa"
            cfg.norm_factor = norm_const
            cfg.scheduling_type = "flooding"

    asced_6_config = channel_code_lib2.Ensemble_config(H, asced_path_configs)

    sim_asced6 = channel_code_lib2.Simulation_Env(H, k, n, "all")

    sim_asced6.use_all_zero_codeword = use_all_zero
    sim_asced6.init(asced_6_config)

    print("start sim")

    sim_asced6.get_error_rates(sim_regime)
    FER_aSCED6 = sim_asced6.error_rates["FER-SNR"]

    print(FER_aSCED6)
    print("aSCED6 finished")


if flag_asced_30:
    parent_folders = [Path("Codes/BCH63_36/aSCED-6"), Path("Codes/BCH63_36/aSCED-24")]

    asced_path_configs = []
    for fld in parent_folders:
        for subdir in fld.iterdir():
            if subdir.is_dir():
                # get the single file inside the subdirectory to setup batch
                file_path = next(subdir.iterdir())
                subcode_ssPCM = np.load(file_path)
                asced_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))
                if not use_all_zero or simulate_affine:
                    offsets = get_affine_offset_structured_PCMs(
                        G_original=G,
                        extended_H_subcode=gf2(subcode_ssPCM),
                        expect_rank=1,
                    )
                    asced_path_configs.append(
                        channel_code_lib2.BP_config(subcode_ssPCM)
                    )
                    asced_path_configs[-1].affine_offset = offsets[0]

    for cfg in asced_path_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "msa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"

    asced_30_config = channel_code_lib2.Ensemble_config(H, asced_path_configs)

    sim_asced30 = channel_code_lib2.Simulation_Env(H, k, n, "all")

    sim_asced30.use_all_zero_codeword = use_all_zero
    sim_asced30.init(asced_30_config)

    print("start sim")

    sim_asced30.get_error_rates(sim_regime)
    FER_aSCED30 = sim_asced30.error_rates["FER-SNR"]

    print(FER_aSCED30)
    print("aSCED30 finished")


if plot_using_tex:
    show_results.save_error_rates(
        (FER_1min, "H1min"),
        (FER_ssPCM2, "ssPCM2"),
        (FER_aSCED6, "aSCED-6"),
        (FER_aSCED30, "aSCED-30"),
        save_name="fig_10.png",
    )
