"""Reproduce Figure 10 results for ASCED on BCH 63_36 codes."""

from time import time

import matplotlib.pyplot as plt
from pathlib import Path
import channel_code_lib2
from Codes.read_AList import read_AList


import numpy as np

import galois

gf2 = galois.GF2


import show_results

from affine_helpers import get_affine_offset_structured_PCMs

use_all_zero = False  # Currently only all-zero since bug in encode of ccsds 256,128


auto_save = True

# stopping after
bool_emulate_stopping = False
target_fraction_coverged_path = 0.5


results_dir = "RESULTS/fig_10"


sim_regime = np.linspace(2, 5, 7)

norm_const = 0.5
max_iter = 20

flag_1min = False  #

flag_ssPCM2 = False  #

flag_asced_6 = False

flag_mbbp_6 = False
flag_mbbp_30 = False

flag_sced_6 = True
flag_sced_6_w_original = True


mbbp_base_dir = Path("Codes/BCH63_36/bch_63_36_sspcm2_mbbp_32_matrices")


flag_asced_30 = False  #


plot_using_tex = False

simulate_affine = True

n, k, H = read_AList("Codes/BCH63_36/BCH_63_36.alist")

G = gf2(H).null_space()


use_all_zero = False  # Currently only all-zero since bug in encode


k, n = G.shape


print(k, n)

if not use_all_zero:
    g_enc_cfg = channel_code_lib2.G_Encoder_config(G, k, n)


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

    sim_msa_1min = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_msa_1min.auto_save = auto_save
    sim_msa_1min.save_dir = results_dir + "/1min"

    if not use_all_zero:
        print("Since only equivalent!")
        sim_msa_1min.all_zero_init(msa_1min_config)
    else:
        sim_msa_1min.all_zero_init(msa_1min_config)

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

    sim_msa_ssPCM2 = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_msa_ssPCM2.auto_save = auto_save
    sim_msa_ssPCM2.save_dir = results_dir + "/ssPCM"
    if not use_all_zero:
        sim_msa_ssPCM2.init(g_enc_cfg, msa_ssPCM2_config, use_all_zero)
    else:
        sim_msa_ssPCM2.all_zero_init(msa_ssPCM2_config)

    sim_msa_ssPCM2.get_error_rates(sim_regime)

    FER_ssPCM2 = sim_msa_ssPCM2.error_rates["FER-SNR"]
    print(FER_ssPCM2)
    print("HssPCM2 finished")


if flag_mbbp_6:
    mbbp_6_paths_configs = []

    for i in range(6):
        file_path = mbbp_base_dir / f"bch_36_63_ssPCM2_e2000_{i+1}.npy"
        if file_path.exists():
            ssPCM = np.load(file_path)
        else:
            raise FileNotFoundError(f"Missing file: {file_path}")

        mbbp_6_paths_configs.append(channel_code_lib2.BP_config(ssPCM))

    for cfg in mbbp_6_paths_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "msa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"
    mbbp_6_config = channel_code_lib2.Ensemble_config(H, mbbp_6_paths_configs)

    sim_mbbp6 = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_mbbp6.auto_save = auto_save
    sim_mbbp6.save_dir = results_dir + "/mbbp6"
    if not use_all_zero:
        sim_mbbp6.init(g_enc_cfg, mbbp_6_config, use_all_zero)
    else:
        sim_mbbp6.all_zero_init(mbbp_6_config)
    sim_mbbp6.get_error_rates(sim_regime)
    FER_mbbp6 = sim_mbbp6.error_rates["FER-SNR"]

    print(FER_mbbp6)
    print("mbbp6 finished")


if flag_sced_6_w_original:
    parent_folders = [Path("Codes/BCH63_36/aSCED-6"), Path("Codes/BCH63_36/aSCED-24")]

    sced_w_o_path_configs = [
        channel_code_lib2.BP_config(
            np.load("Codes/BCH63_36/ssPCM2_20250515_162503_228x139.npy").astype(int)
        )
    ]
    for fld in parent_folders:
        for subdir in fld.iterdir():
            if subdir.is_dir():
                # get the single file inside the subdirectory to setup batch
                file_path = next(subdir.iterdir())
                subcode_ssPCM = np.load(file_path)
                sced_w_o_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))

    sced_w_o_path_configs = sced_w_o_path_configs[:6]  # only take 6 subcode ssPCMs
    for cfg in sced_w_o_path_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "msa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"

    sced_6_w_o_config = channel_code_lib2.Ensemble_config(H, sced_w_o_path_configs)

    sim_sced6_w_o = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_sced6_w_o.auto_save = auto_save
    sim_sced6_w_o.save_dir = results_dir + "/SCED6_w_o"

    if not use_all_zero:
        sim_sced6_w_o.init(g_enc_cfg, sced_6_w_o_config, use_all_zero)
    else:
        print("SCED reqiuires random codewords! so skip")
    print("start sim")

    sim_sced6_w_o.get_error_rates(sim_regime)
    FER_SCED6_w_o = sim_sced6_w_o.error_rates["FER-SNR"]

    print(FER_SCED6_w_o)
    print("SCED6 w o finished")

if flag_sced_6:
    parent_folders = [Path("Codes/BCH63_36/aSCED-6"), Path("Codes/BCH63_36/aSCED-24")]

    sced_path_configs = []
    for fld in parent_folders:
        for subdir in fld.iterdir():
            if subdir.is_dir():
                # get the single file inside the subdirectory to setup batch
                file_path = next(subdir.iterdir())
                subcode_ssPCM = np.load(file_path)
                sced_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))

    sced_path_configs = sced_path_configs[:6]  # only take 6 subcode ssPCMs
    for cfg in sced_path_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "msa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"

    sced_6_config = channel_code_lib2.Ensemble_config(H, sced_path_configs)

    sim_sced6 = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_sced6.auto_save = auto_save
    sim_sced6.save_dir = results_dir + "/SCED6"

    if not use_all_zero:
        sim_sced6.init(g_enc_cfg, sced_6_config, use_all_zero)
    else:
        print("SCED reqiuires random codewords! so skip")
    print("start sim")

    sim_sced6.get_error_rates(sim_regime)
    FER_SCED6 = sim_sced6.error_rates["FER-SNR"]

    print(FER_SCED6)
    print("SCED6 finished")

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

    sim_asced6 = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_asced6.auto_save = auto_save
    sim_asced6.save_dir = results_dir + "/aSCED6"
    if not use_all_zero:
        sim_asced6.init(g_enc_cfg, asced_6_config, use_all_zero)
    else:
        sim_asced6.all_zero_init(asced_6_config)

    print("start sim")

    sim_asced6.get_error_rates(sim_regime)
    FER_aSCED6 = sim_asced6.error_rates["FER-SNR"]

    print(FER_aSCED6)
    print("aSCED6 finished")

if flag_mbbp_30:
    mbbp_30_paths_configs = []

    for i in range(30):
        file_path = mbbp_base_dir / f"bch_36_63_ssPCM2_e2000_{i+1}.npy"
        if file_path.exists():
            ssPCM = np.load(file_path)
        else:
            raise FileNotFoundError(f"Missing file: {file_path}")

        mbbp_30_paths_configs.append(channel_code_lib2.BP_config(ssPCM))

    for cfg in mbbp_30_paths_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "msa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"
    mbbp_30_config = channel_code_lib2.Ensemble_config(H, mbbp_30_paths_configs)

    sim_mbbp30 = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_mbbp30.auto_save = auto_save
    sim_mbbp30.save_dir = results_dir + "/mbbp30"
    if not use_all_zero:
        sim_mbbp30.init(g_enc_cfg, mbbp_30_config, use_all_zero)
    else:
        sim_mbbp30.all_zero_init(mbbp_30_config)
    sim_mbbp30.get_error_rates(sim_regime)
    FER_mbbp30 = sim_mbbp30.error_rates["FER-SNR"]

    print(FER_mbbp30)
    print("mbbp30 finished")

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

    sim_asced30 = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_asced30.auto_save = auto_save
    sim_asced30.save_dir = results_dir + "/aSCED30"

    if not use_all_zero:
        sim_asced30.init(g_enc_cfg, asced_30_config, use_all_zero)
    else:
        sim_asced30.all_zero_init(asced_30_config)
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


# 36 63
# Only equivalent! False
# Since only equivalent!
# Using all-zero codeword assumption
# {2.0: 0.6922374429223744, 2.5: 0.545045045045045, 3.0: 0.3977961432506887, 3.5: 0.25026032627559874, 4.0: 0.14100507885592087}
# H1min finished
# Using random codeword
# {2.0: 0.2951219512195122, 2.5: 0.16874797275381123, 3.0: 0.08242238740708352, 3.5: 0.032531492291843635, 4.0: 0.010745952783487599}
# HssPCM2 finished
# Using random codeword
# start sim
# {2.0: 0.13845917032248134, 2.5: 0.06539864512767066, 3.0: 0.02433279948319718, 3.5: 0.006834806786014397, 4.0: 0.001511431832975676}
# aSCED6 finished
# Using random codeword
# start sim
# {2.0: 0.08405783178827032, 2.5: 0.03256633883837447, 3.0: 0.010391232050243863, 3.5: 0.002534792015098277, 4.0: 0.0004485548599884996}
# aSCED30 finished
