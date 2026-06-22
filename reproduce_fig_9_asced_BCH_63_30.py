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

sim_regime = np.linspace(2, 4.0, 5)

norm_const = 0.5
max_iter = 20

flag_1min = True  # if true simulates AED-11

flag_ssPCM2 = True  # if true simulates spa-32


flag_mbbp_8 = True
flag_mbbp_64 = False

mbbp_base_dir = Path("Codes/BCH63_30/bch_63_30_sspcm2_mbbp_64_matrices")

flag_asced_8 = True
flag_asced_64 = False  # nmsa
flag_asced_spa_64 = False  # spa

asced_base_dir = Path(
    "Codes/BCH63_30/multi_batch_Delta=1/bch_63_30_sspcm2_asced_64_matrices"
)

plot_using_tex = False

simulate_affine = True


n, k, H = read_AList("Codes/BCH63_30/BCH_63_30.alist")

G = gf2(H).null_space()

k, n = G.shape

print(k, n)


if not use_all_zero:
    g_enc_cfg = channel_code_lib2.G_Encoder_config(G, k, n)


if flag_1min:
    _, _, H_1min = read_AList("Codes/BCH63_30/BCH_63_30_1min.alist")

    ## here H1min and H are identical not only equivalent
    print("Equivalent?", np.all(gf2(H_1min) @ G.T == 0))
    msa_1min_config = channel_code_lib2.BP_config(H_1min)

    msa_1min_config.early_stopping = True  # Stop as soon as H@x_hat=0; default is true
    msa_1min_config.max_iterations = max_iter
    msa_1min_config.cn_update_type = "msa"
    msa_1min_config.norm_factor = norm_const
    msa_1min_config.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding

    sim_msa_1min = channel_code_lib2.Simulation_Env( k, n, "all")

    if not use_all_zero:
        sim_msa_1min.init(g_enc_cfg, msa_1min_config, use_all_zero)
    else:
        sim_msa_1min.all_zero_init(msa_1min_config)

    sim_msa_1min.get_error_rates(sim_regime)

    FER_1min = sim_msa_1min.error_rates["FER-SNR"]
    print(FER_1min)
    print("H1min finished")


if flag_ssPCM2:

    H_ssPCM2 = np.load(
        "Codes/BCH63_30/bch_63_30_sspcm2_mbbp_64_matrices/sspcm2_0_mbbp_e2000.npy"
    ).astype(int)
    msa_ssPCM2_config = channel_code_lib2.BP_config(H_ssPCM2)
    msa_ssPCM2_config.use_avns = True

    msa_ssPCM2_config.early_stopping = True
    msa_ssPCM2_config.max_iterations = max_iter
    msa_ssPCM2_config.cn_update_type = "msa"
    msa_ssPCM2_config.norm_factor = norm_const
    msa_ssPCM2_config.scheduling_type = "flooding"  # Scheduling method (flooding, row_layered, column_layered); default is flooding

    sim_msa_ssPCM2 = channel_code_lib2.Simulation_Env( k, n, "all")
    if not use_all_zero:
        sim_msa_ssPCM2.init(g_enc_cfg, msa_ssPCM2_config, use_all_zero)
    else:
        sim_msa_ssPCM2.all_zero_init(msa_ssPCM2_config)

    sim_msa_ssPCM2.get_error_rates(sim_regime)

    FER_ssPCM2 = sim_msa_ssPCM2.error_rates["FER-SNR"]
    print(FER_ssPCM2)
    print("HssPCM2 finished")

if flag_mbbp_8:
    mbbp_8_paths_configs = []

    for i in range(8):
        file_path = mbbp_base_dir / f"sspcm2_{i}_mbbp_e2000.npy"
        if file_path.exists():
            ssPCM = np.load(file_path)
        else:
            raise FileNotFoundError(f"Missing file: {file_path}")

        mbbp_8_paths_configs.append(channel_code_lib2.BP_config(ssPCM))

    for cfg in mbbp_8_paths_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "msa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"
    mbbp_8_config = channel_code_lib2.Ensemble_config(H, mbbp_8_paths_configs)

    sim_mbbp8 = channel_code_lib2.Simulation_Env( k, n, "all")
    if not use_all_zero:
        sim_mbbp8.init(g_enc_cfg, mbbp_8_config, use_all_zero)
    else:
        sim_mbbp8.all_zero_init(mbbp_8_config)
    sim_mbbp8.get_error_rates(sim_regime)
    FER_mbbp8 = sim_mbbp8.error_rates["FER-SNR"]

    print(FER_mbbp8)
    print("mbbp8 finished")


if flag_mbbp_64:
    mbbp_64_paths_configs = []

    for i in range(64):
        file_path = mbbp_base_dir / f"sspcm2_{i}_mbbp_e2000.npy"
        if file_path.exists():
            ssPCM = np.load(file_path)
        else:
            raise FileNotFoundError(f"Missing file: {file_path}")

        mbbp_64_paths_configs.append(channel_code_lib2.BP_config(ssPCM))

    for cfg in mbbp_64_paths_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "msa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"

    mbbp_64_config = channel_code_lib2.Ensemble_config(H, mbbp_64_paths_configs)

    sim_mbbp64 = channel_code_lib2.Simulation_Env( k, n, "all")
    if not use_all_zero:
        sim_mbbp64.init(g_enc_cfg, mbbp_64_config, use_all_zero)
    else:
        sim_mbbp64.all_zero_init(mbbp_64_config)
    sim_mbbp64.get_error_rates(sim_regime)
    FER_mbbp64 = sim_mbbp64.error_rates["FER-SNR"]

    print(FER_mbbp64)
    print("mbbp64 finished")


if flag_asced_8:

    asced8_path_configs = []
    for i in range(4):
        file_path = asced_base_dir / f"sspcm2_{i}.npy"
        if file_path.exists():
            subcode_ssPCM = np.load(file_path)
        else:
            raise FileNotFoundError(f"Missing file: {file_path}")

        subcode_ssPCM = np.load(file_path).astype(int)
        asced8_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))

        if not use_all_zero or simulate_affine:
            offsets = get_affine_offset_structured_PCMs(
                G_original=G,
                extended_H_subcode=gf2(subcode_ssPCM),
                expect_rank=1,
            )
            asced8_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))
            asced8_path_configs[-1].affine_offset = offsets[0]
    for cfg in asced8_path_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "msa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"

    asced_8_config = channel_code_lib2.Ensemble_config(H, asced8_path_configs)

    sim_asced8 = channel_code_lib2.Simulation_Env( k, n, "all")

    if not use_all_zero:
        sim_asced8.init(g_enc_cfg, asced_8_config, use_all_zero)
    else:
        sim_asced8.all_zero_init(asced_8_config)

    sim_asced8.get_error_rates(sim_regime)
    FER_asced8 = sim_asced8.error_rates["FER-SNR"]
    print(FER_asced8)
    print("asced8 finished")


if flag_asced_64:

    asced64_path_configs = []
    for i in range(32):
        file_path = asced_base_dir / f"sspcm2_{i}.npy"
        if file_path.exists():
            subcode_ssPCM = np.load(file_path)
        else:
            raise FileNotFoundError(f"Missing file: {file_path}")

        subcode_ssPCM = np.load(file_path).astype(int)
        asced64_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))

        if not use_all_zero or simulate_affine:
            offsets = get_affine_offset_structured_PCMs(
                G_original=G,
                extended_H_subcode=gf2(subcode_ssPCM),
                expect_rank=1,
            )
            asced64_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))
            asced64_path_configs[-1].affine_offset = offsets[0]
    for cfg in asced64_path_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "msa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"

    asced_64_config = channel_code_lib2.Ensemble_config(H, asced64_path_configs)

    sim_asced64 = channel_code_lib2.Simulation_Env( k, n, "all")

    if not use_all_zero:
        sim_asced64.init(g_enc_cfg, asced_64_config, use_all_zero)
    else:
        sim_asced64.all_zero_init(asced_64_config)
    sim_asced64.get_error_rates(sim_regime)
    FER_asced64 = sim_asced64.error_rates["FER-SNR"]
    print(FER_asced64)
    print("asced64 finished")

if flag_asced_spa_64:

    asced64_spa_path_configs = []
    for i in range(32):
        file_path = asced_base_dir / f"sspcm2_{i}.npy"
        if file_path.exists():
            subcode_ssPCM = np.load(file_path)
        else:
            raise FileNotFoundError(f"Missing file: {file_path}")

        subcode_ssPCM = np.load(file_path).astype(int)
        asced64_spa_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))

        if not use_all_zero or simulate_affine:
            offsets = get_affine_offset_structured_PCMs(
                G_original=G,
                extended_H_subcode=gf2(subcode_ssPCM),
                expect_rank=1,
            )
            asced64_spa_path_configs.append(channel_code_lib2.BP_config(subcode_ssPCM))
            asced64_spa_path_configs[-1].affine_offset = offsets[0]
    for cfg in asced64_spa_path_configs:
        cfg.use_avns = True
        cfg.early_stopping = True
        cfg.max_iterations = max_iter
        cfg.cn_update_type = "spa"
        cfg.norm_factor = norm_const
        cfg.scheduling_type = "flooding"

    asced_64_spa_config = channel_code_lib2.Ensemble_config(H, asced64_spa_path_configs)

    sim_asced64_spa = channel_code_lib2.Simulation_Env( k, n, "all")
    if not use_all_zero:
        sim_asced64_spa.init(g_enc_cfg, asced_64_spa_config, use_all_zero)
    else:
        sim_asced64_spa.all_zero_init(asced_64_spa_config)
    sim_asced64_spa.get_error_rates(sim_regime)
    FER_asced64_spa = sim_asced64_spa.error_rates["FER-SNR"]
    print(FER_asced64_spa)
    print("asced64spa finished")


if plot_using_tex:
    show_results.save_error_rates(
        (FER_1min, "H1min"),
        (FER_ssPCM2, "ssPCM2"),
        (FER_mbbp8, "MBBP-8"),
        # (FER_mbbp64, "MBBP-64"),
        (FER_asced8, "aSCED-8"),
        # (FER_asced64, "aSCED-NMSA-64"),
        # (FER_asced64_spa, "aSCED-NSPA-64"),
        save_name="fig_9.png",
    )
