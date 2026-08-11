import gc
from pathlib import Path

import numpy as np
import galois
import channel_code_lib2

from Codes.read_AList import read_AList
from affine_helpers import get_affine_offset_structured_PCMs


# ============================================================
# User configuration
# ============================================================

USE_ALL_ZERO = False
AUTO_SAVE = True

# Same decoder parameters as the BCH(63,30) reference script
NORM_CONST = 0.5
MAX_ITER = 20

# SNR points to simulate
SIM_REGIME = np.arange(1., 1.1, 0.5)

# ============================================================
# Paths
# ============================================================

CODE_ALIST = r"Codes/BCH63_18/BCH_63_18.alist"
MBBP_MATRIX_DIR = Path(r"Codes/BCH63_18/bch63_18_mbbp_matrices_64/matrices")
ASCED_MATRIX_DIR = Path(r"Codes/BCH63_18/bch63_18_asced_subcode_matrices_52/matrices")
RESULTS_DIR = Path(r"RESULTS/new_fig_12/BCH63_18")

# Ensemble sizes requested
MBBP8_COUNT = 8
MBBP16_COUNT = 16
ASCED8_COUNT = 8
ASCED16_COUNT = 16

# Enable / disable individual simulations
RUN_H = True
RUN_SSPCM = True

RUN_MBBP8 = True
RUN_MBBP16 = True
RUN_ASCSED8 = True
RUN_ASCSED16 = True


gf2 = galois.GF2


def load_npz_matrix(file_path: Path) -> np.ndarray:
    """Load the first array stored in an .npz matrix file."""
    if not file_path.exists():
        raise FileNotFoundError(f"Missing matrix file: {file_path}")

    archive = np.load(file_path)

    if isinstance(archive, np.lib.npyio.NpzFile):
        keys = list(archive.files)
        if not keys:
            archive.close()
            raise ValueError(f"No arrays found in {file_path}")

        matrix = archive[keys[0]]
        archive.close()
    else:
        matrix = archive

    return np.asarray(matrix).astype(int)


def configure_bp(cfg):
    """Use the same BP settings as the reference simulation."""
    cfg.use_avns = True
    cfg.early_stopping = True
    cfg.max_iterations = MAX_ITER
    cfg.cn_update_type = "msa"
    cfg.norm_factor = NORM_CONST
    cfg.scheduling_type = "flooding"
    return cfg


def build_mbbp_configs(matrix_dir: Path, count: int):
    configs = []

    for i in range(count):
        file_path = matrix_dir / f"matrix_{i:03d}.npz"
        matrix = load_npz_matrix(file_path)

        cfg = channel_code_lib2.BP_config(matrix)
        configure_bp(cfg)
        configs.append(cfg)

    return configs


def build_asced_configs(matrix_dir: Path, count: int, G):
    """
    For each aSCED subcode matrix, create:
      1. the normal BP configuration
      2. its affine-offset configuration

    Thus:
      aSCED-8  -> 16 ensemble configurations
      aSCED-16 -> 32 ensemble configurations
    """
    configs = []

    for i in range(count):
        file_path = matrix_dir / f"matrix_{i:03d}.npz"
        subcode_pcm = load_npz_matrix(file_path)

        # Original configuration
        cfg = channel_code_lib2.BP_config(subcode_pcm)
        configs.append(cfg)

        # Affine-offset configuration, following the reference script
        if not USE_ALL_ZERO:
            offsets = get_affine_offset_structured_PCMs(
                G_original=G,
                extended_H_subcode=gf2(subcode_pcm),
                expect_rank=1,
            )

            affine_cfg = channel_code_lib2.BP_config(subcode_pcm)
            affine_cfg.affine_offset = offsets[0]
            configs.append(affine_cfg)

    for cfg in configs:
        configure_bp(cfg)

    return configs


def run_ensemble(name, H, g_enc_cfg, k, n, configs, save_subdir):
    print()
    print("=" * 70)
    print(f"Starting {name}")
    print("=" * 70)

    ensemble_config = channel_code_lib2.Ensemble_config(H, configs)

    sim = channel_code_lib2.Simulation_Env(k, n, "all")
    sim.auto_save = AUTO_SAVE
    sim.save_dir = str(RESULTS_DIR / save_subdir)

    if not USE_ALL_ZERO:
        sim.init(g_enc_cfg, ensemble_config, USE_ALL_ZERO)
    else:
        sim.all_zero_init(ensemble_config)

    sim.get_error_rates(SIM_REGIME)

    fer = sim.error_rates["FER-SNR"]
    print(f"{name} FER:")
    print(fer)
    print(f"{name} finished")

    del sim
    del ensemble_config
    del configs
    gc.collect()

    return fer

def run_H(H, g_enc_cfg, k, n):
    print()
    print("=" * 70)
    print("Starting H")
    print("=" * 70)

    h_config = channel_code_lib2.BP_config(H)

    h_config.use_avns = True
    h_config.early_stopping = True
    h_config.max_iterations = MAX_ITER
    h_config.cn_update_type = "msa"
    h_config.norm_factor = NORM_CONST
    h_config.scheduling_type = "flooding"

    sim_H = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_H.auto_save = AUTO_SAVE
    sim_H.save_dir = str(RESULTS_DIR / "H")

    if not USE_ALL_ZERO:
        sim_H.init(g_enc_cfg, h_config, USE_ALL_ZERO)
    else:
        sim_H.all_zero_init(h_config)

    sim_H.get_error_rates(SIM_REGIME)

    FER_H = sim_H.error_rates["FER-SNR"]

    print("H FER:")
    print(FER_H)
    print("H finished")

    del sim_H
    del h_config
    gc.collect()

    return FER_H


def run_SSPCM(MBBP_MATRIX_DIR, g_enc_cfg, k, n):
    print()
    print("=" * 70)
    print("Starting SSPCM")
    print("=" * 70)

    sspcm_path = MBBP_MATRIX_DIR / "matrix_000.npz"

    if not sspcm_path.exists():
        raise FileNotFoundError(
            f"Missing SSPCM matrix file: {sspcm_path}"
        )

    sspcm_matrix = load_npz_matrix(sspcm_path)

    sspcm_config = channel_code_lib2.BP_config(sspcm_matrix)

    sspcm_config.use_avns = True
    sspcm_config.early_stopping = True
    sspcm_config.max_iterations = MAX_ITER
    sspcm_config.cn_update_type = "msa"
    sspcm_config.norm_factor = NORM_CONST
    sspcm_config.scheduling_type = "flooding"

    sim_sspcm = channel_code_lib2.Simulation_Env(k, n, "all")
    sim_sspcm.auto_save = AUTO_SAVE
    sim_sspcm.save_dir = str(RESULTS_DIR / "SSPCM")

    if not USE_ALL_ZERO:
        sim_sspcm.init(
            g_enc_cfg,
            sspcm_config,
            USE_ALL_ZERO,
        )
    else:
        sim_sspcm.all_zero_init(sspcm_config)

    sim_sspcm.get_error_rates(SIM_REGIME)

    FER_SSPCM = sim_sspcm.error_rates["FER-SNR"]

    print("SSPCM FER:")
    print(FER_SSPCM)
    print("SSPCM finished")

    del sim_sspcm
    del sspcm_config
    del sspcm_matrix
    gc.collect()

    return FER_SSPCM


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Load BCH parity-check matrix and derive G
    # --------------------------------------------------------
    n_alist, k_alist, H = read_AList(CODE_ALIST)

    H = np.asarray(H).astype(int)
    G = gf2(H).null_space()

    k, n = G.shape

    print(f"Code dimensions: BCH({n}, {k})")
    print(f"G shape: {G.shape}")
    print(f"H shape: {H.shape}")

    if n != n_alist:
        print(f"Warning: alist reports n={n_alist}, derived n={n}")

    if k != k_alist:
        print(f"Warning: alist reports k={k_alist}, derived k={k}")

    if not USE_ALL_ZERO:
        g_enc_cfg = channel_code_lib2.G_Encoder_config(G, k, n)
    else:
        g_enc_cfg = None

    # --------------------------------------------------------
    # H
    # --------------------------------------------------------
    if RUN_H:
        run_H(
            H,
            g_enc_cfg,
            k,
            n,
        )

    # --------------------------------------------------------
    # SSPCM
    #
    # Uses matrix_000.npz from the MBBP matrix directory.
    # --------------------------------------------------------
    if RUN_SSPCM:
        run_SSPCM(
            MBBP_MATRIX_DIR,
            g_enc_cfg,
            k,
            n,
        )

    # --------------------------------------------------------
    # MBBP-8
    # --------------------------------------------------------
    if RUN_MBBP8:
        configs = build_mbbp_configs(MBBP_MATRIX_DIR, MBBP8_COUNT)
        run_ensemble(
            "MBBP-8",
            H, g_enc_cfg, k, n,
            configs,
            "MBBP8",
        )

    # --------------------------------------------------------
    # MBBP-16
    # --------------------------------------------------------
    if RUN_MBBP16:
        configs = build_mbbp_configs(MBBP_MATRIX_DIR, MBBP16_COUNT)
        run_ensemble(
            "MBBP-16",
            H, g_enc_cfg, k, n,
            configs,
            "MBBP16",
        )

    # --------------------------------------------------------
    # aSCED-8
    # --------------------------------------------------------
    if RUN_ASCSED8:
        configs = build_asced_configs(
            ASCED_MATRIX_DIR,
            ASCED8_COUNT,
            G,
        )
        run_ensemble(
            "aSCED-8",
            H, g_enc_cfg, k, n,
            configs,
            "aSCED8",
        )

    # --------------------------------------------------------
    # aSCED-16
    # --------------------------------------------------------
    if RUN_ASCSED16:
        configs = build_asced_configs(
            ASCED_MATRIX_DIR,
            ASCED16_COUNT,
            G,
        )
        run_ensemble(
            "aSCED-16",
            H, g_enc_cfg, k, n,
            configs,
            "aSCED16",
        )

    print()
    print("=" * 70)
    print("All requested simulations finished.")
    print("=" * 70)


if __name__ == "__main__":
    main()
