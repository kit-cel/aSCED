# aSCED

This repo introduces how to use Python interface of the channel-code-lib2 C++ repo.
Furthermore, it provides the relevant matrices and scripts for reproducing the results from https://arxiv.org/pdf/2604.06889
(Affine Subcode Ensemble Decoding of Linear Block Codes)

## Setup

conda env create -f environment.yml

conda activate asced

## Install dependencies
uv lock --upgrade

uv sync

## Run exemplary scripts

uv run bp_example.py

uv run avn_bp_example.py

uv run MBBP_example.py

uv run AED_example.py


## General structure

The design is intentionally configuration-driven: instead of constructing encoders and decoders directly, all components are defined via configuration objects and executed inside a unified simulation environment.

The central idea is to separate model definition from execution

1. Simulation Environment

The Simulation_Env object represents a complete experimental setup

code structure (via parity-check matrix H)
block length parameters (k, n)
evaluation settings (e.g., puncturing, shortening, bit comparison scope)
runtime simulation control

It acts as the execution engine for all experiments.

The init of simulation environment takes as input an instance of decoder config (see 3.), which is used to internally construct such a decoder.

2. Encoder

The encoder is a lightweight component constructed from the code definition. It is used to map message vectors to codewords once the code structure is defined.

Encoders are typically initialized once and reused across simulations.

3. Decoder Configuration System

Decoders are not created directly. Instead, they are defined through decoder configuration objects.

A configuration specifies:

decoding algorithm type (e.g., BP-based variants)
algorithm hyperparameters (iterations, scheduling, normalization, etc.)
optional structural modifications (e.g., affine offsets or ensemble structure)

Examples of decoder configuration include 
- BP_config
- Ensemble_config (see 4.)

These configs fully describe how decoding is performed and is used to construct the decoder object within the simulation environment, without requiring manual construction of decoder objects.

4. Ensemble and Advanced Decoding
Ensemble decoding strategies are also expressed purely through configuration objects.

An ensemble configuration consists of two independent components:

a) Each ensemble branch can use a differently configured decoder.
This is specified via a list of decoder configuration objects.

b) Ensemble decoding may include additional pre- and post-processing steps. Since these are independent of the chosen decoder, they are handled via a separate processing configuration object.

An ensemble configuration consists of:

- a list of decoder configurations
- one processing configuration object

This separation allows flexible composition of decoding strategies.

Future convenience constructors may be added, e.g.:

MBBP_config(list of parity-check matrices)
AED_config(single decoder config, list of permutations)

# Practical Usage

All implementation details and parameter choices are demonstrated in the provided example scripts:

bp_example.py
avn_bp_example.py
MBBP_example.py
AED_example.py

These scripts serve as the primary reference the usage of decoder configuration, ensemble construction and simulation workflows.
