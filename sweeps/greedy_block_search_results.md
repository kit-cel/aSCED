# Greedy RL/QC-block selection search -- results

Searches which of the 34 candidate RL/QC blocks (indices into `candidate_rows`,
see `sweeps/greedy_block_search.py`'s header for the exact construction) are
most valuable to add to the aSCED-48 and aSCED-384 ensembles, instead of the
reference script's arbitrary choice of "the next 2 (resp. 4) subsequent
blocks starting at block_offset=0".

**Method**: greedy forward selection. Each round, every not-yet-fixed
candidate block is evaluated *alone* (round 1) or *in combination with the
already-fixed blocks* (round 2+), as its own ensemble with the ensemble's
split pattern (asced48: `[2,4,6,8]`, i.e. sub-blocks of size 2,2,2,2,3,
Delta in {2,3}; asced384: `[5]`, i.e. sub-blocks of size 5,6, Delta in
{5,6}). Cheap evaluation budget: `target_errors=50, max_transmissions=2e5`
(much smaller than the production budget of 200/2e6 used elsewhere -- for
ranking only, not final numbers). FER measured at SNR=1.8dB and SNR=3.3dB.
Candidates are ranked separately at each SNR (1=best/lowest FER), and the two
per-SNR ranks are summed into a combined rank (lower is better) to avoid
being dominated by whichever SNR has larger-magnitude FER. The candidate with
the lowest combined rank is fixed as that round's pick, and the process
repeats with the remaining candidates (now evaluated in combination with all
blocks fixed so far).

The order blocks are picked in is itself the natural priority order for
sequential decoding: the earliest/most-valuable pick should be decoded first.
This order is used directly as the `fixed_sequential` path index order in the
production validation sweep (see below) - i.e. for the "optimized" ensemble,
`fixed_sequential` naturally decodes the most-valuable block's paths first.

## Known issue: blocks 10 and 19 segfault the C++ BP decoder

Blocks 10 and 19 reproducibly crash the decoder (SIGSEGV, no Python
traceback) in every round of both searches (confirmed via a retry - same
crash both times). Root cause (found via direct inspection, not guessed):
these are the only 2 of the 34 candidate blocks where every one of the 11
rows has weight exactly 1 (a "check" connected to a single variable node -
a degenerate, non-informative constraint; every other block has minimum row
weight >= 2). channel-code-lib2's BP/MSA implementation evidently assumes
check-node degree >= 2 somewhere and segfaults on degree 1. This is a C++
library bug, out of scope for this search - worked around here by excluding
blocks 10 and 19 from the candidate pool results (shown as CRASH in the
tables below, never picked as a winner). A separate background task has been
flagged to fix the underlying C++ bug.

## aSCED-48 (split pattern [2,4,6,8], 2 rounds, 67 candidate-evals / 134 SLURM tasks)

Final greedy order (1st pick first): [21, 16]

### Round 1 (fixed before this round: [])

| Candidate | FER@1.8dB | FER@3.3dB | rank@1.8 | rank@3.3 | combined rank |
|---|---|---|---|---|---|
| 21 | 0.094643 | 0.001176 | 1 | 2 | 3 **<- WINNER** |
| 3 | 0.099199 | 0.001195 | 2 | 4 | 6 |
| 31 | 0.107339 | 0.001152 | 5 | 1 | 6 |
| 16 | 0.100351 | 0.001283 | 4 | 7 | 11 |
| 15 | 0.108405 | 0.001453 | 7 | 10 | 17 |
| 17 | 0.108963 | 0.001382 | 8 | 9 | 17 |
| 12 | 0.108040 | 0.001550 | 6 | 12 | 18 |
| 13 | 0.114256 | 0.001276 | 13 | 6 | 19 |
| 23 | 0.112520 | 0.001340 | 12 | 8 | 20 |
| 26 | 0.099437 | 0.001648 | 3 | 18 | 21 |
| 24 | 0.110353 | 0.001516 | 11 | 11 | 22 |
| 28 | 0.125386 | 0.001194 | 20 | 3 | 23 |
| 4 | 0.109659 | 0.001590 | 10 | 15 | 25 |
| 2 | 0.109276 | 0.001676 | 9 | 20 | 29 |
| 8 | 0.121705 | 0.001552 | 17 | 13 | 30 |
| 1 | 0.120468 | 0.001607 | 16 | 16 | 32 |
| 14 | 0.131611 | 0.001268 | 27 | 5 | 32 |
| 7 | 0.117614 | 0.001652 | 15 | 19 | 34 |
| 30 | 0.114286 | 0.001796 | 14 | 23 | 37 |
| 0 | 0.124367 | 0.001744 | 19 | 22 | 41 |
| 33 | 0.128664 | 0.001607 | 26 | 17 | 43 |
| 5 | 0.126649 | 0.001715 | 23 | 21 | 44 |
| 6 | 0.136338 | 0.001567 | 32 | 14 | 46 |
| 22 | 0.125525 | 0.001911 | 21 | 27 | 48 |
| 29 | 0.121729 | 0.002156 | 18 | 30 | 48 |
| 20 | 0.128571 | 0.001818 | 25 | 24 | 49 |
| 32 | 0.127199 | 0.001875 | 24 | 26 | 50 |
| 25 | 0.126005 | 0.002203 | 22 | 31 | 53 |
| 27 | 0.134224 | 0.001824 | 31 | 25 | 56 |
| 9 | 0.132879 | 0.001962 | 28 | 29 | 57 |
| 11 | 0.132884 | 0.001948 | 29 | 28 | 57 |
| 18 | 0.133185 | 0.002364 | 30 | 32 | 62 |
| 10 | CRASH | CRASH | - | - | - |
| 19 | CRASH | CRASH | - | - | - |

### Round 2 (fixed before this round: [21])

| Candidate | FER@1.8dB | FER@3.3dB | rank@1.8 | rank@3.3 | combined rank |
|---|---|---|---|---|---|
| 16 | 0.085714 | 0.000700 | 3 | 1 | 4 **<- WINNER** |
| 17 | 0.089672 | 0.000759 | 6 | 2 | 8 |
| 3 | 0.084837 | 0.000812 | 2 | 9 | 11 |
| 4 | 0.089879 | 0.000782 | 8 | 4 | 12 |
| 14 | 0.091132 | 0.000759 | 9 | 3 | 12 |
| 20 | 0.091989 | 0.000791 | 10 | 5 | 15 |
| 28 | 0.086071 | 0.000842 | 4 | 14 | 18 |
| 11 | 0.089786 | 0.000869 | 7 | 19 | 26 |
| 8 | 0.092308 | 0.000851 | 11 | 17 | 28 |
| 29 | 0.084607 | 0.000932 | 1 | 27 | 28 |
| 22 | 0.095480 | 0.000845 | 15 | 15 | 30 |
| 30 | 0.088769 | 0.000919 | 5 | 25 | 30 |
| 31 | 0.100949 | 0.000810 | 22 | 8 | 30 |
| 33 | 0.093044 | 0.000864 | 12 | 18 | 30 |
| 1 | 0.104485 | 0.000804 | 26 | 6 | 32 |
| 23 | 0.095977 | 0.000846 | 16 | 16 | 32 |
| 24 | 0.103769 | 0.000805 | 25 | 7 | 32 |
| 6 | 0.100899 | 0.000841 | 21 | 12 | 33 |
| 13 | 0.094434 | 0.000885 | 14 | 21 | 35 |
| 15 | 0.105359 | 0.000832 | 28 | 10 | 38 |
| 0 | 0.093910 | 0.000932 | 13 | 28 | 41 |
| 25 | 0.096043 | 0.000909 | 17 | 24 | 41 |
| 26 | 0.098251 | 0.000907 | 20 | 22 | 42 |
| 27 | 0.107952 | 0.000835 | 31 | 11 | 42 |
| 7 | 0.106409 | 0.000841 | 30 | 13 | 43 |
| 12 | 0.096790 | 0.000928 | 18 | 26 | 44 |
| 9 | 0.104528 | 0.000881 | 27 | 20 | 47 |
| 2 | 0.097052 | 0.000965 | 19 | 30 | 49 |
| 5 | 0.106271 | 0.000908 | 29 | 23 | 52 |
| 32 | 0.101322 | 0.000948 | 23 | 29 | 52 |
| 18 | 0.102782 | 0.000998 | 24 | 31 | 55 |
| 10 | CRASH | CRASH | - | - | - |
| 19 | CRASH | CRASH | - | - | - |

## aSCED-384 (split pattern [5], 4 rounds, 130 candidate-evals / 260 SLURM tasks)

Final greedy order (1st pick first): [31, 14, 17, 3]

### Round 1 (fixed before this round: [])

| Candidate | FER@1.8dB | FER@3.3dB | rank@1.8 | rank@3.3 | combined rank |
|---|---|---|---|---|---|
| 31 | 0.068857 | 0.000573 | 1 | 1 | 2 **<- WINNER** |
| 16 | 0.070704 | 0.000601 | 3 | 2 | 5 |
| 28 | 0.078210 | 0.000612 | 5 | 3 | 8 |
| 13 | 0.070365 | 0.000763 | 2 | 9 | 11 |
| 21 | 0.082146 | 0.000614 | 11 | 4 | 15 |
| 1 | 0.076232 | 0.000801 | 4 | 12 | 16 |
| 3 | 0.079452 | 0.000784 | 6 | 10 | 16 |
| 14 | 0.082607 | 0.000621 | 13 | 5 | 18 |
| 23 | 0.085271 | 0.000658 | 15 | 6 | 21 |
| 26 | 0.085586 | 0.000728 | 16 | 7 | 23 |
| 17 | 0.079622 | 0.000870 | 7 | 18 | 25 |
| 7 | 0.081920 | 0.000843 | 10 | 16 | 26 |
| 6 | 0.081091 | 0.000904 | 8 | 19 | 27 |
| 12 | 0.084788 | 0.000831 | 14 | 13 | 27 |
| 4 | 0.089418 | 0.000736 | 21 | 8 | 29 |
| 15 | 0.087459 | 0.000800 | 19 | 11 | 30 |
| 24 | 0.086276 | 0.000833 | 17 | 15 | 32 |
| 27 | 0.081671 | 0.001101 | 9 | 23 | 32 |
| 2 | 0.082589 | 0.001167 | 12 | 24 | 36 |
| 8 | 0.089996 | 0.000832 | 22 | 14 | 36 |
| 5 | 0.086454 | 0.001176 | 18 | 25 | 43 |
| 11 | 0.092522 | 0.001045 | 24 | 20 | 44 |
| 30 | 0.090986 | 0.001048 | 23 | 21 | 44 |
| 0 | 0.093564 | 0.001088 | 25 | 22 | 47 |
| 20 | 0.106577 | 0.000863 | 30 | 17 | 47 |
| 33 | 0.089354 | 0.001223 | 20 | 27 | 47 |
| 22 | 0.101455 | 0.001215 | 26 | 26 | 52 |
| 9 | 0.103044 | 0.001236 | 27 | 28 | 55 |
| 32 | 0.104018 | 0.001659 | 28 | 31 | 59 |
| 18 | 0.108901 | 0.001489 | 32 | 29 | 61 |
| 25 | 0.105600 | 0.001982 | 29 | 32 | 61 |
| 29 | 0.108696 | 0.001515 | 31 | 30 | 61 |
| 10 | CRASH | CRASH | - | - | - |
| 19 | CRASH | CRASH | - | - | - |

### Round 2 (fixed before this round: [31])

| Candidate | FER@1.8dB | FER@3.3dB | rank@1.8 | rank@3.3 | combined rank |
|---|---|---|---|---|---|
| 14 | 0.053652 | 0.000287 | 2 | 1 | 3 **<- WINNER** |
| 16 | 0.055118 | 0.000288 | 4 | 2 | 6 |
| 3 | 0.056445 | 0.000291 | 5 | 4 | 9 |
| 23 | 0.058335 | 0.000289 | 9 | 3 | 12 |
| 17 | 0.056647 | 0.000314 | 6 | 7 | 13 |
| 13 | 0.054430 | 0.000336 | 3 | 12 | 15 |
| 33 | 0.053191 | 0.000387 | 1 | 21 | 22 |
| 12 | 0.057419 | 0.000359 | 7 | 16 | 23 |
| 28 | 0.060159 | 0.000330 | 15 | 9 | 24 |
| 21 | 0.058568 | 0.000351 | 10 | 15 | 25 |
| 6 | 0.064008 | 0.000296 | 21 | 5 | 26 |
| 2 | 0.061352 | 0.000331 | 17 | 10 | 27 |
| 15 | 0.062771 | 0.000330 | 19 | 8 | 27 |
| 1 | 0.065606 | 0.000301 | 23 | 6 | 29 |
| 8 | 0.060800 | 0.000341 | 16 | 13 | 29 |
| 11 | 0.058778 | 0.000379 | 11 | 18 | 29 |
| 26 | 0.061693 | 0.000378 | 18 | 17 | 35 |
| 30 | 0.057748 | 0.000417 | 8 | 27 | 35 |
| 0 | 0.058807 | 0.000412 | 12 | 26 | 38 |
| 5 | 0.059191 | 0.000405 | 13 | 25 | 38 |
| 20 | 0.065924 | 0.000348 | 25 | 14 | 39 |
| 7 | 0.068274 | 0.000335 | 29 | 11 | 40 |
| 32 | 0.063091 | 0.000383 | 20 | 20 | 40 |
| 22 | 0.060052 | 0.000418 | 14 | 28 | 42 |
| 24 | 0.065294 | 0.000388 | 22 | 22 | 44 |
| 4 | 0.066230 | 0.000382 | 28 | 19 | 47 |
| 27 | 0.065944 | 0.000395 | 26 | 23 | 49 |
| 29 | 0.065950 | 0.000397 | 27 | 24 | 51 |
| 25 | 0.065782 | 0.000463 | 24 | 31 | 55 |
| 9 | 0.070159 | 0.000456 | 30 | 30 | 60 |
| 18 | 0.072839 | 0.000419 | 31 | 29 | 60 |
| 10 | CRASH | CRASH | - | - | - |
| 19 | CRASH | CRASH | - | - | - |

### Round 3 (fixed before this round: [31, 14])

| Candidate | FER@1.8dB | FER@3.3dB | rank@1.8 | rank@3.3 | combined rank |
|---|---|---|---|---|---|
| 17 | 0.045326 | 0.000209 | 2 | 2 | 4 **<- WINNER** |
| 2 | 0.049830 | 0.000223 | 5 | 4 | 9 |
| 6 | 0.042618 | 0.000232 | 1 | 10 | 11 |
| 1 | 0.049684 | 0.000235 | 4 | 12 | 16 |
| 23 | 0.052250 | 0.000225 | 11 | 5 | 16 |
| 28 | 0.051263 | 0.000225 | 10 | 6 | 16 |
| 8 | 0.049606 | 0.000245 | 3 | 16 | 19 |
| 12 | 0.049908 | 0.000235 | 7 | 13 | 20 |
| 3 | 0.053333 | 0.000231 | 13 | 8 | 21 |
| 4 | 0.053396 | 0.000230 | 14 | 7 | 21 |
| 13 | 0.055240 | 0.000205 | 22 | 1 | 23 |
| 32 | 0.050805 | 0.000253 | 9 | 19 | 28 |
| 16 | 0.050666 | 0.000259 | 8 | 21 | 29 |
| 21 | 0.049852 | 0.000261 | 6 | 24 | 30 |
| 7 | 0.056671 | 0.000221 | 28 | 3 | 31 |
| 20 | 0.055300 | 0.000232 | 23 | 9 | 32 |
| 30 | 0.053780 | 0.000248 | 15 | 18 | 33 |
| 27 | 0.055596 | 0.000233 | 25 | 11 | 36 |
| 24 | 0.054260 | 0.000257 | 17 | 20 | 37 |
| 5 | 0.053201 | 0.000264 | 12 | 26 | 38 |
| 15 | 0.054254 | 0.000259 | 16 | 22 | 38 |
| 11 | 0.055612 | 0.000247 | 26 | 17 | 43 |
| 22 | 0.055065 | 0.000261 | 20 | 23 | 43 |
| 26 | 0.057554 | 0.000244 | 30 | 14 | 44 |
| 33 | 0.056821 | 0.000245 | 29 | 15 | 44 |
| 29 | 0.054698 | 0.000276 | 18 | 29 | 47 |
| 0 | 0.055227 | 0.000265 | 21 | 27 | 48 |
| 9 | 0.054881 | 0.000289 | 19 | 30 | 49 |
| 18 | 0.056409 | 0.000263 | 27 | 25 | 52 |
| 25 | 0.055556 | 0.000266 | 24 | 28 | 52 |
| 10 | CRASH | CRASH | - | - | - |
| 19 | CRASH | CRASH | - | - | - |

### Round 4 (fixed before this round: [31, 14, 17])

| Candidate | FER@1.8dB | FER@3.3dB | rank@1.8 | rank@3.3 | combined rank |
|---|---|---|---|---|---|
| 3 | 0.042754 | 0.000186 | 1 | 4 | 5 **<- WINNER** |
| 30 | 0.044277 | 0.000181 | 5 | 3 | 8 |
| 13 | 0.044449 | 0.000191 | 6 | 7 | 13 |
| 27 | 0.045187 | 0.000189 | 9 | 5 | 14 |
| 8 | 0.042881 | 0.000197 | 2 | 16 | 18 |
| 20 | 0.044758 | 0.000193 | 8 | 10 | 18 |
| 29 | 0.047054 | 0.000176 | 18 | 1 | 19 |
| 6 | 0.045330 | 0.000193 | 11 | 9 | 20 |
| 7 | 0.046462 | 0.000189 | 15 | 6 | 21 |
| 15 | 0.044631 | 0.000196 | 7 | 14 | 21 |
| 33 | 0.046122 | 0.000191 | 14 | 8 | 22 |
| 12 | 0.044089 | 0.000205 | 3 | 21 | 24 |
| 21 | 0.045887 | 0.000194 | 13 | 12 | 25 |
| 11 | 0.045228 | 0.000198 | 10 | 17 | 27 |
| 28 | 0.048540 | 0.000179 | 25 | 2 | 27 |
| 22 | 0.044177 | 0.000220 | 4 | 26 | 30 |
| 5 | 0.047848 | 0.000193 | 20 | 11 | 31 |
| 1 | 0.046991 | 0.000200 | 16 | 18 | 34 |
| 4 | 0.048154 | 0.000196 | 22 | 15 | 37 |
| 23 | 0.049252 | 0.000196 | 27 | 13 | 40 |
| 25 | 0.047043 | 0.000210 | 17 | 23 | 40 |
| 16 | 0.047873 | 0.000201 | 21 | 20 | 41 |
| 26 | 0.045724 | 0.000233 | 12 | 29 | 41 |
| 24 | 0.047652 | 0.000216 | 19 | 24 | 43 |
| 32 | 0.048347 | 0.000201 | 24 | 19 | 43 |
| 0 | 0.048167 | 0.000219 | 23 | 25 | 48 |
| 2 | 0.051889 | 0.000210 | 29 | 22 | 51 |
| 18 | 0.049092 | 0.000224 | 26 | 27 | 53 |
| 9 | 0.050182 | 0.000226 | 28 | 28 | 56 |
| 10 | CRASH | CRASH | - | - | - |
| 19 | CRASH | CRASH | - | - | - |

## Production validation sweep

Using the winning block sequences above (in greedy-discovered order) and the
same construction logic as `reproduce_fig3_RL_zc11_asced_5G_LDPC_sequential.py`,
a production-quality sweep (`target_errors=200, max_transmissions=2e6`, same
budget as the baseline sweep) was launched via
`reproduce_fig3_RL_zc11_asced_5G_LDPC_greedy.py` /
`sweeps/generate_greedy_production_manifest.py` /
`sweeps/run_greedy_production_array.sbatch`, comparing:

- **optimized** (this search's winning blocks, greedy order) vs.
- **baseline** (the reference script's "next N subsequent blocks" choice,
  already simulated at
  `RESULTS/fig_x_zc11_r4_seq_full_parallel_mpg8_n132/{asced48,asced384}_full_parallel/`
  plus the rest of Sweep A under `RESULTS/fig_x_zc11_r4_seq_<selector>_mpg8_n132/`)

across selectors {full_parallel, fixed_sequential, syndrome_sequential},
SNR {1.0..4.0dB step 0.5}, members_per_group=8,
target_num_converged in {2,6} (ignored by full_parallel), n_simul=132.
Output root: `RESULTS/fig_x_zc11_r4_seq_greedy_<selector>_mpg8_n132/<variant>_<selector>[_mpg8_target<N>]/`
(the `greedy` tag distinguishes it from the baseline's
`fig_x_zc11_r4_seq_<selector>_mpg8_n132` root).

See `sweeps/STATUS.md` for the job ID(s) and completion status of this sweep.
