# Greedy RL/QC-block selection search -- Zc=6 results (Fig. 5)

Extends the Zc=11 greedy RL/QC-block selection search (`sweeps/greedy_block_search_results.md`,
aSCED-48/384 on `C_5G(132,66)`) to the two Zc=6 codes from Fig. 5 of the paper:
`C_5G(78,60)` (37 candidate RL/QC blocks, 6 rows each) and `C_5G(180,60)` (20
candidate blocks). Searches which blocks -- and in what order -- are most
valuable to add to the ensemble, instead of the reference script's
(`reproduce_fig5_scatter_plot_zc6_5G_LDPC.py`) arbitrary choice of "the next L
subsequent blocks starting at block_offset=0".

Dimension sanity check (done before launching anything): for n_simul=78,
`H` is 30x90, `G`'s null space gives k=60, n=90, printed n=90-12=78 -- matches
`C_5G(78,60)`. For n_simul=180, `H` is 132x192, k=60, printed n=180 -- matches
`C_5G(180,60)`. Candidate block counts (`candidate_rows.shape[0]/Zc`):
37 for n=78, 20 for n=180 -- both match the task's stated pool sizes exactly.

## Scouting step (operating-point SNRs for ranking)

Unlike the Zc=11 search (which had a completed production sweep to read
1.8dB/3.3dB off), there was no existing cheap-budget full_parallel baseline
here, so a quick scouting step was run first (`sweeps/scout_zc6_fer_targets.py`):
build the "subsequent blocks" baseline at split3 (Delta=3), L=4 blocks (64
paths, matching the reference script's "64_split3" variant), then use the
reference script's own `search_target_fer()` binary-search helper (cheap
budget: `target_errors=50, max_transmissions=2e5`) to find the SNR achieving
FER=0.1 and FER=1e-3.

| n_simul | code | SNR @ FER=0.1 | measured FER | SNR @ FER=1e-3 | measured FER |
|---|---|---|---|---|---|
| 78  | C_5G(78,60)  | 3.000 dB   | 0.09947 | 5.0605 dB | 0.001015 |
| 180 | C_5G(180,60) | 1.3696 dB  | 0.10055 | 2.9531 dB | 0.000990 |

Both land close to the existing natural-order baseline data at
`/home/pj9034/aSCED/RESULTS/fig_scatter_zc6_n={78,180}/summary/*.dat`
(split=3, numblocks=4 rows: ~5.07dB for n=78, ~2.99dB for n=180 at
target_fer=1e-3) -- consistent, not a coincidence (same L=4/split3 config,
just a cheaper eval budget). These two SNR points are used as the two
ranking operating points for the greedy chain below.

## Method

Identical greedy forward-selection method to the Zc=11 search: each round,
every not-yet-fixed candidate block is evaluated *in combination with the
already-fixed blocks* as its own ensemble, using the split3 (Delta=3)
sub-splitting pattern **only** (a pure block-SELECTION search -- the
`subsplit1..5` within-block variants are NOT re-optimized here, per the
task's explicit scope; they keep their existing natural row-segment order,
applied to whichever block the chain picked first). Cheap evaluation budget:
`target_errors=50, max_transmissions=2e5`. FER measured at the two scouted
SNR points. Candidates ranked separately at each SNR (1=best/lowest FER), two
per-SNR ranks summed into a combined rank (lower is better). Chain depth: 8
rounds for each code (the max L used across the split1/split2/split3 scatter
variants).

**No crashes encountered.** Unlike the Zc=11 search (which hit 2
reproducibly-segfaulting degree-1-check-row blocks out of 34), all candidate
blocks for both Zc=6 codes evaluated cleanly in every round -- `still_missing`
is `[]` for all 16 rounds combined (8 per code). The crash-handling safety net
(missing result after the array goes terminal -> retry once -> exclude if
still missing) was exercised once anyway, incidentally: round 1 of both
searches was first submitted with a 24-hour SLURM time limit that collided
with a scheduled cluster maintenance reservation (`first_tuesday_maint_2026_10`,
2026-10-06 08:00-13:00) and had to be cancelled and resubmitted with a shorter
limit; the orchestrator's built-in "missing result -> retry once" logic
transparently absorbed this as if it were a crash, and both rounds completed
correctly on the retry.

## Results

- **n_simul=78 (C_5G(78,60))** winning blocks, greedy order (1st pick first):
  **`[11, 31, 29, 6, 27, 18, 10, 24]`** (8 rounds; see
  `sweeps/greedy_search_state_zc6_n78.json` for the complete per-round job ID
  list).
- **n_simul=180 (C_5G(180,60))** winning blocks, greedy order (1st pick first):
  **`[7, 0, 6, 3, 5, 19, 14, 2]`** (8 rounds; see
  `sweeps/greedy_search_state_zc6_n180.json` for per-round job IDs).

Only block 6 appears in both lists (n=78's 4th pick, n=180's 3rd pick) --
otherwise the two codes' winners differ, as expected: n=78 uses 37 candidate
blocks vs. n=180's 20, and they're different lifted PCMs with different
column counts, so a block's marginal contribution genuinely differs between
the two codes (same qualitative finding as the Zc=11 search, where asced48
and asced384's winners barely overlapped either).

### Round-by-round FER / rank tables

Generated directly from `sweeps/greedy_search_state_zc6_n{78,180}.json` (not
retyped/re-derived).

## n_simul=78 (C_5G(78,60))

### Round 1 (fixed before this round: [])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 11 | 0.113017 | 0.001155 | 7 | 1 | 8 **<- WINNER** |
| 18 | 0.106707 | 0.001320 | 1 | 8 | 9 |
| 4 | 0.113196 | 0.001239 | 8 | 3 | 11 |
| 31 | 0.114639 | 0.001181 | 9 | 2 | 11 |
| 29 | 0.110660 | 0.001313 | 5 | 7 | 12 |
| 33 | 0.108108 | 0.001345 | 2 | 12 | 14 |
| 34 | 0.115463 | 0.001323 | 10 | 9 | 19 |
| 23 | 0.123069 | 0.001275 | 18 | 5 | 23 |
| 26 | 0.108748 | 0.001517 | 3 | 23 | 26 |
| 30 | 0.120043 | 0.001332 | 16 | 10 | 26 |
| 36 | 0.109539 | 0.001551 | 4 | 25 | 29 |
| 7 | 0.116904 | 0.001386 | 13 | 17 | 30 |
| 16 | 0.111229 | 0.001526 | 6 | 24 | 30 |
| 3 | 0.131854 | 0.001267 | 27 | 4 | 31 |
| 17 | 0.123701 | 0.001347 | 20 | 13 | 33 |
| 5 | 0.135693 | 0.001303 | 30 | 6 | 36 |
| 6 | 0.125070 | 0.001362 | 23 | 15 | 38 |
| 12 | 0.123653 | 0.001451 | 19 | 20 | 39 |
| 22 | 0.116845 | 0.001563 | 12 | 27 | 39 |
| 13 | 0.123870 | 0.001427 | 21 | 19 | 40 |
| 19 | 0.130356 | 0.001358 | 26 | 14 | 40 |
| 27 | 0.124330 | 0.001418 | 22 | 18 | 40 |
| 2 | 0.128460 | 0.001379 | 25 | 16 | 41 |
| 8 | 0.140878 | 0.001336 | 31 | 11 | 42 |
| 10 | 0.115908 | 0.001679 | 11 | 31 | 42 |
| 24 | 0.116931 | 0.001586 | 14 | 28 | 42 |
| 20 | 0.120718 | 0.001561 | 17 | 26 | 43 |
| 15 | 0.118696 | 0.001588 | 15 | 29 | 44 |
| 1 | 0.126230 | 0.001512 | 24 | 22 | 46 |
| 9 | 0.133143 | 0.001466 | 29 | 21 | 50 |
| 0 | 0.131868 | 0.001632 | 28 | 30 | 58 |
| 25 | 0.163804 | 0.001907 | 33 | 33 | 66 |
| 21 | 0.167737 | 0.001872 | 35 | 32 | 67 |
| 35 | 0.163686 | 0.002293 | 32 | 35 | 67 |
| 28 | 0.167568 | 0.002273 | 34 | 34 | 68 |
| 14 | 0.192037 | 0.002322 | 36 | 36 | 72 |
| 32 | 0.199422 | 0.002552 | 37 | 37 | 74 |

### Round 2 (fixed before this round: [11])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 31 | 0.096967 | 0.000891 | 8 | 1 | 9 **<- WINNER** |
| 20 | 0.096293 | 0.000987 | 7 | 6 | 13 |
| 16 | 0.095897 | 0.001031 | 6 | 12 | 18 |
| 9 | 0.093056 | 0.001081 | 4 | 15 | 19 |
| 21 | 0.098271 | 0.001019 | 10 | 11 | 21 |
| 14 | 0.100987 | 0.001002 | 15 | 7 | 22 |
| 30 | 0.095131 | 0.001103 | 5 | 18 | 23 |
| 10 | 0.097934 | 0.001097 | 9 | 17 | 26 |
| 1 | 0.104459 | 0.001004 | 18 | 9 | 27 |
| 5 | 0.108390 | 0.000929 | 26 | 2 | 28 |
| 13 | 0.092671 | 0.001149 | 3 | 27 | 30 |
| 28 | 0.098271 | 0.001117 | 11 | 19 | 30 |
| 2 | 0.109290 | 0.000969 | 28 | 3 | 31 |
| 17 | 0.086233 | 0.001162 | 1 | 30 | 31 |
| 18 | 0.100000 | 0.001119 | 13 | 20 | 33 |
| 29 | 0.109520 | 0.000981 | 29 | 5 | 34 |
| 3 | 0.108302 | 0.001017 | 25 | 10 | 35 |
| 4 | 0.100087 | 0.001124 | 14 | 21 | 35 |
| 8 | 0.106206 | 0.001079 | 21 | 14 | 35 |
| 7 | 0.092357 | 0.001324 | 2 | 35 | 37 |
| 24 | 0.106724 | 0.001095 | 22 | 16 | 38 |
| 34 | 0.113747 | 0.000975 | 35 | 4 | 39 |
| 27 | 0.105120 | 0.001131 | 19 | 23 | 42 |
| 15 | 0.112262 | 0.001070 | 31 | 13 | 44 |
| 35 | 0.117737 | 0.001002 | 36 | 8 | 44 |
| 19 | 0.103035 | 0.001153 | 17 | 28 | 45 |
| 25 | 0.098787 | 0.001234 | 12 | 33 | 45 |
| 0 | 0.107647 | 0.001132 | 24 | 24 | 48 |
| 22 | 0.109141 | 0.001127 | 27 | 22 | 49 |
| 26 | 0.107124 | 0.001148 | 23 | 26 | 49 |
| 36 | 0.105263 | 0.001161 | 20 | 29 | 49 |
| 33 | 0.102505 | 0.001262 | 16 | 34 | 50 |
| 6 | 0.112284 | 0.001135 | 32 | 25 | 57 |
| 12 | 0.112511 | 0.001224 | 33 | 32 | 65 |
| 32 | 0.113226 | 0.001174 | 34 | 31 | 65 |
| 23 | 0.111284 | 0.001335 | 30 | 36 | 66 |

### Round 3 (fixed before this round: [11, 31])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 29 | 0.083219 | 0.000973 | 2 | 4 | 6 **<- WINNER** |
| 13 | 0.086222 | 0.000976 | 7 | 5 | 12 |
| 35 | 0.085233 | 0.001000 | 4 | 12 | 16 |
| 22 | 0.092100 | 0.000938 | 16 | 1 | 17 |
| 15 | 0.084215 | 0.001014 | 3 | 17 | 20 |
| 18 | 0.089347 | 0.000992 | 11 | 9 | 20 |
| 26 | 0.092737 | 0.000962 | 17 | 3 | 20 |
| 27 | 0.089526 | 0.000985 | 12 | 8 | 20 |
| 30 | 0.092095 | 0.000996 | 15 | 10 | 25 |
| 1 | 0.088426 | 0.001017 | 9 | 18 | 27 |
| 8 | 0.085880 | 0.001045 | 6 | 22 | 28 |
| 4 | 0.094421 | 0.000998 | 21 | 11 | 32 |
| 5 | 0.091049 | 0.001026 | 14 | 20 | 34 |
| 19 | 0.093064 | 0.001013 | 18 | 16 | 34 |
| 9 | 0.097412 | 0.000978 | 28 | 7 | 35 |
| 20 | 0.080129 | 0.001149 | 1 | 34 | 35 |
| 24 | 0.094600 | 0.001003 | 22 | 13 | 35 |
| 36 | 0.097772 | 0.000977 | 29 | 6 | 35 |
| 10 | 0.089345 | 0.001052 | 10 | 26 | 36 |
| 32 | 0.106815 | 0.000950 | 34 | 2 | 36 |
| 17 | 0.085558 | 0.001129 | 5 | 33 | 38 |
| 25 | 0.087023 | 0.001126 | 8 | 32 | 40 |
| 2 | 0.096007 | 0.001020 | 24 | 19 | 43 |
| 3 | 0.093108 | 0.001051 | 19 | 25 | 44 |
| 12 | 0.089739 | 0.001114 | 13 | 31 | 44 |
| 16 | 0.098289 | 0.001006 | 30 | 15 | 45 |
| 0 | 0.103226 | 0.001005 | 32 | 14 | 46 |
| 7 | 0.096985 | 0.001047 | 27 | 23 | 50 |
| 33 | 0.094677 | 0.001100 | 23 | 29 | 52 |
| 34 | 0.096245 | 0.001052 | 25 | 27 | 52 |
| 23 | 0.104076 | 0.001028 | 33 | 21 | 54 |
| 6 | 0.100884 | 0.001050 | 31 | 24 | 55 |
| 28 | 0.094174 | 0.001156 | 20 | 35 | 55 |
| 14 | 0.096714 | 0.001104 | 26 | 30 | 56 |
| 21 | 0.107670 | 0.001085 | 35 | 28 | 63 |

### Round 4 (fixed before this round: [11, 31, 29])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 6 | 0.084681 | 0.000916 | 11 | 3 | 14 **<- WINNER** |
| 3 | 0.080830 | 0.000960 | 3 | 12 | 15 |
| 36 | 0.084623 | 0.000923 | 10 | 5 | 15 |
| 14 | 0.084752 | 0.000922 | 12 | 4 | 16 |
| 23 | 0.081038 | 0.000985 | 4 | 15 | 19 |
| 18 | 0.089054 | 0.000906 | 18 | 2 | 20 |
| 19 | 0.087110 | 0.000931 | 14 | 6 | 20 |
| 33 | 0.079638 | 0.000999 | 2 | 21 | 23 |
| 13 | 0.088532 | 0.000946 | 16 | 8 | 24 |
| 34 | 0.081768 | 0.000994 | 6 | 18 | 24 |
| 10 | 0.091481 | 0.000860 | 25 | 1 | 26 |
| 12 | 0.083694 | 0.000995 | 9 | 19 | 28 |
| 7 | 0.081119 | 0.001029 | 5 | 25 | 30 |
| 15 | 0.089988 | 0.000952 | 19 | 11 | 30 |
| 8 | 0.091448 | 0.000943 | 24 | 7 | 31 |
| 20 | 0.075998 | 0.001097 | 1 | 31 | 32 |
| 25 | 0.085693 | 0.001008 | 13 | 23 | 36 |
| 28 | 0.090835 | 0.000960 | 23 | 13 | 36 |
| 16 | 0.087371 | 0.000999 | 15 | 22 | 37 |
| 21 | 0.082626 | 0.001093 | 8 | 29 | 37 |
| 1 | 0.090220 | 0.000988 | 21 | 17 | 38 |
| 0 | 0.082237 | 0.001108 | 7 | 32 | 39 |
| 2 | 0.090198 | 0.000995 | 20 | 20 | 40 |
| 9 | 0.094139 | 0.000951 | 31 | 10 | 41 |
| 32 | 0.095776 | 0.000947 | 32 | 9 | 41 |
| 24 | 0.091736 | 0.000986 | 26 | 16 | 42 |
| 22 | 0.088921 | 0.001034 | 17 | 26 | 43 |
| 27 | 0.096311 | 0.000966 | 33 | 14 | 47 |
| 4 | 0.090507 | 0.001062 | 22 | 28 | 50 |
| 5 | 0.092813 | 0.001020 | 29 | 24 | 53 |
| 17 | 0.093498 | 0.001095 | 30 | 30 | 60 |
| 26 | 0.092712 | 0.001125 | 28 | 33 | 61 |
| 30 | 0.092321 | 0.001134 | 27 | 34 | 61 |
| 35 | 0.098742 | 0.001040 | 34 | 27 | 61 |

### Round 5 (fixed before this round: [11, 31, 29, 6])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 27 | 0.082135 | 0.000879 | 9 | 2 | 11 **<- WINNER** |
| 32 | 0.078032 | 0.000930 | 2 | 9 | 11 |
| 3 | 0.076274 | 0.000943 | 1 | 13 | 14 |
| 23 | 0.085000 | 0.000877 | 13 | 1 | 14 |
| 25 | 0.083020 | 0.000905 | 10 | 4 | 14 |
| 30 | 0.080135 | 0.000940 | 5 | 12 | 17 |
| 5 | 0.081154 | 0.000933 | 7 | 11 | 18 |
| 14 | 0.086631 | 0.000901 | 18 | 3 | 21 |
| 33 | 0.085060 | 0.000912 | 15 | 6 | 21 |
| 17 | 0.078379 | 0.000971 | 3 | 20 | 23 |
| 18 | 0.078517 | 0.000983 | 4 | 22 | 26 |
| 28 | 0.088369 | 0.000911 | 21 | 5 | 26 |
| 2 | 0.083834 | 0.000950 | 12 | 15 | 27 |
| 19 | 0.081153 | 0.001000 | 6 | 24 | 30 |
| 8 | 0.088934 | 0.000931 | 23 | 10 | 33 |
| 20 | 0.094028 | 0.000925 | 29 | 8 | 37 |
| 26 | 0.082100 | 0.001032 | 8 | 30 | 38 |
| 10 | 0.096701 | 0.000916 | 32 | 7 | 39 |
| 22 | 0.083139 | 0.001009 | 11 | 28 | 39 |
| 7 | 0.087846 | 0.000977 | 20 | 21 | 41 |
| 21 | 0.090909 | 0.000958 | 26 | 16 | 42 |
| 16 | 0.089937 | 0.000969 | 25 | 18 | 43 |
| 24 | 0.085898 | 0.001000 | 17 | 26 | 43 |
| 34 | 0.089445 | 0.000969 | 24 | 19 | 43 |
| 9 | 0.090985 | 0.000961 | 27 | 17 | 44 |
| 4 | 0.096539 | 0.000949 | 31 | 14 | 45 |
| 13 | 0.085004 | 0.001137 | 14 | 33 | 47 |
| 15 | 0.085335 | 0.001073 | 16 | 31 | 47 |
| 1 | 0.088676 | 0.001008 | 22 | 27 | 49 |
| 12 | 0.087214 | 0.001122 | 19 | 32 | 51 |
| 0 | 0.092060 | 0.001000 | 28 | 25 | 53 |
| 35 | 0.102710 | 0.000997 | 33 | 23 | 56 |
| 36 | 0.094324 | 0.001010 | 30 | 29 | 59 |

### Round 6 (fixed before this round: [11, 31, 29, 6, 27])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 18 | 0.081228 | 0.000870 | 8 | 2 | 10 **<- WINNER** |
| 33 | 0.076499 | 0.000920 | 3 | 8 | 11 |
| 30 | 0.082022 | 0.000835 | 11 | 1 | 12 |
| 15 | 0.074216 | 0.000938 | 2 | 12 | 14 |
| 17 | 0.073162 | 0.000938 | 1 | 13 | 14 |
| 19 | 0.080256 | 0.000931 | 5 | 9 | 14 |
| 0 | 0.081581 | 0.000920 | 9 | 7 | 16 |
| 9 | 0.080304 | 0.000935 | 6 | 11 | 17 |
| 14 | 0.082366 | 0.000913 | 12 | 5 | 17 |
| 4 | 0.082974 | 0.000899 | 14 | 4 | 18 |
| 24 | 0.084027 | 0.000884 | 17 | 3 | 20 |
| 23 | 0.083210 | 0.000935 | 15 | 10 | 25 |
| 7 | 0.078749 | 0.000980 | 4 | 23 | 27 |
| 28 | 0.081101 | 0.000974 | 7 | 22 | 29 |
| 1 | 0.088905 | 0.000916 | 28 | 6 | 34 |
| 20 | 0.084549 | 0.000945 | 20 | 15 | 35 |
| 16 | 0.083678 | 0.000972 | 16 | 20 | 36 |
| 34 | 0.086074 | 0.000938 | 22 | 14 | 36 |
| 5 | 0.084422 | 0.000969 | 18 | 19 | 37 |
| 32 | 0.081687 | 0.001093 | 10 | 30 | 40 |
| 36 | 0.086325 | 0.000963 | 23 | 17 | 40 |
| 2 | 0.086486 | 0.000969 | 25 | 18 | 43 |
| 25 | 0.082956 | 0.001117 | 13 | 31 | 44 |
| 10 | 0.085984 | 0.001018 | 21 | 24 | 45 |
| 22 | 0.084451 | 0.001058 | 19 | 27 | 46 |
| 35 | 0.090040 | 0.000959 | 31 | 16 | 47 |
| 12 | 0.087520 | 0.000973 | 27 | 21 | 48 |
| 3 | 0.089247 | 0.001048 | 29 | 26 | 55 |
| 21 | 0.086721 | 0.001073 | 26 | 29 | 55 |
| 26 | 0.089798 | 0.001019 | 30 | 25 | 55 |
| 8 | 0.086379 | 0.001129 | 24 | 32 | 56 |
| 13 | 0.092913 | 0.001059 | 32 | 28 | 60 |

### Round 7 (fixed before this round: [11, 31, 29, 6, 27, 18])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 10 | 0.075881 | 0.000913 | 6 | 4 | 10 **<- WINNER** |
| 21 | 0.077677 | 0.000925 | 8 | 6 | 14 |
| 23 | 0.078179 | 0.000913 | 10 | 5 | 15 |
| 30 | 0.074735 | 0.000968 | 5 | 15 | 20 |
| 2 | 0.074475 | 0.000976 | 3 | 18 | 21 |
| 8 | 0.077977 | 0.000946 | 9 | 12 | 21 |
| 33 | 0.083422 | 0.000882 | 20 | 2 | 22 |
| 34 | 0.077421 | 0.000968 | 7 | 16 | 23 |
| 9 | 0.074670 | 0.000993 | 4 | 21 | 25 |
| 1 | 0.081720 | 0.000943 | 15 | 11 | 26 |
| 4 | 0.082687 | 0.000932 | 18 | 9 | 27 |
| 32 | 0.073197 | 0.001023 | 2 | 25 | 27 |
| 36 | 0.082995 | 0.000929 | 19 | 8 | 27 |
| 35 | 0.078270 | 0.000975 | 11 | 17 | 28 |
| 14 | 0.088948 | 0.000909 | 26 | 3 | 29 |
| 19 | 0.081985 | 0.000947 | 16 | 13 | 29 |
| 17 | 0.091995 | 0.000842 | 29 | 1 | 30 |
| 22 | 0.082641 | 0.000960 | 17 | 14 | 31 |
| 0 | 0.072620 | 0.001093 | 1 | 31 | 32 |
| 12 | 0.080214 | 0.000990 | 13 | 20 | 33 |
| 24 | 0.089050 | 0.000928 | 27 | 7 | 34 |
| 20 | 0.078507 | 0.001004 | 12 | 23 | 35 |
| 7 | 0.093391 | 0.000939 | 31 | 10 | 41 |
| 28 | 0.085663 | 0.000988 | 22 | 19 | 41 |
| 13 | 0.081357 | 0.001042 | 14 | 28 | 42 |
| 16 | 0.085640 | 0.000996 | 21 | 22 | 43 |
| 26 | 0.086473 | 0.001034 | 24 | 27 | 51 |
| 3 | 0.090835 | 0.001011 | 28 | 24 | 52 |
| 15 | 0.086350 | 0.001088 | 23 | 30 | 53 |
| 25 | 0.087351 | 0.001064 | 25 | 29 | 54 |
| 5 | 0.092608 | 0.001025 | 30 | 26 | 56 |

### Round 8 (fixed before this round: [11, 31, 29, 6, 27, 18, 10])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 24 | 0.070853 | 0.000885 | 2 | 5 | 7 **<- WINNER** |
| 32 | 0.079151 | 0.000877 | 12 | 3 | 15 |
| 12 | 0.076433 | 0.000909 | 8 | 9 | 17 |
| 17 | 0.078920 | 0.000891 | 11 | 7 | 18 |
| 20 | 0.071957 | 0.000944 | 3 | 15 | 18 |
| 7 | 0.075398 | 0.000943 | 6 | 14 | 20 |
| 13 | 0.070731 | 0.000953 | 1 | 21 | 22 |
| 16 | 0.084264 | 0.000861 | 23 | 2 | 25 |
| 25 | 0.080147 | 0.000909 | 15 | 10 | 25 |
| 26 | 0.074445 | 0.000950 | 5 | 20 | 25 |
| 30 | 0.082991 | 0.000898 | 19 | 8 | 27 |
| 21 | 0.077886 | 0.000945 | 10 | 18 | 28 |
| 9 | 0.086387 | 0.000839 | 28 | 1 | 29 |
| 35 | 0.085737 | 0.000885 | 25 | 4 | 29 |
| 3 | 0.084665 | 0.000887 | 24 | 6 | 30 |
| 34 | 0.075668 | 0.001006 | 7 | 24 | 31 |
| 1 | 0.083815 | 0.000910 | 21 | 11 | 32 |
| 28 | 0.081262 | 0.000944 | 16 | 16 | 32 |
| 5 | 0.079281 | 0.000949 | 14 | 19 | 33 |
| 33 | 0.073359 | 0.001092 | 4 | 29 | 33 |
| 14 | 0.077096 | 0.001076 | 9 | 27 | 36 |
| 36 | 0.079245 | 0.001010 | 13 | 25 | 38 |
| 0 | 0.086254 | 0.000938 | 27 | 12 | 39 |
| 23 | 0.085852 | 0.000942 | 26 | 13 | 39 |
| 22 | 0.082474 | 0.000982 | 18 | 23 | 41 |
| 4 | 0.083272 | 0.000969 | 20 | 22 | 42 |
| 15 | 0.081358 | 0.001063 | 17 | 26 | 43 |
| 19 | 0.089988 | 0.000944 | 29 | 17 | 46 |
| 8 | 0.084034 | 0.001101 | 22 | 30 | 52 |
| 2 | 0.091301 | 0.001077 | 30 | 28 | 58 |
## n_simul=180 (C_5G(180,60))

### Round 1 (fixed before this round: [])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 7 | 0.143633 | 0.002891 | 3 | 3 | 6 **<- WINNER** |
| 9 | 0.146956 | 0.002653 | 6 | 1 | 7 |
| 0 | 0.128058 | 0.003126 | 1 | 9 | 10 |
| 12 | 0.155987 | 0.002954 | 8 | 4 | 12 |
| 5 | 0.145418 | 0.003139 | 4 | 10 | 14 |
| 6 | 0.153210 | 0.003125 | 7 | 8 | 15 |
| 15 | 0.158694 | 0.002959 | 10 | 5 | 15 |
| 17 | 0.157282 | 0.002968 | 9 | 6 | 15 |
| 3 | 0.141230 | 0.003559 | 2 | 14 | 16 |
| 2 | 0.171984 | 0.002879 | 18 | 2 | 20 |
| 1 | 0.167474 | 0.003088 | 17 | 7 | 24 |
| 19 | 0.160625 | 0.003183 | 12 | 12 | 24 |
| 13 | 0.146224 | 0.004032 | 5 | 20 | 25 |
| 10 | 0.164236 | 0.003282 | 13 | 13 | 26 |
| 14 | 0.165746 | 0.003172 | 15 | 11 | 26 |
| 11 | 0.160528 | 0.003768 | 11 | 18 | 29 |
| 16 | 0.165257 | 0.003595 | 14 | 15 | 29 |
| 18 | 0.167183 | 0.003651 | 16 | 16 | 32 |
| 4 | 0.186888 | 0.003718 | 20 | 17 | 37 |
| 8 | 0.180984 | 0.003853 | 19 | 19 | 38 |

### Round 2 (fixed before this round: [7])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 0 | 0.105456 | 0.001281 | 2 | 1 | 3 **<- WINNER** |
| 1 | 0.105658 | 0.001431 | 3 | 4 | 7 |
| 6 | 0.108445 | 0.001505 | 4 | 6 | 10 |
| 12 | 0.115031 | 0.001466 | 7 | 5 | 12 |
| 16 | 0.111826 | 0.001544 | 6 | 7 | 13 |
| 18 | 0.115327 | 0.001584 | 8 | 9 | 17 |
| 3 | 0.100759 | 0.001696 | 1 | 17 | 18 |
| 15 | 0.124027 | 0.001413 | 15 | 3 | 18 |
| 5 | 0.110121 | 0.001668 | 5 | 15 | 20 |
| 9 | 0.118626 | 0.001599 | 10 | 10 | 20 |
| 14 | 0.130688 | 0.001359 | 18 | 2 | 20 |
| 17 | 0.120288 | 0.001576 | 12 | 8 | 20 |
| 2 | 0.118760 | 0.001631 | 11 | 12 | 23 |
| 8 | 0.115942 | 0.001658 | 9 | 14 | 23 |
| 10 | 0.120301 | 0.001629 | 13 | 11 | 24 |
| 11 | 0.124036 | 0.001642 | 16 | 13 | 29 |
| 4 | 0.121109 | 0.001805 | 14 | 19 | 33 |
| 13 | 0.140097 | 0.001686 | 19 | 16 | 35 |
| 19 | 0.130491 | 0.001792 | 17 | 18 | 35 |

### Round 3 (fixed before this round: [7, 0])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 6 | 0.094644 | 0.000878 | 1 | 4 | 5 **<- WINNER** |
| 9 | 0.099955 | 0.000826 | 6 | 1 | 7 |
| 14 | 0.094675 | 0.000889 | 2 | 6 | 8 |
| 3 | 0.099633 | 0.000933 | 5 | 7 | 12 |
| 12 | 0.101368 | 0.000870 | 9 | 3 | 12 |
| 5 | 0.100503 | 0.000938 | 8 | 8 | 16 |
| 17 | 0.098170 | 0.000992 | 4 | 13 | 17 |
| 19 | 0.113112 | 0.000856 | 16 | 2 | 18 |
| 2 | 0.094812 | 0.001046 | 3 | 16 | 19 |
| 13 | 0.102351 | 0.000949 | 10 | 9 | 19 |
| 15 | 0.109686 | 0.000879 | 15 | 5 | 20 |
| 18 | 0.100296 | 0.001020 | 7 | 15 | 22 |
| 1 | 0.103481 | 0.001001 | 11 | 14 | 25 |
| 11 | 0.105947 | 0.000992 | 13 | 12 | 25 |
| 8 | 0.118931 | 0.000965 | 18 | 10 | 28 |
| 16 | 0.113319 | 0.000965 | 17 | 11 | 28 |
| 4 | 0.103543 | 0.001084 | 12 | 17 | 29 |
| 10 | 0.106075 | 0.001085 | 14 | 18 | 32 |

### Round 4 (fixed before this round: [7, 0, 6])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 3 | 0.095780 | 0.000766 | 7 | 1 | 8 **<- WINNER** |
| 5 | 0.091687 | 0.000819 | 3 | 5 | 8 |
| 9 | 0.095983 | 0.000788 | 8 | 3 | 11 |
| 14 | 0.088096 | 0.000890 | 1 | 11 | 12 |
| 18 | 0.092034 | 0.000862 | 4 | 8 | 12 |
| 12 | 0.099635 | 0.000784 | 12 | 2 | 14 |
| 10 | 0.101124 | 0.000801 | 13 | 4 | 17 |
| 11 | 0.097670 | 0.000859 | 10 | 7 | 17 |
| 2 | 0.091203 | 0.000965 | 2 | 16 | 18 |
| 16 | 0.094092 | 0.000945 | 5 | 14 | 19 |
| 1 | 0.097528 | 0.000899 | 9 | 12 | 21 |
| 15 | 0.106101 | 0.000851 | 16 | 6 | 22 |
| 17 | 0.095483 | 0.000968 | 6 | 17 | 23 |
| 8 | 0.105096 | 0.000876 | 15 | 9 | 24 |
| 4 | 0.099516 | 0.000952 | 11 | 15 | 26 |
| 13 | 0.106175 | 0.000881 | 17 | 10 | 27 |
| 19 | 0.101991 | 0.000912 | 14 | 13 | 27 |

### Round 5 (fixed before this round: [7, 0, 6, 3])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 5 | 0.080580 | 0.000609 | 2 | 1 | 3 **<- WINNER** |
| 9 | 0.079118 | 0.000662 | 1 | 4 | 5 |
| 10 | 0.085075 | 0.000672 | 4 | 5 | 9 |
| 14 | 0.094550 | 0.000628 | 10 | 2 | 12 |
| 15 | 0.090458 | 0.000649 | 9 | 3 | 12 |
| 16 | 0.088109 | 0.000677 | 7 | 6 | 13 |
| 12 | 0.083699 | 0.000775 | 3 | 12 | 15 |
| 8 | 0.087633 | 0.000768 | 6 | 11 | 17 |
| 4 | 0.086193 | 0.000777 | 5 | 13 | 18 |
| 17 | 0.089901 | 0.000761 | 8 | 10 | 18 |
| 13 | 0.096542 | 0.000717 | 12 | 8 | 20 |
| 2 | 0.098679 | 0.000685 | 15 | 7 | 22 |
| 11 | 0.099038 | 0.000754 | 16 | 9 | 25 |
| 18 | 0.096837 | 0.000781 | 13 | 14 | 27 |
| 19 | 0.096050 | 0.000808 | 11 | 16 | 27 |
| 1 | 0.097203 | 0.000788 | 14 | 15 | 29 |

### Round 6 (fixed before this round: [7, 0, 6, 3, 5])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 19 | 0.078384 | 0.000536 | 1 | 1 | 2 **<- WINNER** |
| 16 | 0.079204 | 0.000585 | 3 | 5 | 8 |
| 9 | 0.078445 | 0.000617 | 2 | 10 | 12 |
| 10 | 0.079640 | 0.000597 | 4 | 8 | 12 |
| 1 | 0.087921 | 0.000565 | 12 | 2 | 14 |
| 18 | 0.085890 | 0.000590 | 9 | 6 | 15 |
| 2 | 0.083541 | 0.000635 | 6 | 11 | 17 |
| 4 | 0.085798 | 0.000607 | 8 | 9 | 17 |
| 13 | 0.088288 | 0.000579 | 13 | 4 | 17 |
| 17 | 0.093790 | 0.000578 | 15 | 3 | 18 |
| 15 | 0.082474 | 0.000719 | 5 | 15 | 20 |
| 11 | 0.084828 | 0.000679 | 7 | 14 | 21 |
| 14 | 0.090876 | 0.000594 | 14 | 7 | 21 |
| 8 | 0.087448 | 0.000675 | 11 | 12 | 23 |
| 12 | 0.087009 | 0.000678 | 10 | 13 | 23 |

### Round 7 (fixed before this round: [7, 0, 6, 3, 5, 19])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 14 | 0.077797 | 0.000550 | 2 | 1 | 3 **<- WINNER** |
| 12 | 0.076318 | 0.000591 | 1 | 8 | 9 |
| 13 | 0.081151 | 0.000572 | 6 | 3 | 9 |
| 8 | 0.078001 | 0.000598 | 3 | 9 | 12 |
| 16 | 0.087014 | 0.000559 | 10 | 2 | 12 |
| 15 | 0.081481 | 0.000583 | 8 | 5 | 13 |
| 18 | 0.081972 | 0.000578 | 9 | 4 | 13 |
| 4 | 0.081450 | 0.000590 | 7 | 7 | 14 |
| 1 | 0.079526 | 0.000601 | 5 | 11 | 16 |
| 10 | 0.078812 | 0.000632 | 4 | 14 | 18 |
| 2 | 0.089807 | 0.000585 | 14 | 6 | 20 |
| 17 | 0.087240 | 0.000598 | 13 | 10 | 23 |
| 9 | 0.087053 | 0.000616 | 11 | 13 | 24 |
| 11 | 0.087182 | 0.000613 | 12 | 12 | 24 |

### Round 8 (fixed before this round: [7, 0, 6, 3, 5, 19, 14])

| Candidate | FER@SNR1 | FER@SNR2 | rank@1 | rank@2 | combined rank |
|---|---|---|---|---|---|
| 2 | 0.077134 | 0.000493 | 4 | 1 | 5 **<- WINNER** |
| 18 | 0.075524 | 0.000516 | 2 | 4 | 6 |
| 9 | 0.077690 | 0.000509 | 5 | 3 | 8 |
| 4 | 0.081563 | 0.000494 | 10 | 2 | 12 |
| 8 | 0.069550 | 0.000585 | 1 | 11 | 12 |
| 17 | 0.076493 | 0.000581 | 3 | 9 | 12 |
| 11 | 0.079112 | 0.000560 | 7 | 7 | 14 |
| 15 | 0.078381 | 0.000573 | 6 | 8 | 14 |
| 13 | 0.080045 | 0.000555 | 9 | 6 | 15 |
| 1 | 0.079922 | 0.000583 | 8 | 10 | 18 |
| 12 | 0.086337 | 0.000517 | 13 | 5 | 18 |
| 10 | 0.083667 | 0.000600 | 11 | 12 | 23 |
| 16 | 0.086187 | 0.000617 | 12 | 13 | 25 |

## Production validation sweep

Full replacement experiment matrix mirroring
`reproduce_fig5_scatter_plot_zc6_5G_LDPC.py`'s `experiments` list (split1
L=1..8 + 5 subsplit variants, split2 L=1..8, split3 L=1..8, nosplit L=1..3 --
35 variants per code), using `reproduce_fig5_scatter_plot_zc6_5G_LDPC_greedy.py`,
which is identical except block selection uses the first L entries of the
greedy-discovered chain above (in order) in place of the reference's
`block_offset=0, blocks 0..L-1`. Run via the reference script's own
`search_target_fer()` (binary search to target_fer=1e-3, matching Fig. 5's
methodology exactly), `target_errors=1000` (matching the reference's
default), with one **deliberate, documented deviation**: `max_transmissions`
capped at `int(2e6)` (the original script sets no cap at all, which risks an
unbounded run at low-FER operating points -- same convention used elsewhere
in this effort, e.g. the Zc=11 production sweeps).

Manifest: `sweeps/generate_scatter_zc6_greedy_manifest.py` ->
`sweeps/scatter_zc6_greedy_manifest.csv` (70 rows, split by n_simul into
`..._n78.csv` / `..._n180.csv` for separate submission), run via
`sweeps/run_scatter_zc6_greedy_array.sbatch`. Output root:
`RESULTS/fig_scatter_zc6_n={78,180}/greedy/<variant>/` (per-decoder detailed
run output, kept separate from the baseline sweep's identically-named leaf
directories to avoid the collision class documented in STATUS.md's "Fixed a
latent collision bug" note from the Zc=11 effort). Summary file (parallel to,
not overwriting, the existing natural-order one):
`RESULTS/fig_scatter_zc6_n={78,180}/summary/n{78,180}_paths_vs_required_snr_greedy.dat`.

Submitted as jobs 541702 (n=78, 35 tasks) and 541763/541765/541778 (n=180, 35
tasks total, resubmitted across several --time adjustments -- see below).

**Cluster scheduling note**: a scheduled maintenance reservation
(`first_tuesday_maint_2026_10`, 2026-10-06 08:00-13:00 CEST, all compute
nodes) repeatedly blocked newly-submitted array tasks whose requested
walltime didn't fit before 08:00 (`ReqNodeNotAvail, Reserved for
maintenance`); the sbatch `--time` was progressively reduced (24h -> 6h -> 5h
-> 3.5h) across several resubmissions as the available pre-maintenance window
shrank overnight, and two n=180 tasks (`84_split2`, `96_split2`) that were
already RUNNING when their 3.5h limit expired got TIMEOUT'd without finishing
(these, plus whatever is still blocked by the reservation, need resubmitting
with a larger time budget after the window ends at 13:00).

[PRODUCTION_RESULTS_PLACEHOLDER]

See `sweeps/STATUS.md` for job IDs and up-to-date completion status.
