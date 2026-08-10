BCH(63,16) ASCED affine-subcode matrix pool

This ZIP contains all 64 matrices in the cached affine-subcode SSPCM pool used
for the final ASCED-8 and ASCED-16 optimization and 500-observed-error FER
simulation.

Each matrices/matrix_XXX.npz contains:
  matrix                 uint8 extended SSPCM
  punctured              int64 auxiliary/punctured column indices
  affine_representative  uint8 length-63 parent-code coset representative
  appended_check         uint8 length-63 check defining the index-two subcode

The JSON sidecar records the full appended check, support and weight, exact
per-block s and t values, edge/cycle/check-weight metadata, hashes, and final
selection flags.
