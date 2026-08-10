BCH(63,18) ASCED affine-subcode matrix pool

This ZIP contains all 52 matrices in the expanded affine-subcode SSPCM pool
from which the final size-8 and size-16 ensemble curves were selected.

Each matrices/matrix_XXX.npz contains:
  matrix                 uint8 extended SSPCM
  punctured              int64 auxiliary/punctured column indices
  affine_representative  uint8 length-63 parent-code coset representative
  appended_check         uint8 length-63 check defining the index-two subcode

The matching JSON sidecar records the appended check and weight, exact
per-block s and t values and counts, edge/cycle information, hashes, and
whether the matrix was selected for ASCED-8 or ASCED-16.
