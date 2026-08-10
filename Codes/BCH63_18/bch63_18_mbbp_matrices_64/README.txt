BCH(63,18) MBBP matrix bank

This ZIP contains all 64 parent-code SSPCM-II matrices used as the MBBP
candidate bank for the final 500-observed-error FER simulations.

Each matrices/matrix_XXX.npz contains:
  matrix      uint8 parity-check matrix
  punctured   int64 auxiliary/punctured column indices

The matching JSON sidecar records the SSPCM s/t parameters, edge and cycle
counts, provenance, hashes, and whether the matrix was selected for MBBP-8 or
MBBP-16. Parent-code matrices have no appended subcode check, so that field is
null and its weight is zero.
