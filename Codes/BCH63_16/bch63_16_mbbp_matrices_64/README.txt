BCH(63,16) MBBP matrix bank

This ZIP contains all 64 parent-code SSPCMs in the cached pool used by the
final MBBP-8 and MBBP-16 optimization and 500-observed-error FER simulation.

Each matrices/matrix_XXX.npz contains the uint8 parity-check matrix and int64
punctured-column indices. Its JSON sidecar records exact per-block s and t
values, edge/cycle/check-weight metadata, hashes, and final selection flags.

Parent-code matrices have no appended subcode check; the sidecar therefore
records appended_check=null and appended_check_weight=0.
