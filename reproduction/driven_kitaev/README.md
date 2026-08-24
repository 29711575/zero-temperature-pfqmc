# Imaginary-time driven Kitaev chain

`driven_driver` samples an independently HS-decoupled bra and ket contour:

`trial density matrix -> V0 projection -> ket ramp -> measurement -> mirrored bra ramp -> V0 projection`.

Ramp slice `l` uses midpoint `V_l=V0+drive_rate*(l+1/2)*dt`.  Ket slices use
`K/2-even-odd-K/2`; bra slices occur in reversed coupling order and use
`K/2-odd-even-K/2`, the adjoint operator order.  Bra and ket auxiliary fields
are separately allocated Monte Carlo variables.  The legacy helper already
contains the negative onsite contact because the sweep Green matrix has unit
diagonal, so output uses `S_Q=-S_Q_legacy` with no extra contact.  Each bin's
ratio is formed from those full sign-reweighted structure factors.

Run as:

`driven_driver L V0 Vf drive_rate theta_init beta_trial dt delta mu boundary burn measurements seed`

No partial final ramp slice is accepted.

## Production statistics and provenance

Production center values use the ratio of sign-reweighted means. Per-bin
records retain sample_count, sign_sum, signed_S_pi_numerator, and
signed_S_pi_dq_numerator; delete-one-contiguous-bin jackknife supplies the
reported errors. No per-measurement ratio average is used.

Set PFQMC_BIN_RECORDS_PATH to place the per-bin CSV beside the JSON result.
Without it, bin_records_seed_<seed>.csv is written in the current directory.
PFQMC_N_BINS optionally changes the requested bin count from its default of
15. JSON records the schema/convention versions, code identifier,
stabilization/sign-recompute settings, guard and multiprecision state,
operator count, bin path, and numerical diagnostic tolerances.

For `Vf=V0` (`n_ramp=0`), the driver deliberately calls the same static
contour builder as `projector_driver`: all `2*theta_init/dt` physical slices
use the existing projector's `K/2-even-odd-K/2` ordering, with the center after
half the slices.  This guarantees an exact finite-`dt` static regression rather
than replacing the center-right half by a separately adjointed Trotter slice.
