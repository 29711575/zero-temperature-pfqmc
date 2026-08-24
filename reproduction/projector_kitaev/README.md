# Zero-temperature projector PfQMC for the interacting Kitaev chain

The driver samples `Tr[exp(-beta_trial H_T) B_R(theta) B_L(theta)]`.  This
approaches the Gaussian pure-trial-state matrix element when the noninteracting
trial Hamiltonian has a nondegenerate ground state and `beta_trial` is large.
`H_T` uses `V_trial=0`, `delta_trial=delta`, and `mu_trial=mu`; its evolution is
split into slices no larger than `dt` and never participates in HS updates.

Each physical slice is explicitly `K/2 - V_even - V_odd - K/2`.  The Green
matrix is captured after the last `K/2` of the left theta, at the fixed midpoint.
There is no imaginary-time translation averaging or parity projection.

Build after loading oneAPI and setting `EIGEN3_INCLUDE_DIR`, then run:

`projector_driver L theta beta_trial dt V delta mu boundary hs_scheme seed burn measurements threads`

One JSON object is written to stdout. Put outputs in a fresh subdirectory of
`results/`; the driver does not read, remove, or overwrite existing results.

`PfQMC::g` is `2(I+B)^(-1)`, with unit diagonal; the skew Green matrix is
`G=g-I`.  Consequently the legacy `StructureFactorCDW*` helpers already contain
the negative onsite contact `-1/(4L)`, as well as the overall legacy minus sign.
The physical full structure factors are therefore simply minus the legacy
values.  No contact is added at the output layer.  The JSON `onsite_contact`
field is diagnostic only; `S_pi_offsite=S_pi-1/(4L)` is also diagnostic.
`R_cdw=1-S_pi_dq/S_pi` is formed from the full values separately in each bin.
