# Pure-state driven projector stage 1

This directory is isolated from the production driver and default build.  It
implements a pure-state imaginary-time protocol layer over the existing
`PureProjectorFastWalker`; it does not change the static projector, Metropolis
rule, trust thresholds, RNG policy, transported integer Z2 sign, same-proposal
multiprecision fallback, stabilization, or checkpoint implementation.

The physical contour is

`Phi_0 -> U_ket(tau_f,0) -> U_bra(tau_f,0)^dagger -> Phi_0`,

where ket and bra Hubbard--Stratonovich variables occupy separate configuration
slots and are initialized from successive, independent draws of the same
explicit RNG stream.  `U_bra^dagger` is appended as the strict reverse-adjoint
of the canonical bra action factors, including reversal within each
noncommuting checkerboard slice.  Midpoint `V_l` values come only from
`PureImaginaryTimeProtocol`.

The validated Majorana-to-Fock representation reverses the displayed factor
list.  Consequently the canonical ket builder stores physical slices from
late to early imaginary time.  This represents
`U(tau_f,0)=U_{M-1}...U_0`, with
`U_l=K_half V_even(V_l) V_odd(V_l) K_half`.  Do not change the storage loop to
early-to-late without repeating the time-dependent dense enumeration.

Build on the reference HPC environment:

```bash
source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh
export EIGEN3_INCLUDE_DIR=/home/sunxr/software/eigen-3.4.0
export BOOST_INCLUDE_DIR=/home/sunxr/boost_1_70_0
export PFQMC_PFAPACK_DIR=/home/sunxr/new-pfqmc-main/inc/pfapack
reproduction/pure_state_driven_stage1/build.sh
```

Generated validation data must be written below a fresh results directory.
No stage-1 executable writes to production result paths.

The bounded reference campaign is described by `small_campaign.csv` and may
be run with `run_small_campaign.py --jobs 6`.  `analyze_qmc.py` combines signed
observable numerators and the sign denominator with a joint blocked
jackknife; it never averages seed-level ratios.
