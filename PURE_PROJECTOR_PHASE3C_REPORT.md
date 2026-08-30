# Pure-state projector PfQMC Phase 3C trust/initialization report

## Scope and provenance

- Branch: `pure-state-projector-phase3c-trust-init`
- Exact base: `25ad01426ad29d3d4670cfb9e5ddfb05a71841b8`
- Initial implementation commit: `ef89ca8b825b116583a0f94f27db31a56ef85899`
- Explicit-center live-stack correction: `1f7c69b16065f1228999f99e837e90b3e194b259`
- Small-system sequential compatibility audit: `fd35c0262e78dc71284d5133a24c121998980377`
- Parallel mirrored-validation ranges: `79c84c824a18a78640117fb0ebe6e10fda6487b2`
- Center-local two-branch validation stream: `98d94908bf03b6f0a9f8c32d94896ef2cdec4d7a`
  (current code HEAD)
- The base was recorded before edits. `main` and `develop` were not modified.
- Phase 4, driven propagation, condition-aware ratio, left recovery, mutating raw-sign checkpoints, and mutating MP sign checkpoints remain absent/disabled.
- No trust threshold was lowered.

## Implementation

Production RealZ2 initialization now samples only the ket HS protocol and obtains
the bra protocol by applying the contour builder's canonical strict
adjoint/reverse factor order. The combined scalar prefactor must be finite,
nonzero, and positive real. Under these conditions
`<Psi_T|U_ket^dagger U_ket|Psi_T>` is a positive norm, so the initial integer
Z2 is set to `+1` with
`initialization_policy="mirrored_theorem_z2_plus"`. The old sequential
initializer remains available only through the explicit
`--initialization-policy sequential-audit` path, and the driver rejects that
path for `L>6`. It is not the production default.

`PureEndpointRebuildResult` propagates rectangular pure subspaces from both
boundaries using independent thin QR. It does not construct a Green function at
intermediate prefixes. Only the requested cut is tested with a checked overlap
solve; it returns rank/rcond, solve residual, Green residual, log magnitude, and
a structured failure status. No explicit matrix inverse is used.

On a fast trust alarm the walker freezes the same proposal identity,
configuration hash, candidate, and Metropolis uniform. The MP evaluator uses
the canonical physical contour order, high-precision rectangular subspace
propagation, and the local Pfaffian branch formula. It evaluates 160 and 320
digits, then 640 digits only if needed, and trusts a result only when consecutive
precisions agree in sign, magnitude, reality, finiteness, nonzero status, and
endpoint solves. A missing, disabled, malformed, nonfinite, or nonconverged MP
reference is fatal. The original uniform is reused without consuming RNG.

Accepted and rejected MP decisions both rebuild the selected configuration
transactionally. The double live stack is re-anchored at a trustworthy endpoint;
the trusted MP Green may validate a recovery, but it never supplies a periodic
sign correction. Absolute Z2 comes only from the mirrored theorem and subsequent
accepted trusted-ratio transport.

An integration TDD case found that a mirrored center endpoint could be trusted
while construction of an unrelated cut-0 live stack failed. The production
mirrored walker therefore starts at the verified center cut. Sequential audit
still starts at cut 0, preserving the frozen task 4/8 prefix. The stack API now
accepts an explicit initial cut; all thresholds are unchanged.

## Frozen L=6 results

Task 4 reproduces proposal index 146 (attempt 110), pre-configuration hash
`7807527216050905631`, and uniform `0.52928038756354223`. MP160/320 gives
ratio `0.70741914065425615` and accepts, leaving Z2 `+1 -> +1`.

Task 8 reproduces proposal index 553 (attempt 415), pre-configuration hash
`5178673164861548103`, and uniform `0.56577673777691728`. MP160/320 gives
ratio `0.0026573889530685144` and rejects, leaving Z2 `-1 -> -1`.

Both runs continue for another 64 proposals, have zero MP-reference failures,
and retain the expected final configuration hashes. Compiler-dependent alarm
counts are reported rather than used as a physics criterion.

## Mirrored L=12 initialization

Tasks 16, 17, 20, and 21 each passed 128/128 independently seeded mirrored
initializations (512/512 total). Every initial integer Z2 is `+1`; the minimum
center overlap rcond is `0.99999999999999767`, the maximum center solve residual
is `1.1118181626534397e-15`, and the maximum center Green residual is
`3.1613726513817714e-15`.

Each initial state then executed 64 ordinary proposals, alternating ket and bra
factors near the center cut. Every configuration left the mirror-constrained
submanifold. Across 32768 proposals the validation accepted 22700 updates,
observed 3764 final bra/ket mirror mismatches, triggered 16 trusted MP
fallbacks, and had zero MP fallback failures. The maximum endpoint rebuild
Green residual after these proposals was `5.5057463482691647e-13`.

## Regression state

Before the explicit-center stack correction, 63/63 inherited and new checks
passed: Phase 1 (6), Phase 2 core/validation (13), Phase 3A core/validation
(12), Phase 3B I/O/production/ED (16), Phase 3C core/frozen (11), and PBC
structural checks (5). After the correction, Phase 3C core remains 9/9 and the
frozen task 4/8 ratios, decisions, proposal identities, prefix hashes, and final
hashes remain unchanged. The full inherited regression will be repeated from
the final committed tree.

The two L=6 clean-control seeds 706021 and 706022 were also run through the
explicit sequential audit path and through an independently compiled base
`25ad014` driver. Their observables, acceptances, final HS hashes, and final RNG
hashes are exactly equal. Including these compatibility checks gives 65/65
completed checks before the long 128-case and HPC rerun acceptance sets.

The archived full clean-control targets are task 0 HS/RNG hashes
`3859043808499015343`/`5111403547179266565` and task 1 hashes
`15591826770856207301`/`6529542932563746738`. The final HPC compatibility
rows use the sequential audit mode solely to compare against these old hashes;
the six former failures and both OBC controls use the mirrored theorem policy.

Dense configuration maxima retained from Phase 2 are weight error
`8.2295e-15`, Green error `7.6139e-15`, and determinant-identity error
`1.2805e-15`. Phase 3A retains ratio error `1.6007632792608383e-15`, Green
error `1.0332568796446633e-14`, Z2 mismatches 0, and trajectory mismatches 0.
Odd-L PBC rejection, even-L PBC behavior, the PBC boundary interaction-energy
test, OBC/hs1 Z2 `+1`, generic-complex mode, zero-sign JSON, retained-stream
round trip, and `/dev/full` failure all pass.

## Build commands

```bash
source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh
export LC_ALL=C LANG=C OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export EIGEN3_INCLUDE_DIR=/home/sunxr/software/eigen-3.4.0
export BOOST_INCLUDE_DIR=/home/sunxr/boost_1_70_0
export PFQMC_PFAPACK_DIR=/home/sunxr/new-pfqmc-main/inc/pfapack
reproduction/pure_projector_kitaev/build_phase3b.sh
```

## Files

Added:

- `inc/pure_projector_endpoint.h`
- `inc/pure_projector_mp.h`
- `reproduction/pure_projector_kitaev/phase3c_core_test.cpp`
- `reproduction/pure_projector_kitaev/phase3c_frozen_replay.cpp`
- `reproduction/pure_projector_kitaev/phase3c_mirrored_validation.cpp`
- `reproduction/pure_projector_kitaev/phase3c_mirrored_validation.pbs`
- `reproduction/pure_projector_kitaev/phase3c_rerun_manifest.csv`
- `reproduction/pure_projector_kitaev/phase3c_failed_controls_rerun.pbs`
- `reproduction/pure_projector_kitaev/aggregate_phase3c_rerun.py`

Modified:

- `inc/pure_projector_fast.h`
- `inc/pure_projector_stack.h`
- `reproduction/pure_projector_kitaev/pure_projector_driver.cpp`
- `reproduction/pure_projector_kitaev/build_phase3a.sh`
- `reproduction/pure_projector_kitaev/build_phase3b.sh`

Generated validation tables:

- `mirrored_initialization.csv`
- `mp_same_proposal_fallback.csv`
- `endpoint_rebuild_checks.csv`
- `failed_tasks_rerun.csv`
- `regression_summary.csv`

## Pending production acceptance

The four 128-case mirrored checks and the final 10-task HPC array are acceptance
artifacts, not a change to the target distribution. They must complete from the
final committed executable before Phase 3B production benchmark completion can
be declared. Phase 4 must not start from this report.
