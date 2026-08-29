# MP Z2 Trust-Gate Handoff

## Frozen branch

- Worktree: `/home/sunxr/new-pfqmc-fix-mp-z2-trust-gate`
- Branch: `fix-mp-z2-trust-gate`
- Production source tip: `14969585b0aee9b24b7bf221834041f11d1ca0df`
- Validated branch HEAD: `0c584c70ba480956c9fd326fbdfe24b4a6b298f4`
- Base: `1fb85ea5fe7aa19b558fd1586674b01f35ef41a9` (`fix-real-z2-and-raw-checker`)

`1496958` is an ancestor of `0c584c7`. No merge to develop/main was performed.

## Integration commits

Starting from base `1fb85ea`, cherry-pick the cumulative production series in this order:

1. `0c2df511dff146eb15c9cf395c1f54863422ea9e` — trusted, adaptive, canonical MP Z2 result and record-only runtime checkpoints.
2. `89809695c0dea7cb725a0645a0a48bd2f84f378b` — update transported real Z2 for every accepted proposal.
3. `14969585b0aee9b24b7bf221834041f11d1ca0df` — adaptive adjudication for accepted proposals whose ratio phase is indeterminate.

The production source tip is therefore `1496958`; it depends on the preceding two commits. Commit `0c584c70ba480956c9fd326fbdfe24b4a6b298f4` changes only PBS retry scripts and is optional for source integration.

## Executable provenance

- `projector_real_z2_driver`: `62057b5c91d695beddce0ef6f6fefe3ec0f81efe3c79ca74c116b56c1e5996e1`
- `validation_oracle`: `497890799e73d1c602296f0bbed1f4785e7488ed4ff22c8a7941b6e7fd55f8d1`
- `generic_complex_regression`: `0655a92777e6848aaaf8088a4c9ef6b9bd75927ac37e5283feea94732f8a9`
- `driven_driver`: `c2a16a427e195a8e15ae6379289072f9f42cc0f47083432d7b05a0a1c2839dc4`

The complete executable list is retained in `reproduction/projector_kitaev/mp_z2_trust_gate_fix_20260829/executable_sha256.txt`.

## Core files changed

- `inc/pfqmc.h`
- `src/pfqmc.cpp`
- `reproduction/projector_kitaev/projector_mp_z2_oracle.h`
- `reproduction/projector_kitaev/projector_json.h`
- `reproduction/projector_kitaev/real_z2_raw_checker_fix_20260827/projector_mp_z2_oracle.h`
- `reproduction/projector_kitaev/real_z2_raw_checker_fix_20260827/projector_json.h`
- `reproduction/projector_kitaev/real_z2_raw_checker_fix_20260827/projector_real_z2_driver.cpp`
- `reproduction/projector_kitaev/mp_z2_trust_gate_fix_20260829/validation_oracle_driver.cpp` (validation only)

## Required invariants

- Periodic runtime MP checkpoints are record-only: `mp_checkpoint_mutating=false` and `mp_correction_count=0`.
- No periodic MP candidate may overwrite transported `z2_sign`.
- MP results are structured and may be trusted, precision-untrusted, reality-untrusted, condition-untrusted, or unavailable; a bare sign plus nominal success is insufficient.
- Full-contour MP evaluation uses fixed canonical operator ordering and permutation convention.
- Initialization calls MP only when raw double is untrusted. Precision escalates 160 → 320 → 640 and is trusted only after two consecutive precision levels agree and numerical gates pass.
- An accepted proposal with a trustworthy near-real ratio updates Z2 with `sign(Re ratio)`.
- An accepted proposal with an indeterminate ratio may use the dedicated adaptive post-proposal absolute-Z2 adjudication path; this is separate from periodic checkpoint correction and must not change Metropolis decisions.
- Generic-complex mode is unchanged. Complex phase remains diagnostic in real-Z2 projector/driven modes.
- Condition-aware ratio and left recovery remain disabled.
- RNG, HS, fast ratios, accept/reject decisions, Green propagation and Markov trajectory remain unchanged.

## Validation gates

HPC regression job `153817.mgt` exited 0. The following groups pass: generic complex, exact projector, exact finite, OBC+hs1, zero-sign JSON, sign interface, local flips, right-boundary Green, reality symmetry, integration QMC, driven static comparison and driven smoke.

Precision controls:

- For both L18 target configurations, MP160 is untrusted and cyclic-cut dependent.
- MP320 and MP640 agree at `+1` for every tested cut.
- Eight L6 OBC+hs1 configurations agree among dense Fock, MP160, MP320 and transported Z2, all `+1`.
- Three L12 MP-only controls are `+1`.

## Expected replay results

- `L=18, theta=18, V=6, OBC, hs1, seed=2126263300`: initial and all retained physical Z2 values `+1`; average `1.0`; trajectory hash `15327981045417443903`; MP correction count `0`.
- Same cell, seed `2126263301`: initial and all retained physical Z2 values `+1`; average `1.0`; trajectory hash `13371711042269198964`; MP correction count `0`.
- `L=6, theta=12, V=4, PBC, hs0, seed=1926155102`: initial physical Z2 `-1`; average sign `0.3985 ± 0.011385`; trajectory hash `13619898121706317856`; 20/20 center and 20/20 shadow oracle checks match; MP correction count `0`.

The branch is ready for later integration review or a separately authorized archived-data replay. Neither action is part of this freeze.
