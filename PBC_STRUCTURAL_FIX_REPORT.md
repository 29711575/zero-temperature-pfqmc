# PBC structural guards and interaction-energy fix

## Scope and provenance

- HPC working directory: `/home/sunxr/new-pfqmc-fix-pbc-structural-guards`
- Branch: `fix-pbc-structural-guards`
- Frozen base HEAD: `4d79bb7ac2b734bc10667d917b914ce316a39a39`
- Required production commits are ancestors of the base: `0c2df511`,
  `89809695`, and `14969585`.
- No changes were made to `main`, `develop`, or the running proposal-safety
  worktree/job `153850[]`.
- No condition-aware ratio recovery or left recovery was enabled or run.

## Issue #6: odd-L PBC

`SpinlessTvChainUtils`, the configuration class shared by the finite-T,
projector, and driven construction paths, now rejects `boundary == 0` with an
odd `L`.  The exception is raised before any contour or bond operators are
constructed and reports:

`odd-L PBC is unsupported by the two-layer bond decomposition`

The guard is deliberately not applied to OBC.  Even-L PBC follows the prior
construction path unchanged.

Direct checks:

- L=3 PBC shared construction: rejected (PASS).
- L=3 OBC finite-T `Chain_tV`: constructed with bond-layer sizes 1 and 1
  (PASS).
- L=3 PBC projector production entry: exit 2 with the shared guard message
  (PASS).
- L=3 PBC driven production entry: exit 2 with the shared guard message
  (PASS).
- L=6 PBC interaction operators: maximum `||B B_inv - I|| = 0` (PASS).
- L=6 PBC right-boundary Green regression: 180/180 finite boundaries,
  maximum fast/full relative error `2.1837046465355791e-10` (PASS, identical
  to the established regression scale).

## Issue #5: PBC interaction energy

The interaction portion of `SpinlessTvChainUtils::energyFromGreensFunc()` now
uses `L` bonds for PBC and `L-1` bonds for OBC.  Therefore PBC includes the
missing `(L-1,0)` Wick term, while the OBC bond set is unchanged.  The kinetic,
pairing, chemical-potential, structure-factor, ratio, and sign code paths were
not modified.

An independent L=4 deterministic direct-Wick sum was compared with the
interaction-only contribution obtained by subtracting the V=0 helper result
from the V=2 helper result on the same Green matrix:

- PBC: four bonds, relative error `2.341683895620198e-15` (PASS, tolerance
  `1e-13`).
- OBC: three bonds, relative error `8.396944963237188e-16` (PASS, tolerance
  `1e-13`).

No ED energy comparison was added.  The exact direct-Wick test isolates this
helper and its boundary bond without mixing in contour/Trotter or Hamiltonian
convention choices; it is the narrower regression for the reported defect.

## Test execution

Only login-node short tests were used.  External build/results directories:

- `/home/sunxr/pbc_structural_guards_short_build`
- `/home/sunxr/pbc_structural_guards_short_results`

Durable test source:
`reproduction/projector_kitaev/core_regression_tests/pbc_structural_guards_driver.cpp`.
Machine-readable outcomes are summarized in `PBC_STRUCTURAL_REGRESSION_SUMMARY.csv`.

## Observable-path audit

`S_pi`, `S_pi_dq`, `R_CDW`, and average-sign implementations have no diff in
this branch.  The only production-code changes are the shared odd-L PBC
constructor guard and the interaction-energy bond-count bound.
