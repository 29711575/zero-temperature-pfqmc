# 153408 long-contour abort audit

Date: 2026-08-23

## Scope and isolation

This audit used the original `known_physics.csv` parameters.  No HS, Pfaffian,
projector, proposal, RNG, stabilization interval, or sampling logic was changed.
The instrumented source and binaries were kept under
`/home/sunxr/PfQMC-long-contour-audit`; production binaries, raw task directories,
and the running 153407 campaign were not modified.

Two complementary diagnostics were used:

1. an isolated source build checking `zgeqp3`/`zungqr` info, D before reciprocal,
   U/D/T finiteness, and every inverse input in `UDT::onePlusInv` and
   `onePlusInv(UDT&,UDT&)` (the `udtRdivudtL` path), with a light LAPACK `zgecon`
   rcond estimate;
2. the unmodified production `regression_driver` with an `LD_PRELOAD` wrapper
   that only observes `LAPACKE_zgetrf` input and return value, calls the original
   routine unchanged, and prints a backtrace only on the first bad input.

## Exact first failure: task 88

Parameters: `L=26, V=4, theta=39, seed=750049`, OBC, density HS, `dt=0.1`,
`beta_trial=8`, stabilization interval 10, guard off.

The production binary reproducibly fails at:

- phase: burn, second right/left sweep pair (zero-based sweep 1);
- direction: `PfQMC::leftSweep`;
- operator / boundary: `l=10`, boundary 10;
- function: `onePlusInv(UDT&,UDT&)`, i.e. the `udtRdivudtL` scaled core;
- inverse input: `tem2`, dimensions `52 x 52`, column major, `lda=52`;
- state before Eigen inverse: matrix contains NaN/Inf (`finite=0`), while the
  largest finite absolute entry seen by the wrapper is `3.0721709358114948`;
- LAPACK: `LAPACKE_zgetrf` call 11962 returns `info=-4` (invalid fourth argument,
  the matrix pointer/input, due to LAPACKE NaN checking);
- Eigen then asserts `info >= 0`.  The Eigen assertion is therefore downstream,
  not the first numerical failure.

The sweep/operator mapping is exact: the failing call is the 638th left-sweep
two-UDT inverse; an L=26 contour has 319 such calls per left sweep, and call 638
is the last (`l=10`) call of the second left sweep.

## Scaling versus inverse input

The scaling warning precedes the non-finite inverse.  During construction, before
burn sampling, task 88 first overflows Eigen's ordinary Frobenius norm in
`leftInit` at operator/boundary 810.  At that point all matrix elements and U/D/T
elements are still finite, both QR calls return info 0, and the QR diagonal/D
range is approximately

`D_min=1.1040747509507983e-155`, `D_max=9.1109563317260098e154`.

Initialization continues with finite elements but increasingly extreme UDT
scales.  Its full observed range is

`D_min=1.6036493906356989e-207`, `D_max=6.2726840940022158e206`.

No zero, subnormal, NaN, or Inf D and no nonzero `zgeqp3`/`zungqr` info was seen
before the production failure.  The checked scaled inverse cores in the isolated
build were not intrinsically close to singular: minimum estimated rcond was
`3.1916414305300952e-4`.  Thus the evidence is not "a healthy UDT followed by an
independently ill-conditioned inverse".  It is accumulated UDT dynamic-range
growth (ordinary norm overflow already in initialization), followed later by
NaN/Inf formation in the two-UDT scaled-core construction; only then does Eigen
call LU.

## Theta=26 control versus task 88

Control task 48 uses the same `L=26, V=4, seed=750049` and all other settings,
except `theta=26`.  It previously completed all 3000 measurements.  Its audited
initialization has:

| case | first ordinary norm overflow | global D min | global D max | min inverse-core rcond |
|---|---:|---:|---:|---:|
| task 48, theta=26 | none | 1.2874816925492252e-140 | 6.7592223891743904e139 | 3.1916414305300952e-4 |
| task 88, theta=39 | leftInit l=810 | 1.6036493906356989e-207 | 6.2726840940022158e206 | 3.1916414305300952e-4 |

The key distinction is about 67 extra decades at each end of D, while the
minimum scaled-core rcond is essentially unchanged.

## Other failed tasks

Production-binary checks confirm the same failure mechanism:

| task | parameters | first production failure | matrix | finite | zgetrf info |
|---|---|---|---:|---:|---:|
| 89 | L=26, theta=39, seed=750050 | burn sweep 2, leftSweep, l/boundary=120 | 52x52 | no | -4 |
| 92 | L=34, theta=51, seed=750069 | burn sweep 0, leftSweep, l/boundary=870 | 68x68 | no | -4 |

For task 89 the audited initialization range is
`D=[1.7729573770117609e-208, 4.8420327882358735e207]`, with first ordinary norm
overflow at `leftInit l=760`; for task 92 it is
`D=[2.9261818664144133e-267, 3.4899645940519165e266]`, with first overflow at
`leftInit l=1730`.  Task 90 was also checked in the isolated audit and shows the
same pre-failure scaling pattern (`D=[1.0338897554250255e-212,
1.0134449793674043e212]`, first overflow at `leftInit l=860`).

Therefore the sampled failed tasks are consistent with one long-contour UDT
dynamic-range mechanism, not seed/manifest/PBS-specific failures.

## Minimal candidate fix (not implemented)

The smallest principled fix target is the UDT product/scaled-core representation,
not a catch around Eigen inverse.  Candidate implementation:

1. make UDT multiplication scale-safe by separating D mantissas from base-2
   exponents (or log-scales), so `Dl * (Tl*Ur) * Dr` is never materialized with
   both `1e-200` and `1e+200` scales in ordinary doubles;
2. assemble the `udtRdivudtL` core from bounded mantissas/exponents, applying
   row/column powers-of-two scaling before the solve and undoing that scaling
   afterward;
3. replace explicit `.inverse()` with pivoted factorization/solve after the
   representation is made scale-safe, checking LAPACK info and rcond.

A shorter stabilization interval or an `isfinite` catch would only move or hide
the failure, because the problematic scale is accumulated across stabilized UDT
products.  No workaround or candidate fix was applied in this audit.
