# Projector PfQMC core regression tests

This directory consolidates durable, small tests. Generated binaries and
outputs must live outside the repository. The suite covers tiny complete HS
enumeration, the V=0 Gaussian check, same-contour ED, right-boundary Green
comparison, local-update algebra, static-projector reality symmetry,
scale-safe UDT stress, the rank-loss fail-closed guard, L10/task88/task92
smokes, PBC structural/configuration guards, direct interaction-energy Wick
sums, and the optional left-sweep Green recovery.

Provenance: core test sources derive from
`test-core-regression@837552a` and descendants through
`optimize-udt-rank-loss-guard@5b21b110`; integration drivers derive from
`integration-prep@5dcf46b7`. Adaptive MP-ratio production logic is excluded.
