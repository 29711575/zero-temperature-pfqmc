# Pure-state projector PfQMC Phase 2

This directory contains only the slow, rebuild-from-scratch correctness oracle.
It deliberately has no fast Green update, incremental sweep, condition-aware
ratio, left recovery, MP checkpoint, or production-contour machinery.

Build with `EIGEN3_INCLUDE_DIR` and `PFQMC_PFAPACK_DIR` set, then run
`build/phase2_core_test`.  The validation driver added with the implementation
produces the four review CSV files and gates its final complete JSON on all CSV
streams being flushed, checked, closed, and checked successfully.
