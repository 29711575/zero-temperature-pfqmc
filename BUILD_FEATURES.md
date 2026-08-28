# Build features

All supported build paths use C++17 and enable scale-safe UDT by default.

| Entry point | C++ standard | Scale-safe UDT | Notes |
|---|---:|---|---|
| `Makefile` | C++17 | `-DPFQMC_SCALE_SAFE_UDT` | Production `main`; Intel MKL/PFAPACK toolchain |
| `CMakeLists.txt` | C++17 | `PFQMC_SCALE_SAFE_UDT=ON` | Compiler must be selected before `project()`; CMake does not force it afterward |
| `reproduction/projector_kitaev/build.sh` | C++17 | `-DPFQMC_SCALE_SAFE_UDT` | Static-projector driver |

The UDT rank-loss guard is part of the scale-safe production implementation:
45-bit lost-bits fail closed, staged checks, and final normalized adjoint
orthogonality threshold `1e-6`.

Feature defaults reported by projector JSON:

- `scale_safe_udt=true`
- `udt_rank_loss_guard=true`
- `udt_rank_loss_guard_bits=45`
- `udt_orthogonality_gate=1e-6`
- `left_recovery_enabled=false`
- `condition_aware_ratio_enabled=false`

The optional left-sweep recovery must be explicitly enabled. Adaptive MP ratio
and discrete Z2 sign transport are absent from stable main.

Build artifacts belong under ignored build directories and must never be
committed. PFAPACK, Eigen, MKL, and compiler installations are external build
dependencies, not vendored provenance.
