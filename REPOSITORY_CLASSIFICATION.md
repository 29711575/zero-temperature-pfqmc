# Repository classification

## KEEP
- src/**
- inc/pfapack/**
- inc/kitaevChain.h
- inc/pfqmc.h
- inc/spinless_tV.h
- inc/operator.h
- inc/types.h
- inc/skewMatUtils.h
- inc/qr_udt.h
- test/**
- reproduction/projector_kitaev/**
- reproduction/driven_kitaev/**
- driven_fastupdate_check.cpp
- CMakeLists.txt
- Makefile
- main.cpp
- README.md
- LICENSE
- .github/**

## DELETE
- analyze_final_small.py — unrelated historical root analysis/standalone/PBS
- analyze_l6_v4_long.py — unrelated historical root analysis/standalone/PBS
- BUILD_REPORT.md — unrelated historical root analysis/standalone/PBS
- cleanup_manifest.md — unrelated historical root analysis/standalone/PBS
- compare_versions.py — unrelated historical root analysis/standalone/PBS
- final_small_check.pbs — unrelated historical root analysis/standalone/PBS
- l6_v4_long_repeat.pbs — unrelated historical root analysis/standalone/PBS

## Retained due dependency
- inc/square.h, inc/honeycomb.h, inc/singleMajoranaHoneycomb.h remain because root CMake/main build references have not been rewritten.
