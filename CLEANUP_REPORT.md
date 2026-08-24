# Cleanup report

Before bytes: 144744738
After bytes: 144682845
Deleted files: 7
Freed bytes: 61893

## Deleted source/modules
- analyze_final_small.py
- analyze_l6_v4_long.py
- BUILD_REPORT.md
- cleanup_manifest.md
- compare_versions.py
- final_small_check.pbs
- l6_v4_long_repeat.pbs

## Preserved dependency closure
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

## Results
Projector/driven result trees retained; additional campaign deletion deferred because final-report/reference closure has not been mechanically proven.

## Verification
A19/D/E checks retained from previous cleanup.
