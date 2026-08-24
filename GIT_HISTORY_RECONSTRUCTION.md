# Git history reconstruction provenance

The repository history was reconstructed on 2026-08-24 after the working tree
had been maintained without Git metadata.

The pre-scale-safe backup at
`/home/sunxr/PfQMC-main-pre-scale-safe-20260824` is a selective backup, not a
complete tree snapshot. It contains the pre-integration versions of
`inc/qr_udt.h`, `inc/pfqmc.h`, `Makefile`, and
`reproduction/projector_kitaev/build.sh`, plus executable backups. Therefore
the baseline tree is reconstructed from the current trackable source tree with
those four files replaced by their backed-up versions. Runtime data, binaries,
build products, and scheduler logs are intentionally excluded.

The 153450 campaign was treated as immutable during this operation. Its PBS
script, manifest, raw/results, and executable were not modified or rebuilt.
Recorded SHA-256 values at the start of reconstruction:

- `offcritical_sign_scan.pbs`: `ae3da8e202e65f255a7cd6d7e731bb6badf8c41081c663a2e55908988f1b6fc3`
- `manifest.csv`: `d1dc623b23ce635025e6870815b0e5b6177dd7cf4667b48f942907181d8706e4`
- `bin/sign_scan_driver_scale_safe`: `596d8b132eb8a0c0015634bf3c1a95b086c4ac66e78a14f46f08e630ad19e7c0`

The campaign executable path is
`reproduction/projector_kitaev/regression_stress/offcritical_sign_scan/bin/sign_scan_driver_scale_safe`.
