# Retained CSV I/O completion gate fix report

## Scope and isolation

- HPC work directory: `/home/sunxr/new-pfqmc-fix-retained-io-completion-gate`
- Branch: `fix-retained-io-completion-gate`
- Frozen source baseline: `4d79bb7ac2b734bc10667d917b914ce316a39a39`
- Required production ancestry present:
  - `0c2df511dff146eb15c9cf395c1f54863422ea9e`
  - `89809695c0dea7cb725a0645a0a48bd2f84f378b`
  - `14969585b0aee9b24b7bf221834041f11d1ca0df`
- The frozen `/home/sunxr/new-pfqmc-fix-mp-z2-trust-gate` worktree was not modified.
- No changes or merges were made to `main` or `develop`.
- No scheduler job was submitted. All tests ran on the login node with one MKL/OMP thread.
- Build products and raw regression outputs are outside the repository in
  `/home/sunxr/new-pfqmc-fix-retained-io-completion-gate-results`.

## Defect and fix

`reproduction/driven_kitaev/driven_driver.cpp` previously computed and printed
`"status":"complete"` while the retained bin-record and optional time-series
streams were still open. Buffered write, flush, or close errors could therefore
surface only after a successful JSON record and return code.

The fix adds a completion gate before status calculation and JSON output:

1. For every requested retained CSV stream, capture any prior write failure.
2. Call `flush()` and check stream state.
3. Call `close()` and check stream state.
4. Aggregate failures from all requested streams and throw a precise error.
5. Only after all requested streams pass the gate may the driver emit
   `"status":"complete"` and return zero.

The existing top-level exception handler writes the failure to stderr and
returns 2. Since the gate precedes all JSON output, an I/O failure cannot emit a
complete record. Open failures continue to fail immediately with explicit
stderr diagnostics.

The sampling loop, RNG calls, auxiliary-field updates, trajectory formatting,
bin formatting, observable calculation, and physics JSON fields were not
changed.

## Regression coverage

The reusable regression is
`reproduction/driven_kitaev/test_io_completion_gate.py`. It compares binaries
built from the frozen source and the fixed source, and writes
`IO_COMPLETION_GATE_REGRESSION_SUMMARY.csv`.

All 10 checks pass:

- normal single retained stream: return 0 and status complete;
- normal dual retained streams: return 0 and status complete;
- CSV parseability, expected row/column counts, and final newline;
- fixed-seed trajectory byte identity with the frozen binary;
- fixed-seed bin CSV byte identity with the frozen binary;
- physics/provenance JSON identity after excluding only runtime and output path;
- bin-record `/dev/full`: return 2, no complete output, explicit flush error;
- time-series `/dev/full`: return 2, no complete output, explicit flush error;
- unwritable bin-record path: return 2, no complete output, explicit open error;
- unwritable time-series path: return 2, no complete output, explicit open error.

Fixed-seed identity evidence:

- trajectory SHA-256:
  `d0098057730d51e296fb06769671812f6cc5b4abebe6e5729e53d747144b759a`
- bin CSV SHA-256:
  `38dec23be75ffc20d48edc9cc14d7d6120ad8a4681a7b139c01b308fba4e9d5d`
- time-series shape: 31 rows x 8 columns;
- bin-record shape: 16 rows x 5 columns;
- all checked CSV files end with a newline and parse successfully.

## Build and test command

The standard Intel oneAPI/MKL toolchain and Eigen 3.4.0 were used. The final
automated run was:

```bash
source /opt/ohpc/pub/apps/intel/oneapi/setvars.sh
python3 reproduction/driven_kitaev/test_io_completion_gate.py \
  --driver /home/sunxr/new-pfqmc-fix-retained-io-completion-gate-results/bin/driven_driver_fixed \
  --baseline-driver /home/sunxr/new-pfqmc-fix-retained-io-completion-gate-results/bin/driven_driver_baseline \
  --output-dir /home/sunxr/new-pfqmc-fix-retained-io-completion-gate-results/automated \
  --summary IO_COMPLETION_GATE_REGRESSION_SUMMARY.csv
```

The authoritative per-check outcome and diagnostics are in
`IO_COMPLETION_GATE_REGRESSION_SUMMARY.csv`.
