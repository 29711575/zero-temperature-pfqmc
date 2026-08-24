# Static projector regression/stress campaign

This directory is independent of driven dynamics and does not modify PfQMC core algorithms.

Campaigns: complete tiny-HS enumeration, V=0 Gaussian exact checks, matched-seed HS-channel comparison, matched-seed delta-sign symmetry, stabilization/guard stress regression, and an OBC density-channel large-L sanity scan.

All production QMC tasks use one CPU and write atomically from `result.json.tmp` to `result.json`. Existing nonempty results are skipped, so failed/missing array elements can be resubmitted without rerunning completed tasks. `scripts/collect_regression_stress.py` regenerates collected tables, plots, progress, and `summary.md` from existing files only.

The fast/full diagnostics compare the completed-right-sweep live Green and raw sign at boundary 0, ensuring identical HS configuration and boundary. Center observables remain captured at the projector center and are not changed by diagnostics.
