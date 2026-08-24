# PfQMC final results index

Raw historical and high-statistic results are preserved outside the repository at:
/home/sunxr/PfQMC-main_archived_results_20260821/reproduction_driven_kitaev_results/

Key datasets:
- ed_benchmark_v2/: ED-QMC comparison tables and runs.
- ed_diagnostic/: four-seed ED spot validation and blocking diagnostics.
- L6_V4_R2_convergence/: final burn=3000, measurements=10000 convergence data.
- final_small_check_repeat/: clean small validation and fast-update diagnostics.
- Projector/static validation: /home/sunxr/PfQMC-main_archived_results_20260821/reproduction_projector_kitaev_results/.

No raw result files were deleted; historical outputs were moved to the archive.

## Historical/non-production exclusion

Legacy static_note_projector PBC runs with adaptive_guard=true are historical
diagnostics only. They are not reproduced by the current source tree and must
not enter the final paper validation table or production campaigns. The
production guard state is explicit in every current JSON output.
