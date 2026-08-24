# Original rigorous array 153416: collection and QC

Planned tasks: 792. Success: 248. Expected incompatible-grid exclusion (exit_code=2): 336. Other failures: 208. Exit-code counts: {'0': 248, '2': 336, '3': 208}.

All observed other failures have exit_code=3; representative `result.json.tmp` records identify `failure_reason=numerical_diagnostic_exceeds_tolerance`, rather than a scheduler failure.

## Reusable compatible subset

Expected 240 tasks (L=26,34; R=0.5,1.0; Vf=3.8,3.9,4.0,4.1,4.2; 12 seeds): 130/240 complete. QC-flagged successful reuse tasks: 0. Details: `rigorous_153416_reusable_240_qc.csv`.

## QC

Warnings use flags for incomplete measurements, non-finite S_pi/sign, invalid acceptance, imaginary-part thresholds >1e-6, sign corrections, and missing bin records. Expected incompatible-grid tasks are excluded, not numerical failures.

## Preliminary D(V)

`rigorous_153416_preliminary_D.csv` pools each rate separately from sign-reweighted bin numerators/denominators. Rates are treated as independent ensembles; Q and D errors use independent propagation. Rows are marked incomplete unless all 12 seeds exist at each of R=0.5,1,2. No Vc fit is performed here.
