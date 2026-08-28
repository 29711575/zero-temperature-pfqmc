# HS marginal follow-up

Parameters: L=10, V=2, theta=18, beta_trial=8, dt=0.1, PBC; guard OFF; 12 new matched seeds and 10000 measurements per scheme. Common delete-one `(seed, 25-measurement bin)` jackknife units: 4800.

## Observables

| observable | hs0 | hs1 | hs0-hs1 | sigma | ED | hs0-ED sigma | hs1-ED sigma |
|---|---:|---:|---:|---:|---:|---:|---:|
| S_pi | 0.045301785 ± 0.00013 | 0.0456807643 ± 0.00011 | -0.000379 ± 0.00017 | -2.23 | 0.0455136259 | -1.63 | +1.56 |
| S_pi_dq | 0.0353378092 ± 4.2e-05 | 0.0352767484 ± 4e-05 | 6.11e-05 ± 5.8e-05 | +1.05 | 0.0353150237 | +0.54 | -0.97 |
| R_cdw | 0.21994665 ± 0.0024 | 0.227754855 ± 0.0024 | -0.00781 ± 0.0034 | -2.27 | 0.224077999 | -1.69 | +1.54 |

## Diagnostics

* hs0: average sign 0.78965; acceptance 0.745627; max tau_int 0.547 retained measurements; direct comparisons 6000; ±1 mismatches 0; sign corrections 0; max sign imaginary 8.563e-08; max fast/full Green relative error 8.968e-13.
* hs1: average sign 1; acceptance 0.839466; max tau_int 0.909 retained measurements; direct comparisons 6000; ±1 mismatches 0; sign corrections 0; max sign imaginary 3.947e-11; max fast/full Green relative error 9.940e-14.

## Verdict

* Original S_pi_dq 3.43 sigma difference does not reproduce.
* The previous marginal difference is consistent with a statistical multiple-comparison fluctuation.
* No ±1 sign mismatch or material fast/full Green anomaly is present.
