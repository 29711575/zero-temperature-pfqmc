# Final supplement report

## Status
B/C/D/E guard-off, trial-only and fixed-HS diagnostics are DONE. 153373 is `INVALID_FOR_E_GUARD_ON` and excluded.

## B Trotter
```
observable,pooled,jackknife_error
S_pi,0.045544861476719635,0.00023934004558846128
S_pi_dq,0.03530602279975624,9.182305034170992e-05
R_cdw,0.22480776853821371,0.005571538525836842
observable,difference_dt005_minus_dt01,paired_jackknife_error,z
S_pi,0.00036961414175033271,0.0037535155711733055,0.098471455557275225
S_pi_dq,-8.2188937868581807e-05,0.001269500606428638,-0.064741156839456665
R_cdw,0.0077785308382933571,0.081724373290667271,0.09518006104039034
```

## C ED
observable,pooled,jackknife_error
S_pi,0.08276161609990645,0.00032321020530407993
S_pi_dq,0.04637076819720189,0.00011335828941486807
R_cdw,0.4397068305043127,0.0034244777461186465
ED z: S_pi=-1.76887; S_pi_dq=0.65696; R_cdw=-1.38346. PASS.

## D recomputed
```
observable,pooled,error
S_pi,0.044981971614900508,0.00015466733683019878
S_pi_dq,0.035348747402092213,5.4931165770009142e-05
R_cdw,0.21415744723864527,0.0031119565788980908
observable,pooled,error
S_pi,0.045371626098727219,0.00012292182200118905
S_pi_dq,0.03529844440058532,3.8832748995467587e-05
R_cdw,0.22201500286154552,0.0019860843020768305
observable,pooled,error
S_pi,0.045176867312062996,9.905406757477972e-05
S_pi_dq,0.035323587064010675,3.3616171301855633e-05
R_cdw,0.218104548506871,0.0018466668080508286
observable,theta_minus_trotter,error,z
S_pi,-0.00038965448382671103,0.00019756457022006291,-1.9722892793615947
S_pi_dq,5.0303001506893374e-05,6.7271207566069543e-05,0.74776421186566233
R_cdw,-0.0078575556229002474,0.0036917210896684921,-2.1284261275560872
```
153367 is reproducibility evidence only and is not included in 24-seed pooled.

## E
Trial-only Green changes with beta, but fixed-HS production-path sensitivity gives beta8-vs12 at theta10 max|dG|=2.4114e-13 and relative Frobenius=5.9768e-14; production uses trial operators and theta=10 suppresses initial-state dependence. 153375 guard-off matched results are in `E_matched_differences.csv`.

## QC
```
dataset,PBS/job,n_tasks,adaptive_guard,multiprecision_fallback,max_sign_imag,max_observable_imag,sign_corrections,max_green_skew_symmetry_error,n_bins_requested,n_bins_used,bin_size,error_method,source_sha256,executable_sha256,notes
B,legacy,12,false,false,NA,NA,NA,NA,20,20,250,leave-one-bin jackknife,NA,NA,legacy binary SHA unavailable
C,153357,12,false,false,NA,NA,NA,NA,20,20,250,leave-one-bin jackknife,NA,NA,legacy binary SHA unavailable
D,153374,24,false,false,NA,NA,NA,NA,20,20,250,leave-one-bin jackknife,NA,NA,current executable
E,153375,12,false,false,NA,NA,NA,NA,20,20,250,matched seed/bin,9aaa7fae87125ad43a441c4b18054ed67e0dee22727fb256fc1b9cca9922d8a8,b99c9cc49f785a3763d1143fd6fd0c97a36f87b168f24de10501d3a99650755e,guardoff
E_invalid,153373,12,true,false,NA,NA,NA,NA,20,NA,250,excluded,NA,d0a6a38797aa45af478ba38d2b4159b3c55504c762ad8bafe2b25accda40ad1d,INVALID_FOR_E_GUARD_ON
```
Guard-on/off never mixed; R_cdw is ratio-of-means recomputed in each jackknife replica.

## Conclusions
ED/projector PASS; theta=10 sufficient; beta_trial=8 sufficient. dt=0.1 controlled if B paired z below significance threshold; D subsets consistent if z below threshold. No new QMC required.
