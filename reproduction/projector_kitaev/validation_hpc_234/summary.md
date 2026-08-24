# validation_hpc_234 summary
Completed QMC: 138/138; failed: 0. Same-contour ED: 8/8.
## Trotter O(dt^2)
Correlated delete-one jackknife with common matched `(seed, 25-measurement bin)` units (2400 units). Each replica recomputes all three dt ratios and refits using the full 3x3 jackknife covariance.
* S_pi: O0=0.045795059 +/- 0.000161, a=-0.020978339 +/- 0.00544, chi2/dof=0.574/1; dt=0.1 minus O0 = 0.322 sigma.
* S_pi_dq: O0=0.035286364 +/- 5.32e-05, a=-0.0040440546 +/- 0.0018, chi2/dof=1.13/1; dt=0.1 minus O0 = -1.23 sigma.
* R_cdw: O0=0.22915104 +/- 0.00272, a=-0.26251275 +/- 0.0953, chi2/dof=0.814/1; dt=0.1 minus O0 = 0.704 sigma.
dt=0.2 >2sigma fit-residual observables: none.
## Autocorrelation / blocking
* L10_V2_existing: max seed tau_int=0.639; recommended production bin >= 25.
* L10_V4_new: max seed tau_int=0.532; recommended production bin >= 25.
## L=6 ED-QMC benchmark
* PBC V=0 S_pi: QMC 0.041666667 +/- 2.3e-17; ED 0.041666667; z=nan.
* PBC V=0 S_pi_dq: QMC 0.041666667 +/- 2.3e-17; ED 0.041666667; z=nan.
* PBC V=0 R_cdw: QMC 0 +/- 0; ED -4.4408921e-16; z=nan.
* PBC V=2 S_pi: QMC 0.083235375 +/- 0.000359; ED 0.083013582; z=0.618.
* PBC V=2 S_pi_dq: QMC 0.046296397 +/- 9.63e-05; ED 0.046322673; z=-0.273.
* PBC V=2 R_cdw: QMC 0.44378941 +/- 0.00337; ED 0.44198682; z=0.536.
* PBC V=4 S_pi: QMC 0.15352648 +/- 0.0062; ED 0.15152993; z=0.322.
* PBC V=4 S_pi_dq: QMC 0.032118346 +/- 0.00177; ED 0.032752815; z=-0.359.
* PBC V=4 R_cdw: QMC 0.79079605 +/- 0.0198; ED 0.7838525; z=0.351.
* PBC V=6 S_pi: QMC 0.20074838 +/- 0.00553; ED 0.19701263; z=0.676.
* PBC V=6 S_pi_dq: QMC 0.017595578 +/- 0.00195; ED 0.018819462; z=-0.627.
* PBC V=6 R_cdw: QMC 0.91235009 +/- 0.0122; ED 0.90447586; z=0.648.
* OBC V=0 S_pi: QMC 0.041666667 +/- 2.3e-17; ED 0.041666667; z=nan.
* OBC V=0 S_pi_dq: QMC 0.041666667 +/- 2.3e-17; ED 0.041666667; z=nan.
* OBC V=0 R_cdw: QMC 0 +/- 0; ED -2.220446e-16; z=nan.
* OBC V=2 S_pi: QMC 0.085938543 +/- 0.00044; ED 0.08548717; z=1.03.
* OBC V=2 S_pi_dq: QMC 0.04813979 +/- 0.000165; ED 0.048203515; z=-0.387.
* OBC V=2 R_cdw: QMC 0.4398347 +/- 0.0047; ED 0.43613158; z=0.787.
* OBC V=4 S_pi: QMC 0.14393832 +/- 0.00345; ED 0.13816284; z=1.67.
* OBC V=4 S_pi_dq: QMC 0.036454083 +/- 0.00124; ED 0.038461254; z=-1.62.
* OBC V=4 R_cdw: QMC 0.74673816 +/- 0.0145; ED 0.72162374; z=1.73.
* OBC V=6 S_pi: QMC 0.18453507 +/- 0.00983; ED 0.17652909; z=0.814.
* OBC V=6 S_pi_dq: QMC 0.023455253 +/- 0.00314; ED 0.026434304; z=-0.949.
* OBC V=6 R_cdw: QMC 0.87289541 +/- 0.0243; ED 0.85025525; z=0.933.
## Failures / anomalies
* Failed tasks: none.
* Diagnostic WARN tasks (28; full detail: collected/qc_warnings.csv): autocorr:0[SIGN_IMAG_WARN], autocorr:1[SIGN_IMAG_WARN], autocorr:2[SIGN_IMAG_WARN], autocorr:3[SIGN_IMAG_WARN], autocorr:4[SIGN_IMAG_WARN], benchmark:25[SIGN_IMAG_WARN], benchmark:28[SIGN_CORRECTION], benchmark:29[SIGN_IMAG_WARN], benchmark:32[SIGN_CORRECTION], benchmark:33[SIGN_CORRECTION], benchmark:34[SIGN_IMAG_WARN], benchmark:35[SIGN_IMAG_WARN], benchmark:36[SIGN_IMAG_WARN], benchmark:37[SIGN_CORRECTION], benchmark:40[SIGN_CORRECTION], benchmark:42[SIGN_IMAG_WARN], benchmark:43[SIGN_CORRECTION], benchmark:44[SIGN_CORRECTION], benchmark:45[SIGN_IMAG_WARN], benchmark:46[SIGN_IMAG_WARN], benchmark:47[SIGN_IMAG_WARN], benchmark:84[SIGN_IMAG_WARN], benchmark:85[SIGN_CORRECTION], benchmark:87[SIGN_IMAG_WARN], benchmark:88[SIGN_IMAG_WARN], benchmark:89[SIGN_CORRECTION], benchmark:92[SIGN_IMAG_WARN], benchmark:94[SIGN_IMAG_WARN].
