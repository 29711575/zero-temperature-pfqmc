# Projector Convergence Analysis

- Pooling uses per-seed signed numerators reconstructed as `measurements × average_sign × observable`; errors are delete-one-independent-seed jackknife (12 seeds).
- The current driver records `n_bins=15` and per-seed bin errors but does not emit bin records, so a cross-seed pooled jackknife is the available reproducible aggregate.
- `guard_on_diagnostics` is numerical diagnostic only / excluded from production averages.

## Theta
- theta=2.0: S_pi 0.04545785 ± 0.00017182; S_pi_dq 0.03533981 ± 0.00005535; R_cdw 0.22258071 ± 0.00325820; average_sign=0.977900; acceptance=0.747220.
- theta=4.0: S_pi 0.04591181 ± 0.00019206; S_pi_dq 0.03524624 ± 0.00008489; R_cdw 0.23230570 ± 0.00330543; average_sign=0.949333; acceptance=0.746299.
- theta=6.0: S_pi 0.04521075 ± 0.00009639; S_pi_dq 0.03523253 ± 0.00004279; R_cdw 0.22070449 ± 0.00192159; average_sign=0.925000; acceptance=0.745981.
- theta=8.0: S_pi 0.04498197 ± 0.00015116; S_pi_dq 0.03534875 ± 0.00004654; R_cdw 0.21415745 ± 0.00293944; average_sign=0.900933; acceptance=0.745837.
- theta=10.0: S_pi 0.04562226 ± 0.00019285; S_pi_dq 0.03537150 ± 0.00003902; R_cdw 0.22468777 ± 0.00332413; average_sign=0.879867; acceptance=0.745826.
- theta=12.0: S_pi 0.04543440 ± 0.00007507; S_pi_dq 0.03533902 ± 0.00003737; R_cdw 0.22219683 ± 0.00111609; average_sign=0.857667; acceptance=0.745708.

Adjacent theta sigma differences (6→8, 8→10, 10→12):
- 6.0→8.0 S_pi: -1.276 sigma.
- 6.0→8.0 S_pi_dq: +1.838 sigma.
- 6.0→8.0 R_cdw: -1.864 sigma.
- 8.0→10.0 S_pi: +2.613 sigma.
- 8.0→10.0 S_pi_dq: +0.375 sigma.
- 8.0→10.0 R_cdw: +2.373 sigma.
- 10.0→12.0 S_pi: -0.908 sigma.
- 10.0→12.0 S_pi_dq: -0.601 sigma.
- 10.0→12.0 R_cdw: -0.710 sigma.

## Trotter
- dt=0.05: S_pi 0.04577236 ± 0.00021819; S_pi_dq 0.03532353 ± 0.00005795; R_cdw 0.22827825 ± 0.00392228; average_sign=0.867900.
- dt=0.1: S_pi 0.04537163 ± 0.00010043; S_pi_dq 0.03529844 ± 0.00003441; R_cdw 0.22201500 ± 0.00147836; average_sign=0.901567.
- dt=0.2: S_pi 0.04507714 ± 0.00017320; S_pi_dq 0.03511146 ± 0.00003473; R_cdw 0.22108060 ± 0.00248352; average_sign=0.958967.
- The contour uses half-kinetic / interaction / half-kinetic physical slices, i.e. a symmetric second-order splitting; the three-point dt^2 fit is diagnostic only.
- dt^2 diagnostic S_pi: X0=0.04568787 ± 0.00015015.
- dt^2 diagnostic S_pi_dq: X0=0.03534633 ± 0.00004618.
- dt^2 diagnostic R_cdw: X0=0.22634624 ± 0.00267813.

## Sign-imag checks
- seed 990104 (dt=0.2): max_sign_imag=2.512e-07, max_observable_imag=4.994e-15, sign_corrections=0, leave-one-out S_pi shift=-0.083 jackknife sigma.
- seed 990106 (dt=0.05): max_sign_imag=1.602e-06, max_observable_imag=2.175e-14, sign_corrections=0, leave-one-out S_pi shift=-0.045 jackknife sigma.

## Conclusions
- theta=8 is borderline rather than fully established: its 8→10 shifts are +2.613 sigma in S_pi and +2.373 sigma in R_cdw.  Recommend theta=10 for production.
- With beta_trial=8, theta=10 and 12 agree within 1 sigma for all three observables, supporting a plateau from theta=10.
- dt=0.1 is statistically consistent with dt=0.05; retain dt=0.1 for production unless a tighter systematic-error target warrants dt=0.05.
- The two sign-imag threshold seeds show no observable-imag contamination or sign correction; they are not evidence of an observable bias in this data set.
- E can proceed as the next independent trial-density convergence check.
