# Post-cleanup integrity check

| Test | Key result | Status |
|---|---|---|
| Clean build | driven_driver, static_contour_compare, driven_fastupdate_check; mpiicpc -std=c++14 | PASS |
| Static regression V=0 | max operator/Green/observable diff = 0 | PASS |
| Static regression V=4 | max operator/Green/observable diff = 0 | PASS |
| ED spot L=4,Vf=3,R=2 | z(S_pi)=-0.579; z(S_pi_dq)=0.641; z(R_cdw)=-0.601; R is ratio-of-means | PASS |
| Fast update L=4 | max Green=4.64e-14; max ratio=0; post-reset=0 | PASS |
| Fast update L=6 | max Green=1.87e-13; max ratio=0; post-reset=0 | PASS |
| L=8 smoke | finite=True; sign=0.820; acceptance=0.899; max imag sign=2.42e-11; observable=3.1e-15; guard=0 | PASS |
| Retained scripts | see integrity_checks.txt; no archive/backup code reference | PASS |

ED data use two independent seeds (730001, 730002), burn=200, measurements=2000 each.
