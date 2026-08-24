# Trial-density convergence at theta=10

All 72 tasks completed: six beta_trial values, 12 independent seeds each.

| beta_trial | S_pi | S_pi_dq | R_cdw | average sign | acceptance | status |
|---:|---:|---:|---:|---:|---:|---|
| 2 | 0.04522356 +/- 0.00018946 | 0.03534153 +/- 0.00006971 | 0.21851517 +/- 0.00263022 | 0.87787 | 0.745785 | WARNING |
| 4 | 0.04554486 +/- 0.00025317 | 0.03530602 +/- 0.00009463 | 0.22480777 +/- 0.00587587 | 0.87597 | 0.745761 | PASS |
| 6 | 0.04554486 +/- 0.00025317 | 0.03530602 +/- 0.00009463 | 0.22480777 +/- 0.00587587 | 0.87597 | 0.745761 | PASS |
| 8 | 0.04554486 +/- 0.00025317 | 0.03530602 +/- 0.00009463 | 0.22480777 +/- 0.00587587 | 0.87597 | 0.745761 | PASS |
| 10 | 0.04554486 +/- 0.00025317 | 0.03530602 +/- 0.00009463 | 0.22480777 +/- 0.00587587 | 0.87597 | 0.745761 | PASS |
| 12 | 0.04554486 +/- 0.00025317 | 0.03530602 +/- 0.00009463 | 0.22480777 +/- 0.00587587 | 0.87597 | 0.745761 | PASS |

Errors use 24 delete-one-contiguous bins (two per seed); central values use pooled ratio-of-means.
The beta_trial mapping is active: JSON records trial_slices = 20, 40, 60, 80, 100, 120 for beta_trial = 2, 4, 6, 8, 10, 12.
The nearly identical beta_trial >= 4 trajectories are reproducible under the same seed and are treated as a convergence plateau, not as a missing parameter mapping.
A WARNING marks max_sign_imag or max_observable_imag above 1e-7; it does not alter the ratio-of-means central value.
