# MP Z2 Trust-Gate Fix Report

## Result

The bounded validation passes. Runtime periodic MP checkpoints are strictly record-only and never overwrite transported `z2_sign`. The two previously changed OBC+hs1 seeds remain at physical Z2 `+1` for all 4000 retained measurements, while the archived HS/RNG/accept-reject trajectories are reproduced exactly.

The production source change is commit `14969585b0aee9b24b7bf221834041f11d1ca0df` on branch `fix-mp-z2-trust-gate`. The branch HEAD is `0c584c70ba480956c9fd326fbdfe24b4a6b298f4`; the later commit contains portable single-task submission scripts only. No merge to main/develop was performed.

Key executable SHA256 values:

- `projector_real_z2_driver`: `62057b5c91d695beddce0ef6f6fefe3ec0f81efe3c79ca74c116b56c1e5996e1`
- `validation_oracle`: `497890799e73d1c602296f0bbed1f4785e7488ed4ff22c8a7941b6e7fd55f8d1`

Condition-aware ratio and left recovery remained disabled.

## Changed OBC+hs1 seeds

| seed | old mutating physical Z2 average | fixed physical Z2 average | trajectory | MP corrections |
|---:|---:|---:|---|---:|
| 2126263300 | 0.3835 | 1.0 | exact match | 0 |
| 2126263301 | -0.1320 | 1.0 | exact match | 0 |

Both runs start at Z2 `+1` and remain `+1`. Green rebuild maximum relative errors are `8.18e-14` and `1.24e-13`; UDT guard triggers are zero. The earlier changes were caused by runtime MP candidates overwriting the transported state, not by accepted-ratio Z2 updates.

## MP precision and ordering

At 160 digits, each L18 target configuration gives three `+1` and three `-1` results over six cyclic cuts, so MP160 is boundary-dependent and untrusted. At 320 and 640 digits, all six cuts for both configurations give `+1`. Adaptive evaluation escalates through 160/320/640 and returns trusted `+1` only after 320/640 agreement. Thus the converged canonical result is unique and sign-free.

The eight L6 OBC+hs1 configurations agree configuration-by-configuration among dense Fock, MP160, MP320 and transported Z2, all `+1`. Three L12 MP-only controls are also `+1` with transported Z2 `+1`.

## PBC regression target

For `L=6, theta=12, V=4, seed=1926155102`, the initial physical Z2 is `-1`; the fixed average sign is `+0.3985 ± 0.011385`, matching the expected physical result. All 20 measurement-center and 20 shadow-oracle checks agree. Six accepted proposals required the dedicated adaptive ratio adjudication path. Periodic MP checkpoints remained non-mutating, `mp_correction_count=0`, and the trajectory hash exactly matches `13619898121706317856`.

## Regression and trajectory gates

HPC job `153817.mgt` exited 0 in 1m53s. All 12 regression groups pass: generic complex, exact projector/finite enumeration, OBC+hs1, zero-sign JSON, sign interface, local flips, right-boundary Green, reality symmetry, integration QMC, driven static comparison and driven smoke.

All three production replay trajectory gates pass exactly. Runtime MP correction is zero in every replay. Generic-complex behavior remains unchanged.

## Conclusion

- Changed OBC+hs1 seeds restored to all-time Z2 `+1`: **PASS**.
- MP160/320/640 convergence gate: **PASS**; MP160 is correctly rejected for the L18 pathological configurations, and 320/640 agree at `+1`.
- Runtime periodic MP correction strictly zero: **PASS**.
- HS/RNG/accept-reject trajectory unchanged: **PASS**.
- Required regression suite: **PASS**.

The branch is technically ready for a bounded replay of the 150 archived seeds. It has not been merged to develop.
