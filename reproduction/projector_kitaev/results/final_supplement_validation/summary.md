# Final Supplement Validation

## Completed

| Item | Result |
|---|---|
| Production Trotter `dt=0.05` | `153356[]`, 12/12 complete |
| Static projector QMC | `153357[]`, 12/12 complete |
| Static projector PBC ED | `153358` PBS failed because compute nodes lack NumPy; small ED then ran on HPC login with `/usr/bin/python`, NumPy 1.16.6/SciPy 1.2.3, but the independent reference disagrees strongly with QMC and is not accepted |

## Production Trotter comparison

Matched `theta=10`, `beta_trial=8`, `L=10,V=2`, seeds `990201–990212`.

| Observable | `dt=0.1` | `dt=0.05` | sigma difference |
|---|---:|---:|---:|
| `S_pi` | `0.04554486 ± 0.00025317` | `0.04591453 ± 0.00028722` | `0.97` |
| `S_pi_dq` | `0.03530602 ± 0.00009463` | `0.03522343 ± 0.00012780` | `-0.52` |
| `R_cdw` | `0.22480777 ± 0.00587587` | `0.23240449 ± 0.00702870` | `0.83` |

The `dt=0.05` run has average sign `0.83193`, zero sign corrections, maximum sign imaginary part `6.94e-9`, and maximum observable imaginary part `1.55e-12`. The matched Trotter check passes.

## Audit conclusions

- Existing `dt=0.1` data are strictly matched by parameters, seeds, guard settings, observable convention, and source identifier `PfQMC-main_20260821`.
- The theta=8 cross-campaign shift is explained by distinct independent seed sets; no new six-seed replay was submitted.
- The beta-trial 4/6/8/10/12 records have distinct `trial_slices` and result hashes; their near identity is a genuine numerical plateau, not data reuse.
- `supp_guard_on_completion` remains `0/19` and is excluded from guard-off production validation.
- Static projector–ED consistency remains unresolved: the old ED mismatch was traced first to omitted bulk pairing entries. The corrected mechanical mirror returns `S_pi=0.10110146`, `S_pi_dq=0.04485296`, `R_cdw=0.55635695`, while QMC gives `0.08276162`, `0.04637077`, `0.43970683`; z-scores are `-49.6`, `13.9`, and `-32.1`. The Hamiltonian checksum still misses the supplied `E0=-6.464101615` target, so this is not accepted as a physics failure or benchmark PASS.
