# Branch migration map

| Source repository | Source branch / commit | New branch or tag | Purpose | Validation | Production-ready |
|---|---|---|---|---|---|
| `PfQMC-main` | `main@2ae3bb7` | tag `original-pfqmc-main-2ae3bb7` | Original validated scale-safe baseline | A/B/C/D PASS | Historical baseline |
| `PfQMC-integration-staging` | `integration-staging-hygiene@55fee175` | `main`, tag `integration-stable-20260826` | Stable integrated hygiene baseline | Integration hygiene PASS | Yes |
| new repository | stable main | `develop` | Unified future development | Same starting tree as main | Development base |
| `PfQMC-condition-ratio-v3` | `optimize-condition-aware-ratio-v3@243e660` | `feature/condition-aware-ratio-v3` | Optimized adaptive MP ratio v3 | Numerically clean; strict k132 gate incomplete | No |
| `PfQMC-condition-ratio-v3-followup` | dirty overlay on `243e660` | `archive/condition-aware-ratio-v3-followup@f45fec8` | 153577 logging-only observer snapshot | Exit 0; target event not captured | No |
| `PfQMC-main` | visible dirty tree on `d8eebc6` | `archive/condition-aware-ratio-v1@4a31e6f` | 153474/v1 best-effort snapshot | Exact binary/source association unavailable | No |
| `PfQMC-main` | `diagnose-near-zero-sign@557b730` | `archive/near-zero-sign` | Near-zero sign audit | Diagnostic complete | No |
| `PfQMC-main` | `diagnose-sign-green-replay@8490e0c` | `archive/sign-green-replay` | Deterministic sign/Green replay | Diagnostic complete | No |
| `PfQMC-main` | `diagnose-ratio-reality-green-gate@77c7adf` | `archive/ratio-reality-green-gate` | Ratio reality and Green gate audit | Diagnostic complete | No |
| `PfQMC-main` | `fix-green-recovery-prototype@1e7c674` | `archive/green-recovery-prototype` | Fail-closed Green recovery prototype | Prototype validated on saved events | No |
| `PfQMC-main` | `diagnose-mp-ratio-oracle@d8eebc6` | `archive/mp-ratio-oracle` | MP ratio oracle | Oracle audit complete | Diagnostic only |
| `PfQMC-core-regression-tests` | `diagnose-udt-orthogonality@aa0495c` | `archive/udt-orthogonality` | Exponent-QR orthogonality audit | Real failure mode located | No |
| `PfQMC-core-regression-tests` | `fix-udt-mp-fallback-prototype@efe5f4b` | `archive/udt-mp-fallback` | MP UDT fallback prototype | Mixed synthetic recovery | No |
| `PfQMC-core-regression-tests` | `diagnose-udt-householder-feasibility@6be96b7` | `archive/udt-householder-feasibility` | Householder feasibility/rank margin | Audit complete | No |
| `PfQMC-integration-prep` | `integration-prep@5dcf46b` | `archive/integration-prep` | UDT guard plus optional left recovery preparation | Deterministic integration tests PASS | Reference |

Additional remote-tracking refs are retained for provenance. No source branch
was checked out, rewritten, deleted, or modified during migration.
