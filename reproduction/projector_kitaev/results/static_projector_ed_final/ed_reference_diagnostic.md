# Static projector ED reference diagnostic

Parameters: `L=6`, `V=2`, `delta=1`, `mu=0`, PBC, `theta=10`,
`beta_trial=8`, and `dt=0.1`.

## First failed layer and correction

Majorana algebra, kinetic-only reconstruction, and the two single-bond
interaction identities all pass.  The remaining error in the prior static ED
was the even-`L` PBC kinetic boundary reconstruction.  C++ writes
`A[2L-1,L]=-2i` and its antisymmetric partner is `A[L,2L-1]=+2i`.  The old
Python mirror recorded only the former and then summed only `a<b`; because
`2L-1>L`, that PBC bilinear was omitted entirely.

The repaired mirror explicitly completes the antisymmetric matrix and uses the
equivalent checked bilinear forms `1/4 sum_ab gamma_a A_ab gamma_b` and
`1/2 sum_a<b gamma_a A_ab gamma_b`.

## Layered checks

| check | result |
|---|---:|
| `max ||{gamma_a,gamma_b}-2 delta_ab I||` | `0` |
| kinetic bilinear all-pairs vs upper triangle | `0` |
| kinetic hermiticity error | `0` |
| kinetic ground energy (`V=0`) | `-6.0` |
| `max |D_01-M_01|` | `0` |
| `max |D_50-M_50|` | `0` |
| Majorana/full-complex spectrum max difference | `6.22e-15` |

`eig(A) = (-2,-2,-2,-2,-2,-2,2,2,2,2,2,2)` and the first kinetic many-body
levels are `-6, -4, -4, -4, -4, -4, -4, -2` (roundoff omitted).

The raw Majorana and complex-Fock matrices differ by a basis/gauge
representation (`max=1.41421356` at matrix element `(0,3)`), but their spectra
and density observables agree.  The full-Hamiltonian values are
`E0=-6.464101615137757`, `S_pi=0.08333333333333341`,
`S_pi_dq=0.04629629629629638`, and `R_cdw=0.444444444444444`.

## Finite projector contour

Using the validated Hamiltonian with the production-form symmetric contour
`e^(-dt K/2)e^(-dt V E)e^(-dt V O)e^(-dt K/2)`, the PBC groups
`E={(0,1),(2,3),(4,5)}` and `O={(1,2),(3,4),(5,0)}`, trial density
`beta_trial=8`, and projector `theta=10`, the exact finite-contour result is:

| observable | ED | z vs pooled job 153357 |
|---|---:|---:|
| `S_pi` | `0.08301358175` | `+0.6813` |
| `S_pi_dq` | `0.04632267310` | `-0.4390` |
| `R_cdw` | `0.44198681556` | `+0.6278` |

No QMC job was rerun and no PfQMC core or projector-contour source was
modified.  The previous erroneous ED artifacts were retained.
