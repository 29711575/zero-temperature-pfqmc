# CDW observable sign and contact audit

`UDT::onePlusInv` and the two-UDT `onePlusInv` return
`g_code=2(I+B)^(-1)`.  Sweep and captured Green matrices are not shifted.  The
paper's skew matrix is `G=g_code-I`; hence `g_code,aa=1`.

For `i=j`, writing the two Majorana indices as `a,b`, the legacy helper uses

```
tmp = -g_aa g_bb + g_ab g_ab + g_ab g_ba.
```

With `g_aa=g_bb=1` and `g_ba=-g_ab`, the last two terms cancel and `tmp=-1`.
After the helper normalization this supplies `-1/(4L)` at every momentum.
Thus legacy already includes the negative onsite contact.  Together with the
previously audited overall density-sign convention,

```
S_Q_physical = -S_Q_legacy
```

is the complete paper-defined structure factor.  Projector/driven output must
not add another `1/(4L)`.  The reported `onsite_contact=1/(4L)` is diagnostic
only; optional `*_offsite` fields subtract it from the full physical value.
Each bin forms `R=1-S_pi_dq/S_pi` directly from its full sign-reweighted values.
