# first divergence — boundary-corrected audit

Parameters: L=10, V=4, PBC, theta=10, beta_trial=8, dt=0.1, hs_scheme=0, guard OFF, seed=983002.

The prior `GREEN_PRE=0.440544` record was invalid: in `leftSweep()` it compared the Green matrix after `right_propagate` through operator `l=878` with a full-contour rebuild at boundary `879`. The correct post-propagation/local boundary is `l=878`; only the pre-propagation state belongs to `879`.

With the corrected checks, a complete deterministic right sweep followed by a complete left sweep was examined. At each left operator, pre-propagation fast Green was compared with `(l+1)%N`, while post-propagation, local-flip, and stabilization states were compared with `l`.

No event occurred in this sweep. Corrected Green relative errors remain at approximately `1e-13`–`1e-14`; tracked and direct Pfaffian signs agree at every checked stage. The prior claim “Green first at l=878” is withdrawn.

The corrected per-stage trace is `collected/first_divergence.csv`; the original incorrect-boundary trace is retained as `collected/first_divergence_legacy_wrong_boundary.csv`.

Conclusion: the initial apparent split was a boundary-labeling error. This one-sweep audit does not localize the later long-run warning after production burn/measurements. No PfQMC algorithm was changed.
