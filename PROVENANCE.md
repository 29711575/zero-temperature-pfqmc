# Provenance

## Stable lineage

- Original main: `2ae3bb7d81c4c4153b773b81dfe897405a94adab`
- Original tag in this repository: `original-pfqmc-main-2ae3bb7`
- Integration hygiene source: `/home/sunxr/PfQMC-integration-staging`,
  `integration-staging-hygiene@55fee175752a896c12e2acf79d2a252c7db74e2f`
- Stable integration tag: `integration-stable-20260826`
- Git ancestry check: `2ae3bb7` is an ancestor of `55fee175`.

The stable path contains, in order, validated UDT rank-loss guard integration,
zero-average-sign output handling, optional default-off left recovery,
consolidated regression tests, and build/JSON hygiene.

## Protected jobs and executables

| Job | State at inventory | Purpose | Executable SHA-256 |
|---|---|---|---|
| 153474 | running | condition-aware ratio v1 replay | `0f2a6a67f6813611db863ac17695e6d47b0414fb13f501b0d4f4a3096a3c0431` |
| 153573 | finished, exit 0 | optimized v3 replay | `f97549ac6874fb1b8825de1f9e42cc99ada8b1cc38d548023cd8bdd9fb0f1b59` |
| 153577 | finished, exit 0 | v3 logging follow-up | `62c21fd85973babb39b9a717c0fa9a60c03ff292e02efb0a78509430a505b56e` |

153577 is an uncommitted three-file observer overlay on commit `243e660`; its
patch, PBS, manifest, and hashes are stored only on
`archive/condition-aware-ratio-v3-followup`. It did not capture the requested
k132 event and does not clear the v3 integration gate.

153474's executable hash is exact, but the source tree lacks a commit proving
the binary/source mapping. The archive branch records that limitation and must
not be presented as an exact rebuild recipe.

## Reconstruction limits

Committed histories were fetched directly from local read-only repositories.
Untracked raw measurements, scheduler output, binaries, build products, caches,
and third-party trees were not imported. Dirty experimental overlays were
archived only where explicitly documented; absence of a recorded source hash
is never replaced by an inferred commit.
