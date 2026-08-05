# Artifact erratum: computational-window wording

Recorded state boundary: 2026-08-05T13:03:57Z (UTC), immediately after
private artifact mirror commit
`8507bbb172d10a74ad64448f2db259772cf3507f` and before any Spanish
outcome access or validation execution.

## Exact wording correction

The immutable discovery artifact at JSON path
`$.scope.certificate_scope` contains this exact text:

> fixed computational window library and frozen pipeline on source spectra; k acquired windows create 2k digital features

For every interpretation, report, caption, and downstream metadata field, read
that sentence as:

> fixed computational window library and frozen pipeline on source spectra; k computational source-spectrum windows create 2k digital features

The word “acquired” was a metadata wording error. Here, `k` counts exact
computational summaries of fixed intervals in already acquired source
absorbance spectra. It does not count optical filters, detector channels, new
physical measurements, or deployable hardware.

## Why numerical results are unchanged

This correction changes no executable value. The affected string is created
only when the completed result dictionary is serialized. It is not parsed or
used by window integration, feature construction, fold materialization,
scaling, model fitting, candidate selection, threshold selection, nested
audit, negative controls, or the unlabeled-control plausibility check.

The original frozen discovery JSON and SHA sidecar remain byte-for-byte
unchanged. This erratum does not replace, rewrite, or recursively rehash the
scientific artifact.

## Immutable provenance

- Australian-only source run: `31005831329`
- Source job: `92305353385`
- Source code commit: `b012a183484f6525eb52ba4eb4ba40090feda845`
- Source `discover_windows.py` Git blob: `a8f820e21aa51fcbfd2ed0b4e9144e8a51239d74`
- Source `discover_windows.py` SHA-256 recorded inside the artifact:
  `f86232f2c94cf82a4d3b4dd8ae357a13fbf6edf22db2a27bf9b1b527d2eb1fa3`
- Actions artifact ID: `8930796749`
- Artifact ZIP SHA-256:
  `e31aa988444ee07942ade34fbb12c0fc8babdf60016b65402f4447d3abfca99d`
- Frozen discovery JSON SHA-256:
  `41460a6e8e21c6979f8ead43ee383a1691baab37897fb5c7ab29b66b636c531f`
- Fold manifest SHA-256:
  `6ccd3b6fdbeeadb2fa8d8808480f09e9a9eb0fc3d3631cd0f2540121b8011de3`
- Public exact-file promotion commit:
  `ecf9db566c59b99f42e601ebd5113ba558c2d442`
- Private exact-file mirror boundary:
  `8507bbb172d10a74ad64448f2db259772cf3507f`

The artifact itself already records `not_hardware_simulation: true`; the
preregistration, executable docstring, validator, and report scope likewise
exclude hardware simulation and physical-filter claims.

## Spanish lock

At this erratum boundary, the discovery artifact records
`spanish_outcomes_read: false`. No Spanish outcome has been loaded by the
discovery or promotion workflows, and the one-shot validator remains
fail-closed until its reviewed discovery digest is pinned in committed code.
