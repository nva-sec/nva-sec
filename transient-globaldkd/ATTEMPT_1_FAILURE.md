# Spanish validation attempt 1 — technical schema failure

Recorded on 2026-08-05 (UTC), before any retry.

## Immutable failed attempt

- Workflow run: `31009572812`
- Job: `92317956542`
- Trigger commit: `bc45a8092114af4a46d95fc03e8314926f92636d`
- Reviewed boundary commit: `32cb40dc46b1eb4a5d505b6b81db114ba2f23d99`
- Frozen discovery SHA-256:
  `41460a6e8e21c6979f8ead43ee383a1691baab37897fb5c7ab29b66b636c531f`
- Failure site: `validate_spain_locked.py:165`
- Exception: `Spanish PLS axis-label slots differ from audited empty-string schema`
- Workflow artifact count: zero. The always-upload step found no output files.

The fail-closed preflight, pinned source/blob/digest checks, dependency
installation, archive download, and archive MD5 check all passed.

## Outcome-access boundary

Attempt 1 requested and deserialized only the top-level MATLAB variable
`ESP`; it did not request `AUS_Int` or `Global`. SciPy deserialized the
albumin-axis cell as part of `ESP`, and the validator converted that cell to
a local float64 array. Its expected shape and all-finite guard passed.

The exception occurred before the statement that constructs the 30 mg/L
threshold outcome. Therefore attempt 1 did not construct threshold classes,
count positives or negatives, calculate scores, apply either frozen model,
calculate sensitivity, specificity, AUROC, bootstrap intervals, permutation
statistics, calibration, robustness results, or negative-control results.
No Spanish albumin value, subject ID, spectrum, threshold count, score, or
metric was printed or written.

This was not a pristine unopened attempt: the ESP outcome container had been
deserialized. Any subsequent execution is labeled a technical retry, not a
first unblinding.

## Independently confirmed cause

The preregistered schema deepener and redacted probe run `31010114575`
independently confirmed, without printing outcome values, spectra, or IDs,
that:

- `ESP.axisscale` is an object array with shape `(2, 2)`;
- flattened cells 1 and 3 are NumPy Unicode arrays with dtype `<U1`,
  shape `(0,)`, and size zero;
- `ESP.label` is an object array with shape `(2, 2)`;
- flattened label cell 0 is a Unicode string vector with shape `(61,)`; and
- flattened label cells 1, 2, and 3 are Unicode arrays with shape `(0,)`
  and size zero.

The failed guard compared `str(cell)` with `""`. For an empty NumPy array,
`str(cell)` is `"[]"`; the data matched the audited structure and the
assertion encoded it incorrectly.

## Authorized technical correction

The retry may change only the structural schema guard:

- preserve and assert the `(2, 2)` object-container shapes;
- require the empty axis-label cells and empty label cells to be Unicode
  arrays of shape `(0,)` and size zero;
- require label cell 0 to be a Unicode vector of shape `(61,)`; and
- update the corresponding assertion metadata.

The frozen discovery artifact, selected windows, scaler, coefficients,
thresholds, endpoint, gates, resampling, bootstrap, permutation, robustness,
negative controls, and reporting logic remain unchanged.

A retry requires an independently reviewed code/workflow boundary and a new,
exactly-once retry sentinel. The original failed run and original sentinel
must remain in history.
