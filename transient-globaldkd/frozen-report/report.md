# How Many Spectral Windows Are Enough? An Exact Sparse-Window Certificate for Cross-Instrument Screening of Elevated Urinary Albumin Concentration.

## Scope

This secondary computation asks how many fixed source-spectrum windows are enough for the archived 30 mg/L urinary-albumin-concentration screen. Each retained window is a mean-absorbance summary. It is not a physical filter measurement, hardware simulation, uACR assay, or disease diagnosis.

The certificate is bounded to the fixed 24-window library, frozen pipeline, and these archived acquisition settings. Spanish reuse is held out but is not prospective clinical validation.

## Endpoint and provenance

- Archive MD5: 8d5d4ba0e0578d0e39e76e1f1514a841
- Discovery SHA-256: 41460a6e8e21c6979f8ead43ee383a1691baab37897fb5c7ab29b66b636c531f
- Primary Australia: n=155 measured samples, 64 positive and 91 negative at 30 mg/L.
- Excluded controls: 22 missing-albumin spectra, kept unlabeled for a score-distribution plausibility check only.
- Exact/near duplicate pair counts: 0/0.

## Frozen computation

- 24 fixed 64 cm⁻¹ computational windows.
- All 12,950 subsets of one through four windows.
- k raw mean-absorbance summaries plus k row-wise L2-normalized copies, giving 2k digital features.
- Fold-only scaling and balanced L2 logistic regression with fixed seeds, gates, and threshold selection.

## Australian discovery

- Selected k=2 (4 digital features), centers 3360, 3296 cm⁻¹.
- Repeated-CV AUROC 0.950; sensitivity 0.906; specificity 0.868.
- Frozen threshold 0.416377180289.
- Qualifying candidates: 2975.
- Artifact-sourced qualifying counts by k=1/2/3/4: 0/15/282/2678.
- 24-window repeated-CV AUROC: 0.946.
- Nested folds yielding a design: 5/5.
- Nested selection-procedure AUROC: 0.942.

## Unlabeled-control plausibility check

The final primary model was fitted only on 155 measured samples and scored the 22 controls once. No outcome label or performance metric is assigned.
- Median score 0.024; IQR 0.018–0.073.
- Fraction below the frozen threshold 0.955.

## Spanish archived-cohort reuse — schema-only technical retry 1

**Execution disclosure.** The successful Spanish evaluation was technical retry 1 after correcting only the fail-closed MATLAB empty-cell schema assertion. Attempt 1 deserialized the ESP outcome container but stopped before threshold-class construction or any score, metric, or output artifact. It was therefore not a pristine first unblinding. The immutable validation JSON retains its pre-retry phase string; this report supplies the superseding execution label.

- Failed attempt: run 31009572812, job 92317956542.
- Successful technical retry: run 31010663111, job 92321693370, artifact 8932293583.
- Attempt-1 disclosure SHA-256: d84eab971456f7ecc6418b1aa41b3a33f819d27a354d7daadcf78ab55c17100e.
- Validation SHA-256: fb1111b4e0539f854c60bbb11c6211632d3162cb5ccb1505736bed9085a33da4
- N=61 (34 positive, 27 negative); Valencia 46, Madrid 15.
- AUROC 0.980, stratified-bootstrap 95% CI 0.948–1.000.
- Sensitivity 0.941; specificity 0.852; Wilson intervals are in the artifact.
- Sparse-minus-24-window AUROC -0.005.
- The 24-window comparator had the same sensitivity (0.941) but one false positive, versus four for the sparse model.
- The selected 3360 and 3296 cm⁻¹ windows are adjacent closed 64 cm⁻¹ intervals sharing the 3328 cm⁻¹ boundary. They are computational summaries, not evidence for two independent physical channels.
- Prespecified Spanish gates from the frozen validation artifact: PASS.
- Valencia: AUROC 1.000; sensitivity 1.000; specificity 0.840. Undefined single-class quantities are reported explicitly.
- Madrid: AUROC 0.962; sensitivity 0.846; specificity 1.000. Undefined single-class quantities are reported explicitly.
- Within-site label-permutation p=9.999e-05 (0 exceedances in 10000 permutations).
- Exploratory calibration status estimated; intercept -0.596; slope 1.352. Balanced-class-weight scores are not calibrated probabilities.
- Source zero references are left-censored threshold references, not exact concentrations.

### Frozen robustness and negative-control diagnostics

Every robustness perturbation reused the nominal Australian scaler, model, and threshold with no refit, recalibration, or reselection. Common shifts test global spectral-axis tolerance only.

| Perturbation | AUROC | Sensitivity | Specificity | TP/FN/TN/FP |
|---|---:|---:|---:|---:|
| common shift +16 cm inverse | 0.976 | 0.971 | 0.630 | 33/1/17/10 |
| common shift +8 cm inverse | 0.978 | 0.971 | 0.704 | 33/1/19/8 |
| common shift -16 cm inverse | 0.979 | 0.706 | 1.000 | 24/10/27/0 |
| common shift -8 cm inverse | 0.983 | 0.912 | 1.000 | 31/3/27/0 |
| width 48 cm inverse | 0.980 | 0.941 | 0.778 | 32/2/21/6 |
| width 80 cm inverse | 0.980 | 0.941 | 0.926 | 32/2/25/2 |

Negative-control windows were frozen during Australian discovery and were not selected using Spanish outcomes. The atmospheric/instrument block is an artifact diagnostic, not a clean biochemical negative; the clean off-target block is reported separately.

| Negative control | AUROC | Sensitivity | Specificity | TP/FN/TN/FP |
|---|---:|---:|---:|---:|
| atmospheric instrument artifact | 0.946 | 0.971 | 0.741 | 33/1/20/7 |
| clean off target | 0.926 | 0.941 | 0.593 | 32/2/16/11 |

**Operating-point warning.** The ±8/±16 cm⁻¹ common shifts preserved rank discrimination (AUROC 0.976–0.983) but did not preserve the frozen-threshold operating point: specificity fell to 0.704 and 0.630 for +8/+16, while sensitivity fell to 0.706 for −16. This further limits any physical-channel interpretation.

**Interpretive warning.** Both Spanish negative controls showed high discrimination (AUROC above 0.92). That pattern is compatible with broad nuisance, site, or class-correlated spectral signal and materially weakens biochemical-region specificity and mechanistic interpretation. The preregistered external gates still pass, but they do not override this claim ceiling.

## Interpretation and limitations

The permitted result is an exact positive or negative certificate within this fixed computational tile library and pipeline. Selected correlated windows are candidate spectral regions, not causal molecular attribution.

This work does not establish physical-filter behavior, detector or SNR requirements, deployability, clinical readiness, a universal minimum, prospective utility, uACR, or causal kidney-disease diagnosis.

## Reproducibility

- Code SHA-256: f86232f2c94cf82a4d3b4dd8ae357a13fbf6edf22db2a27bf9b1b527d2eb1fa3
- Preregistration SHA-256: 200dde3ff02feac834f6b6b804d33af62226965576cd4acf30388225d0718b68
- Fold manifest SHA-256: 6ccd3b6fdbeeadb2fa8d8808480f09e9a9eb0fc3d3631cd0f2540121b8011de3
- Python/NumPy/SciPy/scikit-learn: 3.11.15 / 2.0.2 / 1.13.1 / 1.5.2

