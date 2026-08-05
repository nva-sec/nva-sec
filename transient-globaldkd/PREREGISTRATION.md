# Preregistration: exact sparse-window cross-instrument albumin-concentration study

Status: **frozen before data inspection or outcome analysis**

Base commit: `c4bf99396ef266b78116f0822d01c74cc93622b3`

Public data source: Navarro-Esteve, Pérez-Guaita, Sánchez-Illana et al., *Dataset for: Towards a Global Model for Diabetic Kidney Disease Screening using ATR-FTIR*, Zenodo record 14762603, version 1.

Archive URL: `https://zenodo.org/records/14762603/files/GlobalDKD.zip?download=1`

Expected archive MD5: `8d5d4ba0e0578d0e39e76e1f1514a841`

## Scientific question

Can a small, explicitly counted set of predefined mid-infrared absorbance-window summaries preserve an Australian-trained screen for urine albumin concentration at the source study's 30 mg/L threshold when transferred, without recalibration, to spectra acquired in Spain on a different instrument?

The intended endpoint is the **laboratory-measured albumin concentration in the submitted urine specimen**, following the source study's standardized protein-extraction procedure. It is not urinary albumin-to-creatinine ratio (uACR), a causal test for diabetic kidney disease (DKD), CKD staging, or a clinical diagnosis.

## Scope of any minimum claim

Any minimum claim will be restricted to:

1. the fixed candidate-filter library defined below;
2. a fixed preprocessing and model family;
3. budgets of one through four passbands; and
4. the datasets and transfer direction stated here.

It will not be described as a universal optical or mathematical minimum. The study emulates rectangular passbands from full ATR-FTIR spectra; it does not validate fabricated hardware.

## Cohorts and lock order

- Discovery/training: Australia, expected 177 biological donors (155 participants with DKD and 22 healthy controls in the source report).
- Confirmatory external validation: Spain, expected 61 biological donors (35 with DKD and 26 controls in the source report), comprising Valencia and Madrid sources.
- Technical replicates and reruns will be collapsed to one row per biological donor before any split.
- Blank spectra are not donors and will be excluded.
- The Spanish albumin outcome will remain unread during schema audit and model discovery.
- A frozen design artifact and its SHA-256 digest will be created from Australian data before the confirmatory Spanish outcome is unlocked.

The source paper reports 64/113 Australian and 34/27 Spanish samples above/below its 30 mg/L concentration threshold. These counts are structural reconciliation checks, not permission to inspect or optimize against Spanish outcomes during discovery.

## Candidate spectral windows

Twenty-four fixed, non-overlapping rectangular windows, each 64 cm^-1 wide, cover only the two modeling regions stated in the source paper.

High-wavenumber centers:

`3616, 3552, 3488, 3424, 3360, 3296, 3232, 3168, 3104, 3040, 2976, 2912 cm^-1`

Fingerprint-region centers:

`1752, 1688, 1624, 1560, 1496, 1432, 1368, 1304, 1240, 1176, 1112, 1048 cm^-1`

Each interval is `[center - 32, center + 32)` cm^-1, except the final upper endpoint may be included without double counting. All 12,950 subsets of one through four windows will be enumerated exactly.

For one selected window, the computational summary is the trapezoidal mean absorbance over that source-spectrum interval. A model uses the selected summary vector and its deterministic L2-normalized copy; this uses no wavelengths outside the counted windows. Full-spectrum SNV, Savitzky-Golay derivatives, or uncounted reference windows are prohibited for the sparse design. This is not a physical-filter measurement or hardware simulation.

## Frozen analysis family

The classifier will be balanced L2 logistic regression with `C=1` and a fixed deterministic seed. Per-feature centering and scaling will be fit on Australian training partitions only. No SMOTE, model-family search, or Spanish recalibration is allowed.

An all-24-passband model using the identical pipeline is the prespecified comparator.

The selected design is the smallest subset satisfying all Australian-only criteria:

- repeated grouped cross-validated AUROC at least 0.90;
- specificity at least 0.80 at an Australian-derived threshold constrained to sensitivity at least 0.90; and
- AUROC no more than 0.03 below the all-24-passband comparator.

Ties will be resolved by higher AUROC and then lexicographic order of filter centers. Nested outer cross-validation will repeat the complete selector to estimate selection-procedure performance.

A biochemical-interpretation claim additionally requires at least one selected window overlapping Amide I or Amide II. Otherwise the design can only be described as an empirical concentration screen.

## External confirmation and success gates

After freezing the Australian design, it will be applied once to all eligible Spanish donors without coefficient, preprocessing, band, or threshold changes.

Primary external outputs:

- AUROC with a 10,000-resample stratified subject bootstrap confidence interval;
- sensitivity and specificity at the frozen Australian threshold;
- paired AUROC difference from the all-24-passband comparator;
- calibration intercept/slope; and
- a 10,000-permutation test of the frozen Spanish score-label association.

A transfer success requires point-estimate AUROC at least 0.90, sensitivity at least 0.85, specificity at least 0.75, and sparse-minus-all-24 AUROC at least -0.10. Confidence intervals will be reported regardless of whether these gates pass.

Valencia-only discrimination and Madrid-stratified performance are mandatory checks. Results will not be pooled silently if site behavior conflicts.

## Robustness and negative controls

Without reselection or recalibration, repeat external scoring after shifting all chosen passband centers by +/-8 and +/-16 cm^-1 and changing widths to 48 and 80 cm^-1. Exact-center-only success will be treated as evidence against physical robustness.

A prespecified negative-control model uses four 64 cm^-1 windows centered at `2368, 2304, 2240, 2176 cm^-1`, a nominally silent region outside the source paper's modeling regions. Unexpected external discrimination by this control is an artifact warning.

Continuous concentration regression, if run, is secondary and uses the already selected outputs. It cannot alter the passbands or primary classifier.

## Data-integrity gates

Before modeling, the audit must establish:

- exact archive checksum;
- archive member inventory;
- MATLAB variable names, shapes, dtypes, and storage version;
- raw versus already transformed spectral representations;
- monotone wavenumber axes with coverage of every candidate and control window;
- donor identifiers and technical-replicate grouping rules;
- exclusion of blanks;
- unique biological-donor counts matching the source report; and
- an unambiguous join to the reference-value table.

Ambiguity, disagreement with published cohort totals, an inability to separate donors from replicates, or evidence that only outcome-transformed features are available is a fail-closed condition.

## Current phase

The first executable phase is schema audit only. It may inspect file names, MATLAB structures, workbook sheet names/dimensions, spectra shapes, identifiers, and donor counts. It must not read albumin values into a model, compute threshold labels, select passbands, fit a classifier, or report performance.

## Dated endpoint-integrity and resampling amendment — 2026-08-05

This amendment was frozen after schema/code inspection and **before any outcome
model was fit**. It corrects the discovery endpoint without changing the
candidate bands, feature construction, model family, subset gates, comparator,
external success gates, or robustness tests above.

### Endpoint populations

- The primary Australian discovery population is the 155 rows of `AUS_Int`
  having a finite laboratory-measured urine albumin concentration. The locked
  structural reconciliation is 64 positive and 91 negative at 30 mg/L.
- The 22 Australian healthy-control means have missing (`NaN`) measured
  albumin. Their source labels are negative, but those labels are not measured
  concentration outcomes. They are excluded from the primary discovery,
  selection, threshold, and nested-audit analyses.
- The 22 missing-albumin control spectra remain unlabeled for computation.
  After discovery, the final primary model fitted only on the 155 measured
  samples is applied to these 22 spectra. Only their score distribution and
  fraction below the frozen threshold are reported as an **unlabeled
  source-convention plausibility check**. No outcome labels are assigned, and
  no sensitivity, specificity, AUROC, or other performance metric is computed
  for these 22 rows.
- The untouched Spanish confirmation population remains 61 donors, with a
  locked structural check of 34 positive and 27 negative at 30 mg/L. Spanish
  target arrays and target-derived labels remain unread during discovery.

Thus the earlier 177-row Australian discovery wording and 64/113 structural
count are superseded for the primary endpoint by 155 measured rows and 64/91.

### Exact deterministic cross-validation

For primary subset evaluation, use
`RepeatedStratifiedKFold(n_splits=5, n_repeats=5, random_state=20260805)`.
Each subject receives exactly one out-of-fold probability in each repeat; the
five probabilities are averaged per subject before computing AUROC,
sensitivity, specificity, and comparator deltas. For each candidate, choose the
highest score threshold whose sensitivity on the aggregated Australian
out-of-fold probabilities is at least 0.90. Threshold ties are resolved by the
highest specificity and then the numerically highest threshold.

The nested selection-procedure audit uses
`StratifiedKFold(n_splits=5, shuffle=True, random_state=20260806)` outside.
For zero-based outer fold `k`, its inner splitter is
`StratifiedKFold(n_splits=5, shuffle=True, random_state=20260807 + k)`.
Within each outer-training set, the same exhaustive selector, gates, threshold
rule, model, scaling discipline, and lexicographic tie rule are run using the
single inner five-fold out-of-fold predictions. The selected pipeline is then
refit on the complete outer-training set and evaluated once on that outer-test
set at its inner-derived threshold. Outer-test predictions are concatenated
once per subject for audit metrics.

All 12,950 nonempty subsets of one through four of the 24 locked bands are
evaluated. The all-24 pipeline remains the comparator and follows the identical
folds, preprocessing discipline, probability aggregation, and threshold rule.

### Exact implementation details

The positive endpoint is albumin concentration greater than or equal to 30 mg/L.
For each locked interval, spectra are linearly interpolated at both interval
boundaries and trapezoidally integrated across the boundary and interior grid
points, then divided by the nominal interval width. The model input
concatenates the selected raw passband means with their per-row L2-normalized
copy (a zero norm maps to zeros).

The fixed estimator is
`StandardScaler()` followed by
`LogisticRegression(C=1, penalty="l2", class_weight="balanced",
solver="liblinear", random_state=20260805, max_iter=5000)`.
Scaling parameters are learned only from each training partition. Candidate
selection follows the original rule: the fewest bands passing every gate, then
higher AUROC, then the lexicographically earliest tuple in the locked band
library order.

The 22 excluded controls are never assigned computational outcome labels and
never enter fitting, cross-validation, selection, or threshold derivation. The
final model fitted only on the 155 measured rows may score them once at the
frozen primary threshold; this is an unlabeled plausibility check, not a
performance analysis.

### Discovery isolation and frozen artifact

Discovery loads only `AUS_Int` from root `Zenodo/dataset.mat` using
`scipy.io.loadmat(..., variable_names=("AUS_Int",), simplify_cells=True)`.
It must not load or dereference `ESP` or `Global`, open
`Reference_Values.xlsx`, or inspect Spanish filenames, targets, class arrays,
or target-derived labels. A fail-closed static/runtime guard records the
requested MAT variable names and prohibited resources.

The discovery result is serialized as canonical UTF-8 JSON using
`json.dumps(result, sort_keys=True, separators=(",", ":"), ensure_ascii=False)`
followed by exactly one newline. Its SHA-256 is computed over those exact bytes
and written to a separate `.sha256` file. Spanish confirmation remains locked
until both canonical files and their digest have been reviewed and committed.

## Dated pre-fit red-team amendment — 2026-08-05

This section was committed before the first discovery fit and supersedes any
inconsistent language above.

### Computational-window scope

The report title is **“How Many Spectral Windows Are Enough? An Exact
Sparse-Window Certificate for Cross-Instrument Screening of Elevated Urinary
Albumin Concentration.”** The counted quantity `k` is the number of acquired
source-spectrum windows. A `k`-window model has `2k` digital features: `k`
mean-absorbance summaries plus their `k` row-wise L2-normalized copies.

Mean absorbance is not a physical-filter measurement:
`mean(-log(T)) != -log(mean(T))` in general, and the source, filter, detector,
and reference responses are unavailable. “Emulated instrument,” “passband
measurement,” hardware minimum, filter design, hardware simulation, and
deployability claims are prohibited. Any certificate is restricted to this
fixed computational tile library, frozen pipeline, source spectra, and
cross-instrument held-out reuse. It is not prospective clinical validation.

### Numerical integration and robustness lock

For center `c` and width `w`, sort wavenumbers ascending, linearly
interpolate absorbance at both closed endpoints `c-w/2` and `c+w/2`, include
all interior samples, and compute `np.trapezoid(...)/(w)`. Nominal windows use
`w=64`; robustness widths divide by their actual width. The endpoint
interpolation avoids the erroneous 62 cm^-1 sampled span produced by a
half-open grid mask.

Set the numerical zero-norm epsilon to `1e-12`. A zero/near-zero L2 norm,
nonfinite feature or score, or logistic-regression convergence warning aborts
the run.

Every robustness perturbation uses the nominal frozen scaler, coefficients,
intercept, and threshold; only shifted/resized window summaries and their L2
copies are recomputed. There is no refit, recalibration, or reselection. Common
center shifts test only tolerance to a global spectral-axis offset.

### Identity, duplication, and fold audit

Australian filename IDs must be unique after stripping whitespace; they are
never model features. Before splitting, scan all 177 common-grid spectra for:

- exact duplicates, defined by element-for-element float64 equality;
- near duplicates after row centering and division by row RMS, defined jointly
  as Pearson correlation at least `0.999999` and normalized RMS distance at
  most `0.0015`.

Only pair counts and hashed pair identifiers are reported. Any exact duplicate
between two measured rows having different threshold labels is a fail-closed
error. No duplicate is silently removed. Every materialized train/test fold is
recorded as exact subject-ID lists, and the canonical fold manifest receives
its own SHA-256 digest.

A full-data search or nested outer fold with no candidate passing every gate
returns **NO DESIGN**. There is no fallback subset. The nested audit reports the
fraction of outer folds yielding a design and selection frequencies for all
selected window tuples.

### Controls and frozen artifact

The original four-window control at centers
`2368, 2304, 2240, 2176 cm^-1` is retained but relabeled the
**atmospheric/instrument-artifact control** because the first two windows are
near atmospheric CO2. A separate cleaner off-target control is fixed, without
outcome-based selection, at centers
`2496, 2560, 2624, 2688 cm^-1`, each 64 cm^-1 wide (union
`2464–2720 cm^-1`). Neither control can select or revise the primary design.

The canonical discovery artifact records the integration definition, epsilon,
library order, `2k` feature mapping, exact IDs/folds and fold-manifest hash,
duplicate counts, package versions, all candidate metrics, nominal scaler,
coefficients, intercept, threshold, comparator, both negative controls, nested
audit, and every assertion. Serialization uses `allow_nan=False`; nonfinite
values abort. The SHA-256 sidecar covers the canonical JSON bytes including
their one terminal newline.

### Locked future Spanish analysis

Spanish confirmation remains prohibited until the discovery JSON and SHA are
reviewed and committed. Its 10,000-resample bootstrap uses seed `20260812` and
resamples within site-by-class strata; sparse and all-window predictions use
the identical resampled indices for paired AUROC differences. Its
10,000-permutation primary test uses seed `20260813`, shuffles labels within
site only, and reports `p=(1+exceed)/(B+1)`.

Valencia receives within-site AUROC, sensitivity, and specificity. Madrid
metrics are reported only when mathematically defined (likely sensitivity).
Sensitivity and specificity receive 95% Wilson binomial intervals. Calibration
intercept/slope are explicitly exploratory because balanced-class-weight
logistic scores are not calibrated probabilities. Spanish below-LOQ or
NaN-to-zero reference entries are left-censored threshold references, not exact
concentrations.

