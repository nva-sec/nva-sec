#!/usr/bin/env python3
"""Locked one-shot Spanish confirmation for the GlobalDKD sparse-window study.

DRAFT ONLY until a human reviews and commits frozen-discovery.json and its
SHA-256. Importing this module does not open the archive. Execution requires an
explicit unlock flag and an exact reviewed discovery digest.

This validates computational absorbance-window summaries. It does not emulate
physical filters, measure uACR, diagnose DKD, or establish deployable hardware.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import platform
import zipfile
from collections import Counter
from pathlib import Path
from typing import Sequence

import numpy as np
import scipy
import discover_windows as discovery_code
from scipy.io import loadmat
from scipy.optimize import minimize
from sklearn.metrics import confusion_matrix, roc_auc_score

from discover_windows import (
    ATMOSPHERIC_CONTROL_CENTERS,
    CENTERS,
    CLEAN_OFF_TARGET_CENTERS,
    CUTOFF_MG_L,
    EXPECTED_MD5,
    MAT_MEMBER,
    WINDOW_WIDTH,
    canonical_bytes,
    integrate_windows,
    md5sum,
    model_matrix,
    sha256_bytes,
)

BOOTSTRAP_REPLICATES = 10000
BOOTSTRAP_SEED = 20260812
PERMUTATIONS = 10000
PERMUTATION_SEED = 20260813
WILSON_Z_95 = 1.959963984540054
ROBUST_CENTER_SHIFTS = (-16, -8, 8, 16)
ROBUST_WIDTHS = (48.0, 80.0)
REQUESTED_VARIABLES = ("ESP",)
FORBIDDEN_VARIABLES = ("AUS_Int", "Global")
# Fail closed until the reviewed discovery artifact is committed.
REVIEWED_DISCOVERY_SHA256: str | None = "41460a6e8e21c6979f8ead43ee383a1691baab37897fb5c7ab29b66b636c531f"


def load_and_verify_discovery(
    json_path: Path,
    sha_path: Path,
    reviewed_digest: str,
) -> tuple[dict, str]:
    payload = json_path.read_bytes()
    observed = sha256_bytes(payload)
    tokens = sha_path.read_text(encoding="utf-8").strip().split()
    if len(tokens) < 1:
        raise RuntimeError("empty frozen-discovery SHA sidecar")
    declared = tokens[0].lower()
    reviewed = reviewed_digest.strip().lower()
    if not (observed == declared == reviewed):
        raise RuntimeError(
            "discovery digest mismatch: "
            f"observed={observed}, declared={declared}, reviewed={reviewed}"
        )
    discovery = json.loads(payload)
    if payload != canonical_bytes(discovery):
        raise RuntimeError("frozen discovery JSON is not canonical")
    if REVIEWED_DISCOVERY_SHA256 is None:
        raise RuntimeError(
            "literal REVIEWED_DISCOVERY_SHA256 is not pinned in committed code"
        )
    pinned = REVIEWED_DISCOVERY_SHA256.strip().lower()
    if reviewed != pinned:
        raise RuntimeError(
            f"reviewed digest does not match committed literal: {reviewed} != {pinned}"
        )
    current_code_sha = sha256_bytes(
        Path(discovery_code.__file__).resolve().read_bytes()
    )
    if current_code_sha != discovery.get("input", {}).get("code_sha256"):
        raise RuntimeError(
            "current discover_windows.py does not match frozen discovery code"
        )
    if discovery.get("schema_version") != "globaldkd-sparse-window-discovery-v1":
        raise RuntimeError("unsupported discovery schema")
    if discovery.get("phase") != (
        "AUSTRALIAN DISCOVERY ONLY; SPANISH OUTCOMES UNREAD"
    ):
        raise RuntimeError("discovery phase/isolation declaration mismatch")
    if discovery.get("status") != "DESIGN":
        raise RuntimeError("discovery returned NO DESIGN; Spanish validation prohibited")
    selected = discovery["discovery"].get("selected_design")
    if not selected:
        raise RuntimeError("missing selected design")
    if discovery["future_spanish_confirmation_lock"].get("status") != "LOCKED":
        raise RuntimeError("discovery artifact did not retain Spanish lock")
    return discovery, observed


def load_spain_after_unlock(archive: Path) -> dict:
    if md5sum(archive) != EXPECTED_MD5:
        raise RuntimeError("archive MD5 mismatch")
    with zipfile.ZipFile(archive) as bundle:
        try:
            payload = bundle.read(MAT_MEMBER)
        except KeyError as error:
            raise RuntimeError(f"missing archive member: {MAT_MEMBER}") from error

    loaded = loadmat(
        io.BytesIO(payload),
        simplify_cells=True,
        variable_names=REQUESTED_VARIABLES,
    )
    public_keys = sorted(key for key in loaded if not key.startswith("__"))
    if public_keys != list(REQUESTED_VARIABLES):
        raise RuntimeError(f"unexpected loaded variables: {public_keys}")
    if any(name in loaded for name in FORBIDDEN_VARIABLES):
        raise RuntimeError(f"forbidden variable deserialized: {public_keys}")

    esp = loaded["ESP"]
    if not isinstance(esp, dict):
        raise RuntimeError(f"ESP did not simplify to dict: {type(esp).__name__}")
    spectra = np.asarray(esp["data"], dtype=np.float64)
    axes = np.asarray(esp["axisscale"], dtype=object).reshape(-1)
    labels = np.asarray(esp["label"], dtype=object).reshape(-1)
    albumin = np.asarray(axes[0], dtype=np.float64).reshape(-1)
    wavenumbers = np.asarray(axes[2], dtype=np.float64).reshape(-1)
    subject_ids = np.char.strip(
        np.asarray(labels[0], dtype=str).reshape(-1)
    ).astype(str)

    if spectra.shape != (61, 1598):
        raise RuntimeError(f"unexpected ESP.data shape: {spectra.shape}")
    if albumin.shape != (61,) or wavenumbers.shape != (1598,):
        raise RuntimeError(
            f"unexpected ESP axes: albumin={albumin.shape}, wn={wavenumbers.shape}"
        )
    if subject_ids.shape != (61,) or len(set(subject_ids.tolist())) != 61:
        raise RuntimeError("Spanish IDs are blank, nonunique, or wrong length")
    if any(not value for value in subject_ids):
        raise RuntimeError("blank Spanish ID")
    if not (
        np.all(np.isfinite(spectra))
        and np.all(np.isfinite(albumin))
        and np.all(np.isfinite(wavenumbers))
    ):
        raise RuntimeError("nonfinite Spanish spectrum, axis, or threshold reference")
    if not (np.all(np.diff(wavenumbers) < 0) or np.all(np.diff(wavenumbers) > 0)):
        raise RuntimeError("Spanish wavenumber axis is not strictly monotone")
    expected_grid = np.arange(3994.0, 798.0, -2.0, dtype=np.float64)
    if not np.array_equal(wavenumbers, expected_grid):
        raise RuntimeError("Spanish axis is not the audited 3994..800 step -2 grid")
    if str(axes[1]) != "" or str(axes[3]) != "":
        raise RuntimeError(
            "Spanish PLS axis-label slots differ from audited empty-string schema"
        )

    outcome = (albumin >= CUTOFF_MG_L).astype(np.uint8)
    if (int(outcome.sum()), int((1 - outcome).sum())) != (34, 27):
        raise RuntimeError(
            f"Spanish threshold count mismatch: "
            f"{int(outcome.sum())}+/{int((1-outcome).sum())}-"
        )
    site_values = []
    for identifier in subject_ids:
        upper = str(identifier).upper()
        if upper.startswith("E_"):
            site_values.append("Valencia")
        elif upper.startswith("DF_"):
            site_values.append("Madrid")
        else:
            raise RuntimeError(f"unknown Spanish site ID pattern: {identifier!r}")
    site = np.asarray(site_values, dtype=object)
    if (
        int(np.sum(site == "Valencia")) != 46
        or int(np.sum(site == "Madrid")) != 15
    ):
        raise RuntimeError(f"ID-derived site counts mismatch: {Counter(site)}")
    expected_order = np.asarray(["Valencia"] * 46 + ["Madrid"] * 15, dtype=object)
    if not np.array_equal(site, expected_order):
        raise RuntimeError("ID-derived sites violate audited Valencia/Madrid block order")
    return {
        "mat_sha256": sha256_bytes(payload),
        "spectra": spectra,
        "albumin": albumin,
        "wavenumbers": wavenumbers,
        "subject_ids": subject_ids,
        "outcome": outcome,
        "site": site,
    }


def score_serialized_model(model: dict, matrix: np.ndarray) -> np.ndarray:
    mean = np.asarray(model["scaler_mean"], dtype=np.float64)
    scale = np.asarray(model["scaler_scale"], dtype=np.float64)
    coefficient = np.asarray(model["coefficient"], dtype=np.float64)
    intercept = float(model["intercept"])
    if matrix.shape[1] != mean.size or not (
        mean.shape == scale.shape == coefficient.shape
    ):
        raise RuntimeError(
            "serialized model dimensionality mismatch: "
            f"X={matrix.shape}, mean={mean.shape}, scale={scale.shape}, "
            f"coef={coefficient.shape}"
        )
    if np.any(scale <= 0) or not np.all(
        np.isfinite(np.concatenate((mean, scale, coefficient, [intercept])))
    ):
        raise RuntimeError("invalid frozen model parameters")
    linear = intercept + ((matrix - mean) / scale) @ coefficient
    scores = np.empty_like(linear)
    positive = linear >= 0
    scores[positive] = 1.0 / (1.0 + np.exp(-linear[positive]))
    exp_linear = np.exp(linear[~positive])
    scores[~positive] = exp_linear / (1.0 + exp_linear)
    if not np.all(np.isfinite(scores)):
        raise RuntimeError("nonfinite frozen score")
    return scores


def wilson_interval(successes: int, total: int) -> list[float] | None:
    if total == 0:
        return None
    proportion = successes / total
    denominator = 1.0 + WILSON_Z_95**2 / total
    center = (
        proportion + WILSON_Z_95**2 / (2.0 * total)
    ) / denominator
    radius = (
        WILSON_Z_95
        * math.sqrt(
            proportion * (1.0 - proportion) / total
            + WILSON_Z_95**2 / (4.0 * total**2)
        )
        / denominator
    )
    return [float(max(0.0, center - radius)), float(min(1.0, center + radius))]


def classification_metrics(
    outcome: np.ndarray,
    scores: np.ndarray,
    threshold: float,
) -> dict:
    predicted = (scores >= threshold).astype(np.uint8)
    tn, fp, fn, tp = confusion_matrix(
        outcome, predicted, labels=[0, 1]
    ).ravel()
    positives = int(tp + fn)
    negatives = int(tn + fp)
    auroc = (
        float(roc_auc_score(outcome, scores))
        if len(np.unique(outcome)) == 2
        else None
    )
    sensitivity = float(tp / positives) if positives else None
    specificity = float(tn / negatives) if negatives else None
    return {
        "n": int(outcome.size),
        "positive": positives,
        "negative": negatives,
        "auroc": auroc,
        "threshold": float(threshold),
        "sensitivity": sensitivity,
        "sensitivity_wilson_95": wilson_interval(int(tp), positives),
        "specificity": specificity,
        "specificity_wilson_95": wilson_interval(int(tn), negatives),
        "tp": int(tp),
        "fn": int(fn),
        "tn": int(tn),
        "fp": int(fp),
    }


def stratified_bootstrap(
    outcome: np.ndarray,
    site: np.ndarray,
    sparse_scores: np.ndarray,
    comparator_scores: np.ndarray,
) -> dict:
    strata = []
    for site_name in ("Valencia", "Madrid"):
        for class_value in (0, 1):
            indices = np.flatnonzero(
                (site == site_name) & (outcome == class_value)
            )
            if indices.size:
                strata.append(indices)
    if sum(indices.size for indices in strata) != outcome.size:
        raise RuntimeError("bootstrap strata do not partition Spain")

    stratum_sizes = {
        f"{site_name}|{class_value}": int(indices.size)
        for site_name in ("Valencia", "Madrid")
        for class_value in (0, 1)
        for indices in [
            np.flatnonzero((site == site_name) & (outcome == class_value))
        ]
    }
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    sparse_auc = np.empty(BOOTSTRAP_REPLICATES, dtype=np.float64)
    comparator_auc = np.empty(BOOTSTRAP_REPLICATES, dtype=np.float64)
    difference = np.empty(BOOTSTRAP_REPLICATES, dtype=np.float64)
    for iteration in range(BOOTSTRAP_REPLICATES):
        sampled = np.concatenate(
            [
                rng.choice(indices, size=indices.size, replace=True)
                for indices in strata
            ]
        )
        sparse_value = roc_auc_score(outcome[sampled], sparse_scores[sampled])
        comparator_value = roc_auc_score(
            outcome[sampled], comparator_scores[sampled]
        )
        sparse_auc[iteration] = sparse_value
        comparator_auc[iteration] = comparator_value
        difference[iteration] = sparse_value - comparator_value

    return {
        "method": "percentile bootstrap within site x class",
        "replicates": BOOTSTRAP_REPLICATES,
        "seed": BOOTSTRAP_SEED,
        "paired_sparse_comparator_indices": True,
        "site_class_stratum_sizes": stratum_sizes,
        "sparse_auroc_ci_95": [
            float(value) for value in np.quantile(sparse_auc, [0.025, 0.975])
        ],
        "all_window_auroc_ci_95": [
            float(value) for value in np.quantile(comparator_auc, [0.025, 0.975])
        ],
        "sparse_minus_all_window_auroc": float(
            roc_auc_score(outcome, sparse_scores)
            - roc_auc_score(outcome, comparator_scores)
        ),
        "sparse_minus_all_window_auroc_ci_95": [
            float(value) for value in np.quantile(difference, [0.025, 0.975])
        ],
    }


def within_site_permutation(
    outcome: np.ndarray,
    site: np.ndarray,
    scores: np.ndarray,
) -> dict:
    observed = float(roc_auc_score(outcome, scores))
    rng = np.random.default_rng(PERMUTATION_SEED)
    exceed = 0
    for _ in range(PERMUTATIONS):
        permuted = outcome.copy()
        for site_name in ("Valencia", "Madrid"):
            indices = np.flatnonzero(site == site_name)
            permuted[indices] = rng.permutation(permuted[indices])
        permuted_auc = float(roc_auc_score(permuted, scores))
        if permuted_auc >= observed:
            exceed += 1
    return {
        "statistic": "AUROC",
        "alternative": "greater",
        "shuffle": "labels within site only",
        "permutations": PERMUTATIONS,
        "seed": PERMUTATION_SEED,
        "observed": observed,
        "exceed": exceed,
        "p_value": float((1 + exceed) / (PERMUTATIONS + 1)),
    }


def exploratory_calibration(
    outcome: np.ndarray,
    scores: np.ndarray,
) -> dict:
    clipped = np.clip(scores, 1e-12, 1.0 - 1e-12)
    logit = np.log(clipped / (1.0 - clipped))

    def objective(parameters: np.ndarray) -> float:
        linear = parameters[0] + parameters[1] * logit
        loss = np.logaddexp(0.0, linear) - outcome * linear
        return float(np.sum(loss))

    fitted = minimize(
        objective,
        x0=np.asarray([0.0, 1.0]),
        method="BFGS",
    )
    if not fitted.success or not np.all(np.isfinite(fitted.x)):
        return {
            "status": "undefined_or_failed",
            "intercept": None,
            "slope": None,
            "exploratory": True,
            "reason": str(fitted.message),
        }
    return {
        "status": "estimated",
        "intercept": float(fitted.x[0]),
        "slope": float(fitted.x[1]),
        "exploratory": True,
        "warning": (
            "balanced-class-weight logistic scores are not calibrated "
            "probabilities"
        ),
    }


def robustness_scores(
    spectra: np.ndarray,
    wavenumbers: np.ndarray,
    selected_centers: Sequence[int],
    model: dict,
) -> dict:
    output = {}
    for shift in ROBUST_CENTER_SHIFTS:
        shifted = [int(center + shift) for center in selected_centers]
        features = integrate_windows(
            spectra, wavenumbers, shifted, WINDOW_WIDTH
        )
        matrix = model_matrix(features, tuple(range(len(shifted))))
        output[f"common_shift_{shift:+d}_cm_inverse"] = {
            "centers_cm_inverse": shifted,
            "width_cm_inverse": WINDOW_WIDTH,
            "scores": [
                float(value)
                for value in score_serialized_model(model, matrix)
            ],
            "interpretation": "global spectral-axis tolerance only",
        }
    for width in ROBUST_WIDTHS:
        features = integrate_windows(
            spectra, wavenumbers, selected_centers, width
        )
        matrix = model_matrix(
            features, tuple(range(len(selected_centers)))
        )
        output[f"width_{int(width)}_cm_inverse"] = {
            "centers_cm_inverse": list(selected_centers),
            "width_cm_inverse": width,
            "scores": [
                float(value)
                for value in score_serialized_model(model, matrix)
            ],
        }
    return output


def validate_no_nonfinite(value: object, path: str = "$") -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            validate_no_nonfinite(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            validate_no_nonfinite(item, f"{path}[{index}]")
    elif isinstance(value, float) and not math.isfinite(value):
        raise RuntimeError(f"nonfinite validation value at {path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("frozen_discovery", type=Path)
    parser.add_argument("frozen_discovery_sha256", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--reviewed-discovery-sha256", required=True)
    parser.add_argument(
        "--unlock-spanish-confirmation",
        action="store_true",
        help="Required explicit one-shot unblinding authorization.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.unlock_spanish_confirmation:
        raise SystemExit(
            "Spanish confirmation remains locked; explicit unlock flag required"
        )
    if args.output_dir.exists():
        raise SystemExit(
            "output directory already exists; refusing before Spanish access"
        )
    if REVIEWED_DISCOVERY_SHA256 is None:
        raise SystemExit(
            "Spanish confirmation remains locked; literal reviewed SHA is unset"
        )
    args.reviewed_discovery_sha256 = (
        args.reviewed_discovery_sha256.strip().lower()
    )
    discovery, discovery_digest = load_and_verify_discovery(
        args.frozen_discovery,
        args.frozen_discovery_sha256,
        args.reviewed_discovery_sha256,
    )

    # This is the first and only point at which the ESP variable is requested.
    spain = load_spain_after_unlock(args.archive)
    spectra = spain["spectra"]
    wavenumbers = spain["wavenumbers"]
    outcome = spain["outcome"]
    site = spain["site"]

    selected = discovery["discovery"]["selected_design"]
    selection = selected["selection"]
    selected_indices = tuple(int(value) for value in selection["indices"])
    selected_centers = [int(CENTERS[index]) for index in selected_indices]
    if selected_centers != selection["centers_cm_inverse"]:
        raise RuntimeError("selected indices/centers mismatch")

    all_features = integrate_windows(
        spectra, wavenumbers, CENTERS, WINDOW_WIDTH
    )
    sparse_matrix = model_matrix(all_features, selected_indices)
    all_matrix = model_matrix(
        all_features, tuple(range(len(CENTERS)))
    )
    sparse_scores = score_serialized_model(
        selected["fitted_model"], sparse_matrix
    )
    all_window = discovery["discovery"]["all_window_comparator"]
    all_scores = score_serialized_model(
        all_window["fitted_model"], all_matrix
    )
    sparse_threshold = float(selected["threshold"])
    all_threshold = float(all_window["threshold"])

    sparse_metrics = classification_metrics(
        outcome, sparse_scores, sparse_threshold
    )
    all_metrics = classification_metrics(
        outcome, all_scores, all_threshold
    )
    bootstrap = stratified_bootstrap(
        outcome, site, sparse_scores, all_scores
    )
    transfer_gates = {
        "auroc": {
            "point": sparse_metrics["auroc"],
            "minimum": 0.90,
            "pass": bool(sparse_metrics["auroc"] >= 0.90),
        },
        "sensitivity": {
            "point": sparse_metrics["sensitivity"],
            "minimum": 0.85,
            "pass": bool(sparse_metrics["sensitivity"] >= 0.85),
        },
        "specificity": {
            "point": sparse_metrics["specificity"],
            "minimum": 0.75,
            "pass": bool(sparse_metrics["specificity"] >= 0.75),
        },
        "sparse_minus_all_window_auroc": {
            "point": bootstrap["sparse_minus_all_window_auroc"],
            "minimum": -0.10,
            "pass": bool(
                bootstrap["sparse_minus_all_window_auroc"] >= -0.10
            ),
        },
    }
    transfer_gates["overall_pass"] = bool(
        all(item["pass"] for item in transfer_gates.values())
    )
    permutation = within_site_permutation(
        outcome, site, sparse_scores
    )

    site_metrics = {}
    for site_name in ("Valencia", "Madrid"):
        mask = site == site_name
        site_metrics[site_name] = classification_metrics(
            outcome[mask], sparse_scores[mask], sparse_threshold
        )

    robustness = robustness_scores(
        spectra,
        wavenumbers,
        selected_centers,
        selected["fitted_model"],
    )
    for item in robustness.values():
        scores = np.asarray(item.pop("scores"), dtype=np.float64)
        item["metrics_at_nominal_frozen_threshold"] = (
            classification_metrics(outcome, scores, sparse_threshold)
        )
        item["score_vector"] = [float(value) for value in scores]
        item["refit"] = False
        item["recalibration"] = False
        item["reselection"] = False

    controls = {}
    control_specs = {
        "atmospheric_instrument_artifact": ATMOSPHERIC_CONTROL_CENTERS,
        "clean_off_target": CLEAN_OFF_TARGET_CENTERS,
    }
    for name, centers in control_specs.items():
        frozen = discovery["negative_controls"][name]
        if list(centers) != frozen["centers_cm_inverse"]:
            raise RuntimeError(f"negative-control centers mismatch: {name}")
        features = integrate_windows(
            spectra, wavenumbers, centers, WINDOW_WIDTH
        )
        matrix = model_matrix(features, tuple(range(len(centers))))
        scores = score_serialized_model(frozen["fitted_model"], matrix)
        controls[name] = {
            "centers_cm_inverse": list(centers),
            "threshold": float(frozen["metrics"]["threshold"]),
            "metrics": classification_metrics(
                outcome, scores, float(frozen["metrics"]["threshold"])
            ),
            "scores": [float(value) for value in scores],
            "refit": False,
        }

    predictions = [
        {
            "subject_id": str(identifier),
            "site": str(site_name),
            "threshold_class": int(label),
            "sparse_score": float(sparse_score),
            "all_window_score": float(all_score),
            "reference_zero_is_left_censored": bool(reference == 0.0),
        }
        for identifier, site_name, label, sparse_score, all_score, reference
        in zip(
            spain["subject_ids"],
            site,
            outcome,
            sparse_scores,
            all_scores,
            spain["albumin"],
            strict=True,
        )
    ]

    result = {
        "schema_version": "globaldkd-sparse-window-validation-v1",
        "phase": "ONE-SHOT SPANISH CONFIRMATION",
        "title": discovery["title"],
        "scope": discovery["scope"],
        "linkage": {
            "frozen_discovery_sha256": discovery_digest,
            "reviewed_discovery_sha256": args.reviewed_discovery_sha256.lower(),
            "selected_centers_cm_inverse": selected_centers,
            "selected_indices": list(selected_indices),
            "sparse_threshold_from_Australia": sparse_threshold,
            "all_window_threshold_from_Australia": all_threshold,
            "no_refit_or_recalibration": True,
        },
        "input": {
            "archive_md5": md5sum(args.archive),
            "mat_member": MAT_MEMBER,
            "mat_member_sha256": spain["mat_sha256"],
            "loadmat_variable_names": list(REQUESTED_VARIABLES),
            "Global_not_loaded": True,
            "ESP_data_already_contains_source_0_8_correction": True,
        },
        "cohort": {
            "n": 61,
            "positive": int(outcome.sum()),
            "negative": int((1 - outcome).sum()),
            "site_counts": {
                "Valencia": int(np.sum(site == "Valencia")),
                "Madrid": int(np.sum(site == "Madrid")),
            },
            "reference_interpretation": (
                "below-LOQ or source NaN-to-zero entries are left-censored "
                "threshold references, not exact concentrations"
            ),
        },
        "primary_sparse": {
            "metrics": sparse_metrics,
            "external_success_gates": transfer_gates,
            "bootstrap": bootstrap,
            "within_site_permutation": permutation,
            "calibration": exploratory_calibration(outcome, sparse_scores),
        },
        "all_window_comparator": {
            "metrics": all_metrics,
        },
        "site_metrics": site_metrics,
        "robustness_nominal_model_only": robustness,
        "negative_controls": controls,
        "predictions": predictions,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "assertions": {
            "reviewed_discovery_digest_exact": True,
            "discovery_status_DESIGN": True,
            "ESP_only": True,
            "shape_61_by_1598": True,
            "axis_grid_3994_to_800_step_minus_2": True,
            "axis_label_slots_match_audited_schema": True,
            "unique_spanish_ids": True,
            "sites_derived_from_audited_ID_prefixes": True,
            "audited_Valencia_then_Madrid_block_order": True,
            "counts_34_27": True,
            "sites_46_15": True,
            "frozen_scaler_model_threshold": True,
            "paired_bootstrap_indices": True,
            "within_site_permutation": True,
            "no_refit_recalibration_reselection": True,
            "all_outputs_finite_or_explicitly_null": True,
        },
    }
    validate_no_nonfinite(result)
    payload = canonical_bytes(result)
    digest = hashlib.sha256(payload).hexdigest()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    result_path = args.output_dir / "frozen-validation.json"
    result_path.write_bytes(payload)
    (args.output_dir / "frozen-validation.sha256").write_text(
        f"{digest}  frozen-validation.json\n",
        encoding="utf-8",
    )
    summary = {
        "phase": result["phase"],
        "frozen_discovery_sha256": discovery_digest,
        "frozen_validation_sha256": digest,
        "selected_centers_cm_inverse": selected_centers,
        "primary_sparse": result["primary_sparse"],
        "all_window_comparator": result["all_window_comparator"],
        "site_metrics": result["site_metrics"],
        "negative_controls": {
            key: value["metrics"] for key, value in controls.items()
        },
        "robustness": {
            key: value["metrics_at_nominal_frozen_threshold"]
            for key, value in robustness.items()
        },
    }
    (args.output_dir / "validation-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
