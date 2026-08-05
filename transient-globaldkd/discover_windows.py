#!/usr/bin/env python3
"""Australian-only exact sparse-window discovery for GlobalDKD.

Scientific scope: fixed computational mean-absorbance window summaries for
cross-instrument urine albumin-concentration screening. This is not uACR, DKD
diagnosis, a physical-filter measurement, hardware simulation, or a universal
minimum claim. The program requests only AUS_Int from the root MATLAB file and
must not deserialize Spanish outcome-bearing variables.
"""
from __future__ import annotations

import hashlib
import io
import itertools
import json
import math
import os
import platform
import sys
import warnings
import zipfile
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import scipy
import sklearn
from scipy.io import loadmat
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold

TITLE = (
    "How Many Spectral Windows Are Enough? An Exact Sparse-Window Certificate "
    "for Cross-Instrument Screening of Elevated Urinary Albumin Concentration."
)
EXPECTED_MD5 = "8d5d4ba0e0578d0e39e76e1f1514a841"
MAT_MEMBER = "Zenodo/dataset.mat"
REQUESTED_VARIABLES = ("AUS_Int",)
FORBIDDEN_VARIABLES = ("ESP", "Global")
CUTOFF_MG_L = 30.0
CENTERS = (
    3616, 3552, 3488, 3424, 3360, 3296, 3232, 3168, 3104, 3040, 2976, 2912,
    1752, 1688, 1624, 1560, 1496, 1432, 1368, 1304, 1240, 1176, 1112, 1048,
)
ATMOSPHERIC_CONTROL_CENTERS = (2368, 2304, 2240, 2176)
CLEAN_OFF_TARGET_CENTERS = (2496, 2560, 2624, 2688)
WINDOW_WIDTH = 64.0
ZERO_EPSILON = 1e-12
NEAR_CORRELATION = 0.999999
NEAR_NORMALIZED_RMS = 0.0015
PRIMARY_SEED = 20260805
OUTER_SEED = 20260806
INNER_SEEDS = (20260807, 20260808, 20260809, 20260810, 20260811)
MIN_SENSITIVITY = 0.90
MIN_SPECIFICITY = 0.80
MIN_AUROC = 0.90
MAX_COMPARATOR_GAP = 0.03
EXPECTED_CANDIDATES = 12950
MODEL_SPEC = {
    "C": 1.0,
    "penalty": "l2",
    "class_weight": "balanced",
    "solver": "liblinear",
    "random_state": PRIMARY_SEED,
    "max_iter": 5000,
}


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_bytes(value: object) -> bytes:
    text = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )
    return (text + "\n").encode("utf-8")


def finite_float(value: float) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"nonfinite result encountered: {result!r}")
    return result


def load_australia_only(archive: Path) -> dict:
    observed_md5 = md5sum(archive)
    if observed_md5 != EXPECTED_MD5:
        raise RuntimeError(f"archive MD5 mismatch: {observed_md5}")

    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
        if MAT_MEMBER not in names:
            raise RuntimeError(f"missing archive member: {MAT_MEMBER}")
        mat_payload = bundle.read(MAT_MEMBER)

    loaded = loadmat(
        io.BytesIO(mat_payload),
        simplify_cells=True,
        variable_names=REQUESTED_VARIABLES,
    )
    public_keys = sorted(key for key in loaded if not key.startswith("__"))
    if public_keys != list(REQUESTED_VARIABLES):
        raise RuntimeError(f"unexpected loaded variables: {public_keys}")
    if any(name in loaded for name in FORBIDDEN_VARIABLES):
        raise RuntimeError("forbidden Spanish-bearing variable was deserialized")

    aus = loaded["AUS_Int"]
    if not isinstance(aus, dict):
        raise RuntimeError(f"AUS_Int did not simplify to dict: {type(aus).__name__}")

    spectra = np.asarray(aus["data"], dtype=np.float64)
    axes = np.asarray(aus["axisscale"], dtype=object).reshape(-1)
    labels = np.asarray(aus["label"], dtype=object).reshape(-1)

    if len(axes) < 7 or len(labels) < 2:
        raise RuntimeError("unexpected PLS dataset cell layout")

    albumin = np.asarray(axes[0], dtype=np.float64).reshape(-1)
    duplicate_albumin = np.asarray(axes[1], dtype=np.float64).reshape(-1)
    wavenumbers = np.asarray(axes[4], dtype=np.float64).reshape(-1)
    subject_ids = np.char.strip(
        np.asarray(labels[0], dtype=str).reshape(-1)
    ).astype(str)

    if spectra.shape != (177, 1598):
        raise RuntimeError(f"unexpected AUS_Int.data shape: {spectra.shape}")
    if albumin.shape != (177,) or wavenumbers.shape != (1598,):
        raise RuntimeError(
            f"unexpected axis shapes: albumin={albumin.shape}, wn={wavenumbers.shape}"
        )
    if subject_ids.shape != (177,):
        raise RuntimeError(f"unexpected subject-ID shape: {subject_ids.shape}")
    if len(set(subject_ids.tolist())) != 177:
        duplicates = sorted(
            key for key, count in Counter(subject_ids.tolist()).items() if count > 1
        )
        raise RuntimeError(f"nonunique Australian IDs: {duplicates}")
    if any(not value for value in subject_ids):
        raise RuntimeError("blank Australian subject ID")
    if not np.array_equal(albumin, duplicate_albumin, equal_nan=True):
        raise RuntimeError("duplicate MATLAB albumin axes disagree")
    if str(axes[2]) != "Albumin" or str(axes[6]) != "Wavenumbers (cm-1)":
        raise RuntimeError(
            f"unexpected axis labels: {axes[2]!r}, {axes[6]!r}"
        )
    if not np.all(np.isfinite(spectra)) or not np.all(np.isfinite(wavenumbers)):
        raise RuntimeError("nonfinite spectrum or wavenumber")
    if not (np.all(np.diff(wavenumbers) < 0) or np.all(np.diff(wavenumbers) > 0)):
        raise RuntimeError("wavenumber axis is not strictly monotone")

    finite = np.isfinite(albumin)
    missing_indices = np.flatnonzero(~finite)
    if not np.array_equal(missing_indices, np.arange(22)):
        raise RuntimeError(
            f"expected first 22 rows to be missing-albumin controls: "
            f"{missing_indices.tolist()}"
        )
    primary_labels = (albumin[finite] >= CUTOFF_MG_L).astype(np.uint8)
    if int(finite.sum()) != 155:
        raise RuntimeError(f"expected 155 measured rows, found {int(finite.sum())}")
    if (int(primary_labels.sum()), int((1 - primary_labels).sum())) != (64, 91):
        raise RuntimeError(
            "Australian threshold count mismatch: "
            f"{int(primary_labels.sum())}+/{int((1-primary_labels).sum())}-"
        )

    return {
        "archive_md5": observed_md5,
        "mat_sha256": sha256_bytes(mat_payload),
        "spectra_all": spectra,
        "wavenumbers": wavenumbers,
        "albumin": albumin,
        "subject_ids_all": subject_ids,
        "finite_mask": finite,
        "labels_primary": primary_labels,
    }


def pair_hash(first: str, second: str) -> str:
    left, right = sorted((first, second))
    return sha256_bytes((left + "\0" + right).encode("utf-8"))


def duplicate_audit(
    spectra: np.ndarray,
    subject_ids: np.ndarray,
    measured_mask: np.ndarray,
    measured_labels: np.ndarray,
) -> dict:
    centered = spectra - spectra.mean(axis=1, keepdims=True)
    row_rms = np.sqrt(np.mean(centered * centered, axis=1))
    if np.any(~np.isfinite(row_rms)) or np.any(row_rms <= ZERO_EPSILON):
        raise RuntimeError("zero/near-zero row RMS in duplicate audit")
    standardized = centered / row_rms[:, None]
    correlations = standardized @ standardized.T / standardized.shape[1]
    correlations = np.clip(correlations, -1.0, 1.0)

    label_by_row: list[int | None] = [None] * spectra.shape[0]
    measured_rows = np.flatnonzero(measured_mask)
    for row, label in zip(measured_rows, measured_labels, strict=True):
        label_by_row[int(row)] = int(label)

    exact_pair_hashes: list[str] = []
    near_pair_hashes: list[str] = []
    cross_label_exact: list[str] = []
    for first in range(spectra.shape[0]):
        for second in range(first + 1, spectra.shape[0]):
            identifier = pair_hash(str(subject_ids[first]), str(subject_ids[second]))
            exact = np.array_equal(spectra[first], spectra[second])
            corr = float(correlations[first, second])
            rms_distance = math.sqrt(max(0.0, 2.0 - 2.0 * corr))
            near = (
                corr >= NEAR_CORRELATION
                and rms_distance <= NEAR_NORMALIZED_RMS
            )
            if exact:
                exact_pair_hashes.append(identifier)
                if (
                    label_by_row[first] is not None
                    and label_by_row[second] is not None
                    and label_by_row[first] != label_by_row[second]
                ):
                    cross_label_exact.append(identifier)
            elif near:
                near_pair_hashes.append(identifier)

    if cross_label_exact:
        raise RuntimeError(
            "cross-label exact duplicate spectra detected: "
            + ",".join(cross_label_exact)
        )
    return {
        "definition": {
            "exact": "element-for-element float64 equality",
            "near": {
                "row_transform": "subtract row mean, divide by row RMS",
                "pearson_correlation_at_least": NEAR_CORRELATION,
                "normalized_rms_distance_at_most": NEAR_NORMALIZED_RMS,
            },
        },
        "exact_pair_count": len(exact_pair_hashes),
        "exact_pair_hashes": sorted(exact_pair_hashes),
        "near_pair_count_excluding_exact": len(near_pair_hashes),
        "near_pair_hashes_excluding_exact": sorted(near_pair_hashes),
        "cross_label_exact_pair_count": 0,
    }


def integrate_windows(
    spectra: np.ndarray,
    wavenumbers: np.ndarray,
    centers: Sequence[int],
    width: float = WINDOW_WIDTH,
) -> np.ndarray:
    if width <= 0:
        raise RuntimeError(f"nonpositive window width: {width}")
    order = np.argsort(wavenumbers)
    axis = np.asarray(wavenumbers[order], dtype=np.float64)
    values = np.asarray(spectra[:, order], dtype=np.float64)
    output = np.empty((spectra.shape[0], len(centers)), dtype=np.float64)

    for column, center in enumerate(centers):
        lower = float(center) - width / 2.0
        upper = float(center) + width / 2.0
        if lower < axis[0] or upper > axis[-1]:
            raise RuntimeError(
                f"window [{lower}, {upper}] outside [{axis[0]}, {axis[-1]}]"
            )
        interior = (axis > lower) & (axis < upper)
        segment_axis = np.concatenate(([lower], axis[interior], [upper]))
        lower_values = np.asarray(
            [np.interp(lower, axis, row) for row in values],
            dtype=np.float64,
        )
        upper_values = np.asarray(
            [np.interp(upper, axis, row) for row in values],
            dtype=np.float64,
        )
        segment_values = np.column_stack(
            (lower_values, values[:, interior], upper_values)
        )
        output[:, column] = (
            np.trapezoid(segment_values, segment_axis, axis=1) / width
        )

    if not np.all(np.isfinite(output)):
        raise RuntimeError("nonfinite integrated window summary")
    return output


def model_matrix(features: np.ndarray, indices: Sequence[int]) -> np.ndarray:
    raw = np.asarray(features[:, indices], dtype=np.float64)
    if raw.ndim == 1:
        raw = raw[:, None]
    norms = np.linalg.norm(raw, axis=1)
    if np.any(~np.isfinite(norms)) or np.any(norms <= ZERO_EPSILON):
        raise RuntimeError("zero/near-zero L2 norm in model feature vector")
    normalized = raw / norms[:, None]
    matrix = np.concatenate((raw, normalized), axis=1)
    if not np.all(np.isfinite(matrix)):
        raise RuntimeError("nonfinite model matrix")
    return matrix


def fit_model(matrix: np.ndarray, labels: np.ndarray) -> dict:
    if not np.all(np.isfinite(matrix)):
        raise RuntimeError("nonfinite training matrix")
    mean = matrix.mean(axis=0)
    variance = np.mean((matrix - mean) ** 2, axis=0)
    scale = np.sqrt(variance)
    scale[scale == 0.0] = 1.0
    standardized = (matrix - mean) / scale
    estimator = LogisticRegression(**MODEL_SPEC)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        estimator.fit(standardized, labels)
    if list(estimator.classes_) != [0, 1]:
        raise RuntimeError(f"unexpected estimator classes: {estimator.classes_}")
    return {
        "mean": mean,
        "scale": scale,
        "coefficient": estimator.coef_.reshape(-1),
        "intercept": float(estimator.intercept_[0]),
        "estimator": estimator,
    }


def predict_model(state: dict, matrix: np.ndarray) -> np.ndarray:
    standardized = (matrix - state["mean"]) / state["scale"]
    probabilities = state["estimator"].predict_proba(standardized)[:, 1]
    if not np.all(np.isfinite(probabilities)):
        raise RuntimeError("nonfinite predicted score")
    return probabilities


def serializable_model(state: dict) -> dict:
    return {
        "feature_order": "raw selected means, then row-L2-normalized selected means",
        "scaler_mean": [finite_float(value) for value in state["mean"]],
        "scaler_scale": [finite_float(value) for value in state["scale"]],
        "coefficient": [finite_float(value) for value in state["coefficient"]],
        "intercept": finite_float(state["intercept"]),
    }


def crossfit_scores(
    matrix: np.ndarray,
    labels: np.ndarray,
    splits: Sequence[tuple[np.ndarray, np.ndarray]],
    expected_multiplicity: int,
) -> np.ndarray:
    totals = np.zeros(labels.size, dtype=np.float64)
    counts = np.zeros(labels.size, dtype=np.int64)
    for train, test in splits:
        state = fit_model(matrix[train], labels[train])
        totals[test] += predict_model(state, matrix[test])
        counts[test] += 1
    if not np.all(counts == expected_multiplicity):
        raise RuntimeError(
            f"OOF multiplicity mismatch: {Counter(counts.tolist())}"
        )
    scores = totals / counts
    if not np.all(np.isfinite(scores)):
        raise RuntimeError("nonfinite aggregated OOF score")
    return scores


def constrained_metrics(labels: np.ndarray, scores: np.ndarray) -> dict:
    positives = np.sort(scores[labels == 1])[::-1]
    required = int(math.ceil(MIN_SENSITIVITY * positives.size))
    if required < 1:
        raise RuntimeError("no positives available for threshold")
    threshold = float(positives[required - 1])
    predicted = (scores >= threshold).astype(np.uint8)
    tn, fp, fn, tp = confusion_matrix(
        labels, predicted, labels=[0, 1]
    ).ravel()
    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)
    return {
        "auroc": finite_float(roc_auc_score(labels, scores)),
        "threshold": finite_float(threshold),
        "sensitivity": finite_float(sensitivity),
        "specificity": finite_float(specificity),
        "tp": int(tp),
        "fn": int(fn),
        "tn": int(tn),
        "fp": int(fp),
        "threshold_rule": (
            "highest aggregated OOF score retaining sensitivity>=0.90; "
            "prediction positive iff score>=threshold"
        ),
        "required_positive_rank": required,
    }


def fixed_threshold_metrics(
    labels: np.ndarray,
    scores: np.ndarray,
    threshold: float,
) -> dict:
    predicted = (scores >= threshold).astype(np.uint8)
    tn, fp, fn, tp = confusion_matrix(
        labels, predicted, labels=[0, 1]
    ).ravel()
    return {
        "auroc": finite_float(roc_auc_score(labels, scores)),
        "threshold": finite_float(threshold),
        "sensitivity": finite_float(tp / (tp + fn)),
        "specificity": finite_float(tn / (tn + fp)),
        "tp": int(tp),
        "fn": int(fn),
        "tn": int(tn),
        "fp": int(fp),
    }


def candidate_library() -> list[tuple[int, ...]]:
    candidates = [
        combination
        for size in range(1, 5)
        for combination in itertools.combinations(range(len(CENTERS)), size)
    ]
    if len(candidates) != EXPECTED_CANDIDATES:
        raise RuntimeError(f"candidate count mismatch: {len(candidates)}")
    if len(set(candidates)) != EXPECTED_CANDIDATES:
        raise RuntimeError("candidate library is not unique")
    return candidates


def candidate_record(
    indices: tuple[int, ...],
    metrics: dict,
    comparator_auroc: float,
) -> dict:
    qualifies = (
        metrics["auroc"] >= MIN_AUROC
        and metrics["specificity"] >= MIN_SPECIFICITY
        and metrics["sensitivity"] >= MIN_SENSITIVITY
        and metrics["auroc"] >= comparator_auroc - MAX_COMPARATOR_GAP
    )
    return {
        "indices": list(indices),
        "centers_cm_inverse": [CENTERS[index] for index in indices],
        "n_windows": len(indices),
        "auroc": metrics["auroc"],
        "threshold": metrics["threshold"],
        "sensitivity": metrics["sensitivity"],
        "specificity": metrics["specificity"],
        "qualifies": bool(qualifies),
    }


def exhaustive_search(
    features: np.ndarray,
    labels: np.ndarray,
    splits: Sequence[tuple[np.ndarray, np.ndarray]],
    multiplicity: int,
    candidates: Sequence[tuple[int, ...]],
    retain_all: bool,
    progress_prefix: str,
) -> dict:
    comparator_matrix = model_matrix(features, tuple(range(len(CENTERS))))
    comparator_scores = crossfit_scores(
        comparator_matrix, labels, splits, multiplicity
    )
    comparator_metrics = constrained_metrics(labels, comparator_scores)

    selected: dict | None = None
    all_records: list[dict] = []
    qualifying_by_size = Counter()
    for position, indices in enumerate(candidates, start=1):
        matrix = model_matrix(features, indices)
        scores = crossfit_scores(matrix, labels, splits, multiplicity)
        metrics = constrained_metrics(labels, scores)
        record = candidate_record(indices, metrics, comparator_metrics["auroc"])
        if retain_all:
            all_records.append(record)
        if record["qualifies"]:
            qualifying_by_size[len(indices)] += 1
            if selected is None:
                selected = {**record, "oof_scores": [
                    finite_float(value) for value in scores
                ]}
            else:
                old_key = (
                    selected["n_windows"],
                    -selected["auroc"],
                    tuple(selected["indices"]),
                )
                new_key = (
                    record["n_windows"],
                    -record["auroc"],
                    tuple(record["indices"]),
                )
                if new_key < old_key:
                    selected = {**record, "oof_scores": [
                        finite_float(value) for value in scores
                    ]}
        if position % 1000 == 0 or position == len(candidates):
            print(
                f"{progress_prefix}: evaluated {position}/{len(candidates)}",
                file=sys.stderr,
                flush=True,
            )

    return {
        "status": "DESIGN" if selected is not None else "NO DESIGN",
        "selected": selected,
        "qualifying_count": int(sum(qualifying_by_size.values())),
        "qualifying_by_size": {
            str(size): int(qualifying_by_size.get(size, 0))
            for size in range(1, 5)
        },
        "comparator": {
            **comparator_metrics,
            "oof_scores": [finite_float(value) for value in comparator_scores],
        },
        "candidate_metrics": all_records if retain_all else None,
    }


def split_manifest(
    name: str,
    splits: Sequence[tuple[np.ndarray, np.ndarray]],
    subject_ids: np.ndarray,
    folds_per_repeat: int,
) -> list[dict]:
    records = []
    for index, (train, test) in enumerate(splits):
        records.append(
            {
                "name": name,
                "split_index": index,
                "repeat_index": index // folds_per_repeat,
                "fold_index": index % folds_per_repeat,
                "train_ids": [str(subject_ids[item]) for item in train],
                "test_ids": [str(subject_ids[item]) for item in test],
            }
        )
    return records


def materialize_primary_splits(
    labels: np.ndarray,
) -> list[tuple[np.ndarray, np.ndarray]]:
    splitter = RepeatedStratifiedKFold(
        n_splits=5,
        n_repeats=5,
        random_state=PRIMARY_SEED,
    )
    return [
        (train.astype(int), test.astype(int))
        for train, test in splitter.split(np.zeros(labels.size), labels)
    ]


def nested_audit(
    features: np.ndarray,
    labels: np.ndarray,
    subject_ids: np.ndarray,
    candidates: Sequence[tuple[int, ...]],
) -> tuple[dict, list[dict]]:
    outer = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=OUTER_SEED,
    )
    outer_splits = [
        (train.astype(int), test.astype(int))
        for train, test in outer.split(np.zeros(labels.size), labels)
    ]
    manifest = split_manifest("nested_outer", outer_splits, subject_ids, 5)

    sparse_scores = np.zeros(labels.size, dtype=np.float64)
    sparse_predicted = np.zeros(labels.size, dtype=np.uint8)
    comparator_scores = np.zeros(labels.size, dtype=np.float64)
    comparator_predicted = np.zeros(labels.size, dtype=np.uint8)
    sparse_seen = np.zeros(labels.size, dtype=np.uint8)
    comparator_seen = np.zeros(labels.size, dtype=np.uint8)
    fold_records = []
    selection_counts = Counter()

    for outer_index, (outer_train, outer_test) in enumerate(outer_splits):
        inner_seed = INNER_SEEDS[outer_index]
        inner_splitter = StratifiedKFold(
            n_splits=5,
            shuffle=True,
            random_state=inner_seed,
        )
        inner_local = [
            (train.astype(int), test.astype(int))
            for train, test in inner_splitter.split(
                np.zeros(outer_train.size), labels[outer_train]
            )
        ]
        inner_global = [
            (outer_train[train], outer_train[test])
            for train, test in inner_local
        ]
        manifest.extend(
            split_manifest(
                f"nested_inner_outer_{outer_index}",
                inner_global,
                subject_ids,
                5,
            )
        )

        search = exhaustive_search(
            features[outer_train],
            labels[outer_train],
            inner_local,
            1,
            candidates,
            retain_all=False,
            progress_prefix=f"outer {outer_index} inner",
        )

        all_indices = tuple(range(len(CENTERS)))
        comparator_inner_matrix = model_matrix(features[outer_train], all_indices)
        comparator_inner_scores = crossfit_scores(
            comparator_inner_matrix, labels[outer_train], inner_local, 1
        )
        comparator_inner_metrics = constrained_metrics(
            labels[outer_train], comparator_inner_scores
        )
        comparator_state = fit_model(
            comparator_inner_matrix, labels[outer_train]
        )
        comparator_test_matrix = model_matrix(features[outer_test], all_indices)
        comparator_fold_scores = predict_model(
            comparator_state, comparator_test_matrix
        )
        comparator_scores[outer_test] = comparator_fold_scores
        comparator_predicted[outer_test] = (
            comparator_fold_scores >= comparator_inner_metrics["threshold"]
        ).astype(np.uint8)
        comparator_seen[outer_test] += 1

        fold_record = {
            "outer_fold": outer_index,
            "inner_seed": inner_seed,
            "status": search["status"],
            "qualifying_count": search["qualifying_count"],
            "qualifying_by_size": search["qualifying_by_size"],
            "selected": None,
            "outer_test_ids": [
                str(subject_ids[index]) for index in outer_test
            ],
            "all_window_inner_threshold": comparator_inner_metrics["threshold"],
        }

        if search["selected"] is not None:
            selected = search["selected"]
            indices = tuple(selected["indices"])
            selection_key = ",".join(str(CENTERS[index]) for index in indices)
            selection_counts[selection_key] += 1
            train_matrix = model_matrix(features[outer_train], indices)
            test_matrix = model_matrix(features[outer_test], indices)
            state = fit_model(train_matrix, labels[outer_train])
            fold_scores = predict_model(state, test_matrix)
            sparse_scores[outer_test] = fold_scores
            sparse_predicted[outer_test] = (
                fold_scores >= selected["threshold"]
            ).astype(np.uint8)
            sparse_seen[outer_test] += 1
            fold_record["selected"] = {
                key: value
                for key, value in selected.items()
                if key != "oof_scores"
            }
        fold_records.append(fold_record)

    if not np.all(comparator_seen == 1):
        raise RuntimeError("nested comparator outer multiplicity mismatch")
    comparator_metrics = fixed_decision_metrics(
        labels, comparator_scores, comparator_predicted
    )
    successful = int(np.sum(sparse_seen == 1))
    selected_folds = int(sum(record["selected"] is not None for record in fold_records))
    if selected_folds == 5:
        sparse_metrics = fixed_decision_metrics(
            labels, sparse_scores, sparse_predicted
        )
        status = "DESIGN_IN_ALL_OUTER_FOLDS"
    else:
        sparse_metrics = None
        status = "NO DESIGN IN ONE OR MORE OUTER FOLDS"

    result = {
        "status": status,
        "outer_fold_count": 5,
        "outer_folds_with_design": selected_folds,
        "qualifying_outer_fold_fraction": finite_float(selected_folds / 5.0),
        "subjects_with_sparse_outer_prediction": successful,
        "selection_frequency": dict(sorted(selection_counts.items())),
        "folds": fold_records,
        "sparse_metrics": sparse_metrics,
        "all_window_metrics": comparator_metrics,
    }
    return result, manifest


def fixed_decision_metrics(
    labels: np.ndarray,
    scores: np.ndarray,
    predicted: np.ndarray,
) -> dict:
    tn, fp, fn, tp = confusion_matrix(
        labels, predicted, labels=[0, 1]
    ).ravel()
    return {
        "auroc": finite_float(roc_auc_score(labels, scores)),
        "sensitivity": finite_float(tp / (tp + fn)),
        "specificity": finite_float(tn / (tn + fp)),
        "tp": int(tp),
        "fn": int(fn),
        "tn": int(tn),
        "fp": int(fp),
        "scores": [finite_float(value) for value in scores],
        "predicted": [int(value) for value in predicted],
    }


def evaluate_fixed_control(
    features: np.ndarray,
    labels: np.ndarray,
    splits: Sequence[tuple[np.ndarray, np.ndarray]],
    centers: Sequence[int],
) -> dict:
    indices = tuple(range(len(centers)))
    matrix = model_matrix(features, indices)
    scores = crossfit_scores(matrix, labels, splits, 5)
    metrics = constrained_metrics(labels, scores)
    state = fit_model(matrix, labels)
    return {
        "centers_cm_inverse": list(centers),
        "metrics": {
            **metrics,
            "oof_scores": [finite_float(value) for value in scores],
        },
        "fitted_model": serializable_model(state),
    }


def score_unlabeled_controls(
    features_all: np.ndarray,
    missing_mask: np.ndarray,
    selected_indices: Sequence[int],
    fitted_state: dict,
    threshold: float,
) -> dict:
    control_matrix = model_matrix(features_all[missing_mask], selected_indices)
    scores = predict_model(fitted_state, control_matrix)
    quantiles = np.quantile(scores, [0.0, 0.25, 0.5, 0.75, 1.0])
    return {
        "name": "unlabeled source-convention plausibility check",
        "n_unlabeled_controls": int(scores.size),
        "outcome_labels_assigned": False,
        "performance_metrics_computed": False,
        "scores": [finite_float(value) for value in scores],
        "score_distribution": {
            "minimum": finite_float(quantiles[0]),
            "q1": finite_float(quantiles[1]),
            "median": finite_float(quantiles[2]),
            "q3": finite_float(quantiles[3]),
            "maximum": finite_float(quantiles[4]),
        },
        "frozen_primary_threshold": finite_float(threshold),
        "fraction_below_threshold": finite_float(np.mean(scores < threshold)),
    }


def code_sha256() -> str:
    return sha256_bytes(Path(__file__).read_bytes())


def main() -> None:
    archive = Path(sys.argv[1] if len(sys.argv) > 1 else "data/GlobalDKD.zip")
    output_dir = Path(
        sys.argv[2] if len(sys.argv) > 2 else "outputs/discovery"
    )
    preregistration = Path(
        sys.argv[3] if len(sys.argv) > 3 else "PREREGISTRATION.md"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    data = load_australia_only(archive)
    spectra_all = data["spectra_all"]
    wavenumbers = data["wavenumbers"]
    finite_mask = data["finite_mask"]
    labels = data["labels_primary"]
    ids_all = data["subject_ids_all"]
    ids_primary = ids_all[finite_mask]
    spectra_primary = spectra_all[finite_mask]

    duplicates = duplicate_audit(
        spectra_all, ids_all, finite_mask, labels
    )
    candidates = candidate_library()
    primary_features_all = integrate_windows(
        spectra_all, wavenumbers, CENTERS, WINDOW_WIDTH
    )
    primary_features = primary_features_all[finite_mask]
    atmospheric_features = integrate_windows(
        spectra_primary,
        wavenumbers,
        ATMOSPHERIC_CONTROL_CENTERS,
        WINDOW_WIDTH,
    )
    clean_control_features = integrate_windows(
        spectra_primary,
        wavenumbers,
        CLEAN_OFF_TARGET_CENTERS,
        WINDOW_WIDTH,
    )

    primary_splits = materialize_primary_splits(labels)
    primary_manifest = split_manifest(
        "primary_repeated_5x5", primary_splits, ids_primary, 5
    )
    primary_search = exhaustive_search(
        primary_features,
        labels,
        primary_splits,
        5,
        candidates,
        retain_all=True,
        progress_prefix="primary",
    )

    nested, nested_manifest = nested_audit(
        primary_features, labels, ids_primary, candidates
    )
    fold_manifest = {
        "primary": primary_manifest,
        "nested": nested_manifest,
    }
    fold_payload = canonical_bytes(fold_manifest)
    fold_digest = sha256_bytes(fold_payload)

    atmospheric_control = evaluate_fixed_control(
        atmospheric_features,
        labels,
        primary_splits,
        ATMOSPHERIC_CONTROL_CENTERS,
    )
    clean_control = evaluate_fixed_control(
        clean_control_features,
        labels,
        primary_splits,
        CLEAN_OFF_TARGET_CENTERS,
    )

    selected_design = None
    unlabeled_check = None
    if primary_search["selected"] is not None:
        selected = primary_search["selected"]
        selected_indices = tuple(selected["indices"])
        selected_matrix = model_matrix(primary_features, selected_indices)
        selected_state = fit_model(selected_matrix, labels)
        selected_design = {
            "selection": selected,
            "fitted_model": serializable_model(selected_state),
            "threshold": selected["threshold"],
        }
        unlabeled_check = score_unlabeled_controls(
            primary_features_all,
            ~finite_mask,
            selected_indices,
            selected_state,
            selected["threshold"],
        )

    all_indices = tuple(range(len(CENTERS)))
    comparator_state = fit_model(
        model_matrix(primary_features, all_indices), labels
    )

    prereg_sha = (
        sha256_bytes(preregistration.read_bytes())
        if preregistration.exists()
        else None
    )
    if prereg_sha is None:
        raise RuntimeError(f"missing preregistration: {preregistration}")

    result = {
        "schema_version": "globaldkd-sparse-window-discovery-v1",
        "phase": "AUSTRALIAN DISCOVERY ONLY; SPANISH OUTCOMES UNREAD",
        "status": primary_search["status"],
        "title": TITLE,
        "scope": {
            "endpoint": (
                "laboratory-measured urinary albumin concentration >=30 mg/L "
                "after the source extraction workflow"
            ),
            "not_uacr": True,
            "not_dkd_diagnosis": True,
            "not_hardware_simulation": True,
            "certificate_scope": (
                "fixed computational window library and frozen pipeline on "
                "source spectra; k acquired windows create 2k digital features"
            ),
        },
        "isolation_guard": {
            "archive_members_read": [MAT_MEMBER],
            "loadmat_variable_names": list(REQUESTED_VARIABLES),
            "loaded_public_variables": list(REQUESTED_VARIABLES),
            "forbidden_variables_not_loaded": list(FORBIDDEN_VARIABLES),
            "AUS_Int_fields_accessed": ["data", "axisscale", "label"],
            "spanish_reference_workbook_opened": False,
            "spanish_outcomes_read": False,
        },
        "input": {
            "archive_md5": data["archive_md5"],
            "mat_member": MAT_MEMBER,
            "mat_member_sha256": data["mat_sha256"],
            "code_sha256": code_sha256(),
            "preregistration_sha256": prereg_sha,
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
            "thread_environment": {
                key: os.environ.get(key)
                for key in (
                    "OMP_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS",
                )
            },
        },
        "data_integrity": {
            "all_australian_rows": 177,
            "measured_primary_rows": int(finite_mask.sum()),
            "unlabeled_missing_albumin_controls": int((~finite_mask).sum()),
            "positive_definition": "finite albumin concentration >=30 mg/L",
            "primary_positive": int(labels.sum()),
            "primary_negative": int((1 - labels).sum()),
            "unique_subject_ids": len(set(ids_all.tolist())),
            "subject_ids_all": [str(value) for value in ids_all],
            "duplicate_audit": duplicates,
        },
        "integration": {
            "definition": (
                "sort ascending; linearly interpolate both closed endpoints; "
                "include interior points; np.trapezoid/actual width"
            ),
            "nominal_width_cm_inverse": WINDOW_WIDTH,
            "zero_norm_epsilon": ZERO_EPSILON,
            "candidate_centers_cm_inverse_in_locked_order": list(CENTERS),
            "atmospheric_instrument_artifact_control_centers": list(
                ATMOSPHERIC_CONTROL_CENTERS
            ),
            "clean_off_target_control_centers": list(
                CLEAN_OFF_TARGET_CENTERS
            ),
            "feature_mapping": (
                "k raw mean-absorbance summaries concatenated with k row-wise "
                "L2-normalized summaries; 2k digital features"
            ),
        },
        "model": {
            "scaler": "training-partition mean and population standard deviation",
            "logistic_regression": MODEL_SPEC,
            "convergence_warnings": "fatal",
        },
        "selection_gates": {
            "candidate_count": len(candidates),
            "subset_sizes": [1, 2, 3, 4],
            "minimum_auroc": MIN_AUROC,
            "minimum_sensitivity": MIN_SENSITIVITY,
            "minimum_specificity": MIN_SPECIFICITY,
            "maximum_auroc_gap_from_all_window_comparator": MAX_COMPARATOR_GAP,
            "tie_order": (
                "fewest windows, higher AUROC, lexicographically earliest "
                "index tuple in locked library order"
            ),
            "no_qualifier_rule": "NO DESIGN; no fallback",
        },
        "cross_validation": {
            "primary": {
                "class": "RepeatedStratifiedKFold",
                "n_splits": 5,
                "n_repeats": 5,
                "seed": PRIMARY_SEED,
                "aggregation": (
                    "arithmetic mean of exactly five OOF probabilities per subject"
                ),
            },
            "nested_outer": {
                "class": "StratifiedKFold",
                "n_splits": 5,
                "shuffle": True,
                "seed": OUTER_SEED,
            },
            "nested_inner": {
                "class": "StratifiedKFold",
                "n_splits": 5,
                "shuffle": True,
                "seeds_by_outer_fold": list(INNER_SEEDS),
            },
            "fold_manifest_sha256": fold_digest,
            "fold_manifest": fold_manifest,
        },
        "discovery": {
            "selected_design": selected_design,
            "qualifying_count": primary_search["qualifying_count"],
            "qualifying_by_size": primary_search["qualifying_by_size"],
            "all_window_comparator": {
                **primary_search["comparator"],
                "fitted_model": serializable_model(comparator_state),
            },
            "candidate_metrics": primary_search["candidate_metrics"],
        },
        "nested_selection_audit": nested,
        "negative_controls": {
            "atmospheric_instrument_artifact": atmospheric_control,
            "clean_off_target": clean_control,
        },
        "unlabeled_control_plausibility_check": unlabeled_check,
        "future_spanish_confirmation_lock": {
            "status": "LOCKED",
            "bootstrap_resamples": 10000,
            "bootstrap_seed": 20260812,
            "bootstrap_strata": "site x class",
            "paired_sparse_comparator_indices": True,
            "permutations": 10000,
            "permutation_seed": 20260813,
            "permutation_rule": "shuffle labels within site",
            "permutation_p": "(1+exceed)/(B+1)",
            "sensitivity_specificity_ci": "95% Wilson",
            "calibration": "exploratory only",
            "below_loq_or_nan_to_zero": (
                "left-censored threshold references, not exact concentrations"
            ),
        },
        "assertions": {
            "archive_md5": True,
            "AUS_Int_only": True,
            "shape_177_by_1598": True,
            "unique_AU_ids": True,
            "first_22_albumin_missing": True,
            "primary_n_155": True,
            "primary_counts_64_91": True,
            "candidate_count_12950": True,
            "no_cross_label_exact_duplicate": True,
            "all_numeric_outputs_finite": True,
        },
    }

    result_payload = canonical_bytes(result)
    result_digest = sha256_bytes(result_payload)
    result_path = output_dir / "frozen-discovery.json"
    result_path.write_bytes(result_payload)
    (output_dir / "frozen-discovery.sha256").write_text(
        f"{result_digest}  {result_path.name}\n",
        encoding="utf-8",
    )
    (output_dir / "fold-manifest.json").write_bytes(fold_payload)
    (output_dir / "fold-manifest.sha256").write_text(
        f"{fold_digest}  fold-manifest.json\n",
        encoding="utf-8",
    )

    selected_summary = None
    if primary_search["selected"] is not None:
        selected_summary = {
            key: value
            for key, value in primary_search["selected"].items()
            if key != "oof_scores"
        }
    summary = {
        "phase": result["phase"],
        "status": result["status"],
        "archive_md5": data["archive_md5"],
        "primary_counts": {
            "n": 155,
            "positive": 64,
            "negative": 91,
            "unlabeled_controls": 22,
        },
        "candidate_count": len(candidates),
        "selected": selected_summary,
        "all_window_comparator": {
            key: value
            for key, value in primary_search["comparator"].items()
            if key != "oof_scores"
        },
        "nested": {
            "status": nested["status"],
            "outer_folds_with_design": nested["outer_folds_with_design"],
            "qualifying_outer_fold_fraction": nested[
                "qualifying_outer_fold_fraction"
            ],
            "selection_frequency": nested["selection_frequency"],
            "sparse_metrics": nested["sparse_metrics"],
            "all_window_metrics": nested["all_window_metrics"],
        },
        "negative_controls": {
            key: value["metrics"]
            for key, value in result["negative_controls"].items()
        },
        "unlabeled_control_plausibility_check": unlabeled_check,
        "fold_manifest_sha256": fold_digest,
        "frozen_discovery_sha256": result_digest,
    }
    (output_dir / "discovery-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
