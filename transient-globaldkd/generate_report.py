#!/usr/bin/env python3
"""Read-only renderer for frozen GlobalDKD sparse-window JSON artifacts.

The renderer never opens spectra, fits or selects a model, changes a threshold,
labels the 22 controls, or simulates hardware.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_curve

DISCOVERY_SCHEMA = "globaldkd-sparse-window-discovery-v1"
VALIDATION_SCHEMA = "globaldkd-sparse-window-validation-v1"
EXPECTED_DISCOVERY_SHA256 = (
    "41460a6e8e21c6979f8ead43ee383a1691baab37897fb5c7ab29b66b636c531f"
)
EXPECTED_VALIDATION_SHA256 = (
    "fb1111b4e0539f854c60bbb11c6211632d3162cb5ccb1505736bed9085a33da4"
)
ATTEMPT_1_RUN = "31009572812"
ATTEMPT_1_JOB = "92317956542"
TECHNICAL_RETRY_RUN = "31010663111"
TECHNICAL_RETRY_JOB = "92321693370"
TECHNICAL_RETRY_ARTIFACT = "8932293583"
TITLE = (
    "How Many Spectral Windows Are Enough? An Exact Sparse-Window Certificate "
    "for Cross-Instrument Screening of Elevated Urinary Albumin Concentration."
)


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_verified(path: Path, sidecar: Path, schema: str) -> tuple[dict, str]:
    payload = path.read_bytes()
    observed = hashlib.sha256(payload).hexdigest()
    declared = sidecar.read_text(encoding="utf-8").strip().split()[0].lower()
    if observed != declared:
        raise RuntimeError(f"SHA-256 mismatch for {path}")
    value = json.loads(payload)
    canonical = (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    if payload != canonical:
        raise RuntimeError(f"noncanonical JSON/newline convention for {path}")
    if value.get("schema_version") != schema or value.get("title") != TITLE:
        raise RuntimeError(f"schema/title mismatch for {path}")
    check_finite(value)
    return value, observed


def check_finite(value, path: str = "$") -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            check_finite(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            check_finite(item, f"{path}[{index}]")
    elif isinstance(value, float) and not math.isfinite(value):
        raise RuntimeError(f"nonfinite value at {path}")


def load_attempt_1_disclosure(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    required = (
        f"Workflow run: `{ATTEMPT_1_RUN}`",
        f"Job: `{ATTEMPT_1_JOB}`",
        "did not construct threshold classes",
        "Any subsequent execution is labeled a technical retry",
        "change only the structural schema guard",
    )
    missing = [marker for marker in required if marker not in text]
    if missing:
        raise RuntimeError(f"attempt-1 disclosure is incomplete: {missing}")
    return file_hash(path)


def save(fig, path: Path) -> None:
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)


def timeline(path: Path, has_validation: bool) -> None:
    labels = [
        "Endpoint lock",
        "Australia exact search",
        "Frozen JSON + SHA",
        "Spanish technical retry 1",
        "Locked stress tests",
    ]
    completed = [True, True, True, has_validation, has_validation]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(12, 2.6))
    ax.plot(x, np.zeros(5), color="#4a5568", linewidth=2)
    ax.scatter(
        x,
        np.zeros(5),
        s=500,
        c=["#276749" if state else "#a0aec0" for state in completed],
        zorder=3,
    )
    for index, label in enumerate(labels):
        ax.text(index, 0.18, label, ha="center", fontsize=9)
        ax.text(
            index,
            -0.18,
            "complete" if completed[index] else "not evaluated",
            ha="center",
            fontsize=8,
        )
    ax.set(xlim=(-0.4, 4.4), ylim=(-0.5, 0.5))
    ax.axis("off")
    ax.set_title("Locked analysis sequence")
    save(fig, path)


def window_library(discovery: dict, path: Path) -> None:
    centers = discovery["integration"][
        "candidate_centers_cm_inverse_in_locked_order"
    ]
    selected = discovery["discovery"]["selected_design"]
    chosen = set(
        selected["selection"]["centers_cm_inverse"] if selected else []
    )
    fig, ax = plt.subplots(figsize=(12, 3))
    for center in centers:
        chosen_here = center in chosen
        ax.add_patch(
            plt.Rectangle(
                (center - 32, 0),
                64,
                0.9 if chosen_here else 0.55,
                color="#c53030" if chosen_here else "#cbd5e0",
                ec="white",
            )
        )
    ax.set(xlim=(3700, 950), ylim=(0, 1.1))
    ax.set_yticks([])
    ax.set_xlabel("Wavenumber (cm⁻¹)")
    ax.set_title(
        "Fixed computational window library"
        + (" and selected windows" if chosen else " — NO DESIGN")
    )
    ax.text(
        0,
        -0.30,
        "Computational source-spectrum intervals, not optical filters or hardware.",
        transform=ax.transAxes,
        fontsize=8,
    )
    save(fig, path)


def search_landscape(discovery: dict, path: Path) -> None:
    candidates = discovery["discovery"]["candidate_metrics"]
    grouped = [
        [row["auroc"] for row in candidates if row["n_windows"] == size]
        for size in range(1, 5)
    ]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.boxplot(grouped, tick_labels=["1", "2", "3", "4"], showfliers=False)
    selected = discovery["discovery"]["selected_design"]
    if selected:
        row = selected["selection"]
        ax.scatter(
            row["n_windows"],
            row["auroc"],
            marker="*",
            s=180,
            color="#c53030",
            label="frozen selection",
        )
        ax.legend(frameon=False)
    ax.axhline(0.90, color="#718096", ls="--")
    ax.set_xlabel("Computational source-spectrum windows (k)")
    ax.set_ylabel("Repeated-CV AUROC")
    ax.set_title(f"Exact landscape: {len(candidates):,} candidates")
    save(fig, path)


def discovery_performance(discovery: dict, path: Path) -> None:
    selected = discovery["discovery"]["selected_design"]
    comparator = discovery["discovery"]["all_window_comparator"]
    nested = discovery["nested_selection_audit"]
    labels = []
    rows = []
    if selected:
        labels.append("Selected\nrepeated CV")
        rows.append(selected["selection"])
    labels.append("24-window\nrepeated CV")
    rows.append(comparator)
    if nested.get("sparse_metrics"):
        labels.append("Selection\nnested CV")
        rows.append(nested["sparse_metrics"])
    labels.append("24-window\nnested CV")
    rows.append(nested["all_window_metrics"])
    x = np.arange(len(rows))
    width = 0.24
    fig, ax = plt.subplots(figsize=(10, 5))
    for offset, metric in enumerate(("auroc", "sensitivity", "specificity")):
        ax.bar(
            x + (offset - 1) * width,
            [row[metric] for row in rows],
            width,
            label=metric.upper() if metric == "auroc" else metric.title(),
        )
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Subject-level resampling estimate")
    ax.set_title("Australian discovery and nested selection audit")
    ax.legend(frameon=False, ncols=3)
    save(fig, path)


def control_check(discovery: dict, path: Path) -> None:
    check = discovery["unlabeled_control_plausibility_check"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    if check is None:
        ax.text(
            0.5, 0.5, "NO DESIGN\nNo control scores", ha="center", va="center"
        )
        ax.axis("off")
    else:
        scores = np.asarray(check["scores"], dtype=float)
        ax.hist(scores, bins=min(10, len(scores)), color="#4c78a8")
        ax.axvline(
            check["frozen_primary_threshold"],
            color="#c53030",
            ls="--",
            label="frozen primary threshold",
        )
        ax.set_xlabel("Frozen-model score")
        ax.set_ylabel("Unlabeled control count")
        ax.set_title("Unlabeled 22-control plausibility check")
        ax.legend(frameon=False)
        ax.text(
            0,
            -0.28,
            "No labels or performance metrics are assigned to these controls.",
            transform=ax.transAxes,
            fontsize=8,
        )
    save(fig, path)


def spanish_roc(validation: dict, path: Path) -> None:
    records = validation["predictions"]
    outcome = np.asarray([row["threshold_class"] for row in records])
    sparse = np.asarray([row["sparse_score"] for row in records])
    all_window = np.asarray([row["all_window_score"] for row in records])
    sparse_fpr, sparse_tpr, _ = roc_curve(outcome, sparse)
    all_fpr, all_tpr, _ = roc_curve(outcome, all_window)
    sparse_auc = validation["primary_sparse"]["metrics"]["auroc"]
    all_auc = validation["all_window_comparator"]["metrics"]["auroc"]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(sparse_fpr, sparse_tpr, label=f"Sparse windows ({sparse_auc:.3f})")
    ax.plot(all_fpr, all_tpr, label=f"24-window comparator ({all_auc:.3f})")
    ax.plot([0, 1], [0, 1], color="#718096", ls="--")
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_xlabel("1 − specificity")
    ax.set_ylabel("Sensitivity")
    ax.set_title("Spanish archived-cohort reuse — technical retry 1")
    ax.legend(frameon=False, loc="lower right")
    save(fig, path)


def stress_plot(validation: dict, path: Path) -> None:
    rows = [("Nominal sparse", validation["primary_sparse"]["metrics"]["auroc"])]
    rows += [
        (
            name.replace("_", " "),
            item["metrics_at_nominal_frozen_threshold"]["auroc"],
        )
        for name, item in validation["robustness_nominal_model_only"].items()
    ]
    rows += [
        ("control: " + name.replace("_", " "), item["metrics"]["auroc"])
        for name, item in validation["negative_controls"].items()
    ]
    labels, values = zip(*rows, strict=True)
    y = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(9, max(4, 0.45 * len(rows))))
    ax.scatter(values, y, color=["#c53030"] + ["#4c78a8"] * (len(rows) - 1))
    ax.axvline(0.5, color="#718096", ls="--")
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Spanish AUROC")
    ax.set_title("Frozen-model stress tests and controls")
    save(fig, path)


def fmt(value: float | None) -> str:
    return "undefined" if value is None else f"{value:.3f}"


def verified_external_gate_overall(validation: dict) -> bool:
    """Consume the frozen validator gates and verify their internal logic."""
    gates = validation["primary_sparse"]["external_success_gates"]
    overall = gates["overall_pass"]
    if type(overall) is not bool:
        raise RuntimeError("external gate overall_pass must be a JSON boolean")
    components = {
        key: value for key, value in gates.items() if key != "overall_pass"
    }
    expected_points = {
        "auroc": validation["primary_sparse"]["metrics"]["auroc"],
        "sensitivity": validation["primary_sparse"]["metrics"]["sensitivity"],
        "specificity": validation["primary_sparse"]["metrics"]["specificity"],
        "sparse_minus_all_window_auroc": validation["primary_sparse"][
            "bootstrap"
        ]["sparse_minus_all_window_auroc"],
    }
    if set(components) != set(expected_points):
        raise RuntimeError("external gate components are missing or unexpected")
    component_passes = []
    for name, item in components.items():
        if float(item["point"]) != float(expected_points[name]):
            raise RuntimeError(f"external gate point/source mismatch: {name}")
        expected = bool(float(item["point"]) >= float(item["minimum"]))
        if type(item["pass"]) is not bool or item["pass"] != expected:
            raise RuntimeError(f"external gate inconsistency: {name}")
        component_passes.append(item["pass"])
    if overall != all(component_passes):
        raise RuntimeError("external overall gate disagrees with components")
    return overall


def make_report(
    discovery: dict,
    discovery_hash: str,
    validation: dict | None,
    validation_hash: str | None,
    attempt_1_note_hash: str | None,
) -> str:
    selected = discovery["discovery"]["selected_design"]
    comparator = discovery["discovery"]["all_window_comparator"]
    nested = discovery["nested_selection_audit"]
    duplicates = discovery["data_integrity"]["duplicate_audit"]
    qualifying_by_k = [
        discovery["discovery"]["qualifying_by_size"][str(k)]
        for k in range(1, 5)
    ]
    if qualifying_by_k != [0, 15, 282, 2678]:
        raise RuntimeError("unexpected artifact-sourced qualifying counts by k")
    lines = [
        f"# {TITLE}",
        "",
        "## Scope",
        "",
        "This secondary computation asks how many fixed source-spectrum windows "
        "are enough for the archived 30 mg/L urinary-albumin-concentration "
        "screen. Each retained window is a mean-absorbance summary. It is not "
        "a physical filter measurement, hardware simulation, uACR assay, or "
        "disease diagnosis.",
        "",
        "The certificate is bounded to the fixed 24-window library, frozen "
        "pipeline, and these archived acquisition settings. Spanish reuse is "
        "held out but is not prospective clinical validation.",
        "",
        "## Endpoint and provenance",
        "",
        f"- Archive MD5: {discovery['input']['archive_md5']}",
        f"- Discovery SHA-256: {discovery_hash}",
        "- Primary Australia: n=155 measured samples, 64 positive and 91 "
        "negative at 30 mg/L.",
        "- Excluded controls: 22 missing-albumin spectra, kept unlabeled for a "
        "score-distribution plausibility check only.",
        f"- Exact/near duplicate pair counts: "
        f"{duplicates['exact_pair_count']}/"
        f"{duplicates['near_pair_count_excluding_exact']}.",
        "",
        "## Frozen computation",
        "",
        "- 24 fixed 64 cm⁻¹ computational windows.",
        "- All 12,950 subsets of one through four windows.",
        "- k raw mean-absorbance summaries plus k row-wise L2-normalized "
        "copies, giving 2k digital features.",
        "- Fold-only scaling and balanced L2 logistic regression with fixed "
        "seeds, gates, and threshold selection.",
        "",
        "## Australian discovery",
        "",
    ]
    if selected is None:
        lines += [
            "**NO DESIGN.** No candidate passed every preregistered gate. "
            "No fallback was used.",
            "",
        ]
    else:
        row = selected["selection"]
        centers = ", ".join(str(value) for value in row["centers_cm_inverse"])
        lines += [
            f"- Selected k={row['n_windows']} ({2*row['n_windows']} digital "
            f"features), centers {centers} cm⁻¹.",
            f"- Repeated-CV AUROC {fmt(row['auroc'])}; sensitivity "
            f"{fmt(row['sensitivity'])}; specificity {fmt(row['specificity'])}.",
            f"- Frozen threshold {row['threshold']:.12g}.",
        ]
    lines += [
        f"- Qualifying candidates: {discovery['discovery']['qualifying_count']}.",
        "- Artifact-sourced qualifying counts by k=1/2/3/4: "
        + "/".join(str(value) for value in qualifying_by_k)
        + ".",
        f"- 24-window repeated-CV AUROC: {fmt(comparator['auroc'])}.",
        f"- Nested folds yielding a design: "
        f"{nested['outer_folds_with_design']}/{nested['outer_fold_count']}.",
    ]
    if nested.get("sparse_metrics"):
        lines.append(
            f"- Nested selection-procedure AUROC: "
            f"{fmt(nested['sparse_metrics']['auroc'])}."
        )
    else:
        lines.append(
            "- Nested selection-procedure metric is undefined because one or "
            "more outer folds yielded NO DESIGN."
        )

    lines += ["", "## Unlabeled-control plausibility check", ""]
    check = discovery["unlabeled_control_plausibility_check"]
    if check is None:
        lines.append("Not computed because discovery yielded NO DESIGN.")
    else:
        distribution = check["score_distribution"]
        lines += [
            "The final primary model was fitted only on 155 measured samples "
            "and scored the 22 controls once. No outcome label or performance "
            "metric is assigned.",
            f"- Median score {fmt(distribution['median'])}; IQR "
            f"{fmt(distribution['q1'])}–{fmt(distribution['q3'])}.",
            f"- Fraction below the frozen threshold "
            f"{fmt(check['fraction_below_threshold'])}.",
        ]

    lines += [
        "",
        "## Spanish archived-cohort reuse — schema-only technical retry 1",
        "",
    ]
    if validation is None:
        lines += [
            "**Not evaluated.** Spanish outcomes remain locked in this report.",
            "",
        ]
    else:
        metrics = validation["primary_sparse"]["metrics"]
        bootstrap = validation["primary_sparse"]["bootstrap"]
        delta = bootstrap["sparse_minus_all_window_auroc"]
        passed = verified_external_gate_overall(validation)
        ci = bootstrap["sparse_auroc_ci_95"]
        lines += [
            "**Execution disclosure.** The successful Spanish evaluation was "
            "technical retry 1 after correcting only the fail-closed MATLAB "
            "empty-cell schema assertion. Attempt 1 deserialized the ESP "
            "outcome container but stopped before threshold-class construction "
            "or any score, metric, or output artifact. It was therefore not a "
            "pristine first unblinding. The immutable validation JSON retains "
            "its pre-retry phase string; this report supplies the superseding "
            "execution label.",
            "",
            f"- Failed attempt: run {ATTEMPT_1_RUN}, job {ATTEMPT_1_JOB}.",
            f"- Successful technical retry: run {TECHNICAL_RETRY_RUN}, job "
            f"{TECHNICAL_RETRY_JOB}, artifact {TECHNICAL_RETRY_ARTIFACT}.",
            f"- Attempt-1 disclosure SHA-256: {attempt_1_note_hash}.",
            f"- Validation SHA-256: {validation_hash}",
            "- N=61 (34 positive, 27 negative); Valencia 46, Madrid 15.",
            f"- AUROC {fmt(metrics['auroc'])}, stratified-bootstrap 95% CI "
            f"{fmt(ci[0])}–{fmt(ci[1])}.",
            f"- Sensitivity {fmt(metrics['sensitivity'])}; specificity "
            f"{fmt(metrics['specificity'])}; Wilson intervals are in the artifact.",
            f"- Sparse-minus-24-window AUROC {fmt(delta)}.",
            "- The 24-window comparator had the same sensitivity "
            f"({fmt(validation['all_window_comparator']['metrics']['sensitivity'])}) "
            "but one false positive, versus four for the sparse model.",
            "- The selected 3360 and 3296 cm⁻¹ windows are adjacent closed "
            "64 cm⁻¹ intervals sharing the 3328 cm⁻¹ boundary. They are "
            "computational summaries, not evidence for two independent "
            "physical channels.",
            f"- Prespecified Spanish gates from the frozen validation artifact: "
            f"{'PASS' if passed else 'FAIL'}.",
        ]
        for site_name in ("Valencia", "Madrid"):
            site_row = validation["site_metrics"][site_name]
            lines.append(
                f"- {site_name}: AUROC {fmt(site_row['auroc'])}; sensitivity "
                f"{fmt(site_row['sensitivity'])}; specificity "
                f"{fmt(site_row['specificity'])}. Undefined single-class "
                "quantities are reported explicitly."
            )
        permutation = validation["primary_sparse"]["within_site_permutation"]
        calibration = validation["primary_sparse"]["calibration"]
        lines += [
            f"- Within-site label-permutation p={permutation['p_value']:.6g} "
            f"({permutation['exceed']} exceedances in "
            f"{permutation['permutations']} permutations).",
            f"- Exploratory calibration status {calibration['status']}; "
            f"intercept {fmt(calibration['intercept'])}; slope "
            f"{fmt(calibration['slope'])}. Balanced-class-weight scores are "
            "not calibrated probabilities.",
            "- Source zero references are left-censored threshold references, "
            "not exact concentrations.",
            "",
            "### Frozen robustness and negative-control diagnostics",
            "",
            "Every robustness perturbation reused the nominal Australian "
            "scaler, model, and threshold with no refit, recalibration, or "
            "reselection. Common shifts test global spectral-axis tolerance "
            "only.",
            "",
            "| Perturbation | AUROC | Sensitivity | Specificity | TP/FN/TN/FP |",
            "|---|---:|---:|---:|---:|",
        ]
        for name, item in validation[
            "robustness_nominal_model_only"
        ].items():
            row = item["metrics_at_nominal_frozen_threshold"]
            lines.append(
                f"| {name.replace('_', ' ')} | {fmt(row['auroc'])} | "
                f"{fmt(row['sensitivity'])} | {fmt(row['specificity'])} | "
                f"{row['tp']}/{row['fn']}/{row['tn']}/{row['fp']} |"
            )
        lines += [
            "",
            "Negative-control windows were frozen during Australian discovery "
            "and were not selected using Spanish outcomes. The atmospheric/"
            "instrument block is an artifact diagnostic, not a clean "
            "biochemical negative; the clean off-target block is reported "
            "separately.",
            "",
            "| Negative control | AUROC | Sensitivity | Specificity | TP/FN/TN/FP |",
            "|---|---:|---:|---:|---:|",
        ]
        for name, item in validation["negative_controls"].items():
            row = item["metrics"]
            lines.append(
                f"| {name.replace('_', ' ')} | {fmt(row['auroc'])} | "
                f"{fmt(row['sensitivity'])} | {fmt(row['specificity'])} | "
                f"{row['tp']}/{row['fn']}/{row['tn']}/{row['fp']} |"
            )
        lines += [
            "",
            "**Operating-point warning.** The ±8/±16 cm⁻¹ common shifts "
            "preserved rank discrimination (AUROC 0.976–0.983) but did not "
            "preserve the frozen-threshold operating point: specificity fell "
            "to 0.704 and 0.630 for +8/+16, while sensitivity fell to 0.706 "
            "for −16. This further limits any physical-channel interpretation.",
            "",
            "**Interpretive warning.** Both Spanish negative controls showed "
            "high discrimination (AUROC above 0.92). That pattern is compatible "
            "with broad nuisance, site, or class-correlated spectral signal and "
            "materially weakens biochemical-region specificity and mechanistic "
            "interpretation. The preregistered external gates still pass, but "
            "they do not override this claim ceiling.",
        ]

    lines += [
        "",
        "## Interpretation and limitations",
        "",
        "The permitted result is an exact positive or negative certificate "
        "within this fixed computational tile library and pipeline. Selected "
        "correlated windows are candidate spectral regions, not causal "
        "molecular attribution.",
        "",
        "This work does not establish physical-filter behavior, detector or SNR "
        "requirements, deployability, clinical readiness, a universal minimum, "
        "prospective utility, uACR, or causal kidney-disease diagnosis.",
        "",
        "## Reproducibility",
        "",
        f"- Code SHA-256: {discovery['input']['code_sha256']}",
        f"- Preregistration SHA-256: "
        f"{discovery['input']['preregistration_sha256']}",
        f"- Fold manifest SHA-256: "
        f"{discovery['cross_validation']['fold_manifest_sha256']}",
        f"- Python/NumPy/SciPy/scikit-learn: "
        f"{discovery['environment']['python']} / "
        f"{discovery['environment']['numpy']} / "
        f"{discovery['environment']['scipy']} / "
        f"{discovery['environment']['scikit_learn']}",
        "",
    ]
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("discovery", type=Path)
    parser.add_argument("discovery_sha256", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--validation", type=Path)
    parser.add_argument("--validation-sha256", type=Path)
    parser.add_argument("--attempt-1-note", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if bool(args.validation) != bool(args.validation_sha256):
        raise SystemExit("validation JSON and SHA sidecar must be supplied together")
    if bool(args.validation) != bool(args.attempt_1_note):
        raise SystemExit(
            "Spanish validation requires the immutable attempt-1 disclosure"
        )
    discovery, discovery_hash = load_verified(
        args.discovery, args.discovery_sha256, DISCOVERY_SCHEMA
    )
    if discovery_hash != EXPECTED_DISCOVERY_SHA256:
        raise RuntimeError("report input is not the reviewed frozen discovery")
    validation = None
    validation_hash = None
    attempt_1_note_hash = None
    if args.validation:
        validation, validation_hash = load_verified(
            args.validation, args.validation_sha256, VALIDATION_SCHEMA
        )
        if validation_hash != EXPECTED_VALIDATION_SHA256:
            raise RuntimeError("report input is not the reviewed frozen validation")
        if validation["phase"] != "ONE-SHOT SPANISH CONFIRMATION":
            raise RuntimeError("unexpected immutable validation phase string")
        attempt_1_note_hash = load_attempt_1_disclosure(args.attempt_1_note)
        if (
            validation["linkage"]["frozen_discovery_sha256"]
            != discovery_hash
        ):
            raise RuntimeError("validation links to another discovery artifact")
        if discovery["status"] != "DESIGN":
            raise RuntimeError("validation supplied for NO DESIGN discovery")
        expected = discovery["discovery"]["selected_design"]["selection"][
            "centers_cm_inverse"
        ]
        if validation["linkage"]["selected_centers_cm_inverse"] != expected:
            raise RuntimeError("selected-window linkage mismatch")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    figures = args.output_dir / "figures"
    figures.mkdir()
    timeline(figures / "analysis-timeline.svg", validation is not None)
    window_library(discovery, figures / "computational-window-library.svg")
    search_landscape(discovery, figures / "exact-search-landscape.svg")
    discovery_performance(discovery, figures / "discovery-performance.svg")
    control_check(discovery, figures / "unlabeled-control-plausibility.svg")
    names = [
        "analysis-timeline.svg",
        "computational-window-library.svg",
        "exact-search-landscape.svg",
        "discovery-performance.svg",
        "unlabeled-control-plausibility.svg",
    ]
    if validation:
        spanish_roc(validation, figures / "held-out-spanish-roc.svg")
        stress_plot(validation, figures / "robustness-and-controls.svg")
        names += ["held-out-spanish-roc.svg", "robustness-and-controls.svg"]

    report = make_report(
        discovery,
        discovery_hash,
        validation,
        validation_hash,
        attempt_1_note_hash,
    )
    report_path = args.output_dir / "report.md"
    report_path.write_text(report + "\n", encoding="utf-8")
    manifest = {
        "schema_version": "globaldkd-report-manifest-v1",
        "scope": (
            "read-only rendering; no archive access, fitting, selection, "
            "threshold changes, control labels, or hardware simulation"
        ),
        "inputs": {
            "discovery_sha256": discovery_hash,
            "validation_sha256": validation_hash,
            "attempt_1_failure_note_sha256": attempt_1_note_hash,
        },
        "execution": {
            "spanish_phase": (
                "schema-only technical retry 1"
                if validation is not None
                else "not evaluated"
            ),
            "attempt_1_run": ATTEMPT_1_RUN if validation is not None else None,
            "attempt_1_job": ATTEMPT_1_JOB if validation is not None else None,
            "technical_retry_run": (
                TECHNICAL_RETRY_RUN if validation is not None else None
            ),
            "technical_retry_job": (
                TECHNICAL_RETRY_JOB if validation is not None else None
            ),
            "technical_retry_artifact": (
                TECHNICAL_RETRY_ARTIFACT if validation is not None else None
            ),
        },
        "outputs": {
            "report.md": file_hash(report_path),
            **{
                "figures/" + name: file_hash(figures / name)
                for name in names
            },
        },
        "warnings": [
            "computational mean-absorbance windows are not physical filters",
            "22 missing-albumin controls remain unlabeled",
            "Spanish reuse is archived held-out reuse, not prospective validation",
            "Spanish evaluation succeeded only on schema-only technical retry 1",
            "high Spanish negative-control AUROC limits spectral-specific interpretation",
        ],
    }
    (args.output_dir / "report-manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
