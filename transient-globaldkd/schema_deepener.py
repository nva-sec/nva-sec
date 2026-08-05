#!/usr/bin/env python3
"""Deep, value-safe schema inspection for the public GlobalDKD MATLAB archive.

Numeric values are never serialized or printed. Only numeric shapes and dtypes
are recorded. String values are retained because MATLAB field/category labels
and source-code identifiers are required to define an unambiguous parser.
"""
from __future__ import annotations

import io
import json
import re
import sys
import tempfile
import zipfile
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
from scipy.io import loadmat

TARGET_MATS = {"Zenodo/Australia/DATASET.mat", "Zenodo/dataset.mat"}
SYMBOLS = re.compile(
    r"AUS_Int|ESP|Global|AlbConcALL|ProCon|wnm|Reference_Values|microal",
    flags=re.IGNORECASE,
)
MAX_NODES = 10000
MAX_OBJECT_ELEMENTS = 32


def extract_mlx_code(payload: bytes) -> str:
    """Return source nodes from a MATLAB live script, excluding output.xml."""
    with zipfile.ZipFile(io.BytesIO(payload)) as live:
        xml_members = [
            name for name in live.namelist()
            if name.lower().endswith("document.xml")
        ]
        blocks = []
        for member in xml_members:
            try:
                root = ElementTree.fromstring(live.read(member))
            except ElementTree.ParseError:
                continue
            preferred = []
            fallback = []
            for node in root.iter():
                tag = node.tag.rsplit("}", 1)[-1].lower()
                text = "".join(node.itertext()).strip()
                if not text:
                    continue
                if tag == "mcode":
                    preferred.append(text)
                elif tag in {"code", "input"}:
                    fallback.append(text)
            candidates = preferred or fallback
            for text in candidates:
                if not blocks or blocks[-1] != text:
                    blocks.append(text)
        return "\n\n% ---- live-script code block ----\n\n".join(blocks)


def append_record(records: list[dict], record: dict) -> bool:
    if len(records) >= MAX_NODES:
        return False
    records.append(record)
    return True


def summarize_value(value, path: str, records: list[dict], depth: int = 0) -> None:
    if len(records) >= MAX_NODES:
        return
    if depth > 30:
        append_record(records, {"path": path, "type": "depth-limit"})
        return

    if isinstance(value, dict):
        keys = sorted(str(key) for key in value)
        append_record(records, {"path": path, "type": "dict", "keys": keys})
        for key in keys:
            summarize_value(value[key], f"{path}.{key}", records, depth + 1)
        return

    fieldnames = getattr(value, "_fieldnames", None)
    if fieldnames:
        names = [str(name) for name in fieldnames]
        append_record(
            records,
            {
                "path": path,
                "type": type(value).__name__,
                "fields": names,
            },
        )
        for name in names:
            summarize_value(getattr(value, name), f"{path}.{name}", records, depth + 1)
        return

    if isinstance(value, (str, bytes, np.str_, np.bytes_)):
        text = value.decode("utf-8", errors="replace") if isinstance(value, bytes) else str(value)
        append_record(records, {"path": path, "type": "string", "value": text})
        return

    if isinstance(value, np.ndarray):
        record = {
            "path": path,
            "type": type(value).__name__,
            "shape": [int(item) for item in value.shape],
            "dtype": str(value.dtype),
            "kind": value.dtype.kind,
        }
        classname = getattr(value, "classname", None)
        if classname is not None:
            record["matlab_classname"] = str(classname)
        if value.dtype.names:
            record["fields"] = list(value.dtype.names)
        append_record(records, record)

        if value.dtype.kind in {"U", "S"}:
            append_record(
                records,
                {
                    "path": f"{path}.__strings__",
                    "type": "string-array-values",
                    "values": [str(item) for item in value.reshape(-1).tolist()],
                },
            )
            return

        if value.dtype.names:
            for name in value.dtype.names:
                summarize_value(value[name], f"{path}.{name}", records, depth + 1)
            return

        if value.dtype.kind == "O":
            flat = value.reshape(-1)
            for index, item in enumerate(flat[:MAX_OBJECT_ELEMENTS]):
                summarize_value(item, f"{path}[{index}]", records, depth + 1)
            if flat.size > MAX_OBJECT_ELEMENTS:
                append_record(
                    records,
                    {
                        "path": f"{path}.__omitted__",
                        "type": "object-elements-omitted",
                        "count": int(flat.size - MAX_OBJECT_ELEMENTS),
                    },
                )
        return

    if isinstance(value, np.void):
        names = list(value.dtype.names or [])
        append_record(
            records,
            {
                "path": path,
                "type": "numpy.void",
                "dtype": str(value.dtype),
                "fields": names,
            },
        )
        for name in names:
            summarize_value(value[name], f"{path}.{name}", records, depth + 1)
        return

    if isinstance(value, (list, tuple)):
        append_record(records, {"path": path, "type": type(value).__name__, "length": len(value)})
        if all(isinstance(item, (str, bytes, np.str_, np.bytes_)) for item in value):
            append_record(
                records,
                {
                    "path": f"{path}.__strings__",
                    "type": "string-list-values",
                    "values": [
                        item.decode("utf-8", errors="replace") if isinstance(item, bytes) else str(item)
                        for item in value
                    ],
                },
            )
        else:
            for index, item in enumerate(value[:MAX_OBJECT_ELEMENTS]):
                summarize_value(item, f"{path}[{index}]", records, depth + 1)
            if len(value) > MAX_OBJECT_ELEMENTS:
                append_record(
                    records,
                    {
                        "path": f"{path}.__omitted__",
                        "type": "sequence-elements-omitted",
                        "count": len(value) - MAX_OBJECT_ELEMENTS,
                    },
                )
        return

    if isinstance(value, (bool, np.bool_)):
        append_record(records, {"path": path, "type": "boolean", "value": bool(value)})
        return

    if isinstance(value, (int, float, complex, np.number)):
        append_record(
            records,
            {
                "path": path,
                "type": "numeric-scalar",
                "dtype": str(np.asarray(value).dtype),
                "value_redacted": True,
            },
        )
        return

    append_record(records, {"path": path, "type": type(value).__name__})


def deep_mat_schema(payload: bytes, archive_path: str) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".mat") as handle:
        handle.write(payload)
        handle.flush()
        try:
            loaded = loadmat(handle.name, simplify_cells=True)
            mode = "simplify_cells=True"
        except Exception as simplify_error:
            loaded = loadmat(handle.name, squeeze_me=True, struct_as_record=False)
            mode = f"fallback squeeze/struct; simplify error={simplify_error!r}"

    records = []
    for key in sorted(key for key in loaded if not key.startswith("__")):
        summarize_value(loaded[key], key, records)
    return {
        "path": archive_path,
        "load_mode": mode,
        "record_count": len(records),
        "truncated": len(records) >= MAX_NODES,
        "records": records,
    }


def symbol_context(code: str, source: str) -> list[str]:
    lines = code.splitlines()
    selected = set()
    for index, line in enumerate(lines):
        if SYMBOLS.search(line):
            for neighbor in range(max(0, index - 1), min(len(lines), index + 2)):
                selected.add(neighbor)
    result = []
    previous = None
    for index in sorted(selected):
        if previous is not None and index > previous + 1:
            result.append("% ...")
        result.append(f"% {source}:{index + 1}\n{lines[index]}")
        previous = index
    return result


def main() -> None:
    archive = Path(sys.argv[1] if len(sys.argv) > 1 else "data/GlobalDKD.zip")
    output_dir = Path(sys.argv[2] if len(sys.argv) > 2 else "outputs/schema-audit")
    output_dir.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(archive) as bundle:
        names = set(bundle.namelist())
        missing = sorted(TARGET_MATS - names)
        if missing:
            raise SystemExit(f"missing required MAT members: {missing}")

        deep = [
            deep_mat_schema(bundle.read(name), name)
            for name in sorted(TARGET_MATS)
        ]

        dataset_mlx = "Zenodo/DATASET.mlx"
        if dataset_mlx not in names:
            raise SystemExit("missing Zenodo/DATASET.mlx")
        dataset_code = extract_mlx_code(bundle.read(dataset_mlx))

        contexts = []
        for name in sorted(
            item for item in names
            if item.lower().endswith(".mlx") and item != dataset_mlx
        ):
            code = extract_mlx_code(bundle.read(name))
            contexts.extend(symbol_context(code, name))

    (output_dir / "deep-mat-schema.json").write_text(
        json.dumps(
            {
                "phase": "schema deepener; numeric values redacted",
                "mat_files": deep,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    (output_dir / "DATASET.complete-code.txt").write_text(dataset_code + "\n", encoding="utf-8")
    (output_dir / "mlx-symbol-context.txt").write_text("\n".join(contexts) + "\n", encoding="utf-8")

    print(
        json.dumps(
            {
                "status": "PASS",
                "mat_files": [
                    {
                        "path": item["path"],
                        "record_count": item["record_count"],
                        "truncated": item["truncated"],
                    }
                    for item in deep
                ],
                "dataset_code_characters": len(dataset_code),
                "symbol_context_lines": len(contexts),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
