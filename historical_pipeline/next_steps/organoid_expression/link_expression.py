"""Join author cell annotations to original integer MTX by sample+barcode.

The author notebook is never run, and processed X/layers are never used as
counts. Aggregate all genes by library and author hepatocyte-like subtype.
These technical replicate libraries support descriptive reference use only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parent
INPUTS = BASE.parent / "organoid_inputs"
SOURCE = BASE / "sources"
MAX_UNCOMPRESSED = 4_000_000_000
OBJECT = SOURCE / "GSE207889_os.h5ad"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_json(name, data):
    (BASE / name).write_text(json.dumps(data, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")


def provenance():
    return {"executed_utc": datetime.now(timezone.utc).isoformat(),
            "command": sys.argv, "interpreter": sys.executable,
            "python_version": sys.version.split()[0],
            "script_sha256": sha256(__file__),
            "matrix_data_source": "Original GEO MTX integer entries; not processed H5AD X or layers",
            "no_statistical_inference": True,
            "no_ingredient_effect_validation": True}


def read_manifest():
    return list(csv.DictReader((BASE.parent / "organoid_sample_manifest.csv").open(encoding="utf-8-sig", newline="")))


def verify_small_input(path):
    audit = json.loads((BASE.parent / "organoid_metadata_audit.json").read_text(encoding="utf-8"))
    record = next(x for x in audit["downloads"] if Path(x["local_file"]).name == path.name)
    if sha256(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
        raise ValueError(f"Small input changed since metadata audit: {path.name}")


def decompress():
    path = SOURCE / "GSE207889_os.h5ad.gz"
    downloaded = json.loads((BASE / "download_log.json").read_text(encoding="utf-8"))
    expected = next(x for x in downloaded["records"] if x["local_file"].endswith(path.name))
    if path.stat().st_size != expected["received_bytes"] or sha256(path) != expected["sha256"]:
        raise ValueError("Compressed H5AD no longer matches recorded download")
    metadata = json.loads((BASE.parent / "organoid_metadata_audit.json").read_text(encoding="utf-8"))
    deposited_md5 = None
    for sheet_rows in metadata["spreadsheet_rows"].values():
        for cells in sheet_rows:
            if "os.h5ad.gz" in cells.values():
                deposited_md5 = next((value for value in cells.values() if re.fullmatch("[0-9a-f]{32}", value)), None)
    with path.open("rb") as handle:
        actual_md5 = hashlib.file_digest(handle, "md5").hexdigest()
    if not deposited_md5 or deposited_md5 != actual_md5:
        raise ValueError(f"H5AD gzip MD5 differs from deposited workbook: {deposited_md5} vs {actual_md5}")
    started = time.monotonic()
    if OBJECT.exists():
        if OBJECT.stat().st_size > MAX_UNCOMPRESSED:
            raise ValueError("Existing decompressed object exceeds budget")
        status = "existing_file_read_without_overwrite"
    else:
        temp = OBJECT.with_suffix(".h5ad.partial")
        if temp.exists():
            raise RuntimeError("Partial decompressed object exists; inspect before resuming")
        count = 0
        with gzip.open(path, "rb") as compressed, temp.open("xb") as output:
            while True:
                block = compressed.read(4*1024*1024)
                if not block:
                    break
                count += len(block)
                if count > MAX_UNCOMPRESSED:
                    raise ValueError("Actual decompressed H5AD exceeds 4 GB; stop and review budget")
                output.write(block)
        temp.rename(OBJECT)
        status = "decompressed"
    report = {**provenance(), "input_file": str(path.relative_to(BASE)),
              "output_file": str(OBJECT.relative_to(BASE)), "status": status,
              "uncompressed_bytes": OBJECT.stat().st_size,
              "deposited_compressed_md5": deposited_md5,
              "actual_compressed_md5": actual_md5,
              "sha256": sha256(OBJECT), "elapsed_seconds": round(time.monotonic()-started, 2),
              "raw_mtx_decompressed_copies_created": False}
    save_json("decompression_log.json", report)
    print(json.dumps(report, ensure_ascii=True), flush=True)


def decode(node):
    import h5py
    import numpy as np
    if isinstance(node, h5py.Dataset):
        if h5py.check_string_dtype(node.dtype) is not None:
            return node.asstr()[...].tolist()
        return node[...].tolist()
    if "codes" in node and "categories" in node:
        categories = decode(node["categories"])
        return [categories[int(x)] if x >= 0 else "" for x in node["codes"][...]]
    raise ValueError(f"Unsupported obs encoding {node.name}")


def inspect_annotations():
    import h5py
    import numpy as np
    manifest = read_manifest()
    aliases = {}
    for row in manifest:
        for alias in [row["author_sample_id"], row["author_short_id_alias"]]:
            aliases[alias] = row
    with h5py.File(OBJECT, "r") as obj:
        obs = obj["obs"]
        index_key = obs.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        fields = ["sample", "condition", "detailed_condition", "celltype_category", "celltype_scTypedb"]
        missing = [name for name in fields if name not in obs]
        if missing:
            raise ValueError(f"Required annotation fields absent: {missing}")
        columns = {name: decode(obs[name]) for name in fields}
        ids = decode(obs[index_key])
        if any(len(x) != len(ids) for x in columns.values()):
            raise ValueError("Annotation column length mismatch")
        schema = {"n_obs": len(ids), "obs_keys": list(obs.keys()),
                  "layers": list(obj["layers"].keys()) if "layers" in obj else [],
                  "raw_present": "raw" in obj, "X_used": False,
                  "sample_values": dict(Counter(str(x) for x in columns["sample"]))}
        kept = []
        keys = set()
        for i, raw_sample in enumerate(columns["sample"]):
            raw_sample = str(raw_sample)
            if raw_sample not in aliases:
                continue
            row = aliases[raw_sample]
            expected_condition = "CTRL" if row["group"] == "CTRL-PA" else "PA"
            if str(columns["condition"][i]) != expected_condition or str(columns["detailed_condition"][i]) != row["group"]:
                raise ValueError(f"Condition labels disagree with manifest: {row['gsm']}")
            raw_id = str(ids[i])
            matches = re.findall(r"(?<![ACGTN])[ACGTN]{16}-[0-9]+", raw_id)
            if len(matches) != 1:
                raise ValueError(f"Cannot uniquely parse original barcode: {raw_id}")
            barcode = matches[0]
            joint_key = (row["author_sample_id"], barcode)
            if joint_key in keys:
                raise ValueError(f"Duplicate sample+barcode: {joint_key}")
            keys.add(joint_key)
            kept.append({"gsm": row["gsm"], "sample": row["author_sample_id"],
                         "group": row["group"], "barcode": barcode, "original_obs_index": raw_id,
                         "raw_sample_label": raw_sample,
                         **{name: str(columns[name][i]) for name in fields if name != "sample"}})
    if set(x["gsm"] for x in kept) != {x["gsm"] for x in manifest}:
        raise ValueError("Not all four libraries are represented in author annotations")
    joins = []
    for row in manifest:
        barcode_file = INPUTS / row["barcode_url"].rsplit("/", 1)[1]
        verify_small_input(barcode_file)
        with gzip.open(barcode_file, "rt", encoding="utf-8") as handle:
            original = [line.strip() for line in handle]
        if len(original) != len(set(original)):
            raise ValueError(f"Duplicate raw barcode in {row['gsm']}")
        subset = [x for x in kept if x["gsm"] == row["gsm"]]
        absent = {x["barcode"] for x in subset} - set(original)
        if absent:
            raise ValueError(f"Processed annotation barcodes missing in raw library: {row['gsm']} {len(absent)}")
        joins.append({"gsm": row["gsm"], "sample": row["author_sample_id"],
                      "raw_barcode_count": len(original), "annotated_cells": len(subset),
                      "unmatched_annotation_barcodes": len(absent),
                      "category_counts": dict(Counter(x["celltype_category"] for x in subset)),
                      "subtype_counts": dict(Counter(x["celltype_scTypedb"] for x in subset)),
                      "condition": dict(Counter(x["condition"] for x in subset)),
                      "detailed_condition": dict(Counter(x["detailed_condition"] for x in subset)),
                      "condition_labels_match_manifest": True,
                      "barcode_source_sha256": sha256(barcode_file)})
    output = BASE / "four_library_annotations.csv.gz"
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(kept[0]))
        writer.writeheader()
        writer.writerows(kept)
    report = {**provenance(), "h5py_version": h5py.__version__, "numpy_version": np.__version__,
              "manifest_sha256": sha256(BASE.parent / "organoid_sample_manifest.csv"),
              "small_input_audit_sha256": sha256(BASE.parent / "organoid_metadata_audit.json"),
              "annotation_object": str(OBJECT.relative_to(BASE)), "schema": schema,
              "sample_join_reports": joins, "annotation_rows_saved": len(kept),
              "annotation_file": output.name, "annotation_file_sha256": sha256(output),
              "barcode_rule": "One 16-base 10x barcode with -gem suffix parsed from original index; always joined with sample identity",
              "post_QC_annotation_scope": "Only cells retained and annotated in author os.h5ad; missing raw cells are not automatically quality failures by our analysis",
              "replication": "Technical culture-well libraries; independent differentiation batch pairing not established"}
    save_json("annotation_join_audit.json", report)
    print(json.dumps({"annotation_rows": len(kept), "libraries": joins}, ensure_ascii=True), flush=True)


def read_features(row):
    filename = INPUTS / row["features_url"].rsplit("/", 1)[1]
    verify_small_input(filename)
    with gzip.open(filename, "rt", encoding="utf-8") as handle:
        values = [line.rstrip("\n").split("\t") for line in handle]
    if any(len(x) != 3 for x in values):
        raise ValueError("Expected gene ID, symbol, feature type")
    if any(x[2] != "Gene Expression" for x in values):
        raise ValueError("Non-gene-expression feature type found; explicit handling required")
    return values


def aggregate(categories):
    import numpy as np
    manifest = read_manifest()
    audit = json.loads((BASE / "annotation_join_audit.json").read_text(encoding="utf-8"))
    if sha256(BASE.parent / "organoid_sample_manifest.csv") != audit["manifest_sha256"]:
        raise ValueError("Sample manifest changed after annotation inspection")
    if sha256(BASE.parent / "organoid_metadata_audit.json") != audit["small_input_audit_sha256"]:
        raise ValueError("Small-file audit changed after annotation inspection")
    annotation_file = BASE / audit["annotation_file"]
    if sha256(annotation_file) != audit["annotation_file_sha256"]:
        raise ValueError("Annotation file hash mismatch")
    with gzip.open(annotation_file, "rt", encoding="utf-8", newline="") as handle:
        annotations = list(csv.DictReader(handle))
    observed = {x["celltype_category"] for x in annotations}
    if not categories or set(categories) - observed:
        raise ValueError(f"Explicit existing category names required; observed: {sorted(observed)}")
    chosen = [x for x in annotations if x["celltype_category"] in categories]
    subtypes = sorted({x["celltype_scTypedb"] for x in chosen})
    if not subtypes or any(not x for x in subtypes):
        raise ValueError("No nonempty selected subtype labels")
    features = read_features(manifest[0])
    if len({x[0] for x in features}) != len(features):
        raise ValueError("Duplicate feature IDs require explicit handling")
    all_columns = []
    count_arrays = []
    matrix_reports = []
    sample_reports = []
    downloaded = json.loads((BASE / "download_log.json").read_text(encoding="utf-8"))
    for row in manifest:
        started = time.monotonic()
        if read_features(row) != features:
            raise ValueError("Feature identity/order differs across libraries")
        barcode_file = INPUTS / row["barcode_url"].rsplit("/", 1)[1]
        verify_small_input(barcode_file)
        with gzip.open(barcode_file, "rt", encoding="utf-8") as handle:
            barcodes = [line.strip() for line in handle]
        subset = [x for x in chosen if x["gsm"] == row["gsm"]]
        lookup = {x["barcode"]: subtypes.index(x["celltype_scTypedb"]) for x in subset}
        if len(lookup) != len(subset):
            raise ValueError("Repeated selected barcode within library")
        subtype_of_cell = np.array([lookup.get(b, -1) for b in barcodes], dtype=np.int16)
        counts = np.zeros((len(features), len(subtypes)), dtype=np.int64)
        matrix_file = SOURCE / row["matrix_url"].rsplit("/", 1)[1]
        expected = next(x for x in downloaded["records"] if x["local_file"].endswith(matrix_file.name))
        actual_sha = sha256(matrix_file)
        if actual_sha != expected["sha256"]:
            raise ValueError("MTX checksum differs from download log")
        entries = 0
        all_umi = 0
        selected_umi = 0
        with gzip.open(matrix_file, "rt", encoding="ascii") as handle:
            header = next(handle).strip()
            if header.lower() != "%%matrixmarket matrix coordinate integer general":
                raise ValueError(f"Unexpected count matrix header: {header}")
            for line in handle:
                if not line.startswith("%"):
                    n_genes, n_cells, declared_nnz = map(int, line.split())
                    break
            if n_genes != len(features) or n_cells != len(barcodes):
                raise ValueError("Matrix dimensions disagree with features/barcodes")
            for line in handle:
                gene, cell, value = map(int, line.split())
                if not (1 <= gene <= n_genes and 1 <= cell <= n_cells and value >= 0):
                    raise ValueError("Invalid count coordinate/value")
                entries += 1
                all_umi += value
                subtype = int(subtype_of_cell[cell-1])
                if subtype >= 0:
                    counts[gene-1, subtype] += value
                    selected_umi += value
        if entries != declared_nnz or int(counts.sum()) != selected_umi or (counts < 0).any():
            raise ValueError("Count aggregation integrity check failed")
        all_hep = counts.sum(axis=1, keepdims=True)
        counts = np.concatenate([counts, all_hep], axis=1)
        labels = subtypes + ["ALL_HEPATOCYTE_LIKE"]
        for j, label in enumerate(labels):
            column = row["gsm"]+"|"+label
            cell_count = len(subset) if label == "ALL_HEPATOCYTE_LIKE" else sum(x["celltype_scTypedb"] == label for x in subset)
            all_columns.append(column)
            sample_reports.append({"count_column": column, "gsm": row["gsm"],
                                   "sample": row["author_sample_id"], "group": row["group"],
                                   "subtype": label, "cells_aggregated": cell_count,
                                   "total_counts": int(counts[:, j].sum()),
                                   "biological_independent_n": "not established; technical library"})
        count_arrays.append(counts)
        matrix_reports.append({"gsm": row["gsm"], "matrix_sha256": actual_sha,
                               "features": n_genes, "raw_cells": n_cells,
                               "declared_nnz": declared_nnz, "read_entries": entries,
                               "all_raw_matrix_counts": all_umi, "selected_hepatocyte_counts": selected_umi,
                               "selected_hepatocyte_cells": len(subset),
                               "elapsed_seconds": round(time.monotonic()-started, 2)})
        print(json.dumps(matrix_reports[-1], ensure_ascii=True), flush=True)
    combined = np.concatenate(count_arrays, axis=1)
    counts_file = BASE / "all_gene_library_subtype_counts.tsv.gz"
    with gzip.open(counts_file, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_id", "gene_symbol", "feature_type"]+all_columns)
        for feature, values in zip(features, combined):
            writer.writerow(feature+values.tolist())
    metadata_file = BASE / "library_subtype_summary.csv"
    with metadata_file.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sample_reports[0]))
        writer.writeheader()
        writer.writerows(sample_reports)
    protocol_file = BASE.parent / "pilot_protocol.json"
    protocol = json.loads(protocol_file.read_text(encoding="utf-8"))
    panel = protocol["fixed_literature_panel"]
    panel_rows = []
    for gene in panel:
        indices = [i for i, feature in enumerate(features) if feature[1] == gene]
        gene_ids = ";".join(features[i][0] for i in indices)
        for j, sample in enumerate(sample_reports):
            count = int(combined[indices, j].sum()) if indices else None
            total = sample["total_counts"]
            no_cells = sample["cells_aggregated"] == 0
            status = "no_selected_cells" if no_cells else ("missing_feature" if not indices else ("measured_zero" if count == 0 else "observed_nonzero"))
            panel_rows.append({**sample, "gene_symbol": gene, "gene_ids": gene_ids,
                               "matched_feature_count": len(indices),
                               "coverage": "present_in_features" if indices else "not_in_feature_annotation",
                               "integer_counts": count,
                               "cpm": count/total*1_000_000 if count is not None and total and not no_cells else None,
                               "zero_vs_missing": status})
    panel_file = BASE / "fixed_panel_counts_cpm.csv"
    with panel_file.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(panel_rows[0]))
        writer.writeheader()
        writer.writerows(panel_rows)
    output_files = [counts_file, metadata_file, panel_file]
    report = {**provenance(), "numpy_version": np.__version__,
              "selected_author_categories": categories, "subtype_labels": subtypes,
              "subtype_label_note": "Original object labels are retained. Author sc_constants.py maps FH2 to display abbreviation cAH (cycling adult hepatocyte-like); do not read raw FH2 as an independently confirmed fetal subtype. All labels are author assignments, not our validation.",
              "all_gene_count_rows": len(features), "count_columns": len(all_columns),
              "libraries": matrix_reports, "panel": panel,
              "panel_basis": "Previously fixed known vitamin-D-associated panel; not selected by organoid results; no expected PA direction imposed",
              "panel_protocol_sha256": sha256(protocol_file),
              "cpm_definition": "Counts summed over matching gene IDs divided by all measured gene counts in the same library/subtype, times 1e6. ALL_HEPATOCYTE_LIKE is a pooled UMI-weighted mixture, not mean expression per cell. No group-level normalization or inferential test.",
              "outputs": [{"file": x.name, "bytes": x.stat().st_size, "sha256": sha256(x)} for x in output_files],
              "limitations": ["Technical replicates; no biological n=2 claim", "Author retained/annotated cells only; no independent cell type reannotation", "Some subtypes have very few cells: HB1 has 4/69/9/128 cells across the four libraries; subtype CPM may be unstable", "Processed annotation labels must be interpreted with 2024 figure correction", "GEO five-day and corrected-paper four-day timing conflict remains", "Disease stimulus reference does not validate an ingredient effect", "Cell composition changes and cell-state differences are not causal or therapeutic evidence"]}
    save_json("aggregation_audit.json", report)
    print(json.dumps({"saved_gene_rows": len(features), "saved_count_columns": len(all_columns), "panel_rows": len(panel_rows)}, ensure_ascii=True), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["decompress", "inspect", "aggregate"])
    parser.add_argument("--category", action="append", default=[])
    args = parser.parse_args()
    if args.mode == "decompress":
        decompress()
    elif args.mode == "inspect":
        inspect_annotations()
    else:
        aggregate(args.category)


if __name__ == "__main__":
    main()
