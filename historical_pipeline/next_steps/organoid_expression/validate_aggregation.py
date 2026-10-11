"""Independent sparse-matrix check of streamed organoid aggregation.

Use scipy.io.mmread on one control and one PA library, then compare every
gene for all subtype columns. This is a computational check, not biological
validation of the organoid model, annotations, or ingredient effects.
"""
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.io import mmread

BASE = Path(__file__).resolve().parent


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    audit = json.loads((BASE / "aggregation_audit.json").read_text(encoding="utf-8"))
    join = json.loads((BASE / "annotation_join_audit.json").read_text(encoding="utf-8"))
    for result in audit["outputs"]:
        assert digest(BASE / result["file"]) == result["sha256"]
    assert audit["script_sha256"] == digest(BASE / "link_expression.py")
    assert join["script_sha256"] == audit["script_sha256"]
    assert all(x["condition_labels_match_manifest"] for x in join["sample_join_reports"])
    counts = pd.read_csv(BASE / "all_gene_library_subtype_counts.tsv.gz", sep="\t")
    summary = pd.read_csv(BASE / "library_subtype_summary.csv")
    panel = pd.read_csv(BASE / "fixed_panel_counts_cpm.csv", keep_default_na=False)
    annotations = pd.read_csv(BASE / "four_library_annotations.csv.gz", keep_default_na=False)
    manifest = pd.read_csv(BASE.parent / "organoid_sample_manifest.csv", keep_default_na=False)
    numeric = counts.iloc[:, 3:].to_numpy()
    assert numeric.dtype.kind in "iu" and (numeric >= 0).all()
    assert len(counts) == 33538 and not counts.gene_id.duplicated().any()
    assert not annotations.duplicated(["sample", "barcode"]).any()
    per_barcode = Counter(annotations.barcode)
    cross_library_reused = sum(n > 1 for n in per_barcode.values())
    independent_checks = []
    for gsm in ["GSM6322577", "GSM6322579"]:
        row = manifest.set_index("gsm").loc[gsm]
        with gzip.open(BASE / "sources" / row.matrix_url.rsplit("/", 1)[1], "rb") as handle:
            sparse = mmread(handle).tocsr()
        with gzip.open(BASE.parent / "organoid_inputs" / row.barcode_url.rsplit("/", 1)[1], "rt") as handle:
            raw_barcodes = [x.strip() for x in handle]
        lookup = {value: i for i, value in enumerate(raw_barcodes)}
        selected = annotations[(annotations.gsm == gsm) & (annotations.celltype_category.isin(audit["selected_author_categories"]))]
        checks = []
        for subtype in audit["subtype_labels"] + ["ALL_HEPATOCYTE_LIKE"]:
            subset = selected if subtype == "ALL_HEPATOCYTE_LIKE" else selected[selected.celltype_scTypedb == subtype]
            indices = [lookup[x] for x in subset.barcode]
            expected = np.asarray(sparse[:, indices].sum(axis=1)).ravel()
            actual = counts[gsm+"|"+subtype].to_numpy()
            assert np.array_equal(actual, expected)
            checks.append({"subtype": subtype, "genes_compared": len(actual), "cells": len(indices), "exact_match": True})
        independent_checks.append({"gsm": gsm, "reader": "scipy.io.mmread to CSR, independent of streamed entry accumulator", "checks": checks})
        del sparse
    for gsm in manifest.gsm:
        subtype_columns = [gsm+"|"+s for s in audit["subtype_labels"]]
        assert np.array_equal(counts[subtype_columns].sum(axis=1).to_numpy(), counts[gsm+"|ALL_HEPATOCYTE_LIKE"].to_numpy())
    for row in summary.itertuples():
        assert int(counts[row.count_column].sum()) == row.total_counts
    for row in panel.itertuples():
        present = counts.gene_symbol == row.gene_symbol
        if int(row.cells_aggregated) == 0:
            assert row.zero_vs_missing == "no_selected_cells" and row.cpm == ""
        elif not present.any():
            assert row.coverage == "not_in_feature_annotation" and row.cpm == ""
        else:
            raw = int(counts.loc[present, row.count_column].sum())
            assert int(row.integer_counts) == raw
            expected_cpm = raw/int(row.total_counts)*1_000_000
            assert np.isclose(float(row.cpm), expected_cpm, rtol=1e-12, atol=1e-12)
            assert row.zero_vs_missing == ("measured_zero" if raw == 0 else "observed_nonzero")
    all_hep = summary[summary.subtype == "ALL_HEPATOCYTE_LIKE"]
    panel_all = panel[panel.subtype == "ALL_HEPATOCYTE_LIKE"]
    report = {"executed_utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
              "interpreter": sys.executable, "script_sha256": digest(Path(__file__)),
              "aggregation_script_sha256": audit["script_sha256"],
              "versions": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
              "status": "passed", "all_output_hashes_match": True,
              "count_shape": [len(counts), numeric.shape[1]],
              "cross_library_reused_barcode_strings": cross_library_reused,
              "joint_key_duplicates": 0,
              "independent_sparse_matrix_checks": independent_checks,
              "all_four_library_subtype_sums_match_pooled_counts": True,
              "all_panel_counts_cpm_zero_missing_checks": True,
              "selected_hepatocyte_cells": int(all_hep.cells_aggregated.sum()),
              "fixed_panel_feature_coverage": panel_all.groupby("gene_symbol").coverage.unique().map(list).to_dict(),
              "panel_zero_status_counts": Counter(panel.zero_vs_missing),
              "biological_inference": "None. Data-connection and arithmetic checks only; technical replicates and disease-reference limitations remain."}
    (BASE / "validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ["status", "count_shape", "cross_library_reused_barcode_strings", "selected_hepatocyte_cells", "panel_zero_status_counts"]}, ensure_ascii=True))


if __name__ == "__main__":
    main()
