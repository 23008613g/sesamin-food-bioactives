"""Validate actual aggregated outputs and plot a descriptive organoid reference.

No inferential statistics, raw-source edits, or model calls. Requires project venv.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "organoid_expression"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qa-scripts", required=True, type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.qa_scripts))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    audit = json.loads((DATA / "aggregation_audit.json").read_text(encoding="utf-8"))
    annotation = json.loads((DATA / "annotation_join_audit.json").read_text(encoding="utf-8"))
    assert audit["script_sha256"] == sha(DATA / "link_expression.py")
    assert audit["no_statistical_inference"] and audit["no_ingredient_effect_validation"]
    for record in audit["outputs"]:
        assert sha(DATA / record["file"]) == record["sha256"], record["file"]
    counts = pd.read_csv(DATA / "all_gene_library_subtype_counts.tsv.gz", sep="\t", keep_default_na=False)
    meta = pd.read_csv(DATA / "library_subtype_summary.csv", keep_default_na=False)
    panel = pd.read_csv(DATA / "fixed_panel_counts_cpm.csv", keep_default_na=False)
    assert counts.shape == (33538, 27) and counts.gene_id.is_unique
    matrix = counts.iloc[:, 3:]
    assert matrix.columns.tolist() == meta.count_column.tolist()
    assert (matrix.dtypes.apply(lambda x: np.issubdtype(x, np.integer))).all()
    assert (matrix.to_numpy() >= 0).all()
    np.testing.assert_array_equal(matrix.sum(axis=0), meta.total_counts)
    gsms = [row["gsm"] for row in audit["libraries"]]
    subtypes = audit["subtype_labels"]
    genes = audit["panel"]
    summaries = []
    for gsm in gsms:
        group = meta[meta.gsm == gsm].set_index("subtype")
        parts = [gsm + "|" + s for s in subtypes]
        total_column = gsm + "|ALL_HEPATOCYTE_LIKE"
        np.testing.assert_array_equal(matrix[parts].sum(axis=1), matrix[total_column])
        assert group.loc[subtypes, "cells_aggregated"].sum() == group.loc["ALL_HEPATOCYTE_LIKE", "cells_aggregated"]
        ann = next(x for x in annotation["sample_join_reports"] if x["gsm"] == gsm)
        assert ann["unmatched_annotation_barcodes"] == 0 and ann["condition_labels_match_manifest"]
        assert ann["category_counts"]["Hepatocytes"] == group.loc["ALL_HEPATOCYTE_LIKE", "cells_aggregated"]
        for subtype in subtypes:
            assert ann["subtype_counts"][subtype] == group.loc[subtype, "cells_aggregated"]
        summaries.append({"gsm": gsm, "group": group.group.iloc[0],
                          "annotated_cells": ann["annotated_cells"],
                          "hepatocyte_cells": int(group.loc["ALL_HEPATOCYTE_LIKE", "cells_aggregated"]),
                          "subtype_cells": {s: int(group.loc[s, "cells_aggregated"]) for s in subtypes}})
    assert len(panel) == 7 * 24
    for row in panel.itertuples(index=False):
        selected = counts[counts.gene_symbol == row.gene_symbol]
        assert len(selected) == row.matched_feature_count == 1
        actual = int(selected[row.count_column].sum())
        assert actual == row.integer_counts
        assert row.total_counts == matrix[row.count_column].sum()
        assert row.cells_aggregated > 0
        assert row.zero_vs_missing == ("measured_zero" if actual == 0 else "observed_nonzero")
        assert np.isclose(row.cpm, actual / row.total_counts * 1e6, rtol=1e-12, atol=1e-12)
    pooled = panel[panel.subtype == "ALL_HEPATOCYTE_LIKE"].pivot(index="gene_symbol", columns="gsm", values="cpm").loc[genes, gsms]
    pooled.to_csv(HERE / "organoid_panel_pooled_cpm.csv")
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 7, "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7})
    fig, axes = plt.subplots(1, 2, figsize=(7.2047244094, 4.7244094488))
    fig.subplots_adjust(left=.085, right=.965, bottom=.36, top=.86, wspace=.55)
    labels = ["CTRL-62", "CTRL-63", "PA-64", "PA-65"]
    palette = ["#567F99", "#C49A6C", "#879F78", "#927CA0", "#BC7A77"]
    bottom = np.zeros(4, dtype=int)
    for subtype, color in zip(subtypes, palette, strict=True):
        values = [r["subtype_cells"][subtype] for r in summaries]
        axes[0].bar(np.arange(4), values, bottom=bottom, width=.66, label=subtype, color=color,
                    edgecolor="white", linewidth=.3)
        bottom += values
    for x, n in enumerate(bottom):
        axes[0].text(x, n + 100, f"{n:,}", ha="center", va="bottom", fontsize=7)
    axes[0].set_ylim(0, 6700)
    axes[0].set_ylabel("Author-annotated hepatocyte-like cells")
    axes[0].set_xticks(range(4), labels, rotation=35, ha="right")
    axes[0].set_title("a  Cells by author subtype", loc="left", fontweight="bold", pad=12)
    axes[0].legend(loc="upper center", bbox_to_anchor=(.5, -.30), ncol=3, frameon=False,
                   handlelength=1.1, columnspacing=.8)
    colored = np.log10(1 + pooled.to_numpy(dtype=float))
    im = axes[1].imshow(colored, cmap="viridis", vmin=0, vmax=max(float(colored.max()), 1), aspect="auto")
    axes[1].set_xticks(range(4), labels, rotation=35, ha="right")
    axes[1].set_yticks(range(7), genes)
    axes[1].set_title("b  Fixed RNA panel (pooled)", loc="left", fontweight="bold", pad=12)
    for i in range(7):
        for j in range(4):
            axes[1].text(j, i, f"{pooled.iloc[i, j]:.2g}", ha="center", va="center", fontsize=6.5,
                         color="black" if colored[i, j] > .62 * im.norm.vmax else "white")
    cax = fig.add_axes([.625, .198, .30, .025])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Color: log10(1 + CPM); numbers: CPM", fontsize=7, labelpad=4)
    cax.tick_params(labelsize=6)
    fig.text(.5, .102, "Four technical libraries from one differentiation experiment; no biological n = 2 inference.", ha="center", fontsize=7)
    fig.text(.5, .063, "CPM pools all selected hepatocyte-like cells; differences can reflect cell mixture and cell state.", ha="center", fontsize=7)
    fig.text(.5, .024, "Disease-stimulus reference only. No food ingredient treatment and no efficacy test.", ha="center", fontsize=7)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=HERE / "organoid_reference.alignment.json",
                                      tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True, exclude_axes=[cax])
    fig.savefig(HERE / "organoid_reference.pdf", facecolor="white")
    fig.savefig(HERE / "organoid_reference.svg", facecolor="white")
    fig.savefig(HERE / "organoid_reference.png", dpi=300, facecolor="white")
    plt.close(fig)
    record = {"executed_utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
              "interpreter": sys.executable, "script_sha256": sha(Path(__file__)),
              "aggregation_audit_sha256": sha(DATA / "aggregation_audit.json"),
              "annotation_audit_sha256": sha(DATA / "annotation_join_audit.json"),
              "checks": {"all_feature_totals_match": True, "subtypes_sum_to_pooled_for_every_gene": True,
                         "annotation_and_aggregated_cells_match": True, "all_168_panel_counts_and_CPM_match": True,
                         "annotation_join_unmatched": 0},
              "libraries": summaries, "all_feature_rows": len(counts), "count_columns": len(matrix.columns),
              "panel_measurements": len(panel), "pooled_panel_cpm": pooled.to_dict(),
              "scope": "Checks aggregated artifacts against annotation totals and all-gene count table, not a second independent scan of all raw MTX entries",
              "interpretation": "Descriptive technical libraries only; pooled CPM is UMI-weighted mixture, not mean single-cell expression; no ingredient/AI validation",
              "outputs": [{"file": name, "sha256": sha(HERE / name)} for name in [
                  "organoid_reference.png", "organoid_reference.pdf", "organoid_reference.svg", "organoid_panel_pooled_cpm.csv"]]}
    (HERE / "organoid_output_validation.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record))


if __name__ == "__main__":
    main()
