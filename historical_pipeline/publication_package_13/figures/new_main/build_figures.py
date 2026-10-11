"""Two additional descriptive figures from unchanged prior-analysis tables."""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import ScaledTranslation
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[3]
SKILLS = Path.home() / ".agents/skills/nature-figure/scripts"
sys.path.insert(0, str(SKILLS))
from audit_panel_alignment import require_matplotlib_panel_alignment

INPUTS = {
    "human": PROJECT / "next_steps/followup_05/bulk_candidate_samples.csv",
    "frozen_facts": PROJECT / "next_steps/manuscript_draft_11/data/facts.json",
    "organoid_subtypes": PROJECT / "next_steps/followup_05/organoid_within_subtype.csv",
    "organoid_umi": PROJECT / "next_steps/followup_05/organoid_library_UMI_fractions.csv",
    "organoid_library_metadata": PROJECT / "next_steps/organoid_expression/library_subtype_summary.csv",
}
INK = "#263443"
DIV = LinearSegmentedColormap.from_list("Directional", ["#216A85", "#FAFAF8", "#B56D32"])
SEQ = LinearSegmentedColormap.from_list("UMIshare", ["#F6FAFB", "#2E758D"])
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7.5, "axes.titlesize": 9, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "text.color": INK,
    "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.edgecolor": "#697581", "axes.linewidth": .6,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(name, rows):
    with (HERE / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(name, payload):
    (HERE / name).write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def plain_heatmap(ax, xlabels, ylabels):
    ax.set_xticks(range(len(xlabels)), xlabels)
    ax.set_yticks(range(len(ylabels)), ylabels)
    ax.tick_params(axis="both", length=0, pad=5)
    ax.set_xticks(np.arange(-.5, len(xlabels), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(ylabels), 1), minor=True)
    ax.tick_params(which="minor", length=0)
    ax.grid(which="minor", color="white", linewidth=.7)
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_letter(ax, letter):
    ax.text(0, 1, letter, transform=ax.transAxes + ScaledTranslation(-28 / 72, 11 / 72, ax.figure.dpi_scale_trans),
            ha="left", va="bottom", fontweight="bold", fontsize=10)


def export(fig, stem, exclude_axes):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=HERE / f"{stem}.alignment.json",
        overlay_svg=HERE / f"{stem}.alignment.svg", tolerance_pt=1.5, gutter_tolerance_pt=1.5,
        require_panel_labels=(len(fig.axes) - len(exclude_axes)) > 1, strict=True, exclude_axes=exclude_axes)
    fig.savefig(HERE / f"{stem}.pdf")
    fig.savefig(HERE / f"{stem}.svg")
    fig.savefig(HERE / f"{stem}.png", dpi=600)
    fig.savefig(HERE / f"{stem}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


human = read_csv(INPUTS["human"])
facts = json.loads(INPUTS["frozen_facts"].read_text(encoding="utf-8"))
frozen = facts["human_candidates_all25_verbatim"]
assert len(frozen) == 25 and len(human) == 225
gene_ids = [r["gene_id"] for r in frozen]
genes = [r["symbol"] for r in frozen]
samples = [f"{group}_{replicate}" for group in ["Control", "Model", "Treat"] for replicate in [1, 2, 3]]
lookup = {(r["gene_id"], r["sample"]): r for r in human}
assert len(lookup) == 225
assert set(lookup) == {(g, s) for g in gene_ids for s in samples}
expected_group = {"Control": "Control", "Model": "OA", "Treat": "OA + sesamin"}
counts = np.array([[float(lookup[g, s]["normalized_count"]) for s in samples] for g in gene_ids])
assert np.isfinite(counts).all() and (counts >= 0).all()
logs = np.log2(counts + 1)
row_centers = logs.mean(axis=1, keepdims=True)
centered = logs - row_centers
assert np.max(np.abs(centered.mean(axis=1))) < 1e-12
assert np.allclose(centered + row_centers, logs, rtol=0, atol=1e-14)
# Independent scalar implementation verifies every derived transform.
scalar = [[math.log2(float(lookup[g, s]["normalized_count"]) + 1) for s in samples] for g in gene_ids]
scalar_centered = np.array([[x - sum(row) / 9 for x in row] for row in scalar])
assert np.allclose(centered, scalar_centered, rtol=0, atol=3e-15)
bound = math.ceil(float(np.abs(centered).max()) * 2) / 2
assert bound >= np.abs(centered).max()
human_source = []
for i, (gene_id, gene) in enumerate(zip(gene_ids, genes)):
    for j, sample in enumerate(samples):
        row = lookup[gene_id, sample]
        assert row["symbol"] == gene and row["group"] == expected_group[sample.split("_")[0]]
        human_source.append({"row_order": i + 1, "column_order": j + 1, "gene_id": gene_id, "symbol": gene,
                             "sample": sample, "source_group": row["group"], "raw_count_verbatim": row["raw_count"],
                             "normalized_count_verbatim": row["normalized_count"], "log2_count_plus1": logs[i,j],
                             "gene_mean_log2_count_plus1": row_centers[i,0], "centered_log2_count_plus1": centered[i,j]})
write_csv("fig2_human_all25_source.csv", human_source)
fig, ax = plt.subplots(figsize=(180 / 25.4, 190 / 25.4))
fig.subplots_adjust(left=.285, right=.865, bottom=.185, top=.85)
im = ax.imshow(centered, cmap=DIV, vmin=-bound, vmax=bound, interpolation="none", aspect="auto")
plain_heatmap(ax, ["1", "2", "3"] * 3, genes)
ax.xaxis.tick_top()
for label in ax.get_yticklabels():
    label.set_fontstyle("italic")
for boundary in [2.5, 5.5]:
    ax.axvline(boundary, color="white", linewidth=2.0)
for x, label in zip([1/6, .5, 5/6], ["Control", "Model", "Treat"]):
    ax.text(x, 1.053, label, transform=ax.transAxes, ha="center", va="bottom", fontsize=8, fontweight="bold")
fig.text(.105, .956, "Complete frozen human candidate set", fontsize=10, fontweight="bold")
fig.text(.105, .927, "25 genes × 9 submitted HepG2 libraries; original candidate order", fontsize=7.5)
cax = fig.add_axes([.89, .35, .024, .33])
colorbar = fig.colorbar(im, cax=cax, ticks=[-bound, 0, bound])
colorbar.outline.set_visible(False)
colorbar.ax.tick_params(length=2, labelsize=7)
colorbar.set_label("Gene-centered log2(count + 1)", fontsize=7.5, labelpad=8)
fig.text(.105, .131, "Each row is centered across its nine libraries after log2(normalized count + 1).", fontsize=7)
fig.text(.105, .102, "Color spans every value without clipping; it does not show absolute expression across genes.", fontsize=7)
fig.text(.105, .064, "Candidates were selected using these group patterns: the heatmap is descriptive,", fontsize=7)
fig.text(.105, .039, "not independent recovery validation or a differential-expression test.", fontsize=7)
export(fig, "fig2_human_all25", exclude_axes=[cax])

subtypes = ["AH", "FH1", "FH2", "HB1", "HB2"]
focal = ["VNN1", "AREG", "JAG2", "CXCL8"]
subtype_rows = read_csv(INPUTS["organoid_subtypes"])
fraction_rows = read_csv(INPUTS["organoid_umi"])
library_rows = read_csv(INPUTS["organoid_library_metadata"])
assert len(subtype_rows) == 125 and len(fraction_rows) == 4
subtype_lookup = {(r["symbol"], r["subtype"]): r for r in subtype_rows}
selected = [subtype_lookup[g, st] for g in focal for st in subtypes]
assert all(r["feature_covered"] == "True" for r in selected)
values = np.array([[float(subtype_lookup[g, st]["log2_CPM_ratio_plus1"]) for st in subtypes] for g in focal])
assert np.isfinite(values).all()
for row in selected:
    recalculated = math.log2((float(row["PA_mean_CPM"]) + 1) / (float(row["control_mean_CPM"]) + 1))
    assert abs(recalculated - float(row["log2_CPM_ratio_plus1"])) < 2e-15
assert np.abs(values).max() < 2.5
libraries = [f"GSM63225{n}" for n in [77, 78, 79, 80]]
fraction_lookup = {r["gsm"]: r for r in fraction_rows}
metadata_lookup = {(r["gsm"], r["subtype"]): r for r in library_rows}
fractions = np.array([[float(fraction_lookup[gsm][st]) for st in subtypes] for gsm in libraries])
assert np.isfinite(fractions).all() and (fractions >= 0).all() and (fractions <= 1).all()
assert np.allclose(fractions.sum(axis=1), 1, atol=1e-14, rtol=0)
percentages = fractions * 100
assert percentages.max() < 55
fraction_source = []
for i, gsm in enumerate(libraries):
    expected = "CTRL-PA" if i < 2 else "PA"
    total = sum(int(metadata_lookup[gsm, st]["total_counts"]) for st in subtypes)
    for j, st in enumerate(subtypes):
        meta = metadata_lookup[gsm, st]
        assert meta["group"] == expected
        assert abs(int(meta["total_counts"]) / total - fractions[i,j]) < 1e-14
        fraction_source.append({"gsm": gsm, "group": meta["group"], "subtype": st,
                                "UMI_fraction_verbatim": fraction_lookup[gsm][st], "UMI_percent": percentages[i,j],
                                "cells_aggregated": meta["cells_aggregated"], "total_counts": meta["total_counts"],
                                "biological_independent_n": meta["biological_independent_n"]})
minimum_cells = {st: min(int(metadata_lookup[gsm, st]["cells_aggregated"]) for gsm in libraries) for st in subtypes}
assert minimum_cells == {"AH": 1430, "FH1": 505, "FH2": 52, "HB1": 4, "HB2": 414}
for row in selected:
    assert int(row["minimum_cells"]) == minimum_cells[row["subtype"]]
    assert (row["any_library_below50_cells"] == "True") == (minimum_cells[row["subtype"]] < 50)
write_csv("fig6_subtype_contrasts_source.csv", selected)
write_csv("fig6_library_UMI_source.csv", fraction_source)
fig, axes = plt.subplots(1, 2, figsize=(180 / 25.4, 111 / 25.4))
fig.subplots_adjust(left=.115, right=.975, bottom=.43, top=.765, wspace=.42)
images = []
for index, (ax, matrix) in enumerate(zip(axes, [values, percentages])):
    im = ax.imshow(matrix, cmap=DIV if index == 0 else SEQ, vmin=-2.5 if index == 0 else 0,
                   vmax=2.5 if index == 0 else 55, interpolation="none", aspect="auto")
    images.append(im)
    plain_heatmap(ax, ["AH", "FH1", "FH2", "HB1*", "HB2"], focal if index == 0 else ["Control 1", "Control 2", "PA 1", "PA 2"])
    for i in range(4):
        for j in range(5):
            value = matrix[i,j]
            text_color = "white" if (abs(value) > 1.6 if index == 0 else value > 35) else INK
            ax.text(j, i, f"{value:+.2f}" if index == 0 else f"{value:.2f}", ha="center", va="center", fontsize=7, color=text_color)
    ax.get_xticklabels()[3].set_fontweight("bold")
    if index == 0:
        for label in ax.get_yticklabels():
            label.set_fontstyle("italic")
    panel_letter(ax, "ab"[index])
    ax.set_title(["Within-subtype PA / control", "Subtype UMI share by library"][index], fontsize=8.5, pad=14)
color_axes = []
for index, ax in enumerate(axes):
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0, .336, pos.width, .022])
    color_axes.append(cax)
    colorbar = fig.colorbar(images[index], cax=cax, orientation="horizontal", ticks=[-2.5, 0, 2.5] if index == 0 else [0, 25, 55])
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(length=2, pad=2, labelsize=7)
    colorbar.set_label("log2[(PA mean CPM + 1) / (control mean CPM + 1)]" if index == 0 else "Share of UMI across the five supplied labels (%)", fontsize=6.6, labelpad=4)
fig.text(.115, .95, "Organoid contrasts vary across subtype labels", fontsize=9.5, fontweight="bold")
fig.text(.115, .902, "Author labels retained; HB1* is sparse (minimum 4 cells per library).", fontsize=7.5)
fig.text(.115, .177, "HB1*: 4 / 69 / 9 / 128 cells in Control 1 / Control 2 / PA 1 / PA 2; two libraries <50 cells.", fontsize=7)
fig.text(.115, .126, "Four technical libraries from one differentiation; no sesamin intervention and no inferential tests.", fontsize=7)
fig.text(.115, .075, "UMI shares describe molecular-count composition, not cell proportions or independent replicates.", fontsize=7)
export(fig, "fig6_organoid_subtypes", exclude_axes=color_axes)

manifest = [{"key": key, "path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size} for key, path in INPUTS.items()]
checks = {
    "human_values": 225, "human_NA": 0, "human_color_limits": [-bound, bound],
    "human_centered_min": float(centered.min()), "human_centered_max": float(centered.max()),
    "human_centered_row_mean_max_abs": float(np.abs(centered.mean(axis=1)).max()),
    "human_scalar_transform_max_abs_difference": float(np.abs(centered - scalar_centered).max()),
    "organoid_contrasts": 20, "organoid_contrast_NA": 0, "organoid_UMI_values": 20,
    "organoid_UMI_fraction_sum_max_abs_error": float(np.abs(fractions.sum(axis=1) - 1).max()),
    "minimum_cells_by_subtype": minimum_cells, "sparse_subtypes": [st for st in subtypes if minimum_cells[st] < 50],
    "HB1_cells": [int(metadata_lookup[gsm, "HB1"]["cells_aggregated"]) for gsm in libraries],
    "subplot_a_color_limits": [-2.5, 2.5], "subplot_b_color_limits_percent": [0, 55],
}
write_json("source_data.json", {"date": "2026-10-09", "inputs": manifest, "checks": checks,
    "human_plot_values": human_source, "organoid_plot_contrasts": selected, "organoid_plot_UMI": fraction_source,
    "environment": {"python": sys.version, "executable": sys.executable, "matplotlib": matplotlib.__version__, "numpy": np.__version__},
    "plotting_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "reanalysis_scope": "Only explicit display transformations and consistency checks on existing tables; no full matrix or sequencing reanalysis, no statistical inference.",
    "randomness": "None; original frozen order, no clustering or random sampling."})
print(json.dumps(checks))
