"""Render three descriptive figures exclusively from archived observed results."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from matplotlib.transforms import ScaledTranslation
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[3]
SKILLS = Path.home() / ".agents/skills/nature-figure/scripts"
sys.path.insert(0, str(SKILLS))
from audit_panel_alignment import require_matplotlib_panel_alignment

SOURCES = {
    "human_library_counts": PROJECT / "next_steps/followup_05/bulk_candidate_samples.csv",
    "mouse_library_tpm": Path("D:/ResearchData/sesamin_mouse_2026-10-08/matrices/candidate_mouse_descriptive.tsv"),
    "archived_facts": PROJECT / "next_steps/manuscript_draft_11/data/facts.json",
    "phh_label_effects": PROJECT / "next_steps/followup_05/phh_donor_effects.csv",
    "organoid_sensitivity": PROJECT / "next_steps/followup_05/organoid_composition_sensitivity.csv",
}
FOCAL = ["VNN1", "AREG", "JAG2"]
BACKGROUND = [*FOCAL, "CXCL8"]
PHH_LABELS = ["F95", "M95", "H20", "H961", "F95P"]
GROUP_COLORS = ["#737B86", "#CA854F", "#3C8793"]
INK = "#263443"
DIVERGE = LinearSegmentedColormap.from_list("Directional", ["#216A85", "#FAFAF8", "#B56D32"])

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7.5, "axes.titlesize": 9, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.edgecolor": "#697581", "axes.linewidth": .65,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False,
    "savefig.facecolor": "white", "figure.facecolor": "white",
    "lines.linewidth": .9,
})


def read_csv(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_json(name, value):
    (HERE / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(name, rows):
    with (HERE / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def panel_letter(ax, letter, x_offset=-27, y_offset=10):
    ax.text(0, 1, letter, transform=ax.transAxes + ScaledTranslation(x_offset / 72, y_offset / 72, ax.figure.dpi_scale_trans),
            ha="left", va="bottom", fontweight="bold", fontsize=10)


def export(fig, stem, exclude_axes=(), **alignment_options):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig, json_out=HERE / f"{stem}.alignment.json",
        overlay_svg=HERE / f"{stem}.alignment.svg",
        tolerance_pt=1.5, gutter_tolerance_pt=1.5,
        require_panel_labels=len(fig.axes) > 1, strict=True,
        exclude_axes=list(exclude_axes), **alignment_options,
    )
    # No tight crop: preserve the explicit 180 mm submission width.
    fig.savefig(HERE / f"{stem}.pdf")
    fig.savefig(HERE / f"{stem}.svg")
    fig.savefig(HERE / f"{stem}.png", dpi=600)
    fig.savefig(HERE / f"{stem}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


human = read_csv(SOURCES["human_library_counts"])
mouse = read_csv(SOURCES["mouse_library_tpm"], "\t")
facts = json.loads(SOURCES["archived_facts"].read_text(encoding="utf-8"))
phh = read_csv(SOURCES["phh_label_effects"])
organoid = read_csv(SOURCES["organoid_sensitivity"])
assert len(human) == 225 and len(phh) == 125 and len(organoid) == 25
assert len(mouse) == 26

# Figure 2: all library observations for the three explicitly selected candidates.
fig, axes = plt.subplots(2, 3, figsize=(180 / 25.4, 132 / 25.4))
fig.subplots_adjust(left=.13, right=.975, bottom=.17, top=.86, hspace=.88, wspace=.48)
fig.text(.13, .957, "Human HepG2 libraries", fontsize=9, fontweight="bold")
fig.text(.13, .925, "Median-ratio-normalized counts", fontsize=7.5)
fig.text(.13, .491, "Mouse libraries (provisional groups)", fontsize=9, fontweight="bold")
fig.text(.13, .460, "Transcripts per million (TPM)", fontsize=7.5)
expression_rows = []
mean_rows = []
for row in range(2):
    for col, gene in enumerate(FOCAL):
        ax = axes[row, col]
        panel_letter(ax, "abcdef"[row * 3 + col], x_offset=-28, y_offset=5)
        ax.set_title(gene, fontstyle="italic", pad=7)
        vals_by_group = []
        if row == 0:
            labels = ["Control", "Model", "Treat"]
            input_groups = ["Control", "OA", "OA + sesamin"]
            for group, display_group in zip(input_groups, labels):
                selected = sorted([r for r in human if r["symbol"] == gene and r["group"] == group], key=lambda r: r["sample"])
                assert len(selected) == 3
                vals = [float(r["normalized_count"]) for r in selected]
                vals_by_group.append(vals)
                expression_rows.extend({"panel": "abcdef"[row * 3 + col], "species": "human", "gene": gene,
                                        "sample": r["sample"], "source_group": group, "display_group": display_group,
                                        "unit": "median_ratio_normalized_count", "value_verbatim": r["normalized_count"]} for r in selected)
        else:
            labels = ["ctrl", "HFD", "sess"]
            selected = [r for r in mouse if r["human_symbol"] == gene]
            assert len(selected) == 1
            selected = selected[0]
            for group_index, group in enumerate(labels):
                accessions = [f"DRR{n}" for n in range(1010532 + 3 * group_index, 1010535 + 3 * group_index)]
                vals = [float(selected[f"{accession}_TPM"]) for accession in accessions]
                vals_by_group.append(vals)
                expression_rows.extend({"panel": "abcdef"[row * 3 + col], "species": "mouse", "gene": gene,
                                        "sample": accession, "source_group": group, "display_group": group,
                                        "unit": "TPM", "value_verbatim": selected[f"{accession}_TPM"]} for accession in accessions)
        for group_index, (label, values) in enumerate(zip(labels, vals_by_group)):
            values = np.asarray(values)
            assert np.all(np.isfinite(values)) and np.all(values >= 0)
            ax.scatter(group_index + np.array([-.14, 0, .14]), values, s=19, facecolors=GROUP_COLORS[group_index],
                       edgecolors="white", linewidths=.3, zorder=4, clip_on=False)
            mean = float(values.mean())
            ax.plot([group_index - .25, group_index + .25], [mean, mean], color=INK, lw=1.0, zorder=3)
            mean_rows.append({"panel": "abcdef"[row * 3 + col], "gene": gene, "species": "human" if row == 0 else "mouse",
                              "group": label, "arithmetic_mean": mean, "library_n": 3})
        maximum = max(max(v) for v in vals_by_group)
        ax.set_ylim(0, 1 if row == 1 and gene == "JAG2" else maximum * 1.22)
        ax.set_xlim(-.48, 2.48)
        ax.set_xticks(range(3), labels)
        ax.tick_params(axis="x", length=0, pad=5)
        ax.tick_params(axis="y", length=3, pad=3)
        ax.set_ylabel("Normalized count" if row == 0 else "TPM", labelpad=5)
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=4, min_n_ticks=3))
        if row == 1 and gene == "JAG2":
            ax.set_yticks([0, .5, 1])
fig.text(.13, .080, "Each point is one submitted library; short lines show means (3 libraries per group).", fontsize=7)
fig.text(.13, .044, "Separate y-scales and units; biological independence is unresolved. No inferential tests.", fontsize=7)
assert len(expression_rows) == 54
write_csv("fig2_source_data.csv", expression_rows)
write_csv("fig2_group_means.csv", mean_rows)
export(fig, "fig2_expression")

# Figure 3: categorical signs, not mixed-unit effect magnitudes.
primary = [r for r in facts["mouse_candidates_all26_verbatim"] if r["primary_comparison_eligible"] == "true"]
assert len(primary) == 15
columns = ["human_OA_vs_ctrl_direction", "human_sesamin_vs_OA_direction", "mouse_HFD_vs_ctrl_direction", "mouse_sess_vs_HFD_direction"]
sign_map = {"down": -1, "unchanged": 0, "up": 1}
symbol_map = {-1: "↓", 0: "0", 1: "↑"}
matrix = np.array([[sign_map[r[c]] for c in columns] for r in primary])
both_n = sum(r["both_cross_species_directions_same"] == "true" for r in primary)
intervention_n = sum(r["cross_species_sesamin_direction"] == "same" for r in primary)
assert both_n == 3 and intervention_n == 8
assert sum(r["all_libraries_TPM_lt1"] == "true" for r in primary) == 5
direction_rows = [{"human_symbol": r["human_symbol"], **{c: r[c] for c in columns},
                   "all_mouse_libraries_below_1_TPM": r["all_libraries_TPM_lt1"],
                   "both_cross_species_directions_same": r["both_cross_species_directions_same"],
                   "cross_species_sesamin_direction": r["cross_species_sesamin_direction"]} for r in primary]
write_csv("fig3_source_data.csv", direction_rows)
fig, ax = plt.subplots(figsize=(180 / 25.4, 158 / 25.4))
fig.subplots_adjust(left=.245, right=.77, top=.81, bottom=.185)
categorical = ListedColormap(["#D5E6EF", "#F0F1F1", "#F2DCC7"])
ax.imshow(matrix, cmap=categorical, norm=BoundaryNorm([-1.5, -.5, .5, 1.5], 3), interpolation="none", aspect="auto")
for i in range(15):
    for j in range(4):
        ax.text(j, i, symbol_map[matrix[i, j]], ha="center", va="center", fontsize=10, color=INK)
ax.set_yticks(range(15), [r["human_symbol"] + (" *" if r["all_libraries_TPM_lt1"] == "true" else "") for r in primary])
for label in ax.get_yticklabels():
    label.set_fontstyle("italic")
ax.set_xticks(range(4), ["Model /\nControl", "Treat /\nModel", "HFD /\nctrl", "sess /\nHFD"])
ax.xaxis.tick_top()
ax.tick_params(axis="both", length=0, pad=6)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
ax.set_yticks(np.arange(-.5, 15, 1), minor=True)
ax.grid(which="minor", color="white", linewidth=1.1)
ax.tick_params(which="minor", length=0)
fig.text(.375, .926, "Human HepG2", ha="center", fontsize=9, fontweight="bold")
fig.text(.638, .926, "Mouse", ha="center", fontsize=9, fontweight="bold")
fig.text(.638, .901, "Provisional groups", ha="center", fontsize=7)
fig.text(.817, .806, "Direction", fontsize=8, fontweight="bold")
fig.text(.817, .766, "↑  Higher", fontsize=8)
fig.text(.817, .730, "↓  Lower", fontsize=8)
fig.text(.817, .694, "0  Unchanged", fontsize=8)
fig.text(.13, .122, "3/15 agree for both contrasts; 8/15 agree for the treatment contrast (including those 3).", fontsize=7.5)
fig.text(.13, .087, "* All 9 mouse libraries <1 TPM (5 genes). All primary one-to-one candidates are shown.", fontsize=7)
fig.text(.13, .052, "Signs only: neither effect magnitude nor statistical significance is encoded.", fontsize=7)
export(fig, "fig3_concordance")

# Figure 4: external disease-only backgrounds, with non-equivalent units separated.
phh_values = np.array([[float(next(r for r in phh if r["symbol"] == gene and r["donor_label"] == label)["log2_FFA_Control_mean_norm_plus1"])
                        for label in PHH_LABELS] for gene in BACKGROUND])
org_values = np.array([[float(next(r for r in organoid if r["symbol"] == gene)[key])
                       for key in ["observed_CPM_log2_ratio_plus1", "fixed_weight_CPM_log2_ratio_plus1"]] for gene in BACKGROUND])
assert np.all(np.isfinite(phh_values)) and np.all(np.isfinite(org_values))
background_rows = []
for gene in BACKGROUND:
    for label in PHH_LABELS:
        record = next(r for r in phh if r["symbol"] == gene and r["donor_label"] == label)
        background_rows.append({"panel": "a", "gene": gene, "label_or_summary": label,
                                "source_field": "log2_FFA_Control_mean_norm_plus1",
                                "value_verbatim": record["log2_FFA_Control_mean_norm_plus1"]})
    record = next(r for r in organoid if r["symbol"] == gene)
    for label, field in [("Original", "observed_CPM_log2_ratio_plus1"), ("Fixed UMI weights", "fixed_weight_CPM_log2_ratio_plus1")]:
        background_rows.append({"panel": "b", "gene": gene, "label_or_summary": label,
                                "source_field": field, "value_verbatim": record[field]})
write_csv("fig4_source_data.csv", background_rows)
fig = plt.figure(figsize=(180 / 25.4, 104 / 25.4))
gs = fig.add_gridspec(1, 2, left=.12, right=.935, bottom=.345, top=.77, width_ratios=[5, 2], wspace=.75)
axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
color_axes = []
for index, (ax, values, columns_display, bound) in enumerate(zip(axes, [phh_values, org_values], [PHH_LABELS, ["Original", "Fixed UMI\nweights"]], [2.0, 2.5])):
    im = ax.imshow(values, vmin=-bound, vmax=bound, cmap=DIVERGE, interpolation="none", aspect="auto")
    for i in range(4):
        for j in range(len(columns_display)):
            ax.text(j, i, f"{values[i,j]:+.2f}", ha="center", va="center", fontsize=7,
                    color="white" if abs(values[i,j]) / bound > .65 else INK)
    ax.set_yticks(range(4), BACKGROUND, fontstyle="italic")
    ax.set_xticks(range(len(columns_display)), columns_display)
    ax.tick_params(axis="both", length=0, pad=5)
    ax.set_xticks(np.arange(-.5, len(columns_display), 1), minor=True)
    ax.set_yticks(np.arange(-.5, 4, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.0)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_letter(ax, "ab"[index], x_offset=-27, y_offset=14)
    ax.set_title(["PHH: FFA / control", "Organoids: PA / control"][index], fontsize=8.5, pad=15)
    pos = ax.get_position()
    cax = fig.add_axes([pos.x0, .255, pos.width, .022])
    color_axes.append(cax)
    colorbar = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[-bound, 0, bound])
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(length=2, pad=2, labelsize=7)
fig.text(.12, .94, "Disease-only backgrounds; no sesamin intervention", fontsize=9, fontweight="bold")
fig.text(.12, .894, "Signs and magnitudes vary with gene, donor label and measurement context.", fontsize=7.5)
fig.text(.12, .178, "a  log2[(FFA mean normalized count + 1) / (control mean normalized count + 1)]", fontsize=7)
fig.text(.12, .136, "b  log2[(PA mean CPM + 1) / (control mean CPM + 1)]; independent color scales", fontsize=7)
fig.text(.12, .083, "PHH: 5 donor labels, F95/F95P independence unresolved.", fontsize=7)
fig.text(.12, .043, "Organoids: 4 technical libraries, 1 differentiation; reweighting is a sensitivity analysis.", fontsize=7)
export(fig, "fig4_background", exclude_axes=color_axes)

source_manifest = []
for key, path in SOURCES.items():
    payload = path.read_bytes()
    source_manifest.append({"key": key, "path": str(path), "sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)})
write_json("source_data.json", {
    "date": "2026-10-09", "source_inputs": source_manifest,
    "figure2_library_values": expression_rows, "figure2_means": mean_rows,
    "figure3_directions": direction_rows, "figure4_background": background_rows,
    "checks": {"figure2_total_library_points": 54, "figure3_primary_genes": 15,
               "figure3_both_direction_n": both_n, "figure3_treatment_direction_n_including_both": intervention_n,
               "figure3_low_expression_genes": 5, "figure4_phh_values": 20, "figure4_organoid_values": 8},
    "environment": {"python": sys.version, "executable": sys.executable, "matplotlib": matplotlib.__version__, "numpy": np.__version__},
    "analysis_scope": "Existing observations copied; arithmetic means only for display. No quantification rerun or inferential tests.",
    "randomness": "None; deterministic horizontal offsets, no random jitter or subsampling.",
    "plotting_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
})
print(json.dumps({"rendered": ["fig2_expression", "fig3_concordance", "fig4_background"], "library_points": len(expression_rows), "primary_genes": len(primary)}))
