"""Render four evidence-bounded Stage 39 scientific figures with matplotlib.

All quantitative inputs are frozen executed results. No data exclusions,
subsampling, inferential recalculation, or changes to upstream result tables.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = ROOT / "figures"
SOURCE = OUT / "source_data"
QA = OUT / "qa"
SKILL = OUT / "qa_tools"
AUDIT_PYTHON = Path(os.environ.get("PFHN_FIGURE_AUDIT_PYTHON", sys.executable))
if not (SKILL / "audit_panel_alignment.py").is_file():
    raise RuntimeError("Bundled figure auditors missing: preserve figures/qa_tools beside the results and scripts directories.")
sys.dont_write_bytecode = True
sys.path.insert(0, str(SKILL))
from audit_panel_alignment import require_matplotlib_panel_alignment

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "legend.frameon": False,
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.spines.right": False, "axes.spines.top": False,
    "axes.linewidth": 0.7, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "text.color": "#172A36", "axes.labelcolor": "#172A36",
    "xtick.color": "#33434D", "ytick.color": "#33434D",
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 300,
})
COLORS = {"control": "#617B92", "sesamin": "#168778", "primary": "#176D70",
          "deposited": "#176D70", "qn10": "#B56F31", "qn20": "#677995"}
RAT_ORDER = [f"GSM{i}" for i in range(1048036, 1048046)]
SHORT = [f"C{i}" for i in range(1, 6)] + [f"S{i}" for i in range(1, 6)]
VARIANTS = ["deposited", "qn10", "qn20"]
DISPLAY = {"deposited": "Deposited", "qn10": "QN10", "qn20": "QN20"}
INPUTS = {}
AUDITS = []


def audit_runtime():
    """Check the explicitly selected auditor interpreter without installation."""
    code = ("import json,sys; import fitz; "
            "print(json.dumps({'python_executable':sys.executable,'python_version':sys.version,"
            "'pymupdf_version':fitz.VersionBind,'mupdf_version':fitz.VersionFitz}))")
    try:
        result = subprocess.run([str(AUDIT_PYTHON), "-c", code], capture_output=True, text=True, encoding="utf-8")
    except OSError as exc:
        raise RuntimeError("Figure audit interpreter is unavailable. Set PFHN_FIGURE_AUDIT_PYTHON to an existing Python with PyMuPDF; no installation was attempted.") from exc
    if result.returncode:
        raise RuntimeError("Figure PDF collision auditing requires PyMuPDF in the audit interpreter. Set PFHN_FIGURE_AUDIT_PYTHON to an existing environment containing it, or prepare dependencies separately. No installation was attempted. Interpreter error: " + result.stderr.strip())
    return json.loads(result.stdout)


def read(path, **kwargs):
    INPUTS[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path, **kwargs)


def source(frame, name):
    frame.to_csv(SOURCE / name, index=False)


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def letter(ax, label):
    # Fixed physical offset maintains common panel-label anchors.
    trans = ax.transAxes + ScaledTranslation(-25 / 72, 17 / 72, ax.figure.dpi_scale_trans)
    ax.text(0, 1, label, transform=trans, fontsize=10, fontweight="bold", ha="left", va="bottom", clip_on=False)


def title(fig, text, subtext):
    fig.text(0.055, 0.965, text, ha="left", va="top", fontsize=11, fontweight="bold")
    fig.text(0.055, 0.91, subtext, ha="left", va="top", fontsize=8, color="#54656E")


def exported(fig, name, axes):
    fig.canvas.draw()
    alignment = require_matplotlib_panel_alignment(
        fig, axes=axes, panel_ids=[chr(97 + i) for i in range(len(axes))],
        row_groups=[{"panels": [chr(97 + i) for i in range(len(axes))]}] if len(axes) > 1 else None,
        json_out=QA / f"{name}.alignment.json", overlay_svg=QA / f"{name}.alignment.svg",
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=len(axes) > 1, strict=True)
    # Fixed canvas size: no tight crop changes the declared physical dimensions.
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.svg")
    fig.savefig(OUT / f"{name}.png", dpi=300)
    plt.close(fig)
    pdf = OUT / f"{name}.pdf"
    t = subprocess.run([str(AUDIT_PYTHON), str(SKILL / "audit_pdf_text.py"), str(pdf), "--min-pt", "5", "--json"], capture_output=True, text=True, encoding="utf-8")
    (QA / f"{name}.text-audit.json").write_text(t.stdout, encoding="utf-8")
    c = subprocess.run([str(AUDIT_PYTHON), str(SKILL / "audit_figure_collisions.py"), str(pdf),
                        "--json-out", str(QA / f"{name}.collision-audit.json"),
                        "--overlay-pdf", str(QA / f"{name}.collision-audit.pdf")],
                       capture_output=True, text=True, encoding="utf-8")
    (QA / f"{name}.collision-audit.txt").write_text(c.stdout + c.stderr, encoding="utf-8")
    result = {"figure": name, "alignment_verdict": alignment.get("verdict"), "text_audit_exit": t.returncode,
              "collision_audit_exit": c.returncode, "width_mm": float(183), "png_dpi": 300}
    AUDITS.append(result)
    print(json.dumps(result), flush=True)


def fig_g1():
    pca = read(RESULTS / "pca_primary10.csv").set_index("sample").loc[RAT_ORDER].reset_index()
    variance_path = RESULTS / "pca_primary10_variance.json"
    INPUTS[str(variance_path.relative_to(ROOT))] = hashlib.sha256(variance_path.read_bytes()).hexdigest()
    explained = json.loads(variance_path.read_text())["explained_variance_ratio"]
    expression = read(RESULTS / "rat_gene_expression_deposited.csv.gz", dtype={"rat_entrez": str})
    vals = expression[RAT_ORDER].to_numpy(dtype=float)
    assert vals.shape == (14604, 10) and np.isfinite(vals).all()
    source(pca.assign(sample_label=SHORT, group=["Control"] * 5 + ["Sesamin"] * 5), "G1a_PCA_all_10_animals.csv")
    source(expression, "G1b_all_14604_gene_distributions.csv")
    q = np.quantile(vals, [0, 0.25, 0.5, 0.75, 1], axis=0)
    source(pd.DataFrame({"sample": RAT_ORDER, "display": SHORT, "minimum": q[0], "q25": q[1], "median": q[2], "q75": q[3], "maximum": q[4]}), "G1b_distribution_summaries.csv")
    fig = plt.figure(figsize=(183 / 25.4, 106 / 25.4))
    a = fig.add_axes([0.10, 0.25, 0.34, 0.53])
    b = fig.add_axes([0.60, 0.25, 0.34, 0.53])
    title(fig, "Whole-transcriptome sample structure", "GSE42693 | 14,604 eligible genes | Five individual rat livers per group")
    for sl, label, color, marker in [(slice(0, 5), "Control", COLORS["control"], "o"), (slice(5, 10), "Sesamin", COLORS["sesamin"], "D")]:
        a.scatter(pca.PC1.iloc[sl], pca.PC2.iloc[sl], s=33, marker=marker, c=color, edgecolors="white", linewidths=0.6, label=label, zorder=3)
    # C1 and S1 are close: labels are intentionally separated vertically.
    offsets = [(3, 6), (4, 4), (4, -10), (4, -11), (-16, 4), (4, -12), (4, 4), (4, -11), (4, -11), (4, 4)]
    for i, row in pca.iterrows():
        a.annotate(SHORT[i], (row.PC1, row.PC2), xytext=offsets[i], textcoords="offset points", fontsize=7.2)
    a.set_xlim(-52, 72); a.set_ylim(-54, 60)
    a.set_xlabel(f"PC1 ({explained[0] * 100:.1f}% variance)")
    a.set_ylabel(f"PC2 ({explained[1] * 100:.1f}% variance)")
    a.set_title("All-gene PCA", pad=9, loc="left")
    a.legend(loc="lower left", bbox_to_anchor=(-0.025, -0.38), ncol=2, handletextpad=0.3, columnspacing=1.2)
    box = b.boxplot([vals[:, i] for i in range(10)], positions=np.arange(10), widths=0.55,
                    patch_artist=True, showfliers=True, whis=(0, 100),
                    medianprops={"color": "white", "linewidth": 1.2},
                    boxprops={"linewidth": 0.65}, whiskerprops={"linewidth": 0.65}, capprops={"linewidth": 0.65})
    for i, patch in enumerate(box["boxes"]):
        patch.set_facecolor(COLORS["control" if i < 5 else "sesamin"])
        patch.set_edgecolor(COLORS["control" if i < 5 else "sesamin"])
    b.set_xticks(np.arange(10), SHORT)
    b.set_ylabel("Gene-level log2(signal + 1)")
    b.set_ylim(0, 17)
    b.set_title("Full gene-signal distributions", pad=9, loc="left")
    b.set_xlabel("Control                      Sesamin", labelpad=7)
    b.tick_params(axis="x", length=0)
    for ax, label in [(a, "a"), (b, "b")]: letter(ax, label)
    exported(fig, "Fig_G1_sample_structure", [a, b])


def fig_g2():
    rows, examples = [], []
    for variant in VARIANTS:
        d = read(RESULTS / f"rat_gene_results_{variant}.csv.gz", dtype={"rat_entrez": str})
        l = read(RESULTS / f"limma/rat_limma_{variant}.csv.gz", dtype={"gene_id": str})
        assert len(d) == len(l) == 14604
        for method, values in [("Welch", d.welch_bh), ("limma", l["adj.P.Val"]), ("Exact labels", d.exact_bh)]:
            rows.append({"normalization": variant, "method": method, "BH_q_less_than_0.05": int((values < 0.05).sum()), "family_size": len(values)})
        e = d.loc[d.rat_entrez.isin(["29142", "29183"])].copy()
        assert len(e) == 2
        e["normalization"] = variant
        examples.append(e)
    counts, ex = pd.DataFrame(rows), pd.concat(examples, ignore_index=True)
    source(counts, "G2a_full_genome_BH_counts.csv")
    source(ex, "G2b_preexisting_examples_complete_results.csv")
    methods = ["Welch", "limma", "Exact labels"]
    matrix = counts.pivot(index="normalization", columns="method", values="BH_q_less_than_0.05").loc[VARIANTS, methods].to_numpy()
    assert np.all(matrix >= 0)
    fig = plt.figure(figsize=(183 / 25.4, 110 / 25.4))
    a = fig.add_axes([0.15, 0.27, 0.29, 0.49])
    b = fig.add_axes([0.72, 0.27, 0.23, 0.49])
    title(fig, "Statistical support depends on the analysis method", "The full 14,604-gene family is retained under every normalization and test")
    cmap = LinearSegmentedColormap.from_list("count_teal", ["#F3F7F6", "#AED0C6", "#276F69"])
    a.imshow(np.log1p(matrix), aspect="auto", cmap=cmap, vmin=0, vmax=float(np.log1p(matrix.max())))
    for i in range(3):
        for j in range(3):
            a.text(j, i, str(matrix[i, j]), ha="center", va="center", fontsize=11,
                   color="white" if matrix[i, j] >= 20 else "#172A36", fontweight="bold")
    a.set_xticks(range(3), ["Welch", "limma", "Exact\nlabels"])
    a.set_yticks(range(3), [DISPLAY[v] for v in VARIANTS])
    a.tick_params(length=0, pad=7)
    a.set_title("Genes with BH q < 0.05", loc="left", pad=12)
    for spine in a.spines.values(): spine.set_visible(False)
    ordered = []
    for gene in ["29142", "29183"]:
        for variant in VARIANTS:
            ordered.append(ex.loc[(ex.rat_entrez == gene) & (ex.normalization == variant)].iloc[0])
    ordered = pd.DataFrame(ordered).reset_index(drop=True)
    y = np.array([6, 5, 4, 2, 1, 0])
    b.axvline(0, color="#B7C2C8", linewidth=0.8, linestyle=(0, (3, 3)), zorder=0)
    for i, row in ordered.iterrows():
        b.errorbar(row.difference, y[i], xerr=np.array([[row.difference - row.welch_ci_low], [row.welch_ci_high - row.difference]]),
                   fmt="o", markersize=4, color=COLORS[row.normalization], ecolor=COLORS[row.normalization], capsize=2, linewidth=1)
    b.set_yticks(y, [f"{r.gene_symbol} · {DISPLAY[r.normalization]}" for _, r in ordered.iterrows()])
    b.tick_params(axis="y", length=0, pad=6)
    b.set_ylim(-0.9, 6.9)
    b.set_xlim(-1.5, 2.65)
    b.set_xticks([-1, 0, 1, 2])
    b.set_xlabel("Mean log2-signal difference\nSesamin − control", labelpad=7)
    b.set_title("Pre-existing examples", loc="left", pad=12)
    for ax, lab in [(a, "a"), (b, "b")]: letter(ax, lab)
    fig.text(0.15, 0.13, "Cell labels give exact counts.\nColor intensity: log(count + 1).", fontsize=7.3, va="top", color="#54656E")
    fig.text(0.56, 0.11, "Points: mean differences; bars: nominal 95% Welch CI.\nVNN1 and AREG were identified before this extension.", fontsize=7.1, va="top", color="#54656E")
    exported(fig, "Fig_G2_method_sensitivity", [a, b])


def fig_g3():
    g = read(RESULTS / "go_results_deposited.csv")
    ranked = g.assign(abs_difference=g.difference.abs()).sort_values(["welch_bh", "abs_difference", "go_id"], ascending=[True, False, True])
    top = ranked.head(15).copy()
    assert len(g) == 1311 and len(top) == 15
    top.insert(0, "display_rank", range(1, 16))
    source(top, "G3_top15_GO_display_ranking.csv")
    source(g, "G3_all1311_GO_complete_results.csv")
    lipid = read(RESULTS / "go_lipid_terms_deposited.csv")
    assert len(lipid) == 31
    source(lipid, "G3_all31_lipid_name_GO_terms_not_selected_for_plot.csv")
    fig = plt.figure(figsize=(183 / 25.4, 203 / 25.4))
    a = fig.add_axes([0.595, 0.14, 0.255, 0.72])
    fig.text(0.045, 0.975, "Coordinated gene-set shifts have limited statistical support", va="top", fontsize=10.5, fontweight="bold")
    fig.text(0.045, 0.94, "Primary deposited signals | 15 sets shown from 1,311 tested GO biological processes", va="top", fontsize=7.8, color="#54656E")
    fig.text(0.045, 0.911, "No set passes BH q < 0.05 under either Welch or exact-label inference", va="top", fontsize=7.8, color="#54656E")
    y = np.arange(len(top))[::-1]
    a.axvline(0, color="#B6C3C8", linewidth=0.8, linestyle=(0, (3, 3)), zorder=0)
    for yi, (_, row) in zip(y, top.iterrows()):
        a.errorbar(row.difference, yi, xerr=np.array([[row.difference - row.welch_ci_low], [row.welch_ci_high - row.difference]]),
                   fmt="o", markersize=4, color=COLORS["primary"], ecolor=COLORS["primary"], capsize=2.2, linewidth=1)
        a.text(1.1, yi, f"{row.welch_bh:.3f}", transform=a.get_yaxis_transform(), ha="left", va="center", fontsize=7.5)
    names = [textwrap.fill(row.go_name, width=43, break_long_words=False) + "\n" + row.go_id for _, row in top.iterrows()]
    a.set_yticks(y, names)
    a.tick_params(axis="y", length=0, pad=9, labelsize=7.2)
    a.set_ylim(-0.7, 14.7)
    a.set_xlim(-1.72, 1.8)
    a.set_xticks([-1.5, 0, 1.5])
    a.set_xlabel("Mean gene-set score difference\nSesamin − control", labelpad=9)
    a.spines["left"].set_visible(False)
    a.text(1.1, 1.035, "Welch\nBH q", transform=a.transAxes, ha="left", va="bottom", fontsize=7.5, fontweight="bold")
    fig.text(0.045, 0.071, "Points: differences between five-rat group means. Bars: nominal 95% Welch CI.\nDisplay order: BH q, then absolute score difference, then GO ID; all tested sets remain in the source data.", fontsize=7.2, color="#54656E", va="top")
    exported(fig, "Fig_G3_GO_score_evidence", [a])


def fig_g4():
    cross = RESULTS / "cross_species"
    d = read(cross / "shared_gene_effects_rat_deposited_human_mr.csv.gz")
    s = read(cross / "gene_comparison_summary.csv")
    assert len(d) == 7912 and np.isfinite(d[["rat_sesamin_minus_control", "human_treat_vs_model"]]).all().all()
    sub = s.loc[s.human_contrast == "treat_vs_model"].copy()
    assert len(sub) == 6
    source(d, "G4a_all7912_shared_gene_effects.csv")
    source(sub, "G4b_all6_normalization_comparisons.csv")
    source(read(cross / "gene_direction_contingency.csv"), "G4_all_direction_contingencies.csv")
    fig = plt.figure(figsize=(183 / 25.4, 118 / 25.4))
    a = fig.add_axes([0.11, 0.31, 0.37, 0.45])
    b = fig.add_axes([0.70, 0.31, 0.245, 0.45])
    title(fig, "Genome-wide effect directions show little cross-study alignment", "All 7,912 shared one-to-one genes | Rat sesamin/control versus HepG2 Treat/Model")
    x, y = d.rat_sesamin_minus_control.to_numpy(), d.human_treat_vs_model.to_numpy()
    hx = a.hexbin(x, y, gridsize=39, mincnt=1, norm=LogNorm(), cmap="viridis", linewidths=0,
                  extent=(float(x.min() - 0.12), float(x.max() + 0.12), float(y.min() - 0.2), float(y.max() + 0.2)))
    assert int(np.sum(hx.get_array())) == 7912
    a.set_xlim(float(x.min() - 0.12), float(x.max() + 0.12))
    a.set_ylim(float(y.min() - 0.2), float(y.max() + 0.2))
    a.set_xlabel("Rat mean log2-signal difference\nSesamin − control", labelpad=6)
    a.set_ylabel("HepG2 log2 mean-count ratio\nTreat / Model", labelpad=6)
    a.set_title("Primary full-gene comparison", loc="left", pad=10)
    cbax = fig.add_axes([0.11, 0.10, 0.37, 0.028])
    cb = fig.colorbar(hx, cax=cbax, orientation="horizontal")
    cb.set_ticks([1, 10, 100])
    cb.set_ticklabels(["1", "10", "100"])
    cb.set_label("Genes per hexagon (log color scale)", fontsize=7, labelpad=3)
    cb.ax.tick_params(labelsize=7)
    order = [(v, norm) for v in VARIANTS for norm in ["mr", "cpm"]]
    positions = np.arange(6)[::-1]
    for pos, (variant, norm) in zip(positions, order):
        row = sub.loc[(sub.rat_normalization == variant) & (sub.human_normalization == norm)].iloc[0]
        b.plot(row.spearman_descriptive, pos, marker="o" if norm == "mr" else "s", markersize=5,
               color=COLORS[variant], linestyle="none")
        b.text(0.97, pos, f"{row.spearman_descriptive:.3f}", transform=b.get_yaxis_transform(), ha="right", va="center", fontsize=7.3)
    b.axvline(0, color="#B6C3C8", linewidth=0.8, linestyle=(0, (3, 3)), zorder=0)
    b.set_xlim(-0.05, 0.035)
    b.set_ylim(-0.7, 5.7)
    b.set_xticks([-0.04, -0.02, 0])
    b.set_yticks(positions, [f"{DISPLAY[v]} / {n.upper()}" for v, n in order])
    b.tick_params(axis="y", length=0, pad=5, labelsize=7.1)
    b.set_xlabel("Descriptive Spearman rho", labelpad=7)
    b.set_title("Normalization sensitivity", loc="left", pad=10)
    primary = sub.loc[(sub.rat_normalization == "deposited") & (sub.human_normalization == "mr")].iloc[0]
    fig.text(0.70, 0.18, f"Primary rho = {primary.spearman_descriptive:.3f}\nSame direction: {int(primary.same_nonzero_direction_n):,} / 7,912 ({100 * primary.same_nonzero_direction_fraction_all:.2f}%)\nNo correlation p value or CI.", fontsize=7.3, va="top", color="#54656E")
    for ax, lab in [(a, "a"), (b, "b")]: letter(ax, lab)
    exported(fig, "Fig_G4_cross_study_alignment", [a, b])


LEGENDS = """Figure G1. Whole-transcriptome sample structure in dietary sesamin-exposed rat liver.
(a) Principal component analysis of the ten deposited control and sesamin arrays using all 14,604 eligible gene summaries. Gene-level log2(signal + 1) values were centered without variance scaling; percentages indicate explained variance. Each point is an individual rat liver. C1-C5 correspond to GSM1048036-GSM1048040 and S1-S5 to GSM1048041-GSM1048045. (b) Distributions of the same 14,604 gene values per array: boxes span the interquartile range, white lines denote medians, and whiskers show the full minimum-to-maximum range. Gene values are features rather than independent biological replicates. No animals or genes were removed for this display. The original investigators selected five of seven animals per diet after excluding terminal-weight extremes. PCA separation is descriptive and is not a group-level hypothesis test.

Figure G2. Gene-level statistical support varies by inferential method and normalization.
(a) Numbers of genes with Benjamini-Hochberg-adjusted q < 0.05 across the same 14,604-gene family under Welch tests, limma empirical-Bayes tests with an intensity trend, and exact-label mean-difference reference tests. Cell labels are exact counts; color intensity is proportional to log(count + 1). Deposited denotes the primary log2(signal + 1) analysis. QN10 and QN20 denote quantile-normalization sensitivity analyses over ten primary arrays and all twenty deposited arrays, respectively, before gene aggregation. Exact-label values enumerate all 252 five-versus-five assignments and require exchangeability, which the reanalysis cannot establish. (b) Rat Vnn1 and Areg, corresponding to the pre-existing human VNN1 and AREG candidates, shown under all normalizations without reselecting candidates. Points are sesamin-minus-control differences in mean log2 signal; bars are nominal 95% Welch intervals, not simultaneous intervals or evidence of genome-wide significance. Every comparison contains five individual rat livers per group. These examples were known before this whole-genome extension; moderated results do not replace the primary and exact-label results.

Figure G3. Coordinated-expression scores for the first 15 GO sets under the frozen display ranking.
All 1,311 platform-annotated GO biological-process sets containing 15-500 eligible genes were tested. The plot displays the first 15 after sorting by primary Welch BH q, descending absolute score difference, and GO identifier. This ranking is a display rule rather than a new testing family. Each gene was standardized across the ten rats with sample SD (ddof = 1), and member gene z scores were averaged equally within each animal. Points show differences between five-animal sesamin and control mean scores; bars show nominal 95% Welch intervals. Listed q values are adjusted over all 1,311 GO identifiers. No primary set passes BH q < 0.05 under either Welch or exact-label inference. These scores describe coordinated expression; they do not measure pathway activation or metabolic flux. Complete results, including all 31 GO terms meeting the frozen lipid-name rule, accompany the figure as source data. GO names/statuses use the frozen 2026-07-26 ontology; gene memberships retain the historical platform annotation.

Figure G4. Full-universe expression effects show little alignment between rat liver and HepG2.
(a) Hexagonal-bin plot of every one of the 7,912 shared genes that passed the fixed human CPM filter, rat eligibility rules, and high-confidence one-to-one Ensembl 116 mapping with unique rat Entrez assignment. The x axis is the primary rat sesamin-minus-control mean log2-signal difference; the y axis is the HepG2 Treat/Model log2 ratio of arithmetic mean median-ratio-normalized counts plus one. These effect measures are not the same physical scale. Color denotes the number of genes in each hexagon on a logarithmic scale; no genes were sampled or selected by effect, significance or prior candidate status. (b) Descriptive Spearman correlations under all three rat normalizations and both human normalizations (MR, median ratio; CPM, counts per million). Correlations have no inferential p values or confidence intervals. The rat comparison contains five individual animals per group; the human comparison contains three deposited libraries per group, with independent culture units unresolved. Species, tumor-derived versus intact-liver context, lipid challenge, exposure and assay differ. Directional agreement is not normal-human validation or proof of treatment rescue.
"""


def main():
    runtime = audit_runtime()
    for path in [OUT, SOURCE, QA]: path.mkdir(parents=True, exist_ok=True)
    contract = {
        "backend": "Python/matplotlib inherited from the established project workflow",
        "journal_context": "Plant Foods for Human Nutrition working manuscript; exploratory supplementary analysis; not a claim of final journal acceptance/production compliance",
        "exports": {"width_mm": 183, "formats": ["PDF", "SVG", "PNG"], "png_dpi": 300, "editable_text": True, "minimum_glyph_pt": 5},
        "data_integrity": "All 10 primary rats, all 14,604 eligible genes for QC; all 7,912 shared genes for cross-study comparison. No sampling, outlier exclusion or significant-gene selection. GO display follows frozen top-15 ranking and complete tables remain available.",
        "figures": {
            "G1": {"claim": "Full-gene variability and per-array signal distributions can be inspected for every primary rat.", "archetype": "quantitative grid", "panels": {"a": "global sample structure, all-gene PCA", "b": "measurement-scale distributions, all genes per sample"}},
            "G2": {"claim": "The strength of gene-level support depends on method and normalization.", "archetype": "quantitative grid", "panels": {"a": "full-family inferential sensitivity", "b": "pre-existing examples with effect sizes and nominal uncertainty"}},
            "G3": {"claim": "None of the primary coordinated-expression sets meets the whole-family BH threshold.", "archetype": "quantitative grid", "panels": {"a": "predefined top-15 display of effect, uncertainty and whole-family multiplicity"}},
            "G4": {"claim": "Whole-universe rat and HepG2 effects have near-zero descriptive rank alignment under all normalization combinations.", "archetype": "quantitative grid", "panels": {"a": "all-gene evidence distribution", "b": "normalization sensitivity of the bounded conclusion"}},
        },
        "qa_contract": "Matplotlib final-axes alignment at 1.5 pt, static source validator, PDF glyph-size audit, PDF collision audit after each render, and manual inspection of every final PNG. Colorbars excluded from unrelated alignment groups. All FAIL findings block delivery.",
    }
    dump(OUT / "figure_contract.json", contract)
    v = subprocess.run([sys.executable, str(SKILL / "validate_figure.py"), str(Path(__file__)), "--json"], capture_output=True, text=True, encoding="utf-8")
    (QA / "source_validation.json").write_text(v.stdout, encoding="utf-8")
    fig_g1(); fig_g2(); fig_g3(); fig_g4()
    (OUT / "figure_legends.txt").write_text(LEGENDS, encoding="utf-8")
    dump(QA / "figure_build.json", {"utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
        "python": sys.executable, "python_version": sys.version, "matplotlib": mpl.__version__, "numpy": np.__version__, "pandas": pd.__version__,
        "audit_runtime": runtime, "audit_interpreter_selection": "PFHN_FIGURE_AUDIT_PYTHON environment variable, otherwise sys.executable",
        "bundled_qa_tool_hashes": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(SKILL.glob("*.py"))},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "input_hashes": INPUTS,
        "source_validation_exit": v.returncode, "figures": AUDITS,
        "visual_inspection": "pending manual view_image inspection; not claimed by script"})
    if v.returncode or any(a["text_audit_exit"] or a["collision_audit_exit"] for a in AUDITS):
        raise SystemExit("Rendered/source QA has findings requiring inspection or repair; see figures/qa.")


if __name__ == "__main__":
    main()
