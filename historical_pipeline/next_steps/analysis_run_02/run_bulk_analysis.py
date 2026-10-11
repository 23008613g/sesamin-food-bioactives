"""Exploratory GSE200765 Welch pilot and data-driven Python figures.

Requires the approved project environment. Raw source files remain read-only.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
STEP = HERE.parent
BASE = STEP.parent
sys.path.insert(0, str(STEP))
from run_descriptive_pilot import read_table

GROUPS = ["vehicle", "FFA", "active_vitamin_D_plus_FFA", "SM_VD2_plus_FFA", "SV_VD3_plus_FFA"]
SHORT = ["V", "F", "A", "S", "D"]
LABELS = ["Vehicle", "FFA", "Active D + FFA", "SM-VD2 + FFA", "SV-VD3 + FFA"]
COLORS = ["#707070", "#B87434", "#4F79A5", "#708F72", "#936990"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export(fig, name, require_matplotlib_panel_alignment, exclude_axes=None):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=HERE / f"{name}.alignment.json", tolerance_pt=1.5,
              gutter_tolerance_pt=1.5, strict=True, exclude_axes=exclude_axes or [])
    fig.savefig(HERE / f"{name}.pdf", facecolor="white")
    fig.savefig(HERE / f"{name}.svg", facecolor="white")
    fig.savefig(HERE / f"{name}.png", dpi=300, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qa-scripts", required=True, type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.qa_scripts))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    contract = json.loads((HERE / "analysis_contract.json").read_text(encoding="utf-8"))
    original = json.loads((STEP / "pilot_protocol.json").read_text(encoding="utf-8"))
    matrix = BASE / "downloads/GSE200765_series_matrix.txt.gz"
    audit = json.loads((BASE / "processed_table_audit.json").read_text(encoding="utf-8"))
    assert digest(matrix) == next(r["sha256"] for r in audit["records"] if r["accession"] == "GSE200765")
    header, rows = read_table(matrix, "!series_matrix_table_begin", "!series_matrix_table_end")
    rows = [r for r in rows if r[0].startswith("TC") and ".hg." in r[0]]
    ids = [r[0] for r in rows]
    samples = header[1:]
    x = np.asarray([list(map(float, r[1:])) for r in rows])
    assert x.shape == (19911, 14) and np.isfinite(x).all()
    manifest = pd.read_csv(BASE / "sample_manifest.csv").set_index("gsm")
    group = [manifest.loc[s, "group"] for s in samples]
    ix = {g: np.flatnonzero(np.asarray(group) == g) for g in GROUPS}
    annotation = pd.read_csv(STEP / "pilot_descriptive/feature_annotation.csv", keep_default_na=False).set_index("feature_id")
    assert annotation.index.tolist() == ids
    contrast_summaries = []
    checks = []
    for a, b in contract["contrasts"]:
        xa, xb = x[:, ix[a]], x[:, ix[b]]
        constant_a = np.ptp(xa, axis=1) == 0
        constant_b = np.ptp(xb, axis=1) == 0
        near_precision_a = np.ptp(xa, axis=1) <= .01 + 1e-12
        near_precision_b = np.ptp(xb, axis=1) <= .01 + 1e-12
        effect = xa.mean(axis=1) - xb.mean(axis=1)
        se = np.sqrt(xa.var(axis=1, ddof=1)/xa.shape[1] + xb.var(axis=1, ddof=1)/xb.shape[1])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            test = stats.ttest_ind(xa, xb, axis=1, equal_var=False, nan_policy="raise")
            ci = test.confidence_interval(confidence_level=.95)
        valid = np.isfinite(test.pvalue) & (se > 0)
        p_for_bh = np.where(valid, test.pvalue, 1.)
        q = stats.false_discovery_control(p_for_bh, method="bh")
        manual_t = np.divide(effect, se, out=np.zeros_like(effect), where=se > 0)
        if valid.any():
            assert np.allclose(manual_t[valid], test.statistic[valid], rtol=1e-10, atol=1e-10)
            p_manual = 2*stats.t.sf(np.abs(manual_t[valid]), test.df[valid])
            assert np.allclose(p_manual, test.pvalue[valid], rtol=1e-10, atol=1e-12)
        tab = pd.DataFrame({"feature_id": ids, "gene_symbols": annotation["gene_symbols"].to_numpy(),
            "mapping_status": annotation["mapping_status"].to_numpy(), "n_treated": xa.shape[1], "n_reference": xb.shape[1],
            "mean_treated": xa.mean(axis=1), "mean_reference": xb.mean(axis=1), "mean_signal_difference": effect,
            "SE_assuming_independent_samples": se, "welch_t": test.statistic, "welch_df": test.df,
            "p_two_sided": test.pvalue, "q_BH_this_contrast": q, "CI95_low_pointwise": ci.low,
            "CI95_high_pointwise": ci.high, "test_defined": valid,
            "constant_treated_on_submitted_scale": constant_a, "constant_reference_on_submitted_scale": constant_b,
            "any_group_range_le_001": near_precision_a | near_precision_b})
        name = a + "_vs_" + b
        tab.to_csv(HERE / (name + ".csv"), index=False)
        contrast_summaries.append({"contrast": name, "n_treated": xa.shape[1], "n_reference": xb.shape[1],
            "tested_clusters": int(valid.sum()), "undefined_clusters": int((~valid).sum()),
            "clusters_q_lt_005": int(np.count_nonzero(valid & (q < .05))),
            "any_constant_group_clusters": int(np.count_nonzero(constant_a | constant_b)),
            "both_constant_group_clusters": int(np.count_nonzero(constant_a & constant_b)),
            "any_group_range_le_001_clusters": int(np.count_nonzero(near_precision_a | near_precision_b)),
            "warnings": sorted({str(w.message) for w in caught})})
        checks.append({"contrast": name, "independent_formula_crosscheck_passed": True})
    mpl_settings = {"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 7, "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7}
    plt.rcParams.update(mpl_settings)
    centered = x.T - x.mean(axis=1)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    scores = u[:, :2] * s[:2]
    fraction = s*s/np.sum(s*s)
    pd.DataFrame({"gsm": samples, "group": group, "PC1": scores[:, 0], "PC2": scores[:, 1]}).to_csv(HERE / "pca_scores.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.2047244094, 4.2519685039))
    fig.subplots_adjust(left=.09, right=.97, bottom=.29, top=.86, wspace=.36)
    bp = axes[0].boxplot(x, tick_labels=[SHORT[GROUPS.index(g)]+str(list(ix[g]).index(i)+1) for i,g in enumerate(group)],
                        showfliers=False, widths=.55, patch_artist=True)
    for box, g in zip(bp["boxes"], group, strict=True):
        box.set_facecolor(COLORS[GROUPS.index(g)])
        box.set_alpha(.65)
    for line in bp["medians"]:
        line.set_color("black")
    axes[0].set_ylabel("Submitted signal")
    axes[0].set_xlabel("Sample")
    axes[0].tick_params(axis="x", labelrotation=90)
    axes[0].set_title("a  Sample distributions", loc="left", fontweight="bold", pad=12)
    for g, label, color in zip(GROUPS, LABELS, COLORS, strict=True):
        indexes = ix[g]
        axes[1].scatter(scores[indexes, 0], scores[indexes, 1], s=27, label=label, color=color,
                        edgecolors="white", linewidths=.5, zorder=3)
    axes[1].set_xlabel(f"PC1 ({fraction[0]:.1%} variance)")
    axes[1].set_ylabel(f"PC2 ({fraction[1]:.1%} variance)")
    axes[1].set_title("b  Overall expression patterns", loc="left", fontweight="bold", pad=12)
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, .075), ncol=3, frameon=False)
    fig.text(.5, .025, "All 14 samples retained; 19,911 human transcript clusters. Descriptive QC only.", ha="center", fontsize=7)
    export(fig, "sample_quality", require_matplotlib_panel_alignment)
    panel = pd.read_csv(STEP / "pilot_descriptive/literature_panel_samples.csv")
    fig, axes = plt.subplots(2, 4, figsize=(7.2047244094, 5.5118110236))
    fig.subplots_adjust(left=.08, right=.975, bottom=.19, top=.90, wspace=.55, hspace=.52)
    for j, symbol in enumerate(original["fixed_literature_panel"]):
        ax = axes.flat[j]
        gene = panel[panel.gene_symbol == symbol]
        assert gene.feature_id.nunique() == 1 and len(gene) == 14
        for k, (g, color) in enumerate(zip(GROUPS, COLORS, strict=True)):
            vals = gene.loc[gene.group == g, "submitted_signal"].to_numpy()
            jitter = np.linspace(-.12, .12, len(vals))
            ax.scatter(k+jitter, vals, color=color, s=17, alpha=.9, linewidths=0, zorder=4)
            ax.errorbar(k, vals.mean(), yerr=vals.std(ddof=1), fmt="_", color="black", capsize=3,
                        markersize=8, elinewidth=.8, linewidth=.8, zorder=3)
        ax.set_xticks(range(5), SHORT)
        ax.set_xlim(-.6, 4.6)
        ax.margins(y=.20)
        ax.set_title(chr(97+j)+"  "+symbol, loc="left", fontweight="bold", pad=11)
        if j % 4 == 0:
            ax.set_ylabel("Submitted signal")
    axes.flat[7].set_axis_off()
    fig.text(.5, .105, "V: Vehicle   F: FFA   A: Active D + FFA   S: SM-VD2 + FFA   D: SV-VD3 + FFA", ha="center", fontsize=7)
    fig.text(.5, .066, "Points: all submitted samples. Black marks: mean and sample SD. n = 3, 3, 2, 3, 3.", ha="center", fontsize=7)
    fig.text(.5, .027, "Fixed literature-related panel; RNA observations do not reproduce protein assays or establish efficacy.", ha="center", fontsize=7)
    export(fig, "literature_panel", require_matplotlib_panel_alignment, exclude_axes=[axes.flat[7]])
    summary = {"executed_at_utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
        "interpreter": sys.executable, "package_versions": {m:importlib.metadata.version(m) for m in ["numpy","scipy","pandas","matplotlib"]},
        "script_sha256": digest(Path(__file__)), "contract_sha256": digest(HERE / "analysis_contract.json"),
        "source_sha256": digest(matrix), "samples": len(samples), "human_clusters": len(ids),
        "contrasts": contrast_summaries, "numerical_checks": checks, "PCA_variance_fraction": fraction[:2].tolist(),
        "inference_limit": "Exploratory unmoderated Welch outputs conditional on unverified independent culture samples and normal errors; no validated biological precision or original-limma reproduction",
        "replicate_independence_verified": False, "AI_calls": 0, "independent_ingredient_validation": False,
        "figures": ["sample_quality", "literature_panel"], "missing_values_removed": 0, "samples_removed": 0}
    (HERE / "analysis_summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
