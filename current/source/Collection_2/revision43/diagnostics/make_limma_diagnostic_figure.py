"""Render proposed Fig S19 from the complete frozen-output diagnostics.

No inference or gene filtering is performed here. The two panels display
empirical signal ranks and exhaustive nonempty hit-intersection patterns.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from PIL import Image

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
VARIANTS = ["deposited", "qn10", "qn20"]
LABELS = {"deposited": "Deposited", "qn10": "QN10", "qn20": "QN20"}
COLORS = {"deposited": "#176D70", "qn10": "#B56F31", "qn20": "#677995"}
WIDTH_MM, HEIGHT_MM = 174, 88


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qa-tools", type=Path, default=Path("D:/ResearchData/pfhn_genome39_2026-10-13/figures/qa_tools"))
    parser.add_argument("--audit-python", type=Path, default=Path("C:/Users/AI/AppData/Local/Programs/Python/Python313/python.exe"))
    args = parser.parse_args()
    assert (args.qa_tools / "audit_panel_alignment.py").is_file()
    sys.path.insert(0, str(args.qa_tools))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    qa = ROOT / "figure_qa"
    qa.mkdir(exist_ok=True)
    frames, input_hashes = {}, {}
    for v in VARIANTS:
        path = ROOT / f"limma_hits_{v}.csv"
        input_hashes[path.name] = sha(path)
        frames[v] = pd.read_csv(path, dtype={"rat_entrez": str})
        assert frames[v].limma_BH_hit.all()
    pattern_path = ROOT / "hit_membership_patterns.csv"
    input_hashes[pattern_path.name] = sha(pattern_path)
    patterns = pd.read_csv(pattern_path)
    displayed = patterns.loc[(patterns.gene_count > 0) & (patterns[["deposited_hit", "qn10_hit", "qn20_hit"]].sum(axis=1) > 0)].copy()
    displayed = displayed.sort_values(["deposited_hit", "qn10_hit", "qn20_hit"], ascending=False)
    assert displayed[["deposited_hit", "qn10_hit", "qn20_hit"]].values.tolist() == [[1, 1, 1], [0, 1, 1]]
    assert displayed.gene_count.tolist() == [24, 17]
    displayed.to_csv(qa / "S19b_nonempty_intersection_source.csv", index=False)
    pd.concat([frames[v][["normalization", "rat_entrez", "gene_symbol", "AveExpr", "signal_empirical_percentile"]] for v in VARIANTS]).to_csv(qa / "S19a_all_hit_percentiles_source.csv", index=False)

    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.8,
        "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7,
        "text.color": "#223744", "axes.labelcolor": "#223744", "xtick.color": "#33434D", "ytick.color": "#33434D",
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4))
    a = fig.add_axes([.105, .22, .40, .63])
    b = fig.add_axes([.74, .22, .225, .63])

    a.plot([0, 100], [0, 100], color="#A4ADB3", linewidth=1, linestyle="--", label="Uniform rank reference", zorder=1)
    styles = {"deposited": ("-", 1.5), "qn10": ("-", 1.4), "qn20": ((0, (3, 2)), 1.2)}
    for v in VARIANTS:
        x = np.sort(frames[v].signal_empirical_percentile.to_numpy())
        y = 100 * np.arange(1, len(x) + 1) / len(x)
        a.step(np.r_[0, x, 100], np.r_[0, y, 100], where="post", color=COLORS[v],
               linestyle=styles[v][0], linewidth=styles[v][1], label=f"{LABELS[v]} (n = {len(x)})", zorder=3)
    a.set_xlim(0, 100); a.set_ylim(0, 100)
    a.set_xticks([0, 25, 50, 75, 100]); a.set_yticks([0, 25, 50, 75, 100])
    a.set_xlabel("Mean-signal percentile within full gene universe", labelpad=7)
    a.set_ylabel("Cumulative percentage of hits", labelpad=6)
    a.set_title("Full hit signal-rank distributions", loc="left", pad=12)
    a.legend(loc="upper left", bbox_to_anchor=(.015, .97), handlelength=1.8, handletextpad=.5)

    b.barh([1, 0], [24, 17], height=.47, color=["#176D70", "#B56F31"])
    b.set_xlim(0, 30); b.set_ylim(-.75, 1.75)
    b.set_xticks([0, 10, 20, 30])
    b.set_yticks([1, 0], ["All three\nnormalizations", "QN10 + QN20\nonly"])
    b.tick_params(axis="y", length=0, pad=8)
    for y, n in [(1, 24), (0, 17)]:
        b.text(n + .7, y, str(n), va="center", fontsize=9, fontweight="bold")
    b.set_xlabel("Genes in each intersection", labelpad=7)
    b.set_title("Hit membership", loc="left", pad=12)
    b.spines["left"].set_visible(False)
    for ax, label in [(a, "a"), (b, "b")]:
        trans = ax.transAxes + ScaledTranslation(-22 / 72, 17 / 72, fig.dpi_scale_trans)
        ax.text(0, 1, label, transform=trans, fontsize=10, fontweight="bold", va="bottom", clip_on=False)
    fig.canvas.draw()
    align = require_matplotlib_panel_alignment(fig, axes=[a, b], panel_ids=["a", "b"], row_groups=[{"panels": ["a", "b"]}],
        json_out=qa / "FigS19.alignment.json", overlay_svg=qa / "FigS19.alignment.svg", tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=True, strict=True)
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(ROOT / f"FigS19.{ext}", dpi=600)
    plt.close(fig)
    subprocess.run([str(args.audit_python), "-c", "import fitz; print(fitz.VersionBind)"], capture_output=True, text=True, check=True)
    audit_commands = [
        [str(args.audit_python), str(args.qa_tools / "audit_pdf_text.py"), str(ROOT / "FigS19.pdf"), "--min-pt", "7", "--json"],
        [str(args.audit_python), str(args.qa_tools / "audit_figure_collisions.py"), str(ROOT / "FigS19.pdf"), "--json-out", str(qa / "FigS19.collision.json"), "--overlay-pdf", str(qa / "FigS19.collision_overlay.pdf")],
    ]
    audits = []
    for i, command in enumerate(audit_commands):
        run = subprocess.run(command, capture_output=True, text=True, encoding="utf-8")
        (qa / ("FigS19.text_audit.json" if i == 0 else "FigS19.collision.log")).write_text(run.stdout + ("\n" + run.stderr if run.stderr else ""), encoding="utf-8")
        audits.append({"command": command, "exit_code": run.returncode})
    with Image.open(ROOT / "FigS19.png") as im:
        raster = {"size_pixels": list(im.size), "dpi": list(im.info.get("dpi", []))}
        assert all(abs(x - 600) < .1 for x in im.info["dpi"])
    for name, digest in input_hashes.items():
        assert sha(ROOT / name) == digest
    legend = (
        "Supplementary Figure S19. Signal-rank distribution and normalization stability of existing limma hits. "
        "a, Empirical cumulative distributions of the mean-signal percentiles of every gene with limma BH q <0.05 in the deposited (24 genes), QN10 (41) and QN20 (41) analyses. "
        "Each gene's percentile is calculated within the complete 14,604-gene mean-log2-signal distribution of the same normalization, using the fraction of gene means less than or equal to its own. "
        "The dashed diagonal is a uniform-rank reference, not a fitted null model; the QN10 and QN20 curves are closely overlapping. "
        "b, All nonempty intersections of the three hit sets. Twenty-four genes occur in all three lists and 17 occur in both QN lists only; no other hit-membership patterns occur. "
        "The union contains 41 genes, all with the same nonzero sesamin-minus-control effect sign across the three analyses. "
        "These post-hoc diagnostics neither refit models nor change the tested gene family. Signal ranks are not detection thresholds. "
        "All three analyses use the same five control and five sesamin animals; overlap is analytical stability, not independent replication. "
        "No new statistical test, confidence interval, expression filter or biological validation is introduced."
    )
    (ROOT / "FigS19_legend.txt").write_text(legend + "\n", encoding="utf-8")
    dump(qa / "figure_execution.json", {"utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
        "python": sys.executable, "matplotlib": mpl.__version__, "input_sha256": input_hashes,
        "script_sha256": sha(__file__), "width_mm": WIDTH_MM, "height_mm": HEIGHT_MM, "png": raster,
        "alignment": align.get("verdict"), "audits": audits, "actual_visual_inspection": "pending",
        "data_or_statistical_changes": False, "output_sha256": {f"FigS19.{ext}": sha(ROOT / f"FigS19.{ext}") for ext in ["pdf", "svg", "png"]}})
    if any(a["exit_code"] for a in audits):
        raise SystemExit("Rendered figure QA requires inspection; see figure_qa.")
    print(json.dumps({"figure": "FigS19", "alignment": align.get("verdict"), "audits_passed": True, "png": raster}))


if __name__ == "__main__":
    main()
