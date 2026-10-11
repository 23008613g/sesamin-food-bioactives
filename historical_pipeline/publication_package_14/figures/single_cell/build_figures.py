"""Six descriptive single-cell figures from real, independently derived tables.

Only author-provided UMAP coordinates are displayed. No embedding is computed.
"""
from __future__ import annotations
import hashlib
import json
import math
import platform
from pathlib import Path
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
import PIL

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[3]
INPUT = PROJECT / "next_steps/single_cell_extension_14/analysis"
SKILLS = Path.home() / ".agents/skills/nature-figure/scripts"
sys.path.insert(0, str(SKILLS))
from audit_panel_alignment import require_matplotlib_panel_alignment

INK = "#263443"
NEUTRAL = "#D9DFE3"
SUBTYPE_ORDER = []  # assigned from the analysis metadata, preserving author labels
PALETTE = ["#276785", "#63A1AF", "#9BCCB8", "#597A55", "#A8B47C",
           "#BB8655", "#D7B578", "#9C739B", "#BC95B3", "#858C96"]
SEQ = LinearSegmentedColormap.from_list("Expression", ["#EDF2F5", "#A6CAD7", "#4C91AA", "#155575", "#07354B"])
FEATURE = LinearSegmentedColormap.from_list("FeatureExpression", ["#DBE1E5", "#97C3D2", "#3D91AC", "#07516B"])
LIBRARY_ORDER = ["GSM6322577", "GSM6322578", "GSM6322579", "GSM6322580"]
SHORT = dict(zip(LIBRARY_ORDER, ["C1", "C2", "P1", "P2"]))
LIBRARY_LABEL = dict(zip(LIBRARY_ORDER, ["Control 1", "Control 2", "PA 1", "PA 2"]))
GENES = ["VNN1", "AREG", "JAG2", "CXCL8"]
DISPLAY = {"Hepatic stellate cells":"HSC", "Smooth muscle cells":"SMC"}
SEED = 20261009  # draw-order only; no data generation or cell downsampling
WIDTH = 180 / 25.4
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


def write_json(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def panel(ax, letter):
    ax.text(0, 1, letter, transform=ax.transAxes + ScaledTranslation(-18 / 72, 12 / 72, ax.figure.dpi_scale_trans),
            ha="left", va="bottom", weight="bold", fontsize=10)


def header(fig, title, subtitle):
    fig.text(.06, .982, title, fontsize=10, weight="bold", va="top")
    fig.text(.06, .951, subtitle, fontsize=7.5, va="top")


def note(fig, text, y=.02):
    fig.text(.06, y, text, fontsize=7, va="bottom", linespacing=1.5)


def export(fig, stem, plot_axes, panel_ids=None):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, axes=plot_axes,
        panel_ids=panel_ids, json_out=HERE / f"{stem}.alignment.json",
        overlay_svg=HERE / f"{stem}.alignment.svg", tolerance_pt=1.5,
        gutter_tolerance_pt=1.5, require_panel_labels=len(plot_axes)>1, strict=True)
    fig.savefig(HERE / f"{stem}.pdf", dpi=600)
    fig.savefig(HERE / f"{stem}.svg", dpi=600)
    fig.savefig(HERE / f"{stem}.png", dpi=600)
    fig.savefig(HERE / f"{stem}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def umap_axes(ax, limits):
    ax.set_xlim(*limits[0]); ax.set_ylim(*limits[1])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Author UMAP 1", labelpad=4)
    ax.set_ylabel("Author UMAP 2", labelpad=4)
    ax.tick_params(length=2.5, width=.6)


def legend_handles(order, colors):
    return [Line2D([0], [0], marker="o", linestyle="none", markersize=4.5,
            markerfacecolor=colors[x], markeredgewidth=0, label=DISPLAY.get(x,x)) for x in order]


def scatter_subtypes(ax, frame, colors, order):
    draw = frame.iloc[np.random.default_rng(SEED).permutation(len(frame))]
    ax.scatter(draw["umap_1"], draw["umap_2"], c=[colors[x] for x in draw["subtype"]],
               s=1.45, linewidths=0, alpha=.82, rasterized=True)


def make_fig7(cells, colors, order, limits):
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 133 / 25.4))
    fig.subplots_adjust(left=.09, right=.98, bottom=.29, top=.82, wspace=.28)
    header(fig, "Single-cell context of the liver organoid samples", "24,034 author-retained cells | Four control/PA technical libraries")
    scatter_subtypes(axes[0], cells, colors, order)
    draw = cells.iloc[np.random.default_rng(SEED).permutation(len(cells))]
    cond_colors = {"Control": "#658DA8", "PA": "#C48A55"}
    axes[1].scatter(draw["umap_1"], draw["umap_2"], c=[cond_colors[x] for x in draw["condition"]],
                    s=1.45, linewidths=0, alpha=.75, rasterized=True)
    for ax, title, letter in zip(axes, ["Author cell subtypes", "Experimental condition"], "ab"):
        umap_axes(ax, limits); ax.set_title(title, pad=12); panel(ax, letter)
    fig.legend(handles=legend_handles(order, colors), loc="lower center", bbox_to_anchor=(.515,.115),
               ncol=5, fontsize=7.5, columnspacing=1.5, handletextpad=.35, labelspacing=.7)
    axes[1].legend(handles=legend_handles(list(cond_colors), cond_colors), loc="upper left",
                   bbox_to_anchor=(.1,-.17), ncol=2, fontsize=7.5, handletextpad=.3)
    note(fig, "Author-provided atlas UMAP; four-library subset, not recomputed. No sesamin intervention.\nPA = palmitate exposure; HSC = hepatic stellate cells; SMC = smooth muscle cells.")
    export(fig, "fig7_sc_umap", list(axes), list("ab"))


def make_fig8(cells, limits):
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH, 181 / 25.4))
    fig.subplots_adjust(left=.09, right=.86, bottom=.15, top=.87, wspace=.31, hspace=.31)
    header(fig, "Candidate expression on the author-provided UMAP", "Four highlighted genes | All 24,034 retained cells | Shared full-range color scale")
    maximum = float(cells[[f"{g}_log1p_cp10k" for g in GENES]].to_numpy().max())
    norm = Normalize(0, maximum)
    for ax, gene, letter in zip(axes.flat, GENES, "abcd"):
        values = cells[f"{gene}_log1p_cp10k"].to_numpy()
        idx = np.argsort(values, kind="stable")
        scatter = ax.scatter(cells["umap_1"].to_numpy()[idx], cells["umap_2"].to_numpy()[idx],
                  c=values[idx], s=1.7, cmap=FEATURE, norm=norm, linewidths=0, rasterized=True)
        umap_axes(ax, limits)
        positive = int((values > 0).sum())
        ax.set_title(f"{gene}  |  {positive:,} cells ({positive/len(cells)*100:.2f}%) detected", fontsize=8, pad=11)
        panel(ax, letter)
    cax = fig.add_axes([.895,.33,.018,.33])
    cb=fig.colorbar(scatter, cax=cax)
    cb.set_label("log1p(CP10k)", fontsize=7.5, labelpad=7)
    cb.ax.tick_params(labelsize=7)
    note(fig, "CP10k uses all-gene raw UMI totals per cell; higher expression is drawn last.\nAll candidates share one scale; no clipping, imputation or per-gene rescaling.\nAuthor-provided coordinates; descriptive control/PA context, without sesamin.")
    export(fig, "fig8_sc_features", list(axes.flat), list("abcd"))


def make_fig9(dots, genes, order, counts):
    fig, ax = plt.subplots(figsize=(WIDTH, 193 / 25.4))
    fig.subplots_adjust(left=.17, right=.78, bottom=.13, top=.89)
    header(fig, "All 25 candidates across author-annotated cell subtypes", "Dot area = detected-cell percentage | Color = mean single-cell log1p(CP10k)")
    vmax = float(dots["mean_log1p_cp10k"].max())
    norm=Normalize(0, vmax)
    for yi, gene in enumerate(genes):
        for xi, subtype in enumerate(order):
            row = dots.loc[(dots.gene==gene) & (dots.subtype==subtype)].iloc[0]
            if pd.isna(row.mean_log1p_cp10k):
                ax.scatter([xi],[yi],marker="x",s=14,c="#919AA1",linewidths=.6)
            elif row.detection_pct > 0:
                ax.scatter([xi],[yi],s=float(row.detection_pct)*.85, c=[float(row.mean_log1p_cp10k)],
                          norm=norm, cmap=SEQ, linewidths=.35, edgecolors="#708A97")
            else:
                ax.scatter([xi],[yi],s=6,facecolors="none",edgecolors="#ADB6BC",linewidths=.5)
    ax.set_xlim(-.65, len(order)-.35); ax.set_ylim(len(genes)-.35,-.65)
    ticks={"Cholangiocytes":"Cholangio-\ncytes", "Ductal cells":"Ductal\ncells", "Fibroblasts":"Fibro-\nblasts", **DISPLAY}
    ax.set_xticks(range(len(order)), [ticks.get(x,x) for x in order], fontsize=7)
    ax.set_yticks(range(len(genes)), genes, fontsize=7.5)
    ax.tick_params(length=0,pad=5)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax.set_xticks(np.arange(-.5,len(order),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(genes),1),minor=True)
    ax.grid(which="minor",color="#EEF1F3",linewidth=.5,zorder=-1)
    ax.tick_params(which="minor",length=0)
    cax=fig.add_axes([.835,.57,.02,.24])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=SEQ),cax=cax)
    cb.set_label("Mean log1p(CP10k)",fontsize=7.5,labelpad=8)
    legend_handles_size=[ax.scatter([],[],s=x*.85,c="#82AEBD",edgecolors="#708A97",linewidths=.35,label=f"{x}%") for x in [10,25,50,100]]
    fig.legend(handles=legend_handles_size,loc="center left",bbox_to_anchor=(.8,.405),title="Detected cells",
               title_fontsize=7.5,fontsize=7.5,labelspacing=1.2,handletextpad=.65)
    fig.text(.82,.265,"×  Missing feature\n\n○  0% detected",fontsize=7,linespacing=1.3)
    note(fig,"All cells pooled across four technical libraries; means include measured zero-count cells.\nMissing candidates are not imputed. Candidate order follows the human screen.\nHSC = hepatic stellate cells; SMC = smooth muscle cells. Subtype cell counts: source table.")
    export(fig,"fig9_sc_dotplot",[ax])


def make_s5(cells,colors,order,limits):
    fig,axes=plt.subplots(2,2,figsize=(WIDTH,185/25.4))
    fig.subplots_adjust(left=.09,right=.98,bottom=.20,top=.87,wspace=.28,hspace=.53)
    header(fig,"Library-resolved coverage of the author-provided UMAP","All retained cells in each technical library | Identical coordinate limits and subtype colors")
    for ax,lib,letter in zip(axes.flat,LIBRARY_ORDER,"abcd"):
        frame=cells[cells.library==lib]
        scatter_subtypes(ax,frame,colors,order);umap_axes(ax,limits)
        ax.set_title(f"{LIBRARY_LABEL[lib]} ({lib})\nn = {len(frame):,} cells",fontsize=8.5,pad=10)
        panel(ax,letter)
    fig.legend(handles=legend_handles(order,colors),loc="lower center",bbox_to_anchor=(.515,.075),ncol=5,
               fontsize=7.5,columnspacing=1.5,handletextpad=.35,labelspacing=.7)
    note(fig,"Author-provided atlas coordinates, subset without recomputation; one differentiation.\nHSC = hepatic stellate cells; SMC = smooth muscle cells. Technical-library coverage only.")
    export(fig,"figS5_sc_library_umap",list(axes.flat),list("abcd"))


def make_s6(cells):
    fig,axes=plt.subplots(1,3,figsize=(WIDTH,115/25.4))
    fig.subplots_adjust(left=.085,right=.98,bottom=.26,top=.80,wspace=.50)
    header(fig,"Raw-count quality metrics in retained cells","All 24,034 author-retained cells | No additional quality filter applied")
    metrics=[("total_counts","Total UMIs (log10)",True),("detected_genes","Detected genes (log10)",True),("pct_mito","Mitochondrial UMIs (%)",False)]
    for ax,(column,label,log),letter in zip(axes,metrics,"abc"):
        values=[cells.loc[cells.library==lib,column].to_numpy() for lib in LIBRARY_ORDER]
        if log: values=[np.log10(x) for x in values]
        violin=ax.violinplot(values,positions=np.arange(4),showmedians=False,showextrema=True,widths=.8,points=150,bw_method="scott")
        for i,body in enumerate(violin["bodies"]):
            body.set_facecolor("#658DA8" if i<2 else "#C48A55");body.set_edgecolor("none");body.set_alpha(.55)
        for key in ["cbars","cmins","cmaxes"]:
            violin[key].set_color("#778590");violin[key].set_linewidth(.65)
        ax.boxplot(values,positions=np.arange(4),widths=.14,showfliers=False,manage_ticks=False,whis=(25,75),
                   medianprops={"color":INK,"linewidth":1.1},boxprops={"color":INK,"linewidth":.6},
                   whiskerprops={"color":INK,"linewidth":.6},capprops={"color":INK,"linewidth":.6})
        ax.set_xticks(range(4),[SHORT[x] for x in LIBRARY_ORDER]);ax.set_ylabel(label,labelpad=5)
        ax.set_xlim(-.55,3.55);ax.margins(y=.08);ax.tick_params(length=2.5);panel(ax,letter)
    note(fig,"C1/C2: control libraries; P1/P2: PA libraries. Violin densities use all cells (Scott bandwidth).\nThin lines span the observed ranges; boxes show the interquartile range and median.\nQC summaries describe retained cells, not independent biological replicates; no P values are calculated.")
    export(fig,"figS6_sc_quality",list(axes),list("abc"))


def make_s7(cells,colors,order):
    broad_order=["Hepatocytes","Cholangiocytes","Fibroblasts","Hepatic stellate cells"]
    broad_colors=dict(zip(broad_order,["#4C8498","#B99868","#9C8CA7","#91A58B"]))
    hep_order=[x for x in ["AH","FH1","FH2","HB1","HB2"] if x in order]
    fig,axes=plt.subplots(1,2,figsize=(WIDTH,139/25.4))
    fig.subplots_adjust(left=.09,right=.98,bottom=.30,top=.80,wspace=.30)
    header(fig,"Cell composition in the four source libraries","Cell-count fractions | Denominators are stated separately for each panel")
    source=[]
    for ax,key,cats,cols,title,letter in [(axes[0],"broad_class",broad_order,broad_colors,"All retained cells","a"),
                                         (axes[1],"subtype",hep_order,colors,"Hepatocyte-like cells only","b")]:
        frame=cells if key=="broad_class" else cells[cells.subtype.isin(hep_order)]
        totals=frame.groupby("library").size().reindex(LIBRARY_ORDER).to_numpy()
        bottom=np.zeros(4)
        for cat in cats:
            n=frame[frame[key]==cat].groupby("library").size().reindex(LIBRARY_ORDER,fill_value=0).to_numpy()
            pct=n/totals*100
            assert np.isfinite(pct).all() and (pct>=0).all()
            ax.bar(np.arange(4),pct,bottom=bottom,color=cols[cat],width=.66,edgecolor="white",linewidth=.5,label=DISPLAY.get(cat,cat))
            for i in range(4):
                source.append({"panel":letter,"library":LIBRARY_ORDER[i],"category":cat,"cells":int(n[i]),"denominator_cells":int(totals[i]),"percent":float(pct[i])})
            bottom+=pct
        assert np.allclose(bottom,100)
        ax.set_ylim(0,100);ax.set_xticks(range(4),[f"{SHORT[l]}\n{n:,}" for l,n in zip(LIBRARY_ORDER,totals)])
        ax.set_ylabel("Retained cells (%)",labelpad=4);ax.set_title(title,pad=12);ax.set_xlabel("Library / cell denominator",labelpad=6)
        ax.tick_params(length=2.5);panel(ax,letter)
        ax.legend(loc="upper left",bbox_to_anchor=(-.1,-.26),ncol=2,fontsize=7.5,columnspacing=1,handlelength=1.2,labelspacing=.65)
    note(fig,"C1/C2: control; P1/P2: PA; HSC = hepatic stellate cells. Fractions use cell counts, not UMI totals.\nTechnical libraries from one differentiation; cell recovery and library effects may contribute.")
    pd.DataFrame(source).to_csv(HERE/"figS7_cell_composition_source.csv",index=False)
    export(fig,"figS7_sc_composition",list(axes),list("ab"))


def main():
    cells=pd.read_csv(INPUT/"cell_level_umap_qc_expression.csv.gz")
    cells=cells.rename(columns={"cell_key":"cell_id","gsm":"library","celltype_scTypedb":"subtype",
        "celltype_category":"broad_class","raw_total_counts":"total_counts","raw_n_genes_by_counts":"detected_genes",
        "raw_pct_mito_MT_prefix":"pct_mito", **{f"{g}__log1p_cp10k":f"{g}_log1p_cp10k" for g in GENES}})
    cells["condition"]=cells["condition"].replace({"CTRL":"Control"})
    dots=pd.read_csv(INPUT/"candidate25_subtype_all10_pooled.csv").rename(columns={
        "symbol":"gene","celltype_scTypedb":"subtype","detected_pct":"detection_pct"})
    metadata=json.loads((INPUT/"analysis_summary.json").read_text(encoding="utf-8"))
    order=["AH","FH1","FH2","HB1","HB2","Cholangiocytes","Ductal cells","Fibroblasts","Hepatic stellate cells","Smooth muscle cells"]
    genes=list(dict.fromkeys(dots.gene))
    assert len(cells)==24034 and cells.cell_id.is_unique
    assert len(genes)==25 and len(order)==10 and len(dots)==250
    assert set(cells.library)==set(LIBRARY_ORDER)
    assert set(cells.subtype)==set(order)
    assert set(cells.condition)=={"Control","PA"}
    assert np.isfinite(cells[["umap_1","umap_2","total_counts","detected_genes","pct_mito"]]).all().all()
    assert (cells.total_counts>0).all() and (cells.detected_genes>0).all()
    assert cells.pct_mito.between(0,100).all()
    assert np.isfinite(cells[[f"{g}_log1p_cp10k" for g in GENES]]).all().all()
    assert (cells[[f"{g}_log1p_cp10k" for g in GENES]]>=0).all().all()
    assert int(cells.subtype.isin(["AH","FH1","FH2","HB1","HB2"]).sum())==18097
    colors=dict(zip(order,PALETTE))
    mins=cells[["umap_1","umap_2"]].min().to_numpy()
    maxs=cells[["umap_1","umap_2"]].max().to_numpy()
    span=maxs-mins
    limits=[(float(lo-.04*s),float(hi+.04*s)) for lo,hi,s in zip(mins,maxs,span)]
    assert (span>0).all()
    copied=[]
    for filename in ["cell_level_umap_qc_expression.csv.gz","candidate25_subtype_all10_pooled.csv","analysis_summary.json"]:
        src=INPUT/filename;dst=HERE/filename
        shutil.copy2(src,dst);assert digest(src)==digest(dst)
        copied.append({"name":filename,"source":str(src),"sha256":digest(src),"bytes":src.stat().st_size})
    # Independently verify all four candidate transformations and all 250 dot entries.
    normalization_residuals={}
    for g in GENES:
        reconstructed=np.log1p(cells[f"{g}__counts"].to_numpy()/cells.total_counts.to_numpy()*10000)
        residual=float(np.max(np.abs(reconstructed-cells[f"{g}_log1p_cp10k"].to_numpy())))
        assert residual<1e-10
        normalization_residuals[g]=residual
    for row in dots.itertuples():
        frame=cells[cells.subtype==row.subtype]
        assert len(frame)==row.n_cells
        if not row.covered:
            assert row.gene=="AOC4P" and pd.isna(row.mean_log1p_cp10k) and pd.isna(row.detection_pct)
        else:
            raw=frame[f"{row.gene}__counts"].to_numpy()
            normed=np.log1p(raw/frame.total_counts.to_numpy()*10000)
            assert np.isclose(normed.mean(),row.mean_log1p_cp10k,rtol=1e-9,atol=1e-12)
            assert np.isclose((raw>0).mean()*100,row.detection_pct,rtol=1e-9,atol=1e-10)
    write_json("plot_parameters.json",{"backend":"Python/Matplotlib","seed":SEED,"seed_scope":"point draw order only",
        "embedding":"author-provided atlas UMAP, four-library subset, not recomputed","cell_count":len(cells),
        "umap_limits":limits,"subtype_order":order,"subtype_colors":colors,"candidate_order":genes,
        "excluded_cells":0,"expression":"log1p(CP10k), all-gene raw UMI denominator","figure_width_mm":180,
        "dpi":600,"rasterized_marks":True,"editable_vector_text":True})
    write_json("source_data.json",{"inputs":copied,"plot_environment":{"executable":sys.executable,"python":platform.python_version(),
        "numpy":np.__version__,"pandas":pd.__version__,"matplotlib":matplotlib.__version__,"Pillow":PIL.__version__}})
    selected=set(sys.argv[1:]) or {"7","8","9","S5","S6","S7"}
    assert selected.issubset({"7","8","9","S5","S6","S7"})
    if "7" in selected: make_fig7(cells,colors,order,limits)
    if "8" in selected: make_fig8(cells,limits)
    if "9" in selected: make_fig9(dots,genes,order,cells.groupby("subtype").size())
    if "S5" in selected: make_s5(cells,colors,order,limits)
    if "S6" in selected: make_s6(cells)
    if "S7" in selected: make_s7(cells,colors,order)
    write_json("numerical_qa.json",{"status":"PASS","all_24034_cells_used":True,"hepatocyte_like_cells":18097,
        "all_25_candidates_preserved":True,"subtypes":10,"dot_entries":250,
        "finite_coordinates":True,"coordinate_transform":"none (equal aspect)","no_color_clipping":True,
        "draw_order_seed":SEED,"no_data_generation":True,"no_embedding_recomputation":True,
        "four_gene_transform_max_residuals":normalization_residuals,
        "dot_summaries_independently_recomputed_from_cells":250,
        "source_files_unchanged":all(digest(INPUT/x["name"])==x["sha256"] for x in copied)})
    print(json.dumps({"status":"rendered","figures":len(selected),"selected":sorted(selected),"cells":len(cells)}))


if __name__=="__main__":
    main()
