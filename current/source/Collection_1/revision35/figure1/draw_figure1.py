"""Original vector evidence map; source semantics match the archived Fig. 1.

Pictograms are qualitative symbols, not measurements or observed cell images.
No raw data, analysis values, or selection rules are modified.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Polygon, Rectangle, PathPatch
from matplotlib.path import Path as MplPath


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("--audit-scripts", type=Path, default=HERE/"qa_tools")
parser.add_argument("--output-dir", type=Path, default=ROOT/"figures"/"main")
args = parser.parse_args()
sys.path.insert(0, str(args.audit_scripts))
from audit_panel_alignment import require_matplotlib_panel_alignment

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.3, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "axes.linewidth": 0.8,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})

W, H = 493.2, 338.4
INK = "#223A49"
MUTED = "#52616B"
TEAL = "#146D75"
TEAL_LIGHT = "#EAF4F4"
TEAL_EDGE = "#BDD6DA"
AMBER = "#95611E"
AMBER_LIGHT = "#FBF5EB"
AMBER_EDGE = "#E3D4BC"
GREY_LIGHT = "#F2F5F7"
GREY_EDGE = "#CFD8DE"

fig = plt.figure(figsize=(6.85, 4.7), dpi=300)
canvas = fig.add_axes([0, 0, 1, 1], label="canvas")
canvas.set(xlim=(0, W), ylim=(0, H)); canvas.axis("off")
axes = []
labels = {}

def text(ax, x, y, value, size=8.3, color=INK, weight="normal", linespacing=1.22):
    return ax.text(x, y, value, fontsize=size, color=color, fontweight=weight,
                   va="center", ha="left", linespacing=linespacing)

def section(y, value, color):
    canvas.add_patch(Rectangle((10, y - 4.8), 3, 9.6, color=color, linewidth=0))
    text(canvas, 20, y, value, 10.4, color, "bold")

def panel(key, x, y, w, h, color, fill, edge, bar=True):
    ax = fig.add_axes([x/W, y/H, w/W, h/H], label=key)
    ax.set(xlim=(0,w), ylim=(0,h)); ax.axis("off")
    ax.add_patch(Rectangle((0.5,0.5),w-1,h-1,facecolor="white",edgecolor=edge,linewidth=.65))
    if bar:
        ax.add_patch(Rectangle((.8,h-31),w-1.6,30,facecolor=fill,edgecolor="none"))
        ax.add_patch(Rectangle((.8,h-2.2),w-1.6,1.6,facecolor=color,edgecolor="none"))
    axes.append(ax); labels[ax]=key
    return ax

def cell(ax, x, y, r, color, nucleus=True):
    ax.add_patch(Circle((x,y),r,facecolor="white",edgecolor=color,linewidth=.85))
    if nucleus:
        ax.add_patch(Circle((x-.7,y+.4),r*.28,facecolor=color,edgecolor="none"))

def dish(ax,x,y,c):
    ax.add_patch(Ellipse((x,y),22,11,facecolor="white",edgecolor=c,linewidth=.85))
    ax.plot([x-11,x-11,x-8,x+8,x+11,x+11],[y,y-3,y-5,y-5,y-3,y],color=c,lw=.85)
    for dx,dy in [(-5,0),(1,2),(5,-1)]:
        cell(ax,x+dx,y+dy,1.6,c,False)

def liver(ax,x,y,c,atlas=False):
    verts=[(x-12,y-2),(x-11,y+6),(x-2,y+8),(x+11,y+5),(x+12,y+1),
           (x+7,y-4),(x+1,y-5),(x-2,y-2),(x-5,y-7),(x-12,y-2)]
    codes=[MplPath.MOVETO]+[MplPath.CURVE3]*8+[MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts,codes),facecolor="white",edgecolor=c,lw=.95))
    ax.plot([x-2,x-1,x+1],[y+5,y,y-4],color=c,lw=.65)
    if atlas:
        for dx,dy in [(-7,0),(-4,3),(3,2),(7,1)]:
            ax.add_patch(Circle((x+dx,y+dy),.85,facecolor=c,edgecolor="none"))

def hepatocyte(ax,x,y,c):
    for dx,dy in [(-5,0),(5,1)]:
        pts=[(x+dx-5,y+dy-3),(x+dx-5,y+dy+3),(x+dx,y+dy+6),
             (x+dx+5,y+dy+3),(x+dx+5,y+dy-3),(x+dx,y+dy-6)]
        ax.add_patch(Polygon(pts,closed=True,facecolor="white",edgecolor=c,lw=.85))
        ax.add_patch(Circle((x+dx,y+dy),1.4,facecolor=c,edgecolor="none"))

def organoid(ax,x,y,c):
    # Symbolic hollow cluster; the number of circles is not a measured cell count.
    for dx,dy in [(0,6),(5.2,3),(5.2,-3),(0,-6),(-5.2,-3),(-5.2,3)]:
        cell(ax,x+dx,y+dy,3.3,c,False)

def mouse(ax,x,y,c):
    ax.add_patch(Ellipse((x,y),15,9,facecolor="white",edgecolor=c,lw=.85))
    ax.add_patch(Ellipse((x+7,y+1),7,6,facecolor="white",edgecolor=c,lw=.85))
    ax.add_patch(Circle((x+6,y+4),2.5,facecolor="white",edgecolor=c,lw=.85))
    ax.add_patch(Circle((x+9,y+1),.6,facecolor=c,edgecolor="none"))
    ax.plot([x-7,x-11,x-13,x-14],[y,y-1,y+1,y+3],color=c,lw=.85)

section(326, "Sesamin-related expression comparisons", TEAL)
pw = 229.6
human=panel("human",10,231,pw,83,TEAL,TEAL_LIGHT,TEAL_EDGE)
rat=panel("rat",253.6,231,pw,83,TEAL,TEAL_LIGHT,TEAL_EDGE)
for ax,title,accession,role,caveat,icon in [
    (human,"Human HepG2 cells","GSE248903  |  9 libraries", "Fixed set of 25 human candidates",
     "Timing and independence unresolved",dish),
    (rat,"Dietary rat liver","GSE42693  |  5 control + 5 sesamin", "Healthy animals; 15-day exposure",
     "Original body-weight selection retained",liver),
]:
    text(ax,12,67,title,10.0,TEAL,"bold"); icon(ax,pw-23,67,TEAL)
    text(ax,12,44,accession,8.4)
    text(ax,12,28,role,8.5)
    text(ax,12,12,caveat,8.1,MUTED)

# Exposure mapping is unresolved: physically separate this neutral strip from
# both the confirmed dietary comparison and the non-sesamin context group.
m=panel("mouse",10,166,473.2,51,MUTED,GREY_LIGHT,GREY_EDGE,bar=False)
m.patches[0].set_facecolor(GREY_LIGHT)
mouse(m,21,34,MUTED)
text(m,42,35,"Mouse source-label comparison",9.4,MUTED,"bold")
text(m,249,35,"CRA013797  |  9 libraries",8.4)
text(m,12,15,"ctrl / HFD / sess: exploration only",8.4)
text(m,249,15,"Exposure mapping remains unconfirmed",8.1,MUTED)

section(148,"Expression and cellular context — no sesamin intervention",AMBER)
cw=(473.2-28)/3
specs=[
    ("phh",10,"Primary human\nhepatocytes","GSE247407","30 libraries; 5 labels",
     "Fatty-acid expression\nbackground","Donor-label relationships\nunresolved",hepatocyte),
    ("organoid",10+cw+14,"Human liver\norganoids","GSE207889","24,034 retained cells",
     "AREG subtypes and\nexpression summaries","One differentiation;\n4 technical libraries",organoid),
    ("atlas",10+2*(cw+14),"Adult human\nliver atlas","GSE115469","8,444 cells; 5 donors",
     "AREG donor and\nimmune-cell context","Author-normalized\nsource expression",liver),
]
for key,x,title,accession,count,role,caveat,icon in specs:
    ax=panel(key,x,18,cw,117,AMBER,AMBER_LIGHT,AMBER_EDGE)
    text(ax,12,100,title,9.2,AMBER,"bold",1.12)
    if key=="atlas": icon(ax,cw-23,100,AMBER,atlas=True)
    else: icon(ax,cw-23,100,AMBER)
    text(ax,12,77,accession,8.4)
    text(ax,12,65,count,8.3)
    text(ax,12,44,role,8.5,INK,linespacing=1.2)
    text(ax,12,16,caveat,8.1,MUTED,linespacing=1.18)

text(canvas,10,7,"Library and cell counts do not establish independent biological replication.",8.0,MUTED)

fig.canvas.draw()
alignment=require_matplotlib_panel_alignment(
    fig, axes=axes, panel_ids=labels,
    row_groups=[["human","rat"],["phh","organoid","atlas"]], column_groups=[],
    json_out=HERE/"Fig1.alignment.json", overlay_svg=HERE/"Fig1.alignment.svg",
    tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True,
)
out=args.output_dir.resolve()
out.mkdir(parents=True,exist_ok=True)
fig.savefig(out/"Fig1.pdf", dpi=600)
fig.savefig(out/"Fig1.svg", dpi=600)
fig.savefig(out/"Fig1.eps", dpi=600)
fig.savefig(out/"Fig1.png", dpi=600)
fig.savefig(out/"Fig1.tiff", dpi=600, pil_kwargs={"compression":"tiff_lzw"})
fig.savefig(HERE/"Fig1.preview.png",dpi=300)

# Text/canvas guard uses final renderer geometry, not character-count guesses.
renderer=fig.canvas.get_renderer()
texts=[]
for ax in [canvas]+axes:
    for artist in ax.texts:
        box=artist.get_window_extent(renderer)
        coords=fig.dpi_scale_trans.inverted().transform([[box.x0,box.y0],[box.x1,box.y1]])*72
        assert coords[0,0]>=0 and coords[0,1]>=0 and coords[1,0]<=W and coords[1,1]<=H
        texts.append({"text":artist.get_text(),"fontsize_pt":artist.get_fontsize(),"bbox_pt":coords.tolist()})
record={"dimensions_in":[6.85,4.7],"dpi":600,"python":sys.version,"matplotlib":mpl.__version__,
        "source_scope":"Schematic dataset roles; no quantitative transforms or new analysis",
        "source_caption":"Current manuscript Fig. 1 legend and Table 1",
        "minimum_source_font_pt":min(t["fontsize_pt"] for t in texts),
        "alignment_verdict":alignment["verdict"], "texts":texts,
        "sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.glob("Fig1.*"))}}
(HERE/"figure1_render_record.json").write_text(json.dumps(record,indent=2),encoding="utf-8")
plt.close(fig)
print(json.dumps({"figure":"Fig1","formats":list(record["sha256"]),"alignment":alignment["verdict"]}))
