"""Render descriptive PFHN figures using frozen real observations only."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import ScaledTranslation
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
DATA = HERE / 'source_data'
QA = HERE / 'qa'
sys.path.insert(0, str(HERE / 'scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
width_mm = 174
WIDTH = width_mm / 25.4
INK = '#223744'
BLUE = '#39748F'
ORANGE = '#B96F3E'
GREY = '#86949C'
LIGHT = '#DCE3E7'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,
    'legend.fontsize':8,'axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,
    'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'figure.facecolor':'white',
    'savefig.facecolor':'white','text.color':INK,'axes.labelcolor':INK,'axes.edgecolor':GREY})

def load(name):
    return pd.read_csv(DATA / name)


def label(ax, letter):
    ax.text(0,1,letter,transform=ax.transAxes+ScaledTranslation(-19/72,12/72,ax.figure.dpi_scale_trans),
            fontsize=10,weight='bold',ha='left',va='bottom')


def export(fig, stem, axes, supplement=False, exemptions=()):
    out = HERE / ('supplement' if supplement else 'main') / stem
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,axes=axes,panel_ids=list('abcdefghijklmnopqrstuvwxyz')[:len(axes)],
        json_out=QA / f'{stem}.alignment.json',overlay_svg=QA / f'{stem}.alignment.svg',
        tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=len(axes)>1,strict=True,exemptions=exemptions)
    fig.savefig(out.with_suffix('.pdf'),dpi=600)
    fig.savefig(out.with_suffix('.svg'),dpi=600)
    fig.savefig(out.with_suffix('.png'),dpi=600)
    fig.savefig(out.with_suffix('.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    with Image.open(out.with_suffix('.png')) as im:
        assert min(im.info['dpi']) > 599 and abs(im.width/600*25.4-width_mm)<.05
        preview = im.copy()
        preview.thumbnail((1400,1700))
        preview.save(QA / f'{stem}.preview.png')
    plt.close(fig)
    print(stem, 'exported at 174 mm and 600 dpi')


def figure3():
    d=load('organoid_umap_AREG_all24034.csv.gz')
    assert len(d)==24034 and d.cell_key.nunique()==24034 and d.celltype_scTypedb.nunique()==10
    assert np.isfinite(d[['umap_1','umap_2','AREG__counts','AREG__log1p_cp10k']]).all().all()
    assert np.allclose(d.AREG__log1p_cp10k,np.log1p(1e4*d.AREG__counts/d.raw_total_counts),atol=1e-10)
    cats=['AH','FH1','FH2','HB1','HB2','Ductal cells','Cholangiocytes','Hepatic stellate cells','Fibroblasts','Smooth muscle cells']
    colors=['#426D90','#5597A9','#94C4CA','#687688','#A5B5C1','#AC7747','#D09154','#876C96','#9F9F69','#6B987E']
    fig,axes=plt.subplots(1,2,figsize=(WIDTH,106/25.4))
    fig.subplots_adjust(left=.07,right=.90,bottom=.31,top=.87,wspace=.18)
    a,b=axes
    for cat,color in zip(cats,colors):
        v=d[d.celltype_scTypedb==cat]
        a.scatter(v.umap_1,v.umap_2,s=.8,c=color,linewidths=0,rasterized=True)
    order=np.argsort(d.AREG__log1p_cp10k.to_numpy(),kind='stable')
    cmap=LinearSegmentedColormap.from_list('AREG',['#E8ECEF','#7CACB3','#225A72'])
    sc=b.scatter(d.umap_1.to_numpy()[order],d.umap_2.to_numpy()[order],s=.8,
                 c=d.AREG__log1p_cp10k.to_numpy()[order],cmap=cmap,vmin=0,
                 vmax=float(d.AREG__log1p_cp10k.max()),linewidths=0,rasterized=True)
    xlo,xhi=d.umap_1.min()-.8,d.umap_1.max()+.8
    ylo,yhi=d.umap_2.min()-.8,d.umap_2.max()+.8
    for ax,letter,head in zip(axes,'ab',['Cell subtypes','AREG expression']):
        ax.set(xlim=(xlo,xhi),ylim=(ylo,yhi),xlabel='UMAP 1',ylabel='UMAP 2')
        ax.set_aspect('equal',adjustable='box');ax.set_xticks([]);ax.set_yticks([])
        ax.set_title(head,pad=13);label(ax,letter)
    cax=fig.add_axes([.925,.375,.018,.35])
    cb=fig.colorbar(sc,cax=cax);cb.set_label('log1p(CP10k)',labelpad=5,fontsize=8)
    handles=[Line2D([0],[0],marker='o',linestyle='',markersize=4,color=color,label=cat) for cat,color in zip(cats,colors)]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.04,.225),ncol=3,columnspacing=1.5,handletextpad=.4,labelspacing=.7)
    counts=d.groupby('celltype_scTypedb',sort=False).agg(n_cells=('cell_key','size'),AREG_detected=('AREG__counts',lambda x:int((x>0).sum())))
    counts.to_csv(DATA/'Figure_03_subtype_counts.csv')
    export(fig,'Fig4',list(axes))


if __name__ == '__main__':
    figure3()
