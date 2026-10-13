"""PFHN revision figures from verified data: Python/matplotlib only."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from PIL import Image

HERE = Path(__file__).resolve().parent
DATA, QA, MAIN = HERE / 'source_data', HERE / 'qa', HERE / 'main'
sys.path.insert(0, str(HERE / 'scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
WIDTH_MM = 174
WIDTH = WIDTH_MM / 25.4
INK, BLUE, ORANGE, GREY = '#223744', '#39748F', '#B96F3E', '#7D8B92'
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
 'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
 'legend.fontsize': 8, 'axes.linewidth': .65, 'axes.spines.top': False, 'axes.spines.right': False,
 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'legend.frameon': False,
 'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'text.color': INK,
 'axes.labelcolor': INK, 'axes.edgecolor': GREY})
CHECKS = []


def check(name, value):
    CHECKS.append({'check': name, 'passed': bool(value)})
    assert value, name


def read(name):
    return pd.read_csv(DATA / name)


def letter(ax, char, dx=-22, dy=12):
    ax.text(0, 1, char, transform=ax.transAxes + ScaledTranslation(dx/72, dy/72, ax.figure.dpi_scale_trans),
            fontsize=10, fontweight='bold', ha='left', va='bottom')


def export(fig, number, axes, exemptions=()):
    stem = f'Fig{number}'
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, axes=axes, panel_ids=list('abcdefghijklmnopqrstuvwxyz')[:len(axes)],
        json_out=QA / f'{stem}.alignment.json', overlay_svg=QA / f'{stem}.alignment.svg',
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=len(axes)>1,
        strict=True, exemptions=exemptions)
    fig.savefig(MAIN / f'{stem}.pdf', dpi=600)
    fig.savefig(MAIN / f'{stem}.svg', dpi=600)
    fig.savefig(MAIN / f'{stem}.eps', dpi=600)
    fig.savefig(MAIN / f'{stem}.png', dpi=600)
    fig.savefig(MAIN / f'{stem}.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    with Image.open(MAIN / f'{stem}.png') as im:
        check(stem + '_600dpi_width', min(im.info['dpi']) > 599 and abs(im.width/600*25.4-WIDTH_MM)<.05)
        preview = im.copy()
        preview.thumbnail((1450, 1700))
        preview.save(QA / f'{stem}.preview.png')
    plt.close(fig)
    print(stem, 'rendered')


def figure1():
    fig, ax = plt.subplots(figsize=(WIDTH, 113/25.4))
    fig.subplots_adjust(left=.03, right=.97, bottom=.03, top=.97)
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    ax.text(.015, .969, 'Expression comparisons', fontsize=9, weight='bold', color=BLUE)
    ax.text(.536, .969, 'Context without sesamin', fontsize=9, weight='bold', color=ORANGE)
    cards = [
      (.005,.903,'Human HepG2 cells','GSE248903 | 9 libraries','Fixed set of 25 human candidates','Timing and independence unresolved',BLUE),
      (.005,.607,'Dietary rat liver','GSE42693 | 5 control + 5 sesamin','Healthy animals; 15-day exposure','Original body-weight selection retained',BLUE),
      (.005,.311,'Mouse source-label comparison','CRA013797 | 9 libraries','ctrl / HFD / sess: exploration only','Exposure mapping remains unconfirmed',GREY),
      (.525,.903,'Primary human hepatocytes','GSE247407 | 30 libraries; 5 labels','Fatty-acid expression background','Donor-label relationships unresolved',ORANGE),
      (.525,.607,'Human liver organoids','GSE207889 | 24,034 retained cells','AREG subtypes and expression summaries','One differentiation; 4 technical libraries',ORANGE),
      (.525,.311,'Adult human liver atlas','GSE115469 | 8,444 cells; 5 donors','AREG donor and immune-cell context','Author-normalized source expression',ORANGE)]
    for x, y, title, accession, role, limit, color in cards:
        ax.add_patch(FancyBboxPatch((x,y-.249),.467,.249,boxstyle='round,pad=.003,rounding_size=.006',
                                   facecolor='#F7F9FA',edgecolor='#C7D2D8',lw=.7))
        for offset, text, size, tone, weight in [(.023,title,8.5,color,'bold'),(.079,accession,8,INK,'normal'),
                                               (.135,role,8,INK,'normal'),(.191,limit,8,'#56666F','normal')]:
            ax.text(x+.012, y-offset, text, fontsize=size, color=tone, weight=weight, va='top')
    export(fig, 1, [ax])


def figure2():
    h = read('human_candidate_samples.csv').rename(columns={'library':'sample','median_ratio_count':'normalized_count'})
    d = read('candidate_direction_order.csv')
    genes = h[['gene_id','symbol']].drop_duplicates().reset_index(drop=True)
    libraries = h['sample'].drop_duplicates().tolist()
    check('Fig2_complete_25_by_9', len(genes)==25 and len(h)==225 and len(libraries)==9)
    check('Fig2_15_orthologues', len(d)==15)
    counts = h.pivot(index='gene_id', columns='sample', values='normalized_count').loc[genes.gene_id,libraries].to_numpy()
    check('Fig2_counts_finite_nonnegative', np.isfinite(counts).all() and (counts>=0).all())
    transformed = np.log2(counts + 1)
    centered = transformed - transformed.mean(axis=1, keepdims=True)
    check('Fig2_heatmap_not_clipped', np.abs(centered).max()<2)
    fields = ['model_control_sign','treat_model_sign','HFD_ctrl_sign','sess_HFD_sign']
    values = d[fields].to_numpy(dtype=int)
    check('Fig2_zero_directions_retained', int((values==0).sum())==5)
    fig = plt.figure(figsize=(WIDTH, 151/25.4))
    grid = fig.add_gridspec(1,2,left=.235,right=.974,bottom=.17,top=.83,width_ratios=[.275,.33],wspace=.134/.3025)
    a, b = fig.add_subplot(grid[0,0]), fig.add_subplot(grid[0,1])
    cmap = LinearSegmentedColormap.from_list('Expression',['#24728C','#F7F7F4','#B86B2D'])
    im = a.pcolormesh(np.arange(10)-.5,np.arange(26)-.5,centered,cmap=cmap,vmin=-2,vmax=2,shading='flat')
    a.set(xlim=(-.5,8.5),ylim=(24.5,-.5))
    a.set_yticks(range(25),genes.symbol.tolist());a.set_xticks(range(9),['1','2','3']*3)
    for x, name in zip([1,4,7],['Control','Model','Treat']):
        a.text(x,1.053,name,transform=a.get_xaxis_transform(),ha='center',va='bottom',fontsize=8)
    for x in [2.5,5.5]: a.axvline(x,color='white',lw=2)
    b.pcolormesh(np.arange(5)-.5,np.arange(16)-.5,values,
        cmap=ListedColormap(['#D5E6EF','#F0F1F1','#F2DCC7']),norm=BoundaryNorm([-1.5,-.5,.5,1.5],3),shading='flat')
    b.set(xlim=(-.5,4.5),ylim=(14.5,-.5))
    chars={-1:'\u2193',0:'0',1:'\u2191'}
    for i in range(15):
        for j in range(4): b.text(j,i,chars[values[i,j]],ha='center',va='center',fontsize=10)
        b.text(4,i,'Y' if d.iloc[i].sess_between else 'N',ha='center',va='center',fontsize=8)
    b.set_yticks(range(15),[r.human_symbol+(' *' if r.all_libraries_TPM_lt1 else '') for r in d.itertuples()])
    b.set_xticks(range(5),['Model /\nControl','Treat /\nModel','HFD /\nctrl','sess /\nHFD','Strict\nbetween'])
    b.text(.5,1.105,'Human',transform=b.get_xaxis_transform(),ha='center',va='bottom',fontsize=8)
    b.text(3,1.105,'Mouse source labels',transform=b.get_xaxis_transform(),ha='center',va='bottom',fontsize=8)
    for ax,nrows,ncols,ch in [(a,25,9,'a'),(b,15,5,'b')]:
        ax.xaxis.tick_top();ax.tick_params(axis='both',length=0,pad=4,labelsize=8)
        ax.set_xticks(np.arange(-.5,ncols,1),minor=True);ax.set_yticks(np.arange(-.5,nrows,1),minor=True)
        ax.grid(which='minor',color='white',lw=.5);ax.tick_params(which='minor',length=0)
        for spine in ax.spines.values(): spine.set_visible(False)
        letter(ax,ch,dy=48)
    cax = fig.add_axes([a.get_position().x0,.080,a.get_position().width,.013])
    cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-2,0,2]);cb.outline.set_visible(False)
    cb.set_label('Gene-centered log2(count + 1)',fontsize=8,labelpad=4)
    cb.ax.tick_params(labelsize=8,length=2)
    fig.text(.60,.087,'\u2191 Higher   \u2193 Lower   0 Unchanged',fontsize=8)
    fig.text(.60,.052,'Y: sess strictly between ctrl and HFD',fontsize=8)
    fig.text(.60,.023,'* All nine mouse libraries <1 TPM',fontsize=8)
    tidy=h.copy();tidy['log2_count_plus1']=np.log2(tidy.normalized_count+1)
    tidy['centered_log2']=tidy.log2_count_plus1-tidy.groupby('gene_id').log2_count_plus1.transform('mean')
    tidy.to_csv(DATA/'Fig2_heatmap_display_values.csv',index=False)
    export(fig,2,[a,b],[{'panels':['a','b'],'checks':['panel-width'],'reason':'A nine-column expression heatmap and a five-column direction/ordering table have intentionally different widths; their rendered top and bottom edges are shared.'}])


def figure3():
    h=read('human_candidate_samples.csv').rename(columns={'library':'sample','median_ratio_count':'normalized_count'}).query("symbol == 'VNN1'")
    m=read('mouse_candidate_samples.csv').rename(columns={'symbol':'human_symbol','library':'run'}).query("human_symbol == 'VNN1'").copy()
    r=read('rat_samples.csv').query("human_symbol == 'VNN1' and method == 'deposited_log2' and primary_comparison")
    check('Fig3_all_28_observations', len(h)==9 and len(m)==9 and len(r)==10)
    specs=[(h,'normalized_count',['Control','Model','Treat'],['Control','Model','Treat'],'Human HepG2','Normalized counts','sample'),
           (m,'TPM',['ctrl','HFD','sess'],['ctrl','HFD','sess'],'Mouse source labels','TPM','run'),
           (r,'aggregate_log2_signal',['control','sesamin'],['Control','Sesamin'],'Dietary rat liver','log2(signal + 1)','sample')]
    fig,axes=plt.subplots(1,3,figsize=(WIDTH,89/25.4))
    fig.subplots_adjust(left=.105,right=.985,bottom=.32,top=.86,wspace=.68)
    display=[]
    for ax,(d,field,order,ticks,head,unit,idcol),ch in zip(axes,specs,'abc'):
        for j,group in enumerate(order):
            sub=d.loc[d.group==group].sort_values(idcol)
            color=[GREY,ORANGE,BLUE][j] if len(order)==3 else [GREY,BLUE][j]
            ax.scatter(j+np.linspace(-.12,.12,len(sub)),sub[field],s=25,c=color,edgecolors='white',lw=.35,zorder=3)
            ax.hlines(sub[field].mean(),j-.24,j+.24,color=INK,lw=1.1,zorder=4)
            display += [{'panel':ch,'observation_id':row[idcol],'group':group,'value':row[field],'unit':unit} for row in sub.to_dict('records')]
        ax.set_xlim(-.5,len(order)-.5);ax.set_xticks(range(len(order)),ticks)
        ax.set_title(head,pad=13);ax.set_ylabel(unit,labelpad=5);letter(ax,ch)
        ax.tick_params(axis='x',length=0,pad=5)
    axes[0].set_ylim(0,max(h.normalized_count)*1.15)
    axes[1].set_ylim(0,max(m.TPM)*1.15)
    rr=read('rat_results.csv').query("human_symbol == 'VNN1' and method == 'deposited_log2'").iloc[0]
    fig.text(.105,.153,'All points retained; short lines show group means',fontsize=8)
    fig.text(.105,.097,'Rat: Welch BH25 = '+format(rr.welch_q_bh25,'.5f')+'; exact-label BH25 = '+format(rr.exact_q_bh25,'.4f'),fontsize=8)
    fig.text(.105,.041,'Mouse exposure mapping and cell-study biological independence remain unconfirmed',fontsize=8)
    pd.DataFrame(display).to_csv(DATA/'Fig3_all_observations.csv',index=False)
    export(fig,3,list(axes))


POPS=['Whole organoid','Hepatocytes','Deposited broad Cholangiocytes','Cholangiocytes + Ductal cells excluding SMC']
POP_SHORT=['Whole organoid','Hepatocyte-like','Broad category','Excluding SMC']
LIBS=['GSM6322577','GSM6322578','GSM6322579','GSM6322580']
LIB_NAMES=['Control 1','Control 2','PA 1','PA 2']
LIB_MARKERS=['o','s','o','s']
LIB_COLORS=[GREY,GREY,BLUE,BLUE]


def figure4():
    sub=read('areg_subtype_by_library.csv')
    pop=read('areg_population_by_library.csv')
    matched=read('areg_matched_scale.csv').query("weighting == 'equal_library'").set_index('population').loc[POPS]
    check('Fig4_all_four_libraries_all_populations',len(pop)==16 and pop.groupby('population').size().eq(4).all())
    check('Fig4_retains_both_negative_contrasts',(matched.PA_minus_control.iloc[2:]<0).all())
    fig,axes=plt.subplots(2,2,figsize=(WIDTH,161/25.4))
    fig.subplots_adjust(left=.21,right=.98,bottom=.16,top=.91,wspace=.64,hspace=.66)
    a,b,c,d=axes.ravel()
    subtype_order=['Cholangiocytes','Ductal cells','Smooth muscle cells']
    for i,cat in enumerate(subtype_order):
        vals=sub.loc[sub.source_subtype_label==cat].set_index('library_accession').loc[LIBS]
        check('Fig4_subtype_four_nonmissing_libraries_'+cat,len(vals)==4 and vals.raw_detection_pct.notna().all())
        for j,gsm in enumerate(LIBS):
            a.scatter(vals.loc[gsm,'raw_detection_pct'],i+(j-1.5)*.13,s=27,marker=LIB_MARKERS[j],c=LIB_COLORS[j],zorder=3)
    a.set(xlim=(0,103),ylim=(2.55,-.55),xlabel='AREG raw-UMI detection (%)')
    a.set_yticks(range(3),['Cholangiocytes','Ductal cells','Smooth muscle\ncells'])
    for i,population in enumerate(POPS[1:]):
        vals=pop.loc[pop.population==population].set_index('library_accession').loc[LIBS]
        for j,gsm in enumerate(LIBS):
            b.scatter(vals.loc[gsm,'raw_detection_pct'],i+(j-1.5)*.13,s=27,marker=LIB_MARKERS[j],c=LIB_COLORS[j],zorder=3)
    b.set(xlim=(0,103),ylim=(2.55,-.55),xlabel='AREG raw-UMI detection (%)')
    b.set_yticks([0,1,2],['Hepatocyte-like','Broad category','Excluding SMC'])
    c.axvline(0,color=GREY,lw=.7,zorder=0)
    for i,(name,row) in enumerate(matched.iterrows()):
        c.scatter(row.PA_minus_control,i,s=33,c=BLUE if row.PA_minus_control>0 else ORANGE,zorder=3)
    c.set(xlim=(-.16,.08),ylim=(3.5,-.5),xlabel='PA minus control mean')
    c.set_yticks(range(4),POP_SHORT);c.set_xticks([-.15,-.10,-.05,0,.05])
    for i,population in enumerate(POPS):
        vals=pop.loc[pop.population==population].set_index('library_accession').loc[LIBS]
        for j,gsm in enumerate(LIBS):
            d.scatter(i+(j-1.5)*.13,vals.loc[gsm,'mean_log1p_CP10k'],s=27,marker=LIB_MARKERS[j],c=LIB_COLORS[j],zorder=3)
    d.set(xlim=(-.5,3.5),ylim=(0,2.85),ylabel='Mean log1p(CP10k)')
    d.set_xticks(range(4),['Whole','Hep.','Broad','No SMC'])
    for ax,ch,head in zip([a,b,c,d],'abcd',['Original subtypes','Population comparison','Matched-scale contrasts','Technical-library means']):
        ax.set_title(head,pad=14);letter(ax,ch)
        ax.tick_params(axis='y',length=0,pad=5)
    handles=[Line2D([0],[0],marker=LIB_MARKERS[j],linestyle='',markersize=5,color=LIB_COLORS[j],label=LIB_NAMES[j]) for j in range(4)]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.064),ncol=4,columnspacing=1.6,handletextpad=.4)
    fig.text(.21,.040,'Four technical libraries from one differentiation; no biological confidence intervals',fontsize=8)
    fig.text(.21,.015,'Contrasts use equal library weights and the same per-cell log1p(CP10k) scale',fontsize=8)
    export(fig,4,list(axes.ravel()))


def figure5():
    paired=read('adult_donor_context.csv')
    cells=read('adult_celltype_donor.csv')
    summary=read('adult_celltype_summary.csv')
    check('Fig5_all_five_donors',len(paired)==5)
    check('Fig5_P1_reversed_retained',paired.iloc[0].delta_mean_source_log2CPM<0)
    check('Fig5_all_15_cell_types',len(summary)==15)
    types=summary.cell_type.tolist()
    names={'Central_venous_LSECs':'Central venous LSECs','Cholangiocytes':'Cholangiocytes','Erythroid_Cells':'Erythroid cells',
      'Hepatic_Stellate_Cells':'Hepatic stellate cells','Hepatocytes':'Hepatocytes','Inflammatory_Macrophage':'Inflammatory macrophages',
      'Mature_B_Cells':'Mature B cells','NK-like_Cells':'NK-like cells','Non-inflammatory_Macrophage':'Non-inflammatory macrophages',
      'Periportal_LSECs':'Periportal LSECs','Plasma_Cells':'Plasma cells','Portal_endothelial_Cells':'Portal endothelial cells',
      'alpha-beta_T_Cells':'Alpha-beta T cells','gamma-delta_T_Cells_1':'Gamma-delta T cells 1','gamma-delta_T_Cells_2':'Gamma-delta T cells 2'}
    fig=plt.figure(figsize=(WIDTH,174/25.4))
    grid=fig.add_gridspec(2,2,left=.365,right=.97,top=.890,bottom=.12,height_ratios=[1,2.3],hspace=.58,wspace=.49)
    a,b,c=fig.add_subplot(grid[0,0]),fig.add_subplot(grid[0,1]),fig.add_subplot(grid[1,:])
    for i,row in paired.iterrows():
        for ax,hval,cval in [(a,row.hep_nonzero_pct,row.chol_nonzero_pct),(b,row.hep_mean_source_log2CPM,row.chol_mean_source_log2CPM)]:
            ax.plot([hval,cval],[i,i],color='#B1BFC6',lw=.8,zorder=1)
            ax.scatter(hval,i,c=BLUE,s=25,zorder=3,marker='o')
            ax.scatter(cval,i,c=ORANGE,s=25,zorder=3,marker='s')
    a.set(xlim=(-.8,17),ylim=(4.65,-.65),xlabel='Cells (%)');a.set_xticks([0,5,10,15])
    b.set(xlim=(-.008,.17),ylim=(4.65,-.65),xlabel='Mean source log2(CPM)');b.set_xticks([0,.05,.10,.15])
    a.set_yticks(range(5),[r.donor.replace('TLH','')+(' *' if not r.both_types_ge20 else '')+' ('+str(r.n_cholangiocytes)+'/'+str(r.n_hepatocytes)+')' for r in paired.itertuples()])
    b.set_yticks(range(5),[r.donor.replace('TLH','') for r in paired.itertuples()])
    a.set_ylabel('Donor (Chol / Hep cells)',labelpad=7)
    a.set_title('Nonzero source values',pad=13);b.set_title('Mean expression',pad=13)
    for i,celltype in enumerate(types):
        v=cells.loc[cells.cell_type==celltype].sort_values('donor')
        mean=summary.loc[summary.cell_type==celltype].iloc[0]
        check('Fig5_available_donor_mean_'+celltype,np.isclose(v.mean_source_log2CPM.mean(),mean.equal_donor_mean_source_log2CPM)
              and np.isclose(v.pct_nonzero_source.mean(),mean.equal_donor_pct_nonzero_source))
        for j,row in enumerate(v.itertuples()):
            c.scatter(row.pct_nonzero_source,i+(j-(len(v)-1)/2)*.11,s=12,c='#A3B1B9',zorder=2)
        color=ORANGE if celltype=='Cholangiocytes' else BLUE if celltype=='Hepatocytes' else INK
        c.scatter(mean.equal_donor_pct_nonzero_source,i,s=23,marker='D',c=color,zorder=3)
        c.text(94,i,str(int(mean.n_donors_present))+'/5',va='center',fontsize=8)
    c.set(xlim=(-2,106),ylim=(14.6,-.6),xlabel='Cells with nonzero source value (%)');c.set_xticks([0,20,40,60,80])
    c.set_yticks(range(15),[names[t] for t in types])
    c.set_title('All cell types: donor values and means',pad=15);c.text(94,-1.22,'Donors',fontsize=8,ha='left')
    for ax,ch in zip([a,b,c],'abc'):ax.tick_params(axis='y',length=0,pad=4);letter(ax,ch)
    handles=[Line2D([0],[0],marker='o',linestyle='',markersize=5,color=BLUE,label='Hepatocytes'),
             Line2D([0],[0],marker='s',linestyle='',markersize=5,color=ORANGE,label='Cholangiocytes')]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.035,.995),ncol=2,handletextpad=.4,columnspacing=1.2)
    fig.text(.10,.044,'* Fewer than 20 cholangiocytes; all five donors retained',fontsize=8)
    fig.text(.10,.018,'Source-normalized nonzero values are not verified raw-UMI detection',fontsize=8)
    export(fig,5,[a,b,c])


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--only',choices=['all','1','2','3','4','5'],default='all');args=parser.parse_args()
    for n,fn in enumerate([figure1,figure2,figure3,figure4,figure5],1):
        if args.only in ['all',str(n)]:fn()
    record={'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DATA.iterdir()) if p.is_file()},
            'checks':CHECKS,'selection':args.only,'width_mm':WIDTH_MM,'minimum_font_pt':8,'dpi':600,
            'environment':{'python':platform.python_version(),'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},
            'output_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(MAIN.iterdir()) if p.is_file()}}
    (QA/f'render_manifest_{args.only}.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('All',len(CHECKS),'executed plotting-source checks passed')
