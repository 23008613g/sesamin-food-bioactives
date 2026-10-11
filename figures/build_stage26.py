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
for output_folder in [QA, HERE / 'main', HERE / 'supplement']:
    output_folder.mkdir(parents=True, exist_ok=True)
from alignment import require_matplotlib_panel_alignment
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

def figure1():
    fig,ax=plt.subplots(figsize=(WIDTH,104/25.4))
    fig.subplots_adjust(left=.03,right=.97,bottom=.03,top=.97)
    ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    ax.text(.012,.964,'Sesamin-related comparisons',fontsize=9,weight='bold',color=BLUE)
    ax.text(.532,.964,'Context without sesamin',fontsize=9,weight='bold',color=ORANGE)
    cards=[
      (.005,.902,'Human HepG2 cells','GSE248903 | 9 libraries','Fixed-candidate treatment comparison','Timing and independence unresolved',BLUE),
      (.005,.605,'Mouse hepatocytes','PRJCA021826 / CRA013797 | 9 libraries','Cross-model direction and sensitivity','Treatment mapping provisional',BLUE),
      (.005,.308,'Independent rat liver','GSE42693 | 5 control + 5 sesamin','Dietary-exposure direction boundary','Healthy liver; selected animals',BLUE),
      (.525,.902,'Primary human hepatocytes','GSE247407 | 30 libraries; 5 labels','Fatty-acid expression background','Donor-label relationships unresolved',ORANGE),
      (.525,.605,'Human liver organoids','GSE207889 | 24,034 cells; 4 libraries','AREG cellular localization','One differentiation; technical libraries',ORANGE),
      (.525,.308,'Independent adult human liver','GSE115469 | 8,444 cells; 5 donors','Donor and immune-cell context','Deposited normalized expression',ORANGE)]
    for x,y,head,accession,role,boundary,color in cards:
        ax.add_patch(FancyBboxPatch((x,y-.248),.467,.248,boxstyle='round,pad=.003,rounding_size=.008',
                facecolor='#F7F9FA',edgecolor='#CAD5DB',lw=.7))
        for dy,text,size,col,weight in [(0.024,head,8.8,color,'bold'),(.080,accession,8,INK,'normal'),
                                       (.135,role,8,INK,'normal'),(.190,boundary,8,'#586972','normal')]:
            ax.text(x+.012,y-dy,text,fontsize=size,color=col,weight=weight,va='top')
    export(fig,'Figure_01',[ax])

def figure2():
    h=load('human_candidate_samples.csv').query("symbol == 'VNN1'").copy()
    m=load('mouse_candidate_samples.csv').query("human_symbol == 'VNN1'").copy()
    r=load('rat_candidate_samples.csv').query("human_symbol == 'VNN1' and method == 'deposited_log2' and group in ['control','sesamin']").copy()
    assert len(h)==9 and len(m)==9 and len(r)==10
    h['display_group']=h['sample'].str.split('_').str[0]
    m['display_group']=m.original_label.str.replace(r'\d+$','',regex=True)
    r['display_group']=r['group']
    specs=[(h,'normalized_count',['Control','Model','Treat'],['Control','OA','OA +\nsesamin'],'Human HepG2','Normalized counts','sample'),
           (m,'TPM',['ctrl','HFD','sess'],['ctrl','HFD','sess'],'Mouse hepatocytes','TPM','run'),
           (r,'aggregate_log2_signal',['control','sesamin'],['Control','Sesamin'],'Rat liver','log2(signal + 1)','sample')]
    fig,axes=plt.subplots(1,3,figsize=(WIDTH,80/25.4))
    fig.subplots_adjust(left=.105,right=.985,bottom=.22,top=.855,wspace=.68)
    tidy=[]
    for ax,(d,field,order,ticks,head,unit,idcol),letter in zip(axes,specs,'abc'):
        for j,g in enumerate(order):
            values=d.loc[d.display_group==g].sort_values(idcol)
            color=[GREY,ORANGE,BLUE][j] if len(order)==3 else [GREY,BLUE][j]
            offsets=np.linspace(-.12,.12,len(values))
            ax.scatter(j+offsets,values[field],s=26,c=color,edgecolors='white',linewidths=.35,zorder=3)
            ax.hlines(values[field].mean(),j-.24,j+.24,color=INK,lw=1.15,zorder=4)
            for row in values.to_dict('records'):
                tidy.append({'panel':letter,'system':head,'observation_id':row[idcol],'group':g,'value':row[field],'unit':unit})
        ax.set_xlim(-.5,len(order)-.5);ax.set_xticks(range(len(order)),ticks)
        ax.set_title(head,pad=13);ax.set_ylabel(unit,labelpad=5);label(ax,letter)
        ax.tick_params(axis='x',length=0,pad=5)
    axes[0].set_ylim(0,120);axes[0].set_yticks([0,40,80,120])
    axes[1].set_ylim(0,14);axes[1].set_yticks([0,4,8,12])
    pd.DataFrame(tidy).to_csv(DATA/'Figure_02_VNN1_all_observations.csv',index=False)
    export(fig,'Figure_02',list(axes))

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
    for ax,letter,head in zip(axes,'ab',['Author cell subtypes','AREG expression']):
        ax.set(xlim=(xlo,xhi),ylim=(ylo,yhi),xlabel='UMAP 1',ylabel='UMAP 2')
        ax.set_aspect('equal',adjustable='box');ax.set_xticks([]);ax.set_yticks([])
        ax.set_title(head,pad=13);label(ax,letter)
    cax=fig.add_axes([.925,.375,.018,.35])
    cb=fig.colorbar(sc,cax=cax);cb.set_label('log1p(CP10k)',labelpad=5,fontsize=8)
    handles=[Line2D([0],[0],marker='o',linestyle='',markersize=4,color=color,label=cat) for cat,color in zip(cats,colors)]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.04,.225),ncol=3,columnspacing=1.5,handletextpad=.4,labelspacing=.7)
    counts=d.groupby('celltype_scTypedb',sort=False).agg(n_cells=('cell_key','size'),AREG_detected=('AREG__counts',lambda x:int((x>0).sum())))
    counts.to_csv(DATA/'Figure_03_subtype_counts.csv')
    export(fig,'Figure_03',list(axes))

def figure4():
    local=load('library_category_candidates_all25.csv').query("symbol == 'AREG' and celltype_category in ['Hepatocytes','Cholangiocytes']")
    matched=load('matched_scale_library_all25.csv').query("symbol == 'AREG'")
    summary=load('matched_scale_summary_all25.csv').query("symbol == 'AREG' and weighting == 'equal_library'")
    libs=['GSM6322577','GSM6322578','GSM6322579','GSM6322580']
    assert len(local)==8 and len(matched)==8 and set(matched.population)=={'Whole organoid','Hepatocytes'}
    assert matched.groupby('population').size().eq(4).all()
    fig,axes=plt.subplots(1,3,figsize=(WIDTH,80/25.4))
    fig.subplots_adjust(left=.12,right=.98,bottom=.255,top=.855,wspace=.66)
    a,b,c=axes
    for i,gsm in enumerate(libs):
        v=local[local.gsm==gsm].set_index('celltype_category')
        h=v.loc['Hepatocytes','detected_pct'];ch=v.loc['Cholangiocytes','detected_pct']
        a.plot([h,ch],[i,i],color='#B7C2C8',lw=1)
        a.scatter(h,i,c=BLUE,s=24,zorder=3);a.scatter(ch,i,c=ORANGE,s=24,zorder=3)
    a.set(yticks=range(4),yticklabels=['Control 1','Control 2','PA 1','PA 2'],ylim=(3.55,-.55),xlim=(0,100),xlabel='Detected cells (%)')
    a.set_title('Raw-UMI detection',pad=13);a.set_xticks([0,50,100]);a.tick_params(axis='y',length=0)
    maximum=float(matched.mean_log1p_cp10k.max())*1.17
    for ax,pop,head in [(b,'Whole organoid','Whole organoid'),(c,'Hepatocytes','Hepatocytes')]:
        d=matched[matched.population==pop]
        groups=[('CTRL','Control',GREY),('PA','PA',ORANGE)]
        for j,(group,display,color) in enumerate(groups):
            values=d[d.condition==group].sort_values('gsm').mean_log1p_cp10k.to_numpy()
            assert len(values)==2
            ax.scatter(j+np.array([-.075,.075]),values,s=28,color=color,edgecolors='white',linewidths=.3,zorder=3)
            ax.hlines(values.mean(),j-.21,j+.21,color=INK,lw=1.1)
        ax.set(xticks=[0,1],xticklabels=['Control','PA'],xlim=(-.5,1.5),ylim=(0,maximum))
        ax.set_title(head,pad=13);ax.tick_params(axis='x',length=0)
        calc=d.groupby('condition').mean_log1p_cp10k.mean()
        expected=summary[summary.population==pop].iloc[0]
        assert np.isclose(calc.PA-calc.CTRL,expected.PA_minus_control,atol=1e-12)
    b.set_ylabel('Mean log1p(CP10k)',labelpad=5)
    for ax,letter in zip(axes,'abc'):label(ax,letter)
    handles=[Line2D([0],[0],marker='o',linestyle='',color=col,markersize=4,label=lab)
             for col,lab in [(BLUE,'Hepatocytes'),(ORANGE,'Cholangiocytes')]]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.035,.09),ncol=2,columnspacing=1.5,handletextpad=.4)
    local.to_csv(DATA/'Figure_04_AREG_detection.csv',index=False)
    matched.to_csv(DATA/'Figure_04_AREG_matched_scale.csv',index=False)
    export(fig,'Figure_04',list(axes))

def figure5():
    paired=load('figure_AREG_paired_donors.csv')
    context=load('figure_AREG_all_cell_context.csv')
    details=load('donor_broad_results_all25.csv')
    details=details[(details.symbol=='AREG')&details.available&(details.n_cells>0)]
    assert len(paired)==5 and paired.donor.nunique()==5 and len(context)==15
    assert (paired.delta_nonzero_percentage_points>0).sum()==4
    assert (paired.delta_mean_source_log2CPM>0).sum()==4
    for row in context.itertuples():
        vals=details[details.cell_type==row.cell_type]
        assert len(vals)==row.n_donors_with_available_expression
        assert np.isclose(vals.pct_nonzero_source.mean(),row.equal_donor_nonzero_pct,atol=1e-8)
        assert np.isclose(vals.mean_source_log2CPM.mean(),row.equal_donor_mean_source_log2CPM,atol=1e-10)
    fig=plt.figure(figsize=(WIDTH,174/25.4))
    grid=fig.add_gridspec(2,2,left=.365,right=.97,top=.890,bottom=.09,height_ratios=[1,2.3],hspace=.58,wspace=.49)
    a=fig.add_subplot(grid[0,0]);b=fig.add_subplot(grid[0,1]);c=fig.add_subplot(grid[1,:])
    handles=[Line2D([0],[0],marker='o',linestyle='',color=col,markersize=4,label=lab)
             for col,lab in [(BLUE,'Hepatocytes'),(ORANGE,'Cholangiocytes')]]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.035,.995),ncol=2,columnspacing=1.4,handletextpad=.4)
    for i,row in enumerate(paired.itertuples()):
        for ax,h,ch in [(a,row.hep_nonzero_pct,row.chol_nonzero_pct),(b,row.hep_mean_source_log2CPM,row.chol_mean_source_log2CPM)]:
            ax.plot([h,ch],[i,i],color='#B7C2C8',lw=1)
            ax.scatter(h,i,color=BLUE,s=23,zorder=3);ax.scatter(ch,i,color=ORANGE,s=23,zorder=3)
    donors=[r.donor.replace('TLH','') for r in paired.itertuples()]
    alabs=[f"{r.donor.replace('TLH','')}{' *' if not r.both_types_ge20 else ''} ({r.n_cholangiocytes}/{r.n_hepatocytes})" for r in paired.itertuples()]
    a.set_yticks(range(5),alabs);b.set_yticks(range(5),donors)
    a.set_ylabel('Donor (Chol / Hep cells)',labelpad=7)
    for ax,head,xlab in [(a,'Nonzero entries','Cells (%)'),(b,'Mean expression','Source log2CPM')]:
        ax.set_title(head,pad=13);ax.set_xlabel(xlab,labelpad=6);ax.set_ylim(4.65,-.65);ax.tick_params(axis='y',length=0)
    a.set_xlim(-.8,17);a.set_xticks([0,5,10,15]);b.set_xlim(-.008,.17);b.set_xticks([0,.05,.1,.15])
    names={'Central_venous_LSECs':'Central venous LSECs','Cholangiocytes':'Cholangiocytes','Erythroid_Cells':'Erythroid cells',
        'Hepatic_Stellate_Cells':'Hepatic stellate cells','Hepatocytes':'Hepatocytes','Inflammatory_Macrophage':'Inflammatory macrophages',
        'Mature_B_Cells':'Mature B cells','NK-like_Cells':'NK-like cells','Non-inflammatory_Macrophage':'Non-inflammatory macrophages',
        'Periportal_LSECs':'Periportal LSECs','Plasma_Cells':'Plasma cells','Portal_endothelial_Cells':'Portal endothelial cells',
        'alpha-beta_T_Cells':'Alpha-beta T cells','gamma-delta_T_Cells_1':'Gamma-delta T cells 1','gamma-delta_T_Cells_2':'Gamma-delta T cells 2'}
    for i,row in enumerate(context.itertuples()):
        vals=details[details.cell_type==row.cell_type].sort_values('donor')
        c.scatter(vals.pct_nonzero_source,i+np.linspace(-.17,.17,len(vals)),s=13,c='#BCC5CA',linewidths=0,zorder=2)
        col=BLUE if row.cell_type=='Hepatocytes' else ORANGE if row.cell_type=='Cholangiocytes' else INK
        c.scatter(row.equal_donor_nonzero_pct,i,marker='D',s=19,color=col,linewidths=0,zorder=3)
        c.text(94,i,f'{row.n_donors_with_cells}/5',fontsize=8,va='center')
    c.set_yticks(range(15),[names[t] for t in context.cell_type]);c.tick_params(axis='y',length=0,pad=5)
    c.set_ylim(14.6,-.6);c.set_xlim(-2,106);c.set_xticks([0,20,40,60,80]);c.set_xlabel('Cells with nonzero source value (%)',labelpad=6)
    c.set_title('All cell types: donor values and means',pad=15);c.text(94,-1.22,'Donors',fontsize=8,ha='left')
    for ax,letter in zip([a,b,c],'abc'):label(ax,letter)
    export(fig,'Figure_05',[a,b,c])

def figureS16():
    summary=load('configuration_effect_summary.csv');candidates=load('frozen_primary_candidates.csv')
    genes=candidates.human_symbol.tolist();cons=load('cross_species_configuration_consensus.csv')
    cons=cons[cons.aggregation=='mean'].set_index('human_symbol').loc[genes]
    ranges=summary[(summary.species=='mouse')&(summary.metric=='TPM')&(summary.aggregation=='mean')&(summary.contrast=='treatment_vs_model')].set_index('human_symbol').loc[genes]
    counts=[];zero=[]
    for species,metric in [('human','median_ratio_reestimated'),('mouse','TPM')]:
        for contrast in ['model_vs_control','treatment_vs_model']:
            x=summary[(summary.species==species)&(summary.metric==metric)&(summary.aggregation=='mean')&(summary.contrast==contrast)].set_index('human_symbol').loc[genes]
            assert x.configurations.eq(64).all()
            counts.append(x.baseline_direction_retained_n.to_numpy(int));zero.append(x.unchanged_n.eq(64).to_numpy())
    counts=np.asarray(counts).T;zero=np.asarray(zero).T
    zero_cells={(genes[i],j) for i,j in zip(*np.where(zero))}
    assert zero_cells=={('NEGR1',2),('NEGR1',3),('BCAS1',2),('BCAS1',3),('NGFR',2)}
    fig,axes=plt.subplots(1,3,figsize=(WIDTH,137/25.4),gridspec_kw={'width_ratios':[1.05,1.52,.73]})
    fig.subplots_adjust(left=.15,right=.98,bottom=.17,top=.91,wspace=.33)
    a,b,c=axes;a.axvline(0,color='#B7C2C8',lw=.6,zorder=0)
    low=dict(zip(candidates.human_symbol,candidates.all_libraries_TPM_lt1))
    for i,(gene,row) in enumerate(ranges.iterrows()):
        col=BLUE if gene=='VNN1' else ORANGE if gene in ['AREG','JAG2'] else GREY
        a.hlines(i,row.min_log2_ratio,row.max_log2_ratio,color=col,lw=1.3)
        a.scatter(row.full_data_log2_ratio,i,s=15,color=col,zorder=3)
        n=cons.loc[gene,'both_model_and_treatment_concordant_n'];c.barh(i,n/4096*100,height=.5,color=col)
    a.set_yticks(range(15),[g+(' *' if low[g] else '') for g in genes])
    a.set_xlim(ranges.min_log2_ratio.min()-.15,ranges.max_log2_ratio.max()+.15)
    a.set_title('Mouse effect range',pad=13);a.set_xlabel('sess / HFD log2 ratio\n0.1 TPM pseudocount',labelpad=7)
    b.imshow(counts,vmin=0,vmax=64,cmap='Blues',aspect='auto',interpolation='nearest',extent=(-.5,3.5,14.5,-.5))
    for i in range(15):
        for j in range(4):
            b.text(j,i,f'{counts[i,j]} (0)' if zero[i,j] else str(counts[i,j]),fontsize=8,ha='center',va='center',color='white' if counts[i,j]>=40 else INK)
    b.set_xticks(range(4),['Human\nModel','Human\nTreat','Mouse\nModel','Mouse\nTreat']);b.tick_params(axis='x',length=0,pad=5)
    b.set_title('Baseline sign retained',pad=13);b.set_xlabel('Configurations out of 64',labelpad=7)
    c.set_title('Both agree',pad=13);c.set_xlim(0,100);c.set_xticks([0,50,100]);c.set_xlabel('Pairs (%)\nOut of 4096',labelpad=7)
    for ax,letter in zip(axes,'abc'):
        ax.set_ylim(14.6,-.6);ax.tick_params(axis='y',length=0);label(ax,letter)
    for ax in [b,c]:ax.set_yticks([])
    for spine in b.spines.values():spine.set_visible(False)
    record=[]
    for i,gene in enumerate(genes):
        for j,name in enumerate(['human_model','human_treat','mouse_model','mouse_treat']):
            record.append({'gene':gene,'contrast':name,'baseline_retained':int(counts[i,j]),'all_zero_direction':bool(zero[i,j]),'display':f'{counts[i,j]} (0)' if zero[i,j] else str(counts[i,j])})
    pd.DataFrame(record).to_csv(DATA/'Figure_S16_sign_counts_zero_explicit.csv',index=False)
    export(fig,'Figure_S16',list(axes),supplement=True,exemptions=[{'panels':list('abc'),'checks':['panel-width'],
        'reason':'Distinct range, four-column count matrix and agreement bar panels require intentional unequal widths; row edges and heights remain identical.'}])

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--only',default='all',choices=['all','1','2','3','4','5','S16']);args=parser.parse_args()
    tasks={'1':figure1,'2':figure2,'3':figure3,'4':figure4,'5':figure5,'S16':figureS16}
    for name,fn in tasks.items():
        if args.only in ['all',name]:fn()
    manifest={'python':platform.python_version(),'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,
              'matplotlib':matplotlib.__version__,'width_mm':width_mm,'export_dpi':600,'minimum_intended_text_pt':8,
              'last_run_selection':args.only,'randomness':'None; deterministic spread offsets only; no sampled observations',
              'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.iterdir() if p.is_file()},
              'figure_outputs':{str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for d in ['main','supplement'] for p in (HERE/d).glob('*') if p.is_file()}}
    (QA/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
