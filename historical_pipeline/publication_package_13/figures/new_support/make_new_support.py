"""Prepare full-candidate descriptive heatmaps from existing result tables."""
from pathlib import Path
import hashlib
import json
import math
import platform
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NEXT = HERE.parents[2]
PROJECT = NEXT.parent
DATA = HERE/'source_data'
DATA.mkdir(exist_ok=True)
SKILL = Path.home()/'.agents/skills/nature-figure/scripts'
sys.path.insert(0,str(SKILL))
from audit_panel_alignment import require_matplotlib_panel_alignment

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],
    'font.size':7.5,'axes.labelsize':7.5,'axes.titlesize':9,
    'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7,
    'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.5,
    'savefig.facecolor':'white'})
WIDTH_IN = 180/25.4
BLUE, ORANGE, GRAY = '#3E6D9C','#C47738','#6C7580'
DONORS = ['F95','M95','H20','H961','F95P']
RUNS = [f'DRR{x}' for x in range(1010532,1010541)]
INPUTS = []

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def register(path):
    path=Path(path)
    INPUTS.append({'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size})
    return path

def export(fig, ax, cax, stem):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,axes=[ax],exclude_axes=[cax],
        json_out=HERE/f'{stem}.alignment.json',overlay_svg=HERE/f'{stem}.alignment.svg',
        tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
    fig.savefig(HERE/f'{stem}.pdf')
    fig.savefig(HERE/f'{stem}.svg')
    fig.savefig(HERE/f'{stem}.png',dpi=600)
    fig.savefig(HERE/f'{stem}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)

def heatmap(fig, matrix, labels, columns, rect, cmap, norm, label):
    ax=fig.add_axes(rect)
    image=ax.imshow(np.ma.masked_invalid(matrix),interpolation='nearest',aspect='auto',cmap=cmap,norm=norm)
    ax.set_xticks(np.arange(len(columns)),columns)
    ax.set_yticks(np.arange(len(labels)),labels)
    ax.tick_params(axis='both',length=0,pad=5)
    ax.set_xticks(np.arange(-.5,len(columns),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(labels),1),minor=True)
    ax.grid(which='minor',color='white',linewidth=.55)
    ax.tick_params(which='minor',bottom=False,left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value=matrix[i,j]
            if not np.isfinite(value):
                text,color='NA','#30363C'
            else:
                rgb=cmap(norm(value))[:3]
                luminance=.2126*rgb[0]+.7152*rgb[1]+.0722*rgb[2]
                color='white' if luminance<.46 else '#17212A'
                text='0' if value==0 else f'{value:.2f}'
            ax.text(j,i,text,ha='center',va='center',fontsize=7,color=color)
    cax=fig.add_axes([.902,rect[1]+rect[3]*.12,.022,rect[3]*.76])
    cbar=fig.colorbar(image,cax=cax)
    cbar.set_label(label,labelpad=6,fontsize=7.5)
    cbar.ax.tick_params(labelsize=7,length=2)
    cbar.outline.set_linewidth(.5)
    return ax,cax,cbar

def fig_s3(values, labels):
    fig=plt.figure(figsize=(WIDTH_IN,185/25.4))
    cmap=LinearSegmentedColormap.from_list('candidate_direction',[BLUE,'#FCFCFC',ORANGE])
    cmap.set_bad('#D8DDE2')
    assert np.nanmin(values)>=-2 and np.nanmax(values)<=2
    ax,cax,bar=heatmap(fig,values,labels,DONORS,[.255,.115,.615,.765],cmap,
        TwoSlopeNorm(vmin=-2,vcenter=0,vmax=2),'log2[(FFA mean + 1) / (Control mean + 1)]')
    bar.set_ticks([-2,-1,0,1,2])
    ax.set_xlabel('Donor label',labelpad=9)
    fig.text(.255,.954,'All 25 candidates: PHH fatty-acid context',fontsize=10,fontweight='bold')
    fig.text(.255,.925,'GSE247407 | No sesamin intervention',fontsize=8,color=GRAY)
    fig.text(.03,.035,'Frozen candidate order; no clustering. 30 libraries, 5 donor labels.\n'
             'F95/F95P relationship and donor-label independence remain unresolved.',fontsize=7.2,color=GRAY,linespacing=1.5)
    export(fig,ax,cax,'figS3_phh_all25')

def fig_s4(values, labels, low):
    fig=plt.figure(figsize=(WIDTH_IN,155/25.4))
    cmap=LinearSegmentedColormap.from_list('candidate_abundance',['#F4F7FA','#91B0CE','#214D78'])
    cmap.set_bad('#D8DDE2')
    assert np.nanmin(values)>=0 and np.nanmax(values)>0
    norm=Normalize(vmin=0,vmax=float(np.nanmax(values)),clip=False)
    rows=[label+(' *' if flag else '') for label,flag in zip(labels,low)]
    columns=['ctrl1','ctrl2','ctrl3','HFD1','HFD2','HFD3','sess1','sess2','sess3']
    ax,cax,bar=heatmap(fig,values,rows,columns,[.165,.23,.705,.625],cmap,norm,'log2(TPM + 1)')
    bar.set_ticks([0,1,2,3,4,5,6])
    for x in [2.5,5.5]:
        ax.axvline(x,color='white',linewidth=2)
    ax.set_xlabel('Submitted mouse library (original label)',labelpad=9)
    fig.text(.165,.954,'All one-to-one candidates: mouse abundance',fontsize=10,fontweight='bold')
    fig.text(.165,.916,'15 primary candidates × 9 libraries | Inferred treatment mapping',fontsize=8,color=GRAY)
    fig.text(.03,.06,'* All nine libraries <1 TPM; asterisks do not denote statistical significance.\n'
             'Color shows relative RNA abundance, not treatment effect; no row scaling.\n'
             'The nine libraries are not assumed to represent nine independent animals.',fontsize=7.2,color=GRAY,linespacing=1.5)
    export(fig,ax,cax,'figS4_mouse_all15')

def main():
    cpath=register(NEXT/'sesamin_exploration_04/priority_all.csv')
    candidates=pd.read_csv(cpath,keep_default_na=False)
    assert len(candidates)==25 and candidates.gene_id.is_unique
    assert sha(cpath)=='d78d674a48448532ecec1e49d7364d447bdc298a202090b93d8a73718517ea78'
    facts=json.loads(register(NEXT/'manuscript_draft_11/data/facts.json').read_text(encoding='utf-8'))
    hpath=register(NEXT/'followup_05/phh_donor_effects.csv')
    human=pd.read_csv(hpath)
    assert len(human)==125
    assert not human.duplicated(['gene_id','donor_label']).any()
    assert set(human.gene_id)==set(candidates.gene_id) and set(human.donor_label)==set(DONORS)
    values=human.pivot(index='gene_id',columns='donor_label',values='log2_FFA_Control_mean_norm_plus1')
    values=values.reindex(index=candidates.gene_id,columns=DONORS)
    # Arithmetic self-check uses stored group means only; normalization is not rerun.
    assert (human.FFA_mean_norm.dropna()>=0).all() and (human.Control_mean_norm.dropna()>=0).all()
    human_check=np.log2((human.FFA_mean_norm.to_numpy()+1)/(human.Control_mean_norm.to_numpy()+1))
    assert np.allclose(human_check,human.log2_FFA_Control_mean_norm_plus1,equal_nan=True,atol=2e-14,rtol=0)
    human_order={gene:i for i,gene in enumerate(candidates.gene_id)}
    donor_order={label:i for i,label in enumerate(DONORS)}
    human['_candidate_order']=human.gene_id.map(human_order)
    human['_donor_order']=human.donor_label.map(donor_order)
    human=human.sort_values(['_candidate_order','_donor_order'])
    human.to_csv(DATA/'figS3_phh_125_values.csv',index=False)
    source_matrix=values.copy()
    source_matrix.insert(0,'symbol',candidates.set_index('gene_id').symbol)
    source_matrix.to_csv(DATA/'figS3_heatmap_matrix.csv')
    mpath=register(Path('D:/ResearchData/sesamin_mouse_2026-10-08/matrices/candidate_mouse_descriptive.tsv'))
    mouse_all=pd.read_csv(mpath,sep='\t')
    assert sha(mpath)==next(x['sha256'] for x in facts['input_sources'] if x['key']=='mouse_candidates')
    mouse=mouse_all.loc[mouse_all.primary_comparison_eligible.eq(True)].copy()
    assert len(mouse)==15 and mouse.original_human_id.is_unique
    assert len(mouse)==facts['counts']['primary_one_to_one_candidates']
    mouse['_candidate_order']=mouse.original_human_id.map(human_order)
    assert mouse['_candidate_order'].notna().all()
    mouse=mouse.sort_values('_candidate_order').reset_index(drop=True)
    tpm=mouse[[run+'_TPM' for run in RUNS]].to_numpy(dtype=float)
    assert tpm.shape==(15,9) and np.isfinite(tpm).all() and (tpm>=0).all()
    transformed=np.log2(tpm+1)
    scalar=np.array([[math.log2(float(v)+1) for v in row] for row in tpm])
    assert np.allclose(transformed,scalar,atol=2e-15,rtol=0)
    # Inverse exp/log roundtrip tolerance scales with representable spacing,
    # rather than one absolute cutoff across 0-97 TPM.
    assert np.all(abs(np.exp2(transformed)-1-tpm)<=8*np.spacing(np.maximum(tpm,1)))
    assert np.all(transformed[tpm==0]==0)
    low=(tpm<1).all(axis=1)
    assert low.sum()==5 and np.array_equal(low,mouse.all_libraries_TPM_lt1.eq(True).to_numpy())
    assert set(mouse.loc[low,'human_symbol'])=={'NEGR1','BCAS1','NGFR','JAG2','RAB3B'}
    mouse.to_csv(DATA/'figS4_primary_mouse_rows.csv',index=False)
    long=[]
    for i,row in mouse.iterrows():
        for j,run in enumerate(RUNS):
            long.append({'original_human_id':row.original_human_id,'human_symbol':row.human_symbol,
                         'mouse_gene_id':row.mouse_gene_id,'mouse_symbol':row.mouse_symbol,'run':run,
                         'original_label':['ctrl','HFD','sess'][j//3]+str(j%3+1),
                         'TPM':float(tpm[i,j]),'log2_TPM_plus1':float(transformed[i,j]),'all9_TPM_lt1':bool(low[i])})
    pd.DataFrame(long).to_csv(DATA/'figS4_mouse_135_values.csv',index=False)
    fig_s3(values.to_numpy(dtype=float),candidates.symbol.tolist())
    fig_s4(transformed,mouse.human_symbol.tolist(),low)
    record={'date':'2026-10-09','status':'figures_rendered_from_existing_results',
            'python':platform.python_version(),'interpreter':sys.executable,
            'packages':{'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},
            'input_fingerprints':INPUTS,'script_sha256':sha(Path(__file__)),
            'PHH':{'cells':int(values.size),'missing':int(values.isna().sum().sum()),'feature_covered':int(human.feature_covered.sum()),
                   'range':[float(np.nanmin(values)),float(np.nanmax(values))],'limits':[-2,2],
                   'arithmetic_max_abs_error':float(np.nanmax(abs(human_check-pd.read_csv(hpath).log2_FFA_Control_mean_norm_plus1))),
                   'row_order':candidates.symbol.tolist(),'column_order':DONORS},
            'mouse':{'all_source_rows':len(mouse_all),'primary_rows':len(mouse),'cells':int(tpm.size),
                     'missing':int(np.isnan(tpm).sum()),'zero_TPM':int((tpm==0).sum()),
                     'TPM_range':[float(tpm.min()),float(tpm.max())],
                     'display_log2_range':[float(transformed.min()),float(transformed.max())],
                     'low_expression_symbols':mouse.loc[low,'human_symbol'].tolist(),
                     'row_order':mouse.human_symbol.tolist(),'run_order':RUNS,
                     'scalar_log2_max_abs_error':float(abs(transformed-scalar).max()),
                     'inverse_transform_max_abs_error':float(abs(np.exp2(transformed)-1-tpm).max()),
                     'excluded_scope':'Two secondary ACAT2 mapping rows and nine noncomparable candidate rows remain in the complete appendix, outside this prespecified primary denominator.'},
            'analysis_not_performed':['new quantification','new normalization','clustering','row standardization','differential expression','p-values'],
            'development_note':'Initial inverse-log self-check used an overstrict fixed 3e-14 absolute tolerance and stopped before any figure rendering. The final check allows eight floating-point spacings at each original value; scalar log2 and original-zero checks also pass. Input data and transformed values were not changed.',
            'source_data_fingerprints':[{'file':str(x.relative_to(HERE)),'sha256':sha(x)} for x in sorted(DATA.glob('*.csv'))]}
    (HERE/'execution_record.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'figures':2,'PHH_cells':125,'mouse_cells':135,'missing':[int(values.isna().sum().sum()),int(np.isnan(tpm).sum())],'low_mouse':int(low.sum())}))

if __name__=='__main__':
    main()
