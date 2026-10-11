"""Frozen-data stage 15 plots: external contrast, omission stability and context."""
from pathlib import Path
import hashlib,json,platform,shutil,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
from matplotlib.lines import Line2D
import PIL

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
SKILLS=Path.home()/'.agents/skills/nature-figure/scripts'
sys.path.insert(0,str(SKILLS))
from audit_panel_alignment import require_matplotlib_panel_alignment
width_mm=180
WIDTH=width_mm/25.4
INK='#263443';BLUE='#3C7F99';ORANGE='#BB7B43';GREY='#8C959E';PALE='#CCD6DC'
INPUTS={
 'rat_results':BASE/'rat_validation/outputs/candidate_results_all25_two_normalizations.csv',
 'rat_samples':BASE/'rat_validation/outputs/candidate_sample_values.csv',
 'rat_exact':BASE/'rat_validation/outputs/posthoc_exact_sensitivity_all25.csv',
 'rat_pca':BASE/'rat_validation/outputs/all20_array_pca.csv',
 'rat_variance':BASE/'rat_validation/outputs/pca_variance.csv',
 'robust_ranges':BASE/'robustness/configuration_effect_summary.csv',
 'robust_consensus':BASE/'robustness/cross_species_configuration_consensus.csv',
 'robust_candidates':BASE/'robustness/frozen_primary_candidates.csv',
 'human_pca':BASE/'robustness/human_pca_scores.csv',
 'human_variance':BASE/'robustness/human_pca_variance.csv',
 'mouse_pca':BASE/'robustness/mouse_pca_scores.csv',
 'mouse_variance':BASE/'robustness/mouse_pca_variance.csv',
 'cell_localization':BASE/'single_cell/cholangiocyte_vs_hepatocyte_per_library_all25.csv',
 'cell_decomposition':BASE/'single_cell/kitagawa_summary_all25.csv',
}
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':7.5,'axes.labelsize':7.5,'axes.titlesize':8.5,'xtick.labelsize':7,'ytick.labelsize':7,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'text.color':INK,
 'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.edgecolor':'#73818C',
 'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'figure.facecolor':'white','savefig.facecolor':'white'})

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,x):(HERE/name).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def header(fig,title,sub):
    fig.text(.045,.98,title,fontsize=10,weight='bold',va='top')
    fig.text(.045,.944,sub,fontsize=7.5,va='top')
def letter(ax,x):
    ax.text(0,1,x,transform=ax.transAxes+ScaledTranslation(-19/72,11/72,ax.figure.dpi_scale_trans),
            fontsize=10,weight='bold',ha='left',va='bottom')
def note(fig,text,y=.025):fig.text(.045,y,text,fontsize=7,va='bottom',linespacing=1.55)
def export(fig,stem,axes):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig,axes=axes,json_out=HERE/(stem+'.alignment.json'),
       overlay_svg=HERE/(stem+'.alignment.svg'),tolerance_pt=1.5,gutter_tolerance_pt=1.5,
       require_panel_labels=len(axes)>1,strict=True)
    fig.savefig(HERE/(stem+'.pdf'),dpi=600)
    fig.savefig(HERE/(stem+'.svg'),dpi=600)
    fig.savefig(HERE/(stem+'.png'),dpi=600)
    fig.savefig(HERE/(stem+'.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)
def qlabel(q):return 'NA' if pd.isna(q) else f'{q:.3g}'

def rat_figure(tables):
    results=tables['rat_results'];r=results[results.method=='deposited_log2'].copy()
    assert len(r)==25 and r.human_symbol.nunique()==25
    valid=r.status.eq('mapped')
    assert int(valid.sum())==15 and r.loc[~valid,'effect_log2_signal'].isna().all()
    assert np.isfinite(r.loc[valid,['effect_log2_signal','ci95_low','ci95_high','q_bh_fixed25']]).all().all()
    assert (r.loc[valid,'ci95_low']<=r.loc[valid,'effect_log2_signal']).all()
    assert (r.loc[valid,'ci95_high']>=r.loc[valid,'effect_log2_signal']).all()
    fig=plt.figure(figsize=(WIDTH,195/25.4))
    gs=fig.add_gridspec(2,2,left=.235,right=.96,bottom=.185,top=.855,wspace=.77,hspace=.85,width_ratios=[1.6,1])
    a=fig.add_subplot(gs[:,0]);b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,1])
    header(fig,'Independent rat dietary study: complete candidate results','Primary contrast: sesamin minus control | Five selected rat livers per group')
    a.axvline(0,color='#BEC7CE',lw=.65,zorder=0)
    lows=r.loc[valid,'ci95_low'];highs=r.loc[valid,'ci95_high']
    limit=max(abs(float(lows.min())),abs(float(highs.max())))*1.15
    a.set_xlim(-limit,limit)
    for i,row in enumerate(r.itertuples()):
        if row.status=='mapped':
            color=ORANGE if row.effect_log2_signal>0 else BLUE
            a.errorbar(row.effect_log2_signal,i,xerr=[[row.effect_log2_signal-row.ci95_low],[row.ci95_high-row.effect_log2_signal]],
                fmt='o',markersize=3.1,color=color,elinewidth=.8,capsize=2,capthick=.7)
        else:a.text(.85,i,'NA',transform=a.get_yaxis_transform(),ha='center',va='center',fontsize=7,color=GREY)
        a.text(1.07,i,qlabel(row.q_bh_fixed25),transform=a.get_yaxis_transform(),ha='left',va='center',fontsize=7)
    a.set_yticks(range(25),r.human_symbol.tolist());a.set_ylim(24.65,-.65)
    a.tick_params(axis='y',length=0,pad=4)
    a.set_title('All 25 fixed candidates',pad=14)
    a.text(1.07,1.018,'Welch\nBH q',transform=a.transAxes,ha='left',va='bottom',fontsize=7)
    a.set_xlabel('Mean log2 signal difference\nNominal 95% Welch CI',labelpad=6)
    letter(a,'a')
    v=tables['rat_samples']
    v=v[(v.method=='deposited_log2')&(v.human_symbol=='VNN1')&(v.group.isin(['control','sesamin']))]
    assert len(v)==10 and v['sample'].nunique()==10 and v.groupby('group').size().eq(5).all()
    offsets=np.array([-.13,-.065,0,.065,.13])
    for i,group,color in [(0,'control',GREY),(1,'sesamin',ORANGE)]:
        values=v[v.group==group].sort_values('sample').aggregate_log2_signal.to_numpy()
        b.scatter(i+offsets,values,s=18,color=color,linewidths=0)
        b.hlines(float(values.mean()),i-.23,i+.23,color=INK,lw=1)
    b.set_xlim(-.5,1.5);b.set_xticks([0,1],['Control','Sesamin']);b.set_ylabel('log2(signal + 1)',labelpad=5)
    b.set_title('VNN1: individual rats',pad=13);letter(b,'b')
    vals=results[results.human_symbol=='VNN1'].set_index('method')
    for i,method,color in [(0,'deposited_log2',ORANGE),(1,'quantile_sensitivity',BLUE)]:
        row=vals.loc[method]
        c.errorbar(row.effect_log2_signal,i,xerr=[[row.effect_log2_signal-row.ci95_low],[row.ci95_high-row.effect_log2_signal]],
                    fmt='o',markersize=4,color=color,elinewidth=.9,capsize=2)
    c.set_yticks([0,1],['Deposited','Quantile']);c.set_ylim(1.65,-.65)
    c.set_xlim(0,2.65);c.axvline(0,color='#BBC5CD',lw=.6)
    c.set_xlabel('Mean log2 signal difference\nNominal 95% Welch CI',labelpad=6)
    c.set_title('Normalization sensitivity',fontsize=8,pad=13);letter(c,'c')
    exact=tables['rat_exact'];ex=exact[(exact.human_symbol=='VNN1')&(exact.method=='deposited_log2')].iloc[0]
    assert np.isclose(ex.permutation_q_bh25,0.198412698413,atol=1e-10)
    note(fig,'VNN1 increases in the rat study, opposite to the original human/mouse treatment directions.\n'
        'VNN1 primary Welch BH q = 0.00257; post-hoc exact-label BH q = 0.1984 (25-candidate family).\n'
        'Quantile sensitivity reuses the same animals. Nominal intervals do not adjust for animal selection.\n'
        'NA = unavailable rat mapping; all 25 candidates are retained.')
    export(fig,'fig15_rat_external',[a,b,c])
    r.to_csv(HERE/'fig15_rat_primary_source.csv',index=False,na_rep='NA')

def robustness_figure(tables):
    f=tables['robust_ranges']
    f=f[(f.species=='mouse')&(f.metric=='TPM')&(f.aggregation=='mean')&(f.contrast=='treatment_vs_model')].copy()
    candidates=tables['robust_candidates'];order=candidates.human_symbol.tolist()
    f=f.set_index('human_symbol').loc[order].reset_index()
    consensus=tables['robust_consensus'];consensus=consensus[consensus.aggregation=='mean'].set_index('human_symbol').loc[order]
    assert len(f)==len(consensus)==15 and (f.configurations==64).all() and (consensus.configuration_pairs==4096).all()
    assert np.isfinite(f[['min_log2_ratio','max_log2_ratio','full_data_log2_ratio']]).all().all()
    assert (f.min_log2_ratio<=f.full_data_log2_ratio+1e-12).all() and (f.max_log2_ratio>=f.full_data_log2_ratio-1e-12).all()
    fig,axes=plt.subplots(1,2,figsize=(WIDTH,161/25.4))
    fig.subplots_adjust(left=.15,right=.91,bottom=.24,top=.84,wspace=.55)
    header(fig,'Candidate sensitivity to the observed library composition','All 15 original primary orthologues | Arithmetic-mean primary analyses only')
    a,b=axes; a.axvline(0,color='#B7C2C9',lw=.7,zorder=0)
    low=dict(zip(candidates.human_symbol,candidates.all_libraries_TPM_lt1))
    for i,row in enumerate(f.itertuples()):
        color=BLUE if row.human_symbol=='VNN1' else ORANGE if row.human_symbol in ['AREG','JAG2'] else GREY
        a.hlines(i,row.min_log2_ratio,row.max_log2_ratio,color=color,lw=1.5)
        a.scatter(row.full_data_log2_ratio,i,s=16,color=color,zorder=3)
        count=int(consensus.loc[row.human_symbol,'both_model_and_treatment_concordant_n'])
        b.barh(i,count/4096*100,color=color,height=.48)
        b.text(104,i,f'{count}/4096',ha='left',va='center',fontsize=7)
    a.set_ylim(14.6,-.6);b.set_ylim(14.6,-.6)
    labels=[g+(' *' if low[g] else '') for g in order]
    a.set_yticks(range(15),labels);b.set_yticks(range(15),labels)
    for ax in axes:ax.tick_params(axis='y',length=0,pad=4)
    a.set_xlim(float(f.min_log2_ratio.min())-.12,float(f.max_log2_ratio.max())+.12)
    b.set_xlim(0,103);b.set_xticks([0,25,50,75,100])
    a.set_xlabel('Mouse sess / HFD log2 ratio\n(0.1 TPM pseudocount)',labelpad=7)
    b.set_xlabel('Configuration pairs (%)',labelpad=7)
    a.set_title('Range across 64 configurations',fontsize=8.5,pad=13)
    b.set_title('Both contrast directions agree',fontsize=8.5,pad=13)
    letter(a,'a');letter(b,'b')
    note(fig,'Left: complete min-max sensitivity ranges; points are full-data means. Ranges are not confidence intervals.\n'
        'Right: fractions of 64 × 64 overlapping configurations, not replication probabilities or independent samples.\n'
        '* All nine mouse libraries below 1 TPM. Original candidate order and zero-direction cases retained.\n'
        'ctrl / HFD / sess are source labels with provisional treatment mapping.')
    export(fig,'fig15_robustness',list(axes))
    f.to_csv(HERE/'fig15_mouse_range_source.csv',index=False,na_rep='NA')
    consensus.reset_index().to_csv(HERE/'fig15_concordance_source.csv',index=False,na_rep='NA')

def cell_context_figure(tables):
    local=tables['cell_localization'];a_data=local[local.symbol=='AREG'].copy()
    lib_order=['GSM6322577','GSM6322578','GSM6322579','GSM6322580']
    a_data=a_data.set_index('gsm').loc[lib_order].reset_index()
    assert len(a_data)==4 and np.isfinite(a_data[['cholangiocyte_detection_pct','hepatocyte_detection_pct']]).all().all()
    df=tables['cell_decomposition'];genes=['VNN1','AREG','JAG2','CXCL8']
    df=df[df.scheme=='equal_library_actual'].set_index('symbol').loc[genes].reset_index()
    assert np.allclose(df.composition_total+df.within_category_total,df.total_PA_minus_control,rtol=0,atol=1e-12)
    fig,axes=plt.subplots(1,2,figsize=(WIDTH,139/25.4))
    fig.subplots_adjust(left=.09,right=.97,bottom=.33,top=.82,wspace=.38)
    header(fig,'Library-resolved localization and descriptive cell-mixture accounting','Organoid control/PA reference | Four technical libraries from one differentiation | No sesamin')
    a,b=axes
    for i,row in enumerate(a_data.itertuples()):
        a.plot([row.hepatocyte_detection_pct,row.cholangiocyte_detection_pct],[i,i],color='#B4BFC6',lw=1.1)
        a.scatter(row.hepatocyte_detection_pct,i,s=30,color=BLUE,zorder=3)
        a.scatter(row.cholangiocyte_detection_pct,i,s=30,color=ORANGE,zorder=3)
    a.set_yticks(range(4),['Control 1','Control 2','PA 1','PA 2']);a.set_ylim(3.6,-.6)
    a.set_xlim(0,100);a.set_xlabel('AREG detection (%)',labelpad=7);a.set_title('AREG in each source library',pad=13)
    a.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=BLUE,markeredgecolor='none',markersize=5,label='Hepatocytes'),
               Line2D([0],[0],marker='o',color='none',markerfacecolor=ORANGE,markeredgecolor='none',markersize=5,label='Cholangiocytes')],
               loc='upper left',bbox_to_anchor=(-.1,-.21),fontsize=7.5,handletextpad=.4)
    b.axvline(0,color='#B5C0C8',lw=.7,zorder=0)
    y=np.arange(4)
    b.barh(y-.16,df.composition_total,height=.25,color=BLUE,label='Composition')
    b.barh(y+.16,df.within_category_total,height=.25,color=ORANGE,label='Within category')
    b.scatter(df.total_PA_minus_control,y,marker='D',s=19,color=INK,zorder=3,label='Total (sum)')
    b.set_yticks(y,genes);b.set_ylim(3.6,-.6)
    b.set_xlim(-.015,.17);b.set_xlabel('PA - control difference\nMean log1p(CP10k)',labelpad=7)
    b.set_title('Equal-library Kitagawa terms',pad=13)
    b.legend(loc='upper left',bbox_to_anchor=(-.06,-.21),fontsize=7.5,handlelength=1.5,ncol=1,labelspacing=.5)
    for ax in axes:ax.tick_params(axis='y',length=0,pad=4)
    letter(a,'a');letter(b,'b')
    note(fig,'Cell categories are the author labels. No cell-level or technical-library hypothesis test is shown.\n'
        'Accounting terms sum to the observed equal-library mean difference; they are not causal contributions.\n'
        'Within-category terms can include finer cell mixtures, depth and technical-library differences.',y=.02)
    export(fig,'fig15_cell_context',list(axes))
    a_data.to_csv(HERE/'fig15_AREG_library_source.csv',index=False)
    df.to_csv(HERE/'fig15_Kitagawa_spotlight_source.csv',index=False)

def pca_figure(tables):
    fig,axes=plt.subplots(1,3,figsize=(WIDTH,128/25.4))
    fig.subplots_adjust(left=.09,right=.98,bottom=.37,top=.80,wspace=.56)
    header(fig,'Global sample structure across three separate datasets','PCA for quality/context assessment | All submitted arrays or libraries retained')
    human=tables['human_pca'];mouse=tables['mouse_pca'];rat=tables['rat_pca'];rat=rat[rat.method=='deposited_log2']
    assert len(human)==len(mouse)==9 and len(rat)==20
    assert human.library.nunique()==9 and mouse.library.nunique()==9 and rat['sample'].nunique()==20
    variances=[tables['human_variance'].set_index('component').explained_variance_ratio.loc[['PC1','PC2']].to_numpy(),
       tables['mouse_variance'].set_index('component').explained_variance_ratio.loc[['PC1','PC2']].to_numpy(),
       tables['rat_variance'][tables['rat_variance'].method=='deposited_log2'][['pc1','pc2']].iloc[0].to_numpy()]
    styles=[{'Control':('Control',GREY,'o'),'Model':('OA',ORANGE,'s'),'Treat':('OA + sesamin',BLUE,'^')},
       {'ctrl':('ctrl',GREY,'o'),'HFD':('HFD',ORANGE,'s'),'sess':('sess',BLUE,'^')},
       {'control':('Control',GREY,'o'),'sesamin':('Sesamin',BLUE,'^'),'episesamin':('Episesamin','#A58B69','s'),'sesamolin':('Sesamolin','#9685AC','D')}]
    titles=['Human: 9 libraries','Mouse: 9 libraries','Rat: all 20 arrays']
    for ax,frame,var,style,title,panel_id in zip(axes,[human,mouse,rat],variances,styles,titles,'abc'):
        assert np.isfinite(frame[['PC1','PC2']]).all().all() and np.isfinite(var).all()
        assert set(frame.group)==set(style)
        for group,(label,color,mark) in style.items():
            sub=frame[frame.group==group]
            ax.scatter(sub.PC1,sub.PC2,s=23,c=color,marker=mark,linewidths=.45,edgecolors='white',label=label)
        ax.margins(x=.19,y=.2);ax.set_xlabel(f'PC1 ({var[0]*100:.1f}%)',labelpad=5)
        ax.set_ylabel(f'PC2 ({var[1]*100:.1f}%)',labelpad=5)
        ax.set_title(title,fontsize=8.5,pad=13);letter(ax,panel_id)
        ax.legend(loc='upper left',bbox_to_anchor=(-.08,-.27),fontsize=7,handletextpad=.3,labelspacing=.55)
    note(fig,'Human: 13,014 filtered genes; mouse: 15,340 filtered genes; rat: 31,099 deposited probes.\n'
        'Each panel uses its own axes and preprocessing. Rat episesamin/sesamolin arrays are not sesamin replicates.\n'
        'PCA positions do not establish treatment attribution, recovery, biological independence or cross-species similarity.',y=.02)
    export(fig,'fig15_global_pca',list(axes))

def main():
    tables={key:pd.read_csv(path) for key,path in INPUTS.items()}
    sources=[];dest=HERE/'source_tables';dest.mkdir(exist_ok=True)
    for key,path in INPUTS.items():
        checksum=sha(path);target=dest/(key+'.csv');shutil.copy2(path,target);assert sha(target)==checksum
        sources.append({'key':key,'source':str(path),'copied_file':str(target.relative_to(HERE)),'sha256':checksum,'rows':len(tables[key])})
    write('source_data.json',{'inputs':sources,'environment':{'python':platform.python_version(),'executable':sys.executable,
         'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__,'Pillow':PIL.__version__},
         'selection':'All candidate families and all submitted samples retained; no inference added for plotting'})
    selected=set(sys.argv[1:]) or {'rat','robustness','cell','pca'}
    assert selected.issubset({'rat','robustness','cell','pca'})
    if 'rat' in selected:rat_figure(tables)
    if 'robustness' in selected:robustness_figure(tables)
    if 'cell' in selected:cell_context_figure(tables)
    if 'pca' in selected:pca_figure(tables)
    assert all(sha(Path(x['source']))==x['sha256'] for x in sources)
    write('numerical_qa.json',{'status':'PASS','rat_candidates':25,'rat_mapped':15,'rat_NA':10,
       'rat_primary_animals':10,'mouse_robustness_candidates':15,'mouse_configurations':64,'cross_species_configuration_pairs':4096,
       'sc_highlight_genes':4,'AREG_technical_libraries':4,'PCA_samples':{'human':9,'mouse':9,'rat':20},
       'all_CI_and_range_values_finite_and_displayed':True,'Kitagawa_sum_identity_checked':True,
       'source_hashes_unchanged':True,'inference_added_for_plotting':False,'selected_render_keys':sorted(selected)})
    print(json.dumps({'rendered':sorted(selected),'sources_unchanged':True}))

if __name__=='__main__':main()
