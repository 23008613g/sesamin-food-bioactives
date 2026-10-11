"""Regenerate S01/S02/S13 from cached figure inputs; display-label changes only."""
from pathlib import Path
import hashlib, json, platform, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from PIL import Image
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
NEXT=ROOT.parent
OUT=ROOT/'figures/supplement'
QA=HERE/'qa'
QA.mkdir(exist_ok=True)
sys.path.insert(0,str(Path.home()/'.agents/skills/nature-figure/scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment
GENES=['VNN1','AREG','JAG2','CXCL8']
WIDTH_IN=WIDTH=180/25.4
MODE='original'
SIGNATURES={}
INPUT_HASHES={}
STEMS={'figS1_quality':'Figure_S01','figS2_sensitivity':'Figure_S02','fig15_global_pca':'Figure_S13'}

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_csv(p):
    p=Path(p); INPUT_HASHES[str(p)]=sha(p); return pd.read_csv(p)
def read_json(p):
    p=Path(p); INPUT_HASHES[str(p)]=sha(p); return json.loads(p.read_text(encoding='utf-8'))
def apply_support_style():
    global BLUE,ORANGE,NEUTRAL,GROUP,MARK
    plt.rcdefaults()
    plt.rcParams.update({
        'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],'font.size':7.5,
        'axes.titlesize':8,'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,
        'legend.fontsize':7,'axes.linewidth':.65,'axes.spines.top':False,'axes.spines.right':False,
        'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'savefig.facecolor':'white'})
    BLUE,ORANGE,NEUTRAL='#3E6D9C','#C47738','#6C7580'
    GROUP={'Control':NEUTRAL,'OA':ORANGE,'OA + sesamin':BLUE}
    MARK={'Control':'o','OA':'s','OA + sesamin':'^'}
def apply_pca_style():
    global INK,BLUE,ORANGE,GREY,PALE
    INK='#263443';BLUE='#3C7F99';ORANGE='#BB7B43';GREY='#8C959E';PALE='#CCD6DC'
    plt.rcdefaults()
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
     'font.size':7.5,'axes.labelsize':7.5,'axes.titlesize':8.5,'xtick.labelsize':7,'ytick.labelsize':7,
     'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'text.color':INK,
     'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,'axes.edgecolor':'#73818C',
     'pdf.fonttype':42,'svg.fonttype':'none','legend.frameon':False,'figure.facecolor':'white','savefig.facecolor':'white'})
def figure_signature(fig):
    result=[]
    for ax in fig.axes:
        result.append({'bounds':list(ax.get_position().bounds),'xlim':list(ax.get_xlim()),'ylim':list(ax.get_ylim()),
            'lines':[{'x':np.asarray(line.get_xdata()).tolist(),'y':np.asarray(line.get_ydata()).tolist(),
                      'color':line.get_color(),'marker':line.get_marker(),'linewidth':line.get_linewidth()} for line in ax.lines],
            'collections':[{'offsets':collection.get_offsets().tolist(),'facecolors':collection.get_facecolors().tolist(),
                            'sizes':collection.get_sizes().tolist()} for collection in ax.collections],
            'xticks':list(ax.get_xticks()),'yticks':list(ax.get_yticks())})
    return result
def export(fig,stem,*positional,**alignment):
    fig.canvas.draw()
    name=STEMS[stem]
    signature=figure_signature(fig)
    if MODE=='original':
        SIGNATURES[name]=signature
    else:
        assert signature==SIGNATURES[name],name+' numerical geometry changed'
        if positional: alignment.update(axes=positional[0],require_panel_labels=len(positional[0])>1)
        require_matplotlib_panel_alignment(fig,json_out=QA/(name+'.alignment.json'),
            overlay_svg=QA/(name+'.alignment.svg'),tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True,**alignment)
        fig.savefig(OUT/(name+'.pdf'),dpi=600)
        fig.savefig(OUT/(name+'.svg'),dpi=600)
        fig.savefig(OUT/(name+'.png'),dpi=600)
        fig.savefig(OUT/(name+'.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'})
        for suffix in ('png','tiff'):
            with Image.open(OUT/(name+'.'+suffix)) as im:
                assert min(im.info['dpi'])>599
                assert im.width>=4251
    plt.close(fig)

def panel_label(ax, letter):
    ax.annotate(letter, xy=(0, 1), xycoords='axes fraction',
                xytext=(-23, 10), textcoords='offset points',
                fontsize=9, fontweight='bold', va='bottom', ha='left')

def figure_s1(summary, pca, qc, mouse):
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, 87/25.4))
    fig.subplots_adjust(left=.09, right=.98, bottom=.20, top=.82, wspace=.40)
    a, b = axes
    for group in ['Control', 'OA', 'OA + sesamin']:
        rows = pca.loc[pca.group == group]
        a.scatter(rows.PC1, rows.PC2, c=GROUP[group], marker=MARK[group], s=30,
                  edgecolors='white', linewidths=.35, label={'Control':'Control','OA':'Model','OA + sesamin':'Treat'}[group], zorder=3)
    a.set_xlim(-23, 16)
    a.set_ylim(-10.5, 11)
    a.set_xlabel(f"PC1 ({100*summary['pca_explained_variance'][0]:.1f}%)")
    a.set_ylabel(f"PC2 ({100*summary['pca_explained_variance'][1]:.1f}%)")
    a.set_title('HepG2 library structure', pad=14)
    a.legend(loc='upper left', ncol=1, handletextpad=.4, borderaxespad=.2, fontsize=7)
    a.text(.01, -.29, '9 libraries; cached PCA coordinates', transform=a.transAxes, fontsize=7, color=NEUTRAL)
    y = np.arange(len(mouse))
    for group, color, marker in [('ctrl', NEUTRAL, 'o'), ('HFD', ORANGE, 's'), ('sess', BLUE, '^')]:
        idx = mouse['original_group'] == group
        b.scatter(mouse.loc[idx, 'pseudoalignment_percent'], y[idx], c=color, marker=marker,
                  s=25, edgecolors='white', linewidths=.35, zorder=3)
    b.set_yticks(y, mouse.original_label)
    b.invert_yaxis()
    b.set_ylim(8.6, -.6)
    b.set_xlim(80, 91)
    b.set_xticks([80, 85, 90])
    b.set_xlabel('Pseudoaligned read pairs (%)')
    b.set_title('Mouse library pseudoalignment', pad=14)
    for i, row in mouse.iterrows():
        b.annotate(f'{row.pseudoalignment_percent:.2f}', (row.pseudoalignment_percent, i),
                   xytext=(5,0), textcoords='offset points', va='center', fontsize=7)
    b.text(.01, -.29, '9 libraries; original sample labels', transform=b.transAxes, fontsize=7, color=NEUTRAL)
    for letter, ax in zip('ab', axes):
        panel_label(ax, letter)
    export(fig, 'figS1_quality', panel_ids=['a','b'], require_panel_labels=True)

def figure_s2(pseudo, omission):
    fig, axes = plt.subplots(2, 4, figsize=(WIDTH_IN, 150/25.4), sharey='row')
    fig.subplots_adjust(left=.09, right=.985, bottom=.16, top=.865, wspace=.33, hspace=.93)
    labels = [('OA_vs_control','Model / Control', ORANGE, 'o'),
              ('sesamin_vs_OA','Treat / Model', BLUE, 's'),
              ('sesamin_vs_control','Treat / Control', NEUTRAL, '^')]
    for j, gene in enumerate(GENES):
        a, b = axes[:,j]
        rows = pseudo.loc[pseudo.symbol == gene].sort_values('pseudocount')
        assert len(rows) == 4 and set(rows.pseudocount) == {0, .5, 1, 5}
        for field, label, color, marker in labels:
            a.plot(range(4), rows[field], color=color, marker=marker, ms=3.7, lw=.9, label=label)
        a.axhline(0, lw=.5, ls='--', color='#A5ADB5', zorder=0)
        a.set_xticks(range(4), ['0', '0.5', '1', '5'])
        a.set_ylim(-1.5, 2.1)
        a.set_xlabel('Pseudocount p')
        a.set_title(gene, fontstyle='italic', pad=12)
        donor_order = ['F95','M95','H20','H961','F95P']
        rows = omission.loc[omission.symbol == gene].set_index('omitted_label').loc[donor_order]
        assert len(rows) == 5
        b.scatter(range(5), rows.median_effect_other4, c=BLUE, s=22, edgecolors='white', linewidths=.35, zorder=3)
        b.axhline(0, lw=.5, ls='--', color='#A5ADB5', zorder=0)
        b.set_xticks(range(5), donor_order, rotation=45, ha='right', rotation_mode='anchor')
        b.set_xlim(-.5, 4.5)
        b.set_ylim(-.65, .55)
        b.set_xlabel('Omitted donor label')
        b.set_title(gene, fontstyle='italic', pad=12)
    axes[0,0].set_ylabel('Descriptive log2 ratio')
    axes[1,0].set_ylabel('Median log2 ratio, other 4 labels')
    for letter, ax in zip('abcdefgh', axes.flat):
        panel_label(ax, letter)
    fig.text(.09, .965, 'HepG2: sensitivity to the additive pseudocount', fontsize=8.5, fontweight='bold')
    handles, texts = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, texts, loc='upper center', bbox_to_anchor=(.54,.944), ncol=3,
               handlelength=1.8, columnspacing=1.7, borderaxespad=0, fontsize=7)
    fig.text(.09, .49, 'PHH: leave-one-label-out FFA / Control summaries (no sesamin)', fontsize=8.5, fontweight='bold')
    fig.text(.09, .035, 'All settings shown; omission subsets overlap and are not independent replicates or confidence intervals.',
             fontsize=7, color=NEUTRAL)
    export(fig, 'figS2_sensitivity', panel_ids=list('abcdefgh'), require_panel_labels=True)

def header(fig,title,sub):
    fig.text(.045,.98,title,fontsize=10,weight='bold',va='top')
    fig.text(.045,.944,sub,fontsize=7.5,va='top')

def letter(ax,x):
    ax.text(0,1,x,transform=ax.transAxes+ScaledTranslation(-19/72,11/72,ax.figure.dpi_scale_trans),
            fontsize=10,weight='bold',ha='left',va='bottom')

def note(fig,text,y=.025):fig.text(.045,y,text,fontsize=7,va='bottom',linespacing=1.55)

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
    styles=[{'Control':('Control',GREY,'o'),'Model':('Model',ORANGE,'s'),'Treat':('Treat',BLUE,'^')},
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
    global MODE
    from shutil import copy2
    before=HERE/'before';before.mkdir(exist_ok=True)
    for name in STEMS.values():
        for ext in ('pdf','svg','png','tiff'):
            p=OUT/(name+'.'+ext)
            if not (before/p.name).exists(): copy2(p,before/p.name)
    support=NEXT/'publication_package_12/figures/support'
    summary=read_json(NEXT/'sesamin_exploration_04/analysis_summary.json')
    pca=read_csv(support/'source_data/figS1_human_library_qc.csv')
    mouse=read_csv(support/'source_data/figS1_mouse_library_qc.csv')
    pseudo=read_csv(support/'source_data/figS2_pseudocounts.csv')
    omission=read_csv(support/'source_data/figS2_label_omission.csv')
    assert len(pseudo)==16 and len(omission)==20
    base=NEXT/'analysis_strengthening_15'
    paths={
      'human_pca':base/'robustness/human_pca_scores.csv','human_variance':base/'robustness/human_pca_variance.csv',
      'mouse_pca':base/'robustness/mouse_pca_scores.csv','mouse_variance':base/'robustness/mouse_pca_variance.csv',
      'rat_pca':base/'rat_validation/outputs/all20_array_pca.csv','rat_variance':base/'rat_validation/outputs/pca_variance.csv'}
    tables={key:read_csv(path) for key,path in paths.items()}
    original={}
    exec(compile((HERE/'extracted_original_functions.py').read_text(encoding='utf-8'),'extracted_original_functions.py','exec'),globals(),original)
    # Helper definitions are identical; only plot functions differ in display labels.
    apply_support_style();MODE='original'
    original['figure_s1'](summary,pca,None,mouse);original['figure_s2'](pseudo,omission)
    apply_pca_style();original['pca_figure'](tables)
    apply_support_style();MODE='corrected'
    figure_s1(summary,pca,None,mouse);figure_s2(pseudo,omission)
    apply_pca_style();pca_figure(tables)
    assert all(sha(Path(path))==checksum for path,checksum in INPUT_HASHES.items())
    record={'scope':'Display labels only; no scientific analysis rerun',
       'environment':{'python':platform.python_version(),'executable':sys.executable,'matplotlib':matplotlib.__version__,
                      'numpy':np.__version__,'pandas':pd.__version__},
       'inputs_sha256_unchanged':INPUT_HASHES,'numeric_artists_and_axes_exactly_equal':True,
       'figure_signatures':SIGNATURES,'pseudocount_rows':16,'omission_rows':20,
       'S01_human_libraries':9,'S01_mouse_libraries':9,'S13_human_libraries':9,'S13_mouse_libraries':9,'S13_rat_arrays':20,
       'outputs_sha256':{p.name:sha(p) for name in STEMS.values() for p in OUT.glob(name+'.*')},
       'source_functions_sha256':sha(HERE/'extracted_original_functions.py'),'generator_sha256':sha(Path(__file__))}
    (HERE/'numeric_preservation.json').write_text(json.dumps(record,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n',encoding='utf-8')
    print(json.dumps({'numeric_artists_and_axes_exactly_equal':True,'figures':list(STEMS.values())}))
if __name__=='__main__':main()
