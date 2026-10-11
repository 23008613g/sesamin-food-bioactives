"""Extract only three existing figure functions, correcting display labels.

No original generator is imported or executed in full; no analyses rerun.
"""
from pathlib import Path
import ast
import difflib

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
NEXT=ROOT.parent
SUPPORT=NEXT/'publication_package_12/figures/support/make_support_figures.py'
PCA=NEXT/'analysis_strengthening_15/figures/build_figures.py'

def extract(path,names):
    text=path.read_text(encoding='utf-8')
    tree=ast.parse(text)
    found={node.name:ast.get_source_segment(text,node) for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in names}
    assert set(found)==set(names)
    return '\n\n'.join(found[name] for name in names)+'\n'

support=extract(SUPPORT,['panel_label','figure_s1','figure_s2'])
pca=extract(PCA,['header','letter','note','pca_figure'])
old=support+'\n'+pca
new=old
for a,b,count in [
    ("label=group, zorder=3", "label={'Control':'Control','OA':'Model','OA + sesamin':'Treat'}[group], zorder=3",1),
    ("'OA / Control'", "'Model / Control'",1),
    ("'OA + Sesamin / OA'", "'Treat / Model'",1),
    ("'OA + Sesamin / Control'", "'Treat / Control'",1),
    ("'Model':('OA',ORANGE,'s'),'Treat':('OA + sesamin',BLUE,'^')", "'Model':('Model',ORANGE,'s'),'Treat':('Treat',BLUE,'^')",1),
]:
    assert new.count(a)==count,(a,new.count(a))
    new=new.replace(a,b)
(HERE/'extracted_original_functions.py').write_text(old,encoding='utf-8')
(HERE/'display_label_functions.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='original functions',tofile='display-label correction')),encoding='utf-8')
prefix='''"""Regenerate S01/S02/S13 from cached figure inputs; display-label changes only."""
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

'''
suffix='''
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
    (HERE/'numeric_preservation.json').write_text(json.dumps(record,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\\n',encoding='utf-8')
    print(json.dumps({'numeric_artists_and_axes_exactly_equal':True,'figures':list(STEMS.values())}))
if __name__=='__main__':main()
'''
(HERE/'revise_figure_labels.py').write_text(prefix+new+suffix,encoding='utf-8')
print(HERE/'revise_figure_labels.py')
