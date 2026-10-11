"""Eight-main-figure layout using frozen measurements and cached PCA coordinates.

Default: redraw only new Figures 2/3. --only all also regenerates the six inherited
artworks from the portable stage26 implementation. Delivered inherited files are
byte-exact copies; regeneration is optional and not required to reproduce inputs.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
import build_stage26 as inherited

HERE=Path(__file__).resolve().parent
DATA=HERE/'source_data'
QA=HERE/'qa'
WIDTH=174/25.4
INK,BLUE,ORANGE,GREY='#223744','#39748F','#B96F3E','#86949C'
CHECKS=[]

def check(name, condition):
    CHECKS.append({'check':name,'pass':bool(condition)})
    if not condition:raise AssertionError(name)

def read(name):return pd.read_csv(DATA/f'{name}.csv')

def letter(ax, value, dy=12):
    ax.text(0,1,value,transform=ax.transAxes+ScaledTranslation(-19/72,dy/72,ax.figure.dpi_scale_trans),
            fontsize=10,fontweight='bold',ha='left',va='bottom')

_export=inherited.export
RENAME={'Figure_01':'Figure_01','Figure_02':'Figure_04','Figure_S16':'Figure_05',
        'Figure_03':'Figure_06','Figure_04':'Figure_07','Figure_05':'Figure_08'}
def inherited_export(fig,stem,axes,supplement=False,exemptions=()):
    return _export(fig,RENAME[stem],axes,supplement=False,exemptions=exemptions)
inherited.export=inherited_export

def pca():
    frames=[read('human_pca'),read('mouse_pca'),read('rat_pca').query("method == 'deposited_log2'")]
    variances=[read('human_variance').set_index('component').loc[['PC1','PC2'],'explained_variance_ratio'].to_numpy(),
               read('mouse_variance').set_index('component').loc[['PC1','PC2'],'explained_variance_ratio'].to_numpy(),
               read('rat_variance').query("method == 'deposited_log2'")[['pc1','pc2']].iloc[0].to_numpy()]
    styles=[{'Control':('Control',GREY,'o'),'Model':('Model',ORANGE,'s'),'Treat':('Treat',BLUE,'^')},
            {'ctrl':('ctrl',GREY,'o'),'HFD':('HFD',ORANGE,'s'),'sess':('sess',BLUE,'^')},
            {'control':('Control',GREY,'o'),'sesamin':('Sesamin',BLUE,'^'),
             'episesamin':('Episesamin','#A58B69','s'),'sesamolin':('Sesamolin','#9685AC','D')}]
    fig,axes=plt.subplots(1,3,figsize=(WIDTH,80/25.4))
    fig.subplots_adjust(left=.105,right=.985,bottom=.365,top=.84,wspace=.69)
    tidy=[]
    for ax,frame,var,style,title,panel,n in zip(axes,frames,variances,styles,
            ['Human HepG2','Mouse hepatocytes','Rat liver'],'abc',[9,9,20]):
        check(f'PCA {panel}: all {n} cached observations',len(frame)==n)
        check(f'PCA {panel}: finite coordinates',np.isfinite(frame[['PC1','PC2']]).all().all())
        check(f'PCA {panel}: complete group set',set(frame.group)==set(style))
        for group,(name,color,mark) in style.items():
            subset=frame.loc[frame.group==group]
            check(f'PCA {panel}/{group}: complete group',len(subset)==(5 if panel=='c' else 3))
            ax.scatter(subset.PC1,subset.PC2,s=23,c=color,marker=mark,edgecolors='white',linewidths=.4,label=name)
            for _,row in subset.iterrows():
                tidy.append({'panel':panel,'dataset':title,'observation':row.get('library',row.get('sample')),
                             'source_group':group,'display_group':name,'PC1':row.PC1,'PC2':row.PC2,
                             'PC1_explained_fraction':var[0],'PC2_explained_fraction':var[1]})
        ax.margins(x=.19,y=.2)
        ax.set_xlabel(f'PC1 ({var[0]*100:.1f}%)',labelpad=4)
        ax.set_ylabel(f'PC2 ({var[1]*100:.1f}%)',labelpad=4)
        ax.set_title(title,pad=12,fontsize=9)
        ax.xaxis.set_major_locator(plt.MaxNLocator(3))
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        letter(ax,panel)
        ax.legend(loc='upper left',bbox_to_anchor=(-.12,-.28),fontsize=8,
                  handletextpad=.3,labelspacing=.32,borderaxespad=0)
    pd.DataFrame(tidy).to_csv(DATA/'Figure_02_PCA_all38_observations.csv',index=False)
    _export(fig,'Figure_02',list(axes))

def candidates():
    human=read('human_all25_heatmap').sort_values(['row_order','column_order'])
    direct=read('primary15_directions')
    check('Human heatmap contains all 225 values',len(human)==225)
    check('Human heatmap has 25 candidates and 9 libraries',human.gene_id.nunique()==25 and human['sample'].nunique()==9)
    check('Human heatmap contains every row/column exactly once',not human.duplicated(['row_order','column_order']).any())
    check('Normalized heatmap counts are finite and nonnegative',
          np.isfinite(human.normalized_count_verbatim).all() and (human.normalized_count_verbatim >= 0).all())
    expected=np.log2(human.normalized_count_verbatim.to_numpy()+1)
    check('Existing log2 transform matches copied normalized counts',np.allclose(expected,human.log2_count_plus1,rtol=0,atol=1e-12))
    centered=human.log2_count_plus1-human.groupby('gene_id').log2_count_plus1.transform('mean')
    check('Existing gene centering matches all nine libraries',np.allclose(centered,human.centered_log2_count_plus1,rtol=0,atol=1e-12))
    matrix=human.pivot(index='row_order',columns='column_order',values='centered_log2_count_plus1').to_numpy()
    check('No heatmap values missing or color-clipped',np.isfinite(matrix).all() and np.max(np.abs(matrix))<2)
    fields=['human_OA_vs_ctrl_direction','human_sesamin_vs_OA_direction',
            'mouse_HFD_vs_ctrl_direction','mouse_sess_vs_HFD_direction']
    signs={'down':-1,'unchanged':0,'up':1}
    values=np.array([[signs[row[col]] for col in fields] for _,row in direct.iterrows()])
    check('Direction panel contains all 15 eligible candidates',len(direct)==15)
    check('Exactly five observed zero direction cells retained',int((values==0).sum())==5)
    check('Exactly five low-TPM rows retained',int(direct.all_mouse_libraries_below_1_TPM.sum())==5)
    check('Three both-direction and eight treatment-direction agreements retained',
          int(direct.both_cross_species_directions_same.sum())==3 and int((direct.cross_species_sesamin_direction=='same').sum())==8)
    fig=plt.figure(figsize=(WIDTH,160/25.4))
    grid=fig.add_gridspec(1,2,left=.225,right=.97,bottom=.16,top=.84,
                         width_ratios=[.31,.245],wspace=.19/.2775)
    a=fig.add_subplot(grid[0,0])
    b=fig.add_subplot(grid[0,1])
    cmap=LinearSegmentedColormap.from_list('source_continuous',['#24728C','#F7F7F4','#B86B2D'])
    im=a.pcolormesh(np.arange(10)-.5,np.arange(26)-.5,matrix,
                    cmap=cmap,vmin=-2,vmax=2,shading='flat',rasterized=False)
    a.set(xlim=(-.5,8.5),ylim=(24.5,-.5))
    gene_names=human.drop_duplicates('row_order').symbol.tolist()
    a.set_yticks(range(25),gene_names)
    a.set_xticks(range(9),['1','2','3']*3)
    for pos,name in zip([1,4,7],['Control','Model','Treat']):
        a.text(pos,1.053,name,transform=a.get_xaxis_transform(),ha='center',va='bottom',fontsize=8)
    for x in [2.5,5.5]:a.axvline(x,color='white',lw=2)
    b.pcolormesh(np.arange(5)-.5,np.arange(16)-.5,values,
                 cmap=ListedColormap(['#D5E6EF','#F0F1F1','#F2DCC7']),
                 norm=BoundaryNorm([-1.5,-.5,.5,1.5],3),shading='flat',rasterized=False)
    b.set(xlim=(-.5,3.5),ylim=(14.5,-.5))
    chars={-1:'\u2193',0:'0',1:'\u2191'}
    for i in range(15):
        for j in range(4):b.text(j,i,chars[values[i,j]],ha='center',va='center',fontsize=10,color=INK)
    b.set_yticks(range(15),[row.human_symbol+(' *' if row.all_mouse_libraries_below_1_TPM else '') for _,row in direct.iterrows()])
    b.set_xticks(range(4),['Model /\nControl','Treat /\nModel','HFD /\nctrl','sess /\nHFD'])
    b.text(.5,1.09,'Human',transform=b.get_xaxis_transform(),ha='center',va='bottom',fontsize=8)
    b.text(2.5,1.09,'Mouse',transform=b.get_xaxis_transform(),ha='center',va='bottom',fontsize=8)
    for ax,nrows,ncols,head,panel in [(a,25,9,'Human: all 25','a'),(b,15,4,'Primary one-to-one','b')]:
        ax.xaxis.tick_top()
        ax.tick_params(axis='both',length=0,pad=4,labelsize=8)
        ax.set_xticks(np.arange(-.5,ncols,1),minor=True)
        ax.set_yticks(np.arange(-.5,nrows,1),minor=True)
        ax.grid(which='minor',color='white',linewidth=.5)
        ax.tick_params(which='minor',length=0)
        for spine in ax.spines.values():spine.set_visible(False)
        ax.text(.5,1.16,head,transform=ax.transAxes,ha='center',va='bottom',fontsize=9)
        letter(ax,panel,dy=44)
    cax=fig.add_axes([.225,.080,.31,.015])
    cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-2,0,2])
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')
    cb.outline.set_visible(False)
    cb.set_label('Gene-centered log2(count + 1)',fontsize=8,labelpad=4)
    cb.ax.tick_params(labelsize=8,length=2)
    for x,txt in [(.61,'\u2191 Higher'),(.745,'\u2193 Lower'),(.87,'0 Unchanged')]:
        fig.text(x,.080,txt,fontsize=8,va='center')
    fig.text(.61,.045,'* All 9 mouse libraries <1 TPM',fontsize=8,va='center')
    _export(fig,'Figure_03',[a,b],exemptions=[{'panels':['a','b'],'checks':['panel-width'],
        'reason':'The 9-column continuous human heatmap and 4-column categorical direction panel intentionally use different widths; plot top/bottom remain identical.'}])

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--only',choices=['new','all','2','3'],default='new')
    args=parser.parse_args()
    if args.only in ['all','new','2']:pca()
    if args.only in ['all','new','3']:candidates()
    if args.only=='all':
        for fn in (inherited.figure1,inherited.figure2,inherited.figureS16,inherited.figure3,inherited.figure4,inherited.figure5):fn()
    result={'scope':'Display-only expansion; no new analysis or inference','checks':CHECKS,
            'environment':{'python':platform.python_version(),'executable':sys.executable,
                           'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},
            'width_mm':174,'heights_mm':{'1':104,'2':80,'3':160,'4':80,'5':137,'6':106,'7':80,'8':174},
            'dpi':600,'minimum_intended_font_pt':8,'last_run_selection':args.only,
            'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.iterdir() if p.is_file()},
            'output_sha256':{str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'main').glob('*') if p.is_file()}}
    (QA/'render_manifest.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('All',len(CHECKS),'executed source/transform checks passed')
