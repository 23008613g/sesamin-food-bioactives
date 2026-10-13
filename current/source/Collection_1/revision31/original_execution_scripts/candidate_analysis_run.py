"""Recompute candidates from deposited counts and archived mouse gene TPM.

Original inputs are read-only. No statistical independence is assumed.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, itertools, json, platform, shutil, sys
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NEXT = HERE.parents[1]
INPUT = HERE / 'inputs'
INPUT.mkdir(exist_ok=True)
sys.stdout.reconfigure(encoding='utf-8')
def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()
def js(name, value):
    (HERE/name).write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
def factors(a):
    ref=a[(a>0).all(axis=1)]
    assert len(ref)>100
    sf=np.median(ref/np.exp(np.log(ref).mean(axis=1))[:,None],axis=0)
    return sf/np.exp(np.log(sf).mean())
def effects(a, groups):
    means=np.column_stack([a[:,g].mean(axis=1) for g in groups])
    c,m,t=means.T
    d=np.log2((m+1)/(c+1)); e=np.log2((t+1)/(m+1))
    between=(t>np.minimum(c,m)) & (t<np.maximum(c,m))
    return means,d,e,between

started=datetime.now(timezone.utc).isoformat()
freeze=HERE/'protocol_freeze.json'
if freeze.exists():
    assert json.loads(freeze.read_text())['sha256']==sha(HERE/'protocol.json')
else:
    js('protocol_freeze.json',{'utc':started,'sha256':sha(HERE/'protocol.json'),'before_new_expression_calculation':True,'prior_results_known':True})
sources={
 'human_counts.tsv.gz': NEXT/'candidate_screen_03/sesamin_inputs/GSE248903_total_count.anno.txt.gz',
 'frozen_human25.csv':NEXT/'sesamin_exploration_04/priority_all.csv',
 'frozen_mouse15.csv':NEXT/'analysis_strengthening_15/robustness/frozen_primary_candidates.csv',
 'mouse_mappings26.tsv':Path('D:/ResearchData/sesamin_mouse_2026-10-08/matrices/candidate_mouse_descriptive.tsv'),
 'mouse_gene_tpm.tsv.gz':Path('D:/ResearchData/sesamin_mouse_2026-10-08/matrices/gene_tpm.tsv.gz')}
for name,p in sources.items():
    if not (INPUT/name).exists(): shutil.copy2(p,INPUT/name)
    assert sha(p)==sha(INPUT/name)
js('input_manifest.json',[{'file':n,'original_path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for n,p in sources.items()])
hn=[f'{g}_{i}' for g in ['Control','Model','Treat'] for i in range(1,4)]
mn=[f'DRR{i}' for i in range(1010532,1010541)]
h=pd.read_csv(INPUT/'human_counts.tsv.gz',sep='\t',usecols=['name','GeneName']+hn)
x=h[hn].to_numpy(float); assert x.shape==(60605,9) and h.name.is_unique
assert np.isfinite(x).all() and (x>=0).all() and (x%1==0).all()
cpm=x/x.sum(axis=0)*1e6; keep=(cpm>=1).sum(axis=1)>=3
sf=factors(x); norm=x/sf; groups=[list(range(i,i+3)) for i in [0,3,6]]
means,d,e,between=effects(norm,groups)
_,cd,ce,cb=effects(cpm,groups)
agree=(np.sign(d)==np.sign(cd))&(np.sign(e)==np.sign(ce))&cb
spec=json.loads((HERE/'protocol.json').read_text()); thresholds=spec['adjacent_thresholds']
loo={t:np.zeros(len(h),int) for t in thresholds}
for drop in itertools.product(*groups):
    inds=[j for j in range(9) if j not in drop]
    sub=x[:,inds]; _,di,ei,bi=effects(sub/factors(sub),[list(range(i,i+2)) for i in [0,2,4]])
    signmatch=(np.sign(di)==np.sign(d))&(np.sign(ei)==np.sign(e))&bi
    for t in thresholds: loo[t]+=signmatch&(abs(di)>=np.log2(t))&(abs(ei)>=np.log2(t))
old=pd.read_csv(INPUT/'frozen_human25.csv'); hi=h.set_index('name').index.get_indexer(old.gene_id)
assert (hi>=0).all() and len(hi)==25
table=pd.DataFrame({'gene_id':h.name,'symbol':h.GeneName,'cpm_filter':keep,'control_mean':means[:,0], 'model_mean':means[:,1], 'treat_mean':means[:,2], 'model_control_log2':d,'treat_model_log2':e,'model_control_sign':np.sign(d).astype(int),'treat_model_sign':np.sign(e).astype(int),'treat_between':between,'cpm_agreement':agree})
summary=[]
for t in thresholds:
    base=keep&(abs(d)>=np.log2(t))&(abs(e)>=np.log2(t))&(d*e<0)&between
    selected=base&agree&(loo[t]>=22)
    table[f'loo_pass_{t}']=loo[t]; table[f'selected_{t}']=selected
    summary.append({'fold_threshold':t,'initial_n':int(base.sum()),'cpm_agree_n':int((base&agree).sum()),'selected_n':int(selected.sum()),'frozen25_retained':int(selected[hi].sum()),**{g:bool(selected[np.flatnonzero(h.GeneName.eq(g))[0]]) for g in spec['named_retention']}})
    if t==1.5: assert set(h.loc[selected,'name'])==set(old.gene_id)
table.to_csv(HERE/'human_all_gene_recalculation.csv.gz',index=False)
table.iloc[hi].to_csv(HERE/'human_frozen25_direction_order.csv',index=False)
pd.DataFrame(summary).to_csv(HERE/'adjacent_threshold_summary.csv',index=False)
samples=[]
for i in hi:
    for j,n in enumerate(hn): samples.append({'gene_id':h.name[i],'symbol':h.GeneName[i],'library':n,'group':['Control','Model','Treat'][j//3],'raw_count':int(x[i,j]),'median_ratio_count':norm[i,j],'CPM':cpm[i,j]})
pd.DataFrame(samples).to_csv(HERE/'human_candidate_samples.csv',index=False)
pd.DataFrame({'library':hn,'median_ratio_factor':sf}).to_csv(HERE/'human_size_factors.csv',index=False)
m=pd.read_csv(INPUT/'mouse_gene_tpm.tsv.gz',sep='\t'); primary=pd.read_csv(INPUT/'frozen_mouse15.csv')
mi=m.set_index('gene_id').index.get_indexer(primary.mouse_gene_id); assert (mi>=0).all()
a=m.loc[mi,mn].to_numpy(float); assert np.isfinite(a).all() and (a>=0).all()
mm=np.column_stack([a[:,g].mean(axis=1) for g in groups]); c,md,t=mm.T
mouse=primary.copy(); mouse['ctrl_mean_TPM']=c;mouse['HFD_mean_TPM']=md;mouse['sess_mean_TPM']=t
mouse['HFD_ctrl_log2']=np.log2((md+.1)/(c+.1));mouse['sess_HFD_log2']=np.log2((t+.1)/(md+.1))
mouse['HFD_ctrl_sign']=np.sign(md-c).astype(int);mouse['sess_HFD_sign']=np.sign(t-md).astype(int)
mouse['sess_between']= (t>np.minimum(c,md)) & (t<np.maximum(c,md));mouse['all_libraries_TPM_lt1']=(a<1).all(axis=1)
joined=mouse.merge(table.iloc[hi],left_on='original_human_id',right_on='gene_id',validate='one_to_one')
joined['two_signs_concordant']=(joined.HFD_ctrl_sign==joined.model_control_sign)&(joined.sess_HFD_sign==joined.treat_model_sign)&joined.HFD_ctrl_sign.ne(0)&joined.sess_HFD_sign.ne(0)
joined['intervention_label_sign_concordant']=(joined.sess_HFD_sign==joined.treat_model_sign)&joined.sess_HFD_sign.ne(0)
joined['mouse_exposure_mapping']='unconfirmed: label-level exploratory comparison'
joined.to_csv(HERE/'candidate_direction_order.csv',index=False)
ms=[]
for k,row in primary.iterrows():
    for j,n in enumerate(mn): ms.append({'symbol':row.human_symbol,'mouse_gene_id':row.mouse_gene_id,'library':n,'group':['ctrl','HFD','sess'][j//3],'TPM':a[k,j]})
pd.DataFrame(ms).to_csv(HERE/'mouse_candidate_samples.csv',index=False)
full=pd.read_csv(INPUT/'mouse_mappings26.tsv',sep='\t'); full.to_csv(HERE/'mouse_all_mappings26.csv',index=False)
checks={'human_old_new_means_max_abs':float(np.max(abs(means[hi]-old[['mean_norm_control','mean_norm_OA','mean_norm_sesamin_OA']].to_numpy()))),'human_fixed25_identical':True,'primary_mouse_n':len(primary),'both_signs':joined.loc[joined.two_signs_concordant,'human_symbol'].tolist(),'intervention_sign_n':int(joined.intervention_label_sign_concordant.sum()),'strict_mouse_between_n':int(joined.sess_between.sum())}
assert checks['human_old_new_means_max_abs']<1e-8
js('run_log.json',{'started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),'interpreter':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,'script_sha256':sha(__file__),'checks':checks,'scope':'human deposited-count recomputation; mouse archived gene-TPM recomputation, no new read quantification','threshold_results':summary})
print(json.dumps({'checks':checks,'thresholds':summary},indent=2))
