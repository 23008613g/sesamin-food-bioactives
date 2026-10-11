"""Frozen candidate analysis of an independent rat-liver sesamin experiment."""
from pathlib import Path
from datetime import datetime, timezone
import csv, gzip, hashlib, io, json, platform, sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats

S=Path(__file__).resolve().parent;N=S.parents[1];O=S/'outputs';O.mkdir(exist_ok=True)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,j): p.write_text(json.dumps(j,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def out(df,name): df.to_csv(O/name,index=False,na_rep='NA',float_format='%.12g')
matrix=N/'evidence_closure_06/intervention_access/sources/GSE42693_series_matrix.txt.gz'
annotation=S/'sources/GPL1355.annot.gz'
candidatefile=N/'publication_package_14/data/human_candidates_all25.csv'
orthfile=S/'sources/human_rat_orthology_ensembl116.tsv'
xreffile=S/'sources/rat_entrez_ensembl116.tsv'
inputs=[matrix,annotation,candidatefile,orthfile,xreffile,S/'protocol.md']
fingerprints=[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in inputs]
text=gzip.open(matrix,'rt',encoding='utf-8').read()
body=text.split('!series_matrix_table_begin\n')[1].split('!series_matrix_table_end')[0]
m=pd.read_csv(io.StringIO(body),sep='\t',index_col=0)
assert m.shape==(31099,20) and m.index.is_unique and m.columns.is_unique
assert np.isfinite(m.to_numpy()).all() and (m.to_numpy()>=0).all()
expected=[f'GSM{i}' for i in range(1048036,1048056)];assert m.columns.tolist()==expected
annotext=gzip.open(annotation,'rt',encoding='utf-8').read()
annbody=annotext.split('!platform_table_begin\n')[1].split('!platform_table_end')[0]
ann=pd.read_csv(io.StringIO(annbody),sep='\t',dtype=str,keep_default_na=False)
assert len(ann)==31099 and ann['ID'].is_unique and set(ann.ID)==set(m.index)
ann['entrez_set']=ann['Gene ID'].map(lambda x:tuple(sorted(set(y.strip() for y in x.split('///') if y.strip() and y.strip()!='---'))))
ann['symbol_set']=ann['Gene symbol'].map(lambda x:tuple(sorted(set(y.strip() for y in x.split('///') if y.strip() and y.strip()!='---'))))
orth=pd.read_csv(orthfile,sep='\t',dtype=str,keep_default_na=False)
orth.columns=['human_id','current_human_symbol','rat_id','rat_symbol','orthology_type','orthology_confidence','rat_percent_id','human_percent_id']
xref=pd.read_csv(xreffile,sep='\t',dtype=str,keep_default_na=False);xref.columns=['rat_id','rat_symbol','rat_entrez']
candidates=pd.read_csv(candidatefile,dtype={'gene_id':str})
mapping=[];excluded=[]
for _,c in candidates.iterrows():
    hs=orth.loc[orth.human_id==c.gene_id]
    targets=hs.loc[hs.rat_id!='']
    r={'human_id':c.gene_id,'human_symbol':c.symbol,'rat_id':'','rat_symbol':'','rat_entrez':'','orthology_type':'','orthology_confidence':'','probe_ids':'','probe_n':0,'status':''}
    if len(targets)!=1 or targets.iloc[0].orthology_type!='ortholog_one2one':
        r['status']='no_rat_orthologue' if len(targets)==0 else 'non_one_to_one_orthology'
        r['orthology_type']=';'.join(sorted(set(targets.orthology_type)))
    else:
        t=targets.iloc[0];r.update(t[['rat_id','rat_symbol','orthology_type','orthology_confidence']].to_dict())
        ids=sorted(set(xref.loc[xref.rat_id==r['rat_id'],'rat_entrez'])-{''})
        if len(ids)!=1:r['status']='ambiguous_or_missing_rat_entrez'
        else:
            eid=ids[0];r['rat_entrez']=eid
            hit=ann.loc[ann.entrez_set.map(lambda x:eid in x)]
            clean=hit.loc[hit.entrez_set.map(len)==1]
            for _,p in hit.loc[hit.entrez_set.map(len)>1].iterrows():
                excluded.append({'human_id':c.gene_id,'probe_id':p.ID,'entrez_ids':'///'.join(p.entrez_set),'reason':'multiple_gene_IDs'})
            r['probe_ids']=';'.join(sorted(clean.ID));r['probe_n']=len(clean)
            r['status']='mapped' if len(clean) else 'no_unambiguous_probe'
    mapping.append(r)
mapdf=pd.DataFrame(mapping);out(mapdf,'candidate_mapping_all25.csv')
out(pd.DataFrame(excluded,columns=['human_id','probe_id','entrez_ids','reason']),'excluded_ambiguous_probes.csv')

# Freeze mapping before deriving candidate contrasts.
freeze=S/'mapping_freeze.json'
mr={'created_utc':datetime.now(timezone.utc).isoformat(),'mapping_sha256':sha(O/'candidate_mapping_all25.csv'),'protocol_sha256':sha(S/'protocol.md'),'mapped_candidates':int((mapdf.status=='mapped').sum()),'effects_inspected_when_frozen':False}
if freeze.exists(): assert json.loads(freeze.read_text(encoding='utf-8'))['mapping_sha256']==mr['mapping_sha256']
else: write(freeze,mr)

log=np.log2(m.to_numpy()+1.)
def quantile_normalize(a):
    target=np.sort(a,axis=0).mean(axis=1);result=np.empty_like(a)
    for j in range(a.shape[1]):
        order=np.argsort(a[:,j],kind='stable');values=a[order,j]
        starts=np.r_[0,np.flatnonzero(np.diff(values)!=0)+1];ends=np.r_[starts[1:],len(values)]
        for start,end in zip(starts,ends):result[order[start:end],j]=target[start:end].mean()
    return result
qn=quantile_normalize(log)
def compare(v):
    a=v[:5];b=v[5:10];effect=float(b.mean()-a.mean())
    va=a.var(ddof=1)/5;vb=b.var(ddof=1)/5;se=np.sqrt(va+vb)
    if se==0:return {'effect_log2_signal':effect,'p_value':np.nan,'ci95_low':np.nan,'ci95_high':np.nan,'welch_df':np.nan}
    df=float((va+vb)**2/(va**2/4+vb**2/4));delta=stats.t.ppf(.975,df)*se
    test=stats.ttest_ind(b,a,equal_var=False)
    assert np.isclose(test.statistic,effect/se,rtol=1e-11,atol=1e-11)
    return {'effect_log2_signal':effect,'p_value':float(test.pvalue),'ci95_low':float(effect-delta),'ci95_high':float(effect+delta),'welch_df':df}
def bh(p):
    p=np.nan_to_num(np.array(p,float),nan=1.);order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    r=np.empty_like(q);r[order]=np.minimum(q,1);return r

results=[];samples=[];probes=[]
for method,arr in [('deposited_log2',log),('quantile_sensitivity',qn)]:
    for _,r in mapdf.iterrows():
        base=r.to_dict();base['method']=method
        if r.status!='mapped':
            results.append({**base,'effect_log2_signal':np.nan,'p_value':np.nan,'ci95_low':np.nan,'ci95_high':np.nan});continue
        pr=r.probe_ids.split(';');ix=m.index.get_indexer(pr);assert (ix>=0).all()
        vals=np.median(arr[ix,:],axis=0);comp=compare(vals)
        probeeffects=[float(arr[i,5:10].mean()-arr[i,:5].mean()) for i in ix]
        results.append({**base,**comp,'control_mean_log2':float(vals[:5].mean()),'sesamin_mean_log2':float(vals[5:10].mean()),'control_sd_log2':float(vals[:5].std(ddof=1)),'sesamin_sd_log2':float(vals[5:10].std(ddof=1)),'median_group_difference':float(np.median(vals[5:10])-np.median(vals[:5])),'probe_effect_min':min(probeeffects),'probe_effect_max':max(probeeffects),'probe_directions_mixed':min(probeeffects)<0<max(probeeffects),'min_primary_raw_signal':float(m.iloc[ix,:10].to_numpy().min()),'max_primary_raw_signal':float(m.iloc[ix,:10].to_numpy().max())})
        for sample,v in zip(m.columns,vals):
            samples.append({'human_id':r.human_id,'human_symbol':r.human_symbol,'method':method,'sample':sample,'group':['control','sesamin','episesamin','sesamolin'][expected.index(sample)//5],'aggregate_log2_signal':float(v),'primary_comparison':expected.index(sample)<10})
        for pid,i,eff in zip(pr,ix,probeeffects):
            ap=ann.loc[ann.ID==pid].iloc[0]
            probes.append({'human_id':r.human_id,'human_symbol':r.human_symbol,'rat_symbol':r.rat_symbol,'probe_id':pid,'deposited_symbols':ap['Gene symbol'],'entrez_ids':ap['Gene ID'],'method':method,'effect_log2_signal':eff,**{sample:float(arr[i,j]) for j,sample in enumerate(m.columns)}})
res=pd.DataFrame(results)
for method in res.method.unique():
    ix=res.index[res.method==method];assert len(ix)==25
    res.loc[ix,'q_bh_fixed25']=bh(res.loc[ix,'p_value'])
    res.loc[ix[res.loc[ix,'p_value'].isna()],'q_bh_fixed25']=np.nan
res['fdr05']=res.q_bh_fixed25<.05
out(res,'candidate_results_all25_two_normalizations.csv');out(pd.DataFrame(samples),'candidate_sample_values.csv');out(pd.DataFrame(probes),'candidate_probe_values.csv')
qc=[];pca=[];variance=[]
for method,arr in [('deposited_log2',log),('quantile_sensitivity',qn)]:
    centered=arr-arr.mean(axis=1,keepdims=True)
    u,d,vt=np.linalg.svd(centered.T,full_matrices=False);scores=u*d;explained=d*d/(d*d).sum()
    corr=np.corrcoef(arr.T);out(pd.DataFrame(corr,columns=m.columns).assign(sample=m.columns).set_index('sample').reset_index(),f'correlation_{method}.csv')
    variance.append({'method':method,'pc1':float(explained[0]),'pc2':float(explained[1]),'n_probes':len(m)})
    for j,sample in enumerate(m.columns):
        qc.append({'sample':sample,'method':method,'group':['control','sesamin','episesamin','sesamolin'][j//5],'raw_min':m.iloc[:,j].min(),'raw_max':m.iloc[:,j].max(),'log2_q25':np.quantile(arr[:,j],.25),'log2_median':np.median(arr[:,j]),'log2_q75':np.quantile(arr[:,j],.75),'log2_mean':arr[:,j].mean(),'log2_sd':arr[:,j].std(ddof=1)})
        pca.append({'sample':sample,'method':method,'group':['control','sesamin','episesamin','sesamolin'][j//5],'PC1':scores[j,0],'PC2':scores[j,1],'primary_comparison':j<10})
out(pd.DataFrame(qc),'all20_array_quality.csv');out(pd.DataFrame(pca),'all20_array_pca.csv');out(pd.DataFrame(variance),'pca_variance.csv')
human=candidates.set_index('gene_id');mouse=pd.read_csv(N/'publication_package_14/data/mouse_candidates_all26.csv',keep_default_na=False)
cross=[]
for _,r in res.loc[res.method=='deposited_log2'].iterrows():
    h=human.loc[r.human_id];mr=mouse.loc[(mouse.original_human_id==r.human_id)&(mouse.primary_comparison_eligible.astype(str).str.lower()=='true')]
    hv=float(h.log2_sesamin_vs_OA);rv=r.effect_log2_signal
    mv=float(mr.iloc[0].log2_sess_vs_HFD_pc0p1) if len(mr)==1 else np.nan
    q=res.loc[(res.human_id==r.human_id)&(res.method=='quantile_sensitivity')].iloc[0]
    cross.append({'human_id':r.human_id,'human_symbol':r.human_symbol,'rat_mapping_status':r.status,'human_intervention_sign':int(np.sign(hv)),'mouse_intervention_sign':int(np.sign(mv)) if np.isfinite(mv) else np.nan,'rat_primary_effect':rv,'rat_q_bh_fixed25':r.q_bh_fixed25,'rat_quantile_effect':q.effect_log2_signal,'rat_quantile_q_bh_fixed25':q.q_bh_fixed25,'rat_primary_sign':int(np.sign(rv)) if np.isfinite(rv) else np.nan,'rat_human_direction_agrees':bool(np.sign(rv)==np.sign(hv)) if np.isfinite(rv) else None,'rat_mouse_direction_agrees':bool(np.sign(rv)==np.sign(mv)) if np.isfinite(rv) and np.isfinite(mv) else None,'rat_normalization_direction_stable':bool(np.sign(rv)==np.sign(q.effect_log2_signal)) if np.isfinite(rv) else None})
out(pd.DataFrame(cross),'cross_model_all25.csv')
assert all(sha(p['path'])==p['sha256'] for p in fingerprints)
summary={'status':'EXECUTED_EXPLORATORY_EXTERNAL_ANALYSIS','executed_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(__file__),'inputs':fingerprints,'protocol_sha256':sha(S/'protocol.md'),'environment':{'python':sys.executable,'python_version':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'platform':platform.platform()},'dimensions':{'probes':len(m),'arrays':20,'primary_control_animals':5,'primary_sesamin_animals':5,'candidate_family':25,'mapped_candidates':int((mapdf.status=='mapped').sum())},'mapping_status_counts':mapdf.status.value_counts().to_dict(),'primary_fdr05_symbols':res.loc[(res.method=='deposited_log2')&res.fdr05,'human_symbol'].tolist(),'sensitivity_fdr05_symbols':res.loc[(res.method=='quantile_sensitivity')&res.fdr05,'human_symbol'].tolist(),'limits':['Original study omitted 2 of 7 rats per group by sacrifice body weight; all deposited arrays retained.','Healthy rat whole liver, dietary 15-day exposure, no fatty-acid injury model or sesamin rescue test.','Processed globally scaled signals; no CEL reprocessing or array detection calls.','Welch intervals and tests rely on a small-sample log-scale approximation.','Probe annotations are deposited 2016-era annotations; no sequence-level remapping.','Differences between genes/species/platforms are not pooled effect sizes.','Mapping and hypotheses fixed before candidate effect calculation.'],'output_files':[{'file':p.relative_to(S).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(O.glob('*.csv'))]}
write(S/'analysis_summary.json',summary)
print(json.dumps({k:summary[k] for k in ['status','dimensions','mapping_status_counts','primary_fdr05_symbols','sensitivity_fdr05_symbols']},ensure_ascii=False))
print(res.loc[res.human_symbol.isin(['VNN1','AREG','JAG2','PALLD','SOCS2']),['human_symbol','method','status','probe_n','effect_log2_signal','ci95_low','ci95_high','p_value','q_bh_fixed25']].to_string(index=False))
