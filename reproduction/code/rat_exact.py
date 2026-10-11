"""Post hoc, all-candidate sensitivity to finite-sample inference assumptions."""
from pathlib import Path
from itertools import combinations
from math import comb
from datetime import datetime,timezone
import hashlib,json,sys
import numpy as np
import pandas as pd

S=Path(sys.argv[1]).resolve();O=S/'outputs'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
values=O/'candidate_sample_values.csv';results=O/'candidate_results_all25_two_normalizations.csv'
fingerprints={str(p):sha(p) for p in [values,results,S/'posthoc_sensitivity_plan.md']}
v=pd.read_csv(values);r=pd.read_csv(results);alloc=np.array(list(combinations(range(10),5)));assert len(alloc)==comb(10,5)==252
output=[]
for _,row in r.iterrows():
    base={'human_id':row.human_id,'human_symbol':row.human_symbol,'method':row.method,'mapping_status':row.status,'welch_p':row.p_value,'welch_q_bh25':row.q_bh_fixed25,'permutation_p':np.nan,'exceedances':np.nan,'allocations':252,'complete_group_separation':None}
    if row.status=='mapped':
        x=v.loc[(v.human_id==row.human_id)&(v.method==row.method)&v.primary_comparison].sort_values('sample').aggregate_log2_signal.to_numpy()
        assert len(x)==10
        observed=x[5:].mean()-x[:5].mean();summed=x[alloc].sum(axis=1);d=summed/5-(x.sum()-summed)/5
        k=int(np.sum(np.abs(d)>=abs(observed)-1e-12))
        assert k%2==0 and k>=2
        base.update({'observed_mean_difference':observed,'permutation_p':k/252,'exceedances':k,'complete_group_separation':bool(x[:5].max()<x[5:].min() or x[5:].max()<x[:5].min())})
    output.append(base)
df=pd.DataFrame(output)
for method in df.method.unique():
    idx=df.index[df.method==method];p=df.loc[idx,'permutation_p'].fillna(1).to_numpy();order=np.argsort(p);q=np.minimum.accumulate((p[order]*25/np.arange(1,26))[::-1])[::-1];z=np.empty_like(q);z[order]=np.clip(q,0,1)
    df.loc[idx,'permutation_q_bh25']=z;df.loc[idx[df.loc[idx,'permutation_p'].isna()],'permutation_q_bh25']=np.nan
df.to_csv(O/'posthoc_exact_sensitivity_all25.csv',index=False,na_rep='NA',float_format='%.12g')
assert all(sha(p)==h for p,h in fingerprints.items())
summary={'status':'POSTHOC_SENSITIVITY_EXECUTED','created_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(__file__),'input_sha256':fingerprints,'output_sha256':sha(O/'posthoc_exact_sensitivity_all25.csv'),'candidate_family':25,'allocations':252,'new_primary_analysis':False,'fdr05_symbols':df.loc[df.permutation_q_bh25<.05,['human_symbol','method']].to_dict('records'),'scope':'Exchangeability reference, not a verified randomized experiment; selected animals only.'}
(S/'exact_sensitivity_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(df.loc[df.human_symbol=='VNN1'].to_string(index=False))
