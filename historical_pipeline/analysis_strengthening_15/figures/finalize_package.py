"""Freeze figure captions, independent numeric cross-checks, and file checksums."""
from pathlib import Path
import hashlib,json,re,sys
import numpy as np
import pandas as pd
from scipy.stats import t

HERE=Path(__file__).resolve().parent
STEMS=['fig15_rat_external','fig15_robustness','fig15_cell_context','fig15_global_pca']
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(name,value):(HERE/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def csv(key):return pd.read_csv(HERE/'source_tables'/(key+'.csv'))

captions={}
for chunk in (HERE/'figlegends.md').read_text(encoding='utf-8').split('\n## ')[1:]:
    key,body=chunk.split('\n',1);captions[key.strip()]=body.strip()
assert list(captions)==STEMS
write('figlegends.json',captions)

checks={}
rat=csv('rat_results');samples=csv('rat_samples')
assert rat.loc[rat.status.eq('mapped'),'probe_n'].eq(1).all()
max_mean_error=max_ci_error=0.
for row in rat[rat.status.eq('mapped')].itertuples():
    sample=samples[(samples.human_symbol==row.human_symbol)&(samples.method==row.method)]
    c=sample.loc[sample.group.eq('control'),'aggregate_log2_signal'].to_numpy()
    s=sample.loc[sample.group.eq('sesamin'),'aggregate_log2_signal'].to_numpy()
    assert len(c)==len(s)==5
    effect=s.mean()-c.mean();v1=s.var(ddof=1)/5;v0=c.var(ddof=1)/5
    df=(v1+v0)**2/(v1*v1/4+v0*v0/4)
    half=t.ppf(.975,df)*np.sqrt(v1+v0)
    max_mean_error=max(max_mean_error,abs(effect-row.effect_log2_signal))
    max_ci_error=max(max_ci_error,abs(effect-half-row.ci95_low),abs(effect+half-row.ci95_high))
assert max_mean_error<1e-9 and max_ci_error<1e-8
checks['rat']={'mapped_rows_checked_across_two_normalizations':30,'animals_per_group':5,
 'one_eligible_probe_for_each_mapped_candidate':True,
 'max_effect_reconstruction_absolute_error':max_mean_error,'max_CI_reconstruction_absolute_error':max_ci_error,
 'purpose':'Arithmetic verification of the archived Welch CI; no new inferential result generated'}
exact=csv('rat_exact');primary_exact=exact[exact.method.eq('deposited_log2')]
assert len(primary_exact)==25 and not (primary_exact.permutation_q_bh25<.05).any()
vnn1=primary_exact[primary_exact.human_symbol.eq('VNN1')].iloc[0]
assert vnn1.allocations==252 and vnn1.exceedances==2
checks['rat'].update({'VNN1_exact_p':float(vnn1.permutation_p),'VNN1_exact_q_BH25':float(vnn1.permutation_q_bh25),
 'no_exact_FDR_below_005':True})

consensus=pd.read_csv(HERE/'fig15_concordance_source.csv')
assert len(consensus)==15 and (consensus.configuration_pairs==4096).all()
count=consensus.both_model_and_treatment_concordant_n
assert ((count>=0)&(count<=4096)).all()
assert np.allclose(count/4096,consensus.both_concordant_configuration_share,atol=1e-12,rtol=0)
ranges=pd.read_csv(HERE/'fig15_mouse_range_source.csv')
assert ranges.configurations.eq(64).all()
assert (ranges.min_log2_ratio<=ranges.full_data_log2_ratio+1e-12).all()
assert (ranges.max_log2_ratio>=ranges.full_data_log2_ratio-1e-12).all()
checks['robustness']={'retained_candidates':15,'range_contains_full_data_estimate':True,
 'integer_numerators_and_4096_denominator_verified':True,'no_CI_interpretation':True}

local=pd.read_csv(HERE/'fig15_AREG_library_source.csv')
assert len(local)==4
for column in ['cholangiocyte_detection_pct','hepatocyte_detection_pct']:
    assert local[column].between(0,100).all()
assert (local.cholangiocyte_detection_pct>local.hepatocyte_detection_pct).all()
decomp=pd.read_csv(HERE/'fig15_Kitagawa_spotlight_source.csv')
res=(decomp.composition_total+decomp.within_category_total-decomp.total_PA_minus_control).abs().max()
assert res<1e-12
checks['single_cell']={'AREG_cholangiocyte_detection_higher_in_each_library':True,
 'library_n':4,'Kitagawa_max_identity_error':float(res),'gene_n':4,'no_new_P_values_or_CI':True}

pca={}
for key,n in [('human_pca',9),('mouse_pca',9),('rat_pca',20)]:
    frame=csv(key)
    if key=='rat_pca':frame=frame[frame.method.eq('deposited_log2')]
    assert len(frame)==n and np.isfinite(frame[['PC1','PC2']]).all().all()
    assert abs(frame[['PC1','PC2']].mean()).max()<1e-8
    pca[key]={'sample_n':n,'source_group_counts':frame.groupby('group').size().to_dict(),
      'max_PC_center_residual':float(abs(frame[['PC1','PC2']].mean()).max())}
checks['pca']=pca
sources=json.loads((HERE/'source_data.json').read_text(encoding='utf-8'))
assert all(sha(Path(item['source']))==item['sha256']==sha(HERE/item['copied_file']) for item in sources['inputs'])
checks['input_source_and_copy_checksums_unchanged']=True
checks['status']='PASS'
checks['caption_word_counts']={key:len(re.findall(r'\S+',value)) for key,value in captions.items()}
write('independent_numerical_checks.json',checks)
print(json.dumps({'status':'PASS','caption_word_counts':checks['caption_word_counts'],'max_rat_CI_error':max_ci_error}))

if '--freeze' in sys.argv:
    visual=json.loads((HERE/'visual_review.json').read_text(encoding='utf-8'))
    assert visual['status']=='PASS'
    for item in visual['figures']:
        assert sha(HERE/item['png'])==item['png_sha256']
        assert sha(HERE/item['pdf'])==item['pdf_sha256']
    files=[]
    for path in sorted(HERE.rglob('*')):
        if path.is_file() and path.name!='artifact_manifest.json':
            files.append({'path':str(path.relative_to(HERE)).replace('\\','/'),'bytes':path.stat().st_size,'sha256':sha(path)})
    write('artifact_manifest.json',{'scope':'Stage15 four scientific figures only','files':files,'source_tables_read_only':True})
