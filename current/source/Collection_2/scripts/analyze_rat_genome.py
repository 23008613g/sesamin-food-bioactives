"""Execute the frozen exploratory GSE42693 whole-genome protocol.

Original data are read-only. This script writes only inside its stage directory.
No HepG2 candidates are used to define rat genes or biological-process sets.
"""
from __future__ import annotations
from pathlib import Path
from datetime import datetime, timezone
from itertools import combinations
from collections import defaultdict
import csv, gzip, hashlib, io, json, platform, sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
PROJECT = Path('C:/Users/AI/OneDrive/Desktop/公司文章/research_planning/hepatic_lipotoxicity_2026-10-07')
NEXT = PROJECT/'next_steps'
OUT = ROOT/'results'
INPUTS = {
    'rat_matrix': NEXT/'evidence_closure_06/intervention_access/sources/GSE42693_series_matrix.txt.gz',
    'platform_annotation': NEXT/'analysis_strengthening_15/rat_validation/sources/GPL1355.annot.gz',
    'go_ontology': ROOT/'sources/go/go-basic.obo',
    'protocol': ROOT/'PROTOCOL.txt',
    'protocol_addendum_1': ROOT/'PROTOCOL_ADDENDUM_1.txt',
}
SAMPLES = [f'GSM{i}' for i in range(1048036,1048056)]
VARIANTS = ['deposited','qn10','qn20']

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')

def write_csv(name, frame):
    frame.to_csv(OUT/name,index=False,na_rep='NA',float_format='%.15g',compression='infer')

def read_table(path, kind):
    text = gzip.open(path,'rt',encoding='utf8').read()
    data = text.split(f'!{kind}_table_begin\n',1)[1].split(f'!{kind}_table_end',1)[0]
    return text,pd.read_csv(io.StringIO(data),sep='\t',dtype=str,keep_default_na=False)

def split_ids(value):
    return sorted(set(x.strip() for x in value.split('///') if x.strip() and x.strip()!='---'))

def quantile_normalize(array):
    target = np.sort(array,axis=0).mean(axis=1)
    output = np.empty_like(array)
    for j in range(array.shape[1]):
        _, inverse, counts = np.unique(array[:,j],return_inverse=True,return_counts=True)
        starts = np.r_[0,np.cumsum(counts)[:-1]]
        output[:,j] = (np.add.reduceat(target,starts)/counts)[inverse]
    return output

def fdr(pvalues, by=False):
    p = np.asarray(pvalues,dtype=float)
    p = np.where(np.isfinite(p),p,1.)
    assert ((p>=0)&(p<=1)).all()
    n = len(p)
    order = np.argsort(p,kind='stable')
    factor = np.sum(1/np.arange(1,n+1)) if by else 1.
    adjusted = np.minimum.accumulate((p[order]*n*factor/np.arange(1,n+1))[::-1])[::-1]
    answer = np.empty(n)
    answer[order] = np.clip(adjusted,0,1)
    return answer

def allocation_matrix():
    allocations = np.array(list(combinations(range(10),5)),dtype=int)
    assert allocations.shape==(252,5)
    weights = np.full((252,10),-1/5.)
    for k,indices in enumerate(allocations): weights[k,indices] = 1/5.
    assert len(np.unique(weights,axis=0))==252
    assert np.isclose(weights.sum(axis=1),0).all()
    return allocations,weights

def infer_rows(values, weights):
    """Values rows are endpoints; ten columns are biological animals."""
    assert values.ndim==2 and values.shape[1]==10 and np.isfinite(values).all()
    control,treat = values[:,:5],values[:,5:]
    ma,mb = control.mean(axis=1),treat.mean(axis=1)
    va,vb = control.var(axis=1,ddof=1),treat.var(axis=1,ddof=1)
    diff = mb-ma
    se2 = va/5+vb/5
    estimable = se2>0
    se = np.sqrt(se2)
    df = np.full(len(values),np.nan)
    df[estimable] = se2[estimable]**2/((va[estimable]/5)**2/4+(vb[estimable]/5)**2/4)
    t = np.full(len(values),np.nan)
    t[estimable] = diff[estimable]/se[estimable]
    p = np.full(len(values),np.nan)
    p[estimable] = 2*stats.t.sf(np.abs(t[estimable]),df[estimable])
    half = stats.t.ppf(.975,df)*se
    permutations = values@weights.T
    tol = 1e-12*np.maximum(1,np.abs(diff))
    exact_p = (np.abs(permutations)>=np.abs(diff[:,None])-tol[:,None]).mean(axis=1)
    assert (exact_p>=2/252-1e-15).all() and (exact_p<=1).all()
    return pd.DataFrame({
        'mean_control':ma,'sd_control':np.sqrt(va),'n_control':5,
        'mean_sesamin':mb,'sd_sesamin':np.sqrt(vb),'n_sesamin':5,
        'difference':diff,'welch_se':se,'welch_t':t,'welch_df':df,
        'welch_ci_low':diff-half,'welch_ci_high':diff+half,
        'welch_p':p,'welch_estimable':estimable,
        'welch_bh':fdr(p),'welch_by':fdr(p,True),
        'exact_p':exact_p,'exact_bh':fdr(exact_p),'exact_by':fdr(exact_p,True),
        'median_control':np.median(control,axis=1),'median_sesamin':np.median(treat,axis=1),
        'min_control':control.min(axis=1),'max_control':control.max(axis=1),
        'min_sesamin':treat.min(axis=1),'max_sesamin':treat.max(axis=1),
        'mean_sample_signal':values.mean(axis=1),
    })

def parse_ontology(path):
    terms = {}
    for block in path.read_text(encoding='utf8').split('\n[Term]\n')[1:]:
        block = block.split('\n[Typedef]',1)[0]
        fields = defaultdict(list)
        for line in block.splitlines():
            if ': ' in line:
                key,value = line.split(': ',1)
                fields[key].append(value)
        if 'id' not in fields: continue
        identity = fields['id'][0]
        terms[identity] = {'name':fields['name'][0],
            'namespace':fields.get('namespace',[''])[0],
            'obsolete':fields.get('is_obsolete',['false'])[0]=='true',
            'alt_ids':fields.get('alt_id',[])}
    aliases = defaultdict(set)
    for identity,term in terms.items():
        aliases[identity].add(identity)
        for alt in term['alt_ids']: aliases[alt].add(identity)
    return terms,aliases

def build_go(ann, gene_ids):
    terms,aliases = parse_ontology(INPUTS['go_ontology'])
    raw_members = defaultdict(set)
    for row in ann.itertuples(index=False):
        if row.rat_entrez not in gene_ids: continue
        for term in split_ids(row.go_process_ids): raw_members[term].add(row.rat_entrez)
    members = defaultdict(set)
    audits = []
    for source_id,genes in sorted(raw_members.items()):
        matched = aliases.get(source_id,set())
        canonical,name,status = '', '', ''
        if not matched: status='unresolved_id'
        elif len(matched)!=1: status='ambiguous_alias'
        else:
            canonical = next(iter(matched));term = terms[canonical];name=term['name']
            if term['obsolete']: status='obsolete_excluded'
            elif term['namespace']!='biological_process': status='not_biological_process'
            else:
                status='retained_id' if source_id==canonical else 'unique_alt_id_resolved'
                members[canonical].update(genes)
        audits.append({'source_go_id':source_id,'canonical_go_id':canonical,'go_name':name,
                       'status':status,'n_eligible_genes_for_source_id':len(genes)})
    write_csv('go_term_status_audit.csv',pd.DataFrame(audits))
    universe,rows = [],[]
    for term,genes in sorted(members.items()):
        selected = 15<=len(genes)<=500
        signature = hashlib.sha256(','.join(sorted(genes,key=int)).encode()).hexdigest()
        universe.append({'go_id':term,'go_name':terms[term]['name'],'n_genes':len(genes),
                         'tested':selected,'member_set_sha256':signature})
        if selected:
            rows.extend({'go_id':term,'go_name':terms[term]['name'],'rat_entrez':g}
                        for g in sorted(genes,key=int))
    definition = pd.DataFrame(universe)
    definition['identical_membership_go_id_count'] = definition.groupby('member_set_sha256')['go_id'].transform('size')
    write_csv('go_gene_set_universe.csv',definition)
    membership = pd.DataFrame(rows)
    write_csv('go_membership.csv.gz',membership)
    return definition.loc[definition.tested].reset_index(drop=True), membership

def basic_tests():
    p=np.array([.001,.04,.02,1.,.6,.9,.5,.05])
    assert np.allclose(fdr(p),stats.false_discovery_control(p,method='bh'))
    assert np.allclose(fdr(p,True),stats.false_discovery_control(p,method='by'))
    tied=np.array([[1.,3.],[1.,1.],[3.,2.]])
    assert np.array_equal(quantile_normalize(tied),np.array([[1.25,3.],[1.25,1.],[3.,1.5]]))
    _,w=allocation_matrix()
    fixture=np.array([[0.,1.,2.,3.,4.,7.,8.,9.,10.,11.], [2.]*10])
    result=infer_rows(fixture,w)
    assert np.isclose(result.loc[0,'exact_p'],2/252)
    assert np.isclose(result.loc[0,'welch_p'],stats.ttest_ind(fixture[0,5:],fixture[0,:5],equal_var=False).pvalue)
    assert result.loc[1,'exact_p']==1 and result.loc[1,'welch_bh']==1
    return ['BH and BY match scipy','tie-aware QN fixture','all allocations',
            'Welch matches scipy fixture','constant endpoint retained conservatively']

def main():
    started=datetime.now(timezone.utc).isoformat()
    OUT.mkdir(exist_ok=True);(ROOT/'review').mkdir(exist_ok=True)
    assert not (ROOT/'analysis_summary.json').exists(), 'Refusing to silently rerun completed analysis'
    freeze=json.loads((ROOT/'protocol_freeze.json').read_text(encoding='utf8'))
    assert sha(INPUTS['protocol'])==freeze['protocol_sha256']
    input_manifest={role:{'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size,'read_only':True}
                    for role,path in INPUTS.items()}
    write_json(ROOT/'input_manifest.json',input_manifest)
    software_tests=basic_tests()
    text,matrix=read_table(INPUTS['rat_matrix'],'series_matrix')
    assert matrix.shape==(31099,21) and matrix.columns[1:].tolist()==SAMPLES
    assert matrix.ID_REF.is_unique
    raw=matrix[SAMPLES].to_numpy(dtype=float)
    assert np.isfinite(raw).all() and (raw>=0).all()
    metadata=defaultdict(list)
    for line in text.splitlines():
        if line.startswith('!Sample_'):
            row=next(csv.reader([line],delimiter='\t'));metadata[row[0]].append(row[1:])
    assert 'global scaling' in metadata['!Sample_data_processing'][0][0]
    sample_info=[]
    for j,sample in enumerate(SAMPLES):
        sample_info.append({'sample':sample,'group':['control','sesamin','episesamin','sesamolin'][j//5],
            'source_title':metadata['!Sample_title'][0][j],'primary_comparison':j<10,
            'log_probe_q01':float(np.quantile(np.log2(raw[:,j]+1),.01)),
            'log_probe_median':float(np.median(np.log2(raw[:,j]+1))),
            'log_probe_q99':float(np.quantile(np.log2(raw[:,j]+1),.99)),
            'animal_selection':'original seven minus highest and lowest terminal body weight',
            'current_exclusion':'none'})
    write_csv('sample_metadata_and_qc.csv',pd.DataFrame(sample_info))
    _,annotation=read_table(INPUTS['platform_annotation'],'platform')
    assert annotation.ID.is_unique and set(annotation.ID)==set(matrix.ID_REF)
    annotation=annotation.set_index('ID').loc[matrix.ID_REF].reset_index()
    entrez_sets=annotation['Gene ID'].map(split_ids)
    annotation['rat_entrez']=[ids[0] if len(ids)==1 and ids[0].isdigit() else '' for ids in entrez_sets]
    annotation['n_entrez_ids']=entrez_sets.map(len)
    write_csv('probe_mapping_audit.csv.gz',annotation[['ID','Gene ID','Gene symbol','rat_entrez','n_entrez_ids']])
    eligible=annotation.rat_entrez!=''
    group_indices={gene:np.asarray(indices,dtype=int) for gene,indices in
                   annotation.loc[eligible].groupby('rat_entrez',sort=False).groups.items()}
    gene_ids=sorted(group_indices,key=int)
    log_raw=np.log2(raw+1)
    probe_values={'deposited':log_raw[:,:10],
                  'qn10':quantile_normalize(log_raw[:,:10]),
                  'qn20':quantile_normalize(log_raw)[:,:10]}
    values={variant:np.array([np.median(probes[group_indices[g]],axis=0) for g in gene_ids])
            for variant,probes in probe_values.items()}
    retained=np.isfinite(values['deposited']).all(axis=1)&(values['deposited'].var(axis=1,ddof=1)>0)
    write_csv('rat_gene_universe.csv',pd.DataFrame({'rat_entrez':gene_ids,'included':retained,
              'reason':np.where(retained,'unique Entrez probes; finite positive total variance','zero/nonfinite total variance')}))
    gene_ids=np.array(gene_ids)[retained].tolist()
    values={variant:arr[retained] for variant,arr in values.items()}
    annotation_rows=[]
    for gene in gene_ids:
        probes=annotation.iloc[group_indices[gene]]
        symbols=sorted(set(s for value in probes['Gene symbol'] for s in split_ids(value)))
        annotation_rows.append({'rat_entrez':gene,'gene_symbol':';'.join(symbols),
            'n_symbols':len(symbols),'n_probes':len(probes),'probe_ids':';'.join(probes.ID),
            'gene_titles':';'.join(sorted(set(probes['Gene title'])-{''}))})
    gene_annotation=pd.DataFrame(annotation_rows)
    write_csv('rat_gene_annotation.csv',gene_annotation)
    allocations,weights=allocation_matrix()
    write_csv('exact_allocations.csv',pd.DataFrame(allocations,columns=[f'treatment_index_{i+1}' for i in range(5)]))
    all_summary={}
    for variant,arr in values.items():
        write_csv(f'rat_gene_expression_{variant}.csv.gz',pd.DataFrame(arr,index=pd.Index(gene_ids,name='rat_entrez'),columns=SAMPLES[:10]).reset_index())
        result=pd.concat([gene_annotation[['rat_entrez','gene_symbol','n_probes']],infer_rows(arr,weights)],axis=1)
        result['variant']=variant
        write_csv(f'rat_gene_results_{variant}.csv.gz',result)
        all_summary[variant]={'n_genes':len(result),'welch_bh_lt_005':int((result.welch_bh<.05).sum()),
            'welch_by_lt_005':int((result.welch_by<.05).sum()),'exact_bh_lt_005':int((result.exact_bh<.05).sum()),
            'exact_by_lt_005':int((result.exact_by<.05).sum()),'nonestimable_welch':int((~result.welch_estimable).sum())}
        print(json.dumps({'phase':'gene_results','variant':variant,**all_summary[variant]}),flush=True)
    go_ann=annotation.rename(columns={'GO:Process ID':'go_process_ids'})
    go_defs,membership=build_go(go_ann,set(gene_ids))
    row_index={gene:i for i,gene in enumerate(gene_ids)}
    go_indices={term:np.array([row_index[g] for g in group.rat_entrez])
                for term,group in membership.groupby('go_id')}
    go_summaries={}
    for variant,arr in values.items():
        sd=arr.std(axis=1,ddof=1)
        z=np.divide(arr-arr.mean(axis=1,keepdims=True),sd[:,None],out=np.zeros_like(arr),where=sd[:,None]>0)
        scores=np.array([z[go_indices[term]].mean(axis=0) for term in go_defs.go_id])
        result=pd.concat([go_defs,infer_rows(scores,weights)],axis=1)
        result['variant']=variant
        result['lipid_name_subset']=result.go_name.str.contains('lipid|fatty acid|cholesterol|sterol|triglyceride|triacylglycerol',case=False,regex=True)
        write_csv(f'go_results_{variant}.csv',result)
        write_csv(f'go_sample_scores_{variant}.csv.gz',pd.DataFrame(scores,index=pd.Index(go_defs.go_id,name='go_id'),columns=SAMPLES[:10]).reset_index())
        write_csv(f'go_lipid_terms_{variant}.csv',result.loc[result.lipid_name_subset])
        go_summaries[variant]={'n_gene_sets':len(result),'welch_bh_lt_005':int((result.welch_bh<.05).sum()),
            'exact_bh_lt_005':int((result.exact_bh<.05).sum()),'exact_by_lt_005':int((result.exact_by<.05).sum()),
            'n_lipid_name_sets':int(result.lipid_name_subset.sum())}
        print(json.dumps({'phase':'go_results','variant':variant,**go_summaries[variant]}),flush=True)
    # Label-independent PCA of primary arrays; all twenty arrays retained for QC.
    for scope,arr in [('primary10',values['deposited']),('all20_probes',log_raw)]:
        centred=arr.T-arr.T.mean(axis=0)
        u,s,_=np.linalg.svd(centred,full_matrices=False)
        pcs=u[:,:3]*s[:3]
        evr=s**2/np.sum(s**2)
        frame=pd.DataFrame(pcs,columns=['PC1','PC2','PC3']);frame.insert(0,'sample',SAMPLES[:len(frame)])
        write_csv(f'pca_{scope}.csv',frame)
        write_json(OUT/f'pca_{scope}_variance.json',{'explained_variance_ratio':evr.tolist(),'gene_scaled':False})
    # Pre-existing examples are extracted only after the whole-genome calculations.
    all_deposited=pd.read_csv(OUT/'rat_gene_results_deposited.csv.gz',dtype={'rat_entrez':str})
    examples=all_deposited.loc[all_deposited.gene_symbol.str.split(';').map(lambda x:any(s.lower() in ['vnn1','areg'] for s in x))]
    write_csv('preexisting_examples_not_new_hits.csv',examples)
    for path_meta in input_manifest.values(): assert sha(path_meta['path'])==path_meta['sha256']
    summary={'started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),
        'protocol_sha256':freeze['protocol_sha256'],'script_sha256':sha(__file__),
        'python':sys.executable,'python_version':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,
        'n_probes':len(matrix),'eligible_single_entrez_probes':int(eligible.sum()),
        'n_genes_before_variance_filter':len(retained),'n_genes':len(gene_ids),
        'excluded_zero_nonfinite_genes':int((~retained).sum()),'gene_results':all_summary,'go_results':go_summaries,
        'software_checks_passed':software_tests,'input_hashes_unchanged':True,
        'analysis_status':'executed; pending independent computational and scientific review',
        'scientific_limitations':['post-hoc exploratory analysis','five selected healthy rats per group',
            'exchangeability unverified','historical gene-term associations','no normal-human sesamin validation']}
    write_json(ROOT/'analysis_summary.json',summary)
    write_json(ROOT/'output_manifest.json',[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size}
               for p in sorted(OUT.rglob('*')) if p.is_file()])
    print(json.dumps(summary,ensure_ascii=False),flush=True)

if __name__=='__main__':
    main()
