"""Independent numerical audit. No production functions are imported.

Writes only review/independent_validation.json. Raw input files are read-only.
"""
from pathlib import Path
import gzip, io, json, hashlib, itertools, sys, platform
import numpy as np
import pandas as pd
from scipy import stats
import scipy

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
manifest = json.loads((ROOT/'input_manifest.json').read_text(encoding='utf8'))
checks = []

def check(name, condition, **details):
    record = {'check':name, 'passed':bool(condition), **details}
    checks.append(record)
    if not condition:
        raise AssertionError(record)

def close(name, observed, expected, atol=1e-10, rtol=1e-8):
    a,b = np.asarray(observed,dtype=float),np.asarray(expected,dtype=float)
    delta = float(np.nanmax(np.abs(a-b)))
    check(name, np.allclose(a,b,atol=atol,rtol=rtol,equal_nan=True),max_absolute_difference=delta)

def raw_table(path, begin, end):
    lines = gzip.open(path,'rt',encoding='utf8').read().splitlines()
    i,j = lines.index(begin), lines.index(end)
    return pd.read_csv(io.StringIO('\n'.join(lines[i+1:j])),sep='\t',dtype=str,keep_default_na=False),lines

raw, lines = raw_table(manifest['rat_matrix']['path'],'!series_matrix_table_begin','!series_matrix_table_end')
ann,_ = raw_table(manifest['platform_annotation']['path'],'!platform_table_begin','!platform_table_end')
sample = [f'GSM{i}' for i in range(1048036,1048056)]
check('raw accession order', raw.columns.tolist()==['ID_REF']+sample)
check('raw probe uniqueness', raw.ID_REF.is_unique and ann.ID.is_unique)
raw = raw.set_index('ID_REF')
ann = ann.set_index('ID')
titles = next(x for x in lines if x.startswith('!Sample_title\t')).split('\t')[1:]
check('raw group titles', all('control diet' in x for x in titles[:5]) and all('liver_sesamin_' in x for x in titles[5:10]))
ids = ann['Gene ID'].map(lambda x: sorted(set(y.strip() for y in x.split('///') if y.strip() not in ('','---'))))
ids = ids.map(lambda x:x[0] if len(x)==1 and x[0].isdigit() else None).dropna()
check('independent single-Entrez probe count',len(ids)==21221,count=len(ids))
alllog = np.log2(raw.loc[:,sample].astype(float)+1)

def independent_qn(frame):
    arr = frame.to_numpy()
    target = np.sort(arr,axis=0).mean(axis=1)
    result = np.zeros_like(arr)
    for j in range(arr.shape[1]):
        order = np.argsort(arr[:,j],kind='stable')
        sorted_values = arr[order,j]
        boundaries = np.r_[0,1+np.flatnonzero(np.diff(sorted_values)!=0),len(order)]
        for start,stop in zip(boundaries[:-1],boundaries[1:]):
            result[order[start:stop],j] = target[start:stop].mean()
    return pd.DataFrame(result,index=frame.index,columns=frame.columns)

probe_variants = {'deposited':alllog.iloc[:,:10],
                  'qn10':independent_qn(alllog.iloc[:,:10]),
                  'qn20':independent_qn(alllog).iloc[:,:10]}
gene_matrices = {}
gene_ids = sorted(set(ids),key=int)
for variant, frame in probe_variants.items():
    eligible = frame.loc[ids.index].copy()
    eligible['independent_entrez'] = ids
    aggregated = eligible.groupby('independent_entrez').median().loc[gene_ids]
    gene_matrices[variant] = aggregated
    saved = pd.read_csv(OUT/f'rat_gene_expression_{variant}.csv.gz',dtype={'rat_entrez':str}).set_index('rat_entrez')
    check(f'{variant}: gene universe/order', saved.index.tolist()==gene_ids and list(saved.columns)==sample[:10])
    close(f'{variant}: all probe aggregation and normalization', saved,aggregated)
check('independent no zero-variance exclusion', (gene_matrices['deposited'].var(axis=1,ddof=1)>0).all())

subsets = list(itertools.combinations(range(10),5))
check('252 allocations and complements',len(subsets)==252 and len(set(subsets))==252 and all(tuple(i for i in range(10) if i not in a) in subsets for a in subsets))

def recompute(frame):
    x = np.asarray(frame)
    t = stats.ttest_ind(x[:,5:],x[:,:5],axis=1,equal_var=False)
    d = x[:,5:].mean(axis=1)-x[:,:5].mean(axis=1)
    ci = t.confidence_interval(.95)
    exceed = np.zeros(len(x),dtype=int)
    threshold = np.abs(d)-1e-12*np.maximum(1,np.abs(d))
    for a in subsets:
        b = tuple(i for i in range(10) if i not in a)
        shuffled = x[:,a].mean(axis=1)-x[:,b].mean(axis=1)
        exceed += np.abs(shuffled)>=threshold
    p = t.pvalue
    p_clean = np.where(np.isfinite(p),p,1)
    ep = exceed/len(subsets)
    return {
        'difference':d, 'mean_control':x[:,:5].mean(axis=1), 'mean_sesamin':x[:,5:].mean(axis=1),
        'sd_control':x[:,:5].std(axis=1,ddof=1), 'sd_sesamin':x[:,5:].std(axis=1,ddof=1),
        'welch_t':t.statistic, 'welch_df':t.df,'welch_p':p,
        'welch_ci_low':ci.low,'welch_ci_high':ci.high,
        'welch_bh':stats.false_discovery_control(p_clean,method='bh'),
        'welch_by':stats.false_discovery_control(p_clean,method='by'),
        'exact_p':ep,'exact_bh':stats.false_discovery_control(ep,method='bh'),
        'exact_by':stats.false_discovery_control(ep,method='by'),
    }

summaries, selected = {}, {}
for variant, mat in gene_matrices.items():
    result = recompute(mat)
    saved = pd.read_csv(OUT/f'rat_gene_results_{variant}.csv.gz',dtype={'rat_entrez':str}).set_index('rat_entrez')
    for column, value in result.items():
        close(f'{variant}: all-gene {column}', saved[column], value)
    summaries[variant] = {'n_genes':len(mat),'welch_bh_lt_005':int((result['welch_bh']<.05).sum()),
                          'exact_bh_lt_005':int((result['exact_bh']<.05).sum()),
                          'min_exact_p':float(result['exact_p'].min())}
    if variant=='deposited':
        examples = saved[(saved.welch_bh<.05)|saved.gene_symbol.fillna('').str.fullmatch('Vnn1|Areg',case=False)]
        for eid,row in examples.iterrows():
            selected[eid] = {'symbol':row.gene_symbol,'expression':mat.loc[eid].tolist(),
                             **{k:float(row[k]) for k in ('difference','welch_ci_low','welch_ci_high','welch_p','welch_bh','exact_p','exact_bh')}}

# Reconstruct the historical platform annotation and ontology status independently.
obo = Path(manifest['go_ontology']['path']).read_text(encoding='utf8').splitlines()
terms, current, item = {},None,{}
def finish():
    if current=='Term' and 'id' in item:
        terms[item['id'][0]]=dict(item)
for line in obo+['[end]']:
    if line.startswith('[') and line.endswith(']'):
        finish(); current=line[1:-1]; item={}
    elif ': ' in line:
        key,value=line.split(': ',1)
        item.setdefault(key,[]).append(value)
aliases={}
for term,data in terms.items():
    for name in [term]+data.get('alt_id',[]):
        aliases.setdefault(name,set()).add(term)
members={}
for probe,eid in ids.items():
    for source in ann.loc[probe,'GO:Process ID'].split('///'):
        candidates=aliases.get(source.strip(),set())
        if len(candidates)!=1: continue
        term=next(iter(candidates)); data=terms[term]
        if data.get('namespace')!=['biological_process'] or data.get('is_obsolete')==['true']: continue
        members.setdefault(term,set()).add(eid)
members={term:genes for term,genes in members.items() if 15<=len(genes)<=500}
go_ids=sorted(members)
saved_members=pd.read_csv(OUT/'go_membership.csv.gz',dtype={'rat_entrez':str})
saved_dict={term:set(group.rat_entrez) for term,group in saved_members.groupby('go_id')}
check('all tested GO membership from raw annotation/OBO',members==saved_dict,n_go_sets=len(members))
go_summary, go_spots = {},{}
for variant,mat in gene_matrices.items():
    z=pd.DataFrame(stats.zscore(mat,axis=1,ddof=1),index=mat.index,columns=mat.columns)
    score=pd.DataFrame([z.loc[sorted(members[t],key=int)].mean(axis=0) for t in go_ids],index=go_ids)
    saved_score=pd.read_csv(OUT/f'go_sample_scores_{variant}.csv.gz').set_index('go_id')
    close(f'{variant}: all GO sample scores',saved_score.loc[go_ids],score)
    result=recompute(score)
    saved=pd.read_csv(OUT/f'go_results_{variant}.csv').set_index('go_id').loc[go_ids]
    for column,value in result.items():
        close(f'{variant}: all-GO {column}',saved[column],value)
    go_summary[variant]={'n_sets':len(go_ids),'welch_bh_lt_005':int((result['welch_bh']<.05).sum()),
                         'exact_bh_lt_005':int((result['exact_bh']<.05).sum())}
    if variant=='deposited':
        for term in ['GO:0006629','GO:0006631']:
            if term in saved.index:
                row=saved.loc[term]
                go_spots[term]={'name':row.go_name,'n_genes':int(row.n_genes),'scores':score.loc[term].tolist(),
                               **{k:float(row[k]) for k in ('difference','welch_p','welch_bh','exact_p','exact_bh')}}

for role, meta in manifest.items():
    h=hashlib.sha256(Path(meta['path']).read_bytes()).hexdigest()
    check(f'input unchanged: {role}',h==meta['sha256'])
report={'status':'passed','scope':'Independent full rat gene and GO numeric recomputation from read-only deposited matrix and platform/ontology inputs; no production statistical function imports.',
        'not_validated':['Original animal randomization','Unselected animals absent from archive','Raw CEL processing','Human mapping/cross-species script','limma empirical Bayes hyperparameter implementation'],
        'software':{'python':sys.executable,'version':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
        'checks_passed':len(checks),'checks':checks,'gene_summary':summaries,'go_summary':go_summary,
        'primary_gene_spots':selected,'primary_lipid_GO_spots':go_spots,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'claim_boundary':'These checks validate reproduction of the specified exploratory calculation; they do not establish biological truth, normal-human sesamin response, or absence of lipid-metabolic effects.'}
(ROOT/'review/independent_validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')
print(json.dumps({'status':report['status'],'checks_passed':len(checks),'gene_summary':summaries,'go_summary':go_summary,'primary_gene_spots':selected,'primary_lipid_GO_spots':go_spots},ensure_ascii=False,indent=2))
