"""Portable numerical reproduction from supplied cell extracts, not raw extraction."""
from pathlib import Path
from datetime import datetime, timezone
import json, hashlib, platform
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent
org=pd.read_csv(HERE/'organoid_areg_cells.csv.gz')
adult=pd.read_csv(HERE/'adult_areg_cells.csv.gz')
checks=[]
def check(name,value):
    checks.append({'name':name,'passed':bool(value)})
    assert value,name
check('organoid cells 24034',len(org)==24034)
check('all measured features',org.feature_count_in_T.eq(33538).all())
check('unique library barcode',not org.duplicated(['library_accession','barcode']).any())
check('raw counts to transform',np.allclose(np.log1p(1e4*org.AREG_raw_UMI/org.all_feature_raw_UMI),org.AREG_log1p_CP10k,atol=1e-12))
populations={
 'Whole organoid':np.ones(len(org),bool),
 'Hepatocytes':org.source_broad_label.eq('Hepatocytes'),
 'Broad Cholangiocytes':org.source_broad_label.eq('Cholangiocytes'),
 'Cholangiocytes + Ductal cells (excluding SMC)':org.source_subtype_label.isin(['Cholangiocytes','Ductal cells'])}
original=pd.read_csv(HERE/'areg_population_by_library.csv')
names=original.population.unique().tolist()
rows=[]
for (our_name,mask),actual_name in zip(populations.items(),names):
    # Exact declared population order is checked by source membership below.
    subset=org.loc[mask]
    for lib,d in subset.groupby('library_accession'):
        saved=original.loc[(original.population==actual_name)&(original.library_accession==lib)].iloc[0]
        check(f'{actual_name}/{lib}/cells',len(d)==saved.n_cells)
        check(f'{actual_name}/{lib}/mean',np.isclose(d.AREG_log1p_CP10k.mean(),saved.mean_log1p_CP10k,atol=1e-12))
        check(f'{actual_name}/{lib}/nonzero',int(d.AREG_raw_UMI.gt(0).sum())==saved.n_raw_nonzero)
        rows.append({'population':actual_name,'library':lib,'condition':d.analysis_condition.iloc[0],'n':len(d),'mean_log1p_CP10k':d.AREG_log1p_CP10k.mean()})
check('adult cells 8444',len(adult)==8444)
saved=pd.read_csv(HERE/'areg_adult_donor_context.csv')
for row in saved.itertuples():
    d=adult.loc[adult.donor.eq(row.donor)]
    c=d.loc[d.analysis_broad_label.eq('Cholangiocytes'),'AREG_source_log2CPM']
    h=d.loc[d.analysis_broad_label.eq('Hepatocytes'),'AREG_source_log2CPM']
    check(row.donor+'/coverage',len(c)==row.n_cholangiocytes and len(h)==row.n_hepatocytes)
    check(row.donor+'/nonzero_delta',np.isclose(100*(c.gt(0).mean()-h.gt(0).mean()),row.delta_nonzero_percentage_points,atol=1e-12))
    check(row.donor+'/source_mean_delta',np.isclose(c.mean()-h.mean(),row.delta_mean_source_log2CPM,atol=1e-12))
out=HERE/'portable_check_results';out.mkdir(exist_ok=True)
pd.DataFrame(rows).to_csv(out/'organoid_library_recalculation.csv',index=False)
(out/'verification.json').write_text(json.dumps({'utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),'scope':'processed per-cell extracts; not a fresh raw-matrix or author-H5AD extraction','checks':checks,'all_passed':all(c['passed'] for c in checks)},indent=2)+'\n')
print(len(checks),'processed-input checks passed')
