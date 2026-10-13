"""Independent focused output audit, including CSV-module AREG re-read."""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',type=Path,default=Path('inputs'))
    ap.add_argument('--output-dir',type=Path,default=Path('results'))
    args=ap.parse_args()
    here=Path(__file__).resolve().parent
    inp,out=args.input_dir,args.output_dir
    checks={}
    def ck(k,v):
        checks[k]=bool(v)
        assert v,k
    raw={}
    with gzip.open(inp/'GSE115469_Data.csv.gz','rt',encoding='utf-8',newline='') as f:
        reader=csv.reader(f)
        cells=next(reader)[1:]
        for row in reader:
            if row[0]=='AREG':
                raw=dict(zip(cells,map(float,row[1:])))
                break
    ck('AREG_independently_reread_from_deposited_csv',len(raw)==8444)
    with gzip.open(inp/'GSE115469_CellClusterType.txt.gz','rt',encoding='utf-8',newline='') as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    paired=pd.read_csv(out/'paired_cholangiocyte_minus_hepatocyte_all25.csv')
    apairs=paired.loc[paired.symbol.eq('AREG')]
    scalar=[]
    max_error=0.
    for donor in [f'P{i}TLH' for i in range(1,6)]:
        observed=apairs.loc[apairs.donor.eq(donor)].iloc[0]
        entry={'donor':donor}
        for typ,key in [('Cholangiocytes','chol'),('Hepatocyte_','hep')]:
            values=[raw[r['CellName']] for r in rows if r['Sample']==donor and (r['CellType']==typ or (typ=='Hepatocyte_' and r['CellType'].startswith(typ)))]
            entry[key+'_n']=len(values)
            entry[key+'_nonzero']=sum(v>0 for v in values)
            entry[key+'_pct']=100*entry[key+'_nonzero']/len(values)
            entry[key+'_mean']=math.fsum(values)/len(values)
            max_error=max(max_error,abs(entry[key+'_pct']-observed[key+'_nonzero_pct']),abs(entry[key+'_mean']-observed[key+'_mean_source_log2CPM']))
        entry['delta_pct']=entry['chol_pct']-entry['hep_pct']
        entry['delta_mean']=entry['chol_mean']-entry['hep_mean']
        scalar.append(entry)
    ck('all10_source_AREG_aggregates_agree_with_saved_results',max_error<1e-9)
    ck('P1_negative_is_preserved',scalar[0]['chol_nonzero']==0 and scalar[0]['chol_n']==8 and scalar[0]['hep_nonzero']==1 and scalar[0]['hep_n']==96 and scalar[0]['delta_pct']<0 and scalar[0]['delta_mean']<0)
    ck('AREG_four_positive_one_negative',sum(r['delta_pct']>0 for r in scalar)==4 and sum(r['delta_mean']>0 for r in scalar)==4)
    ck('pooled_AREG_counts_13_of119_and26_of3501',sum(r['chol_nonzero'] for r in scalar)==13 and sum(r['chol_n'] for r in scalar)==119 and sum(r['hep_nonzero'] for r in scalar)==26 and sum(r['hep_n'] for r in scalar)==3501)
    for metric in ['pct','mean']:
        mean_delta=math.fsum(r['delta_'+metric] for r in scalar)/5
        difference_of_means=math.fsum(r['chol_'+metric] for r in scalar)/5-math.fsum(r['hep_'+metric] for r in scalar)/5
        ck('equal_donor_identity_'+metric,abs(mean_delta-difference_of_means)<1e-12)
    ck('coverage_sensitivity_exactly_P2_P3_P5',[r['donor'] for r in scalar if min(r['chol_n'],r['hep_n'])>=20]==['P2TLH','P3TLH','P5TLH'])
    protocol=(here/'protocol.md').read_text('utf-8')
    freeze=json.loads((here/'protocol_freeze.json').read_text('utf-8'))
    scale=json.loads((here/'source_scale_decision_freeze.json').read_text('utf-8'))
    ck('original_protocol_hash_unchanged',hashlib.sha256((here/'protocol.md').read_bytes()).hexdigest()==freeze['files']['protocol.md'])
    ck('sensitivity_was_prespecified','at least 20 cells in both types' in protocol)
    ck('source_scale_decision_hash_unchanged',hashlib.sha256((here/'source_scale_decision.md').read_bytes()).hexdigest()==scale['sha256'])
    first_execution=json.loads((out/'execution_record.json').read_text('utf-8'))
    ck('both_freezes_predate_execution',freeze['frozen_utc']<scale['frozen_utc']<first_execution['started_utc'])
    candidates=pd.read_csv(inp/'human_candidates_all25.csv')
    cov=pd.read_csv(out/'candidate_feature_coverage_all25.csv')
    ck('all25_ids_order_preserved',cov.gene_id.tolist()==candidates.gene_id.tolist() and len(cov)==25)
    ck('all25_pairs_have_all5_donors',paired.groupby('gene_id').size().eq(5).all() and set(paired.gene_id)==set(candidates.gene_id))
    missing=set(cov.loc[cov.mapping_status.ne('matched'),'gene_id'])
    ck('six_missing_candidates_preserved',len(missing)==6 and cov.mapping_status.eq('matched').sum()==19)
    ck('missing_candidates_paired_values_NA',paired.loc[paired.gene_id.isin(missing),['chol_nonzero_pct','hep_nonzero_pct','delta_nonzero_percentage_points','chol_mean_source_log2CPM','hep_mean_source_log2CPM','delta_mean_source_log2CPM']].isna().all().all())
    for label,ntypes in [('author',20),('broad',15)]:
        table=pd.read_csv(out/f'donor_{label}_results_all25.csv')
        ck(f'{label}_complete_25x5x{ntypes}',len(table)==25*5*ntypes and not table.duplicated(['gene_id','donor','cell_type']).any())
        ck(f'{label}_absent_groups_NA',table.loc[table.n_cells.eq(0),['pct_nonzero_source','mean_source_log2CPM']].isna().all().all())
        ck(f'{label}_missing_features_NA',table.loc[table.gene_id.isin(missing),['pct_nonzero_source','mean_source_log2CPM']].isna().all().all())
    context=pd.read_csv(out/'all_celltype_context_all25.csv')
    ck('all25_by15_context_retained',len(context)==375 and context.gene_id.nunique()==25)
    nk=context.loc[context.symbol.eq('AREG')&context.cell_type.eq('NK-like_Cells')].iloc[0]
    ck('NK_missing_donor_denominator_is3',nk.n_donors_with_cells==3 and nk.n_donors_with_available_expression==3)
    result={'verified_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'all_passed':all(checks.values()),'checks':checks,
            'max_abs_saved_csv_vs_independent_AREG_difference':max_error,'AREG_independent_source_aggregates':scalar,
            'original_protocol_frozen_utc':freeze['frozen_utc'],'source_scale_decision_frozen_utc':scale['frozen_utc'],
            'scope':'Independent CSV-reader AREG checks, equal-donor arithmetic, frozen sensitivity and all25 NA/grid preservation. No inferential or raw-read validation.'}
    (out/'focused_self_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({'all_passed':result['all_passed'],'n_checks':len(checks),'max_error':max_error}))

if __name__=='__main__':
    main()
