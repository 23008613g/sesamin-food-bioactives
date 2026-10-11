"""Bounded independent checks: derived tables, original embedding, one raw library.

Read-only outside this review directory. This is arithmetic/source-link checking,
not a new normalization analysis, embedding, annotation or biological validation.
"""
from pathlib import Path
import gzip
import hashlib
import json
import math
import sys

import h5py
import numpy as np
import pandas as pd
from scipy.io import mmread

HERE=Path(__file__).resolve().parent
ANALYSIS=HERE.parent/'analysis'
NEXT=HERE.parents[1]
ORIG=NEXT/'organoid_expression'
CHECKS=[]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check(name,passed,details=None):
    CHECKS.append({'check':name,'passed':bool(passed),'details':details})
    if not passed:
        raise AssertionError(name)

def decode(node):
    if isinstance(node,h5py.Dataset):
        return node.asstr()[...] if h5py.check_string_dtype(node.dtype) else node[...]
    labels=decode(node['categories'])
    return np.asarray([labels[i] if i>=0 else None for i in node['codes'][...]],dtype=object)

def compare_groups(cells,coverage,file,groupcols,subset=None,genes=None):
    actual=pd.read_csv(ANALYSIS/file)
    check(file+'_unique_rows',not actual.duplicated(groupcols+['gene_id']).any())
    if subset is None:
        subset=np.ones(len(cells),dtype=bool)
    if genes is None:
        genes=coverage.symbol.tolist()
    selected=cells.loc[subset]
    groups=selected.groupby(groupcols,dropna=False,sort=False)
    check(file+'_expected_rows',len(actual)==len(groups)*len(genes),{'rows':len(actual)})
    largest=0.0
    for keys,rows in groups:
        keys=keys if isinstance(keys,tuple) else (keys,)
        match=np.ones(len(actual),dtype=bool)
        for k,value in zip(groupcols,keys):
            match &= actual[k].eq(value).to_numpy()
        got=actual.loc[match].set_index('symbol')
        for gene in genes:
            item=got.loc[gene]
            covered=bool(coverage.set_index('symbol').at[gene,'covered'])
            checkname=f'{file}:{keys}:{gene}'
            if not (bool(item.covered)==covered and int(item.n_cells)==len(rows)):
                raise AssertionError(checkname+' coverage/n')
            fields=['n_detected_cells','detected_pct','mean_raw_counts','mean_cp10k','mean_log1p_cp10k','sum_raw_counts']
            if not covered:
                if not item[fields].isna().all():
                    raise AssertionError(checkname+' missing-gene values are not NA')
                continue
            raw=rows[gene+'__counts'].to_numpy(dtype=float)
            cp=raw/rows.raw_total_counts.to_numpy()*10000
            expected=np.asarray([(raw>0).sum(),(raw>0).mean()*100,raw.mean(),cp.mean(),np.log1p(cp).mean(),raw.sum()])
            stored=item[fields].to_numpy(dtype=float)
            if not np.allclose(stored,expected,rtol=6e-12,atol=1e-10):
                raise AssertionError(checkname+' aggregation mismatch')
            largest=max(largest,float(abs(stored-expected).max()))
    check(file+'_all_aggregates_recomputed',True,{'largest_absolute_rounding_difference':largest})

def main():
    summary=json.loads((ANALYSIS/'analysis_summary.json').read_text(encoding='utf-8'))
    cells=pd.read_csv(ANALYSIS/'cell_level_umap_qc_expression.csv.gz')
    coverage=pd.read_csv(ANALYSIS/'candidate_gene_coverage.csv')
    check('cell_rows_24034',len(cells)==24034)
    check('unique_library_plus_barcode',not cells.duplicated(['gsm','barcode']).any() and cells.cell_key.is_unique)
    check('compound_key_explicit',(cells.cell_key==cells.gsm+'|'+cells.barcode).all())
    repeated=cells.loc[cells.barcode.duplicated(False),'barcode'].nunique()
    check('cross_library_barcode_repeats_retained',repeated==57 and cells.barcode.duplicated().sum()==57,{'repeated_barcode_strings':int(repeated),'excess_rows':int(cells.barcode.duplicated().sum())})
    check('unique_original_observations',cells.original_obs_index.is_unique)
    old_annotations=pd.read_csv(ORIG/'four_library_annotations.csv.gz',dtype=str)
    joined=cells.merge(old_annotations,on=['gsm','barcode'],suffixes=('_new','_old'),validate='one_to_one')
    check('prior_annotation_join_full_coverage',len(joined)==len(cells)==len(old_annotations))
    for key in ['original_obs_index','raw_sample_label','condition','detailed_condition','celltype_category','celltype_scTypedb']:
        check('annotation_'+key,joined[key+'_new'].astype(str).eq(joined[key+'_old'].astype(str)).all())
    counts=cells.groupby('gsm').size().to_dict()
    check('library_counts',counts=={'GSM6322577':3093,'GSM6322578':5510,'GSM6322579':7399,'GSM6322580':8032},counts)
    hep=cells.celltype_category.eq('Hepatocytes')
    check('hepatocyte_cells',int(hep.sum())==18097,cells.loc[hep].groupby('gsm').size().to_dict())
    check('four_categories_ten_subtypes',cells.celltype_category.nunique()==4 and cells.celltype_scTypedb.nunique()==10)
    hb=cells.loc[cells.celltype_scTypedb.eq('HB1')].groupby('gsm').size().reindex(sorted(counts)).tolist()
    check('HB1_sparse_cell_counts',hb==[4,69,9,128],hb)
    check('raw_total_positive',(cells.raw_total_counts>0).all())
    check('raw_QC_integer_nonnegative',all(np.all(cells[k]>=0) and np.all(cells[k]%1==0) for k in ['raw_total_counts','raw_n_genes_by_counts','raw_mito_counts']))
    mito=100*cells.raw_mito_counts/cells.raw_total_counts
    check('MT_percentage_full_raw_denominator',np.allclose(mito,cells.raw_pct_mito_MT_prefix,rtol=6e-12,atol=1e-10))
    check('MT_proportion_range',cells.raw_pct_mito_MT_prefix.between(0,100).all())
    # Compare the author field only after documenting its fraction convention.
    author_fraction=cells.raw_mito_counts/cells.author_total_counts
    check('author_mito_fraction_definition',np.allclose(author_fraction,cells.author_pct_mito,rtol=1e-6,atol=1e-8),{'max_fraction_difference':float(abs(author_fraction-cells.author_pct_mito).max())})
    check('raw_vs_author_mito_percentage_difference',float(abs(cells.raw_pct_mito_MT_prefix-100*cells.author_pct_mito).max())<.306,{'maximum_percentage_point_difference':float(abs(cells.raw_pct_mito_MT_prefix-100*cells.author_pct_mito).max())})
    expectedmissing=coverage.loc[~coverage.covered,['gene_id','symbol']].to_dict('records')
    check('candidate_coverage',len(coverage)==25 and coverage.covered.sum()==24 and expectedmissing==[{'gene_id':'ENSG00000260105','symbol':'AOC4P'}])
    max_normalization=0.0
    for gene in coverage.symbol:
        raw=cells[gene+'__counts'].to_numpy()
        stored=cells[gene+'__log1p_cp10k'].to_numpy()
        covered=coverage.set_index('symbol').at[gene,'covered']
        if not covered:
            check('AOC4P_per_cell_count_and_normalization_are_NA',np.isnan(raw).all() and np.isnan(stored).all())
            continue
        check(gene+'_integer_counts',np.isfinite(raw).all() and (raw>=0).all() and (raw%1==0).all())
        expected=np.log1p(raw/cells.raw_total_counts.to_numpy()*10000)
        check(gene+'_log1p_cp10k',np.allclose(expected,stored,rtol=6e-12,atol=1e-10))
        check(gene+'_zero_preserved',np.all(stored[raw==0]==0))
        max_normalization=max(max_normalization,float(abs(expected-stored).max()))
    # The h5ad is inspected only for observation fields and existing coordinates.
    h5path=ORIG/'sources/GSE207889_os.h5ad'
    hstat=h5path.stat()
    with h5py.File(h5path,'r') as h5:
        obs=h5['obs']; all_ids=decode(obs[obs.attrs['_index']])
        check('atlas_observation_count',len(all_ids)==82467 and len(set(all_ids))==82467)
        source_pos=pd.Index(all_ids).get_indexer(cells.original_obs_index)
        check('all_source_obs_resolved',(source_pos>=0).all() and len(np.unique(source_pos))==len(cells))
        check('source_obs_positions',np.array_equal(source_pos,cells.author_obs_row.to_numpy()))
        coords=h5['obsm/X_umap'][...]
        check('original_embedding_shape',coords.shape==(82467,2) and np.isfinite(coords).all())
        e=abs(coords[source_pos,:]-cells[['umap_1','umap_2']].to_numpy())
        check('all_UMAP_coordinates_match_original',float(e.max())<1e-9,{'maximum_serialization_difference':float(e.max())})
        for key in ['condition','detailed_condition','celltype_category','celltype_scTypedb']:
            check('original_h5_'+key,np.array_equal(decode(obs[key])[source_pos].astype(str),cells[key].astype(str).to_numpy()))
        check('original_h5_sample',np.array_equal(decode(obs['sample'])[source_pos].astype(str),cells.raw_sample_label.astype(str).to_numpy()))
        if 'barcode' in obs:
            original_barcodes=decode(obs['barcode'])[source_pos].astype(str)
            check('original_h5_barcode_is_full_obs_index',np.array_equal(original_barcodes,cells.original_obs_index.astype(str).to_numpy()))
            parsed=pd.Series(original_barcodes).str.extract(r'^([ACGT]{16}-[0-9]+)',expand=False)
            check('original_h5_barcode_parsed_prefix',parsed.notna().all() and np.array_equal(parsed,cells.barcode.astype(str).to_numpy()))
    check('h5_size_mtime_unchanged',h5path.stat().st_size==hstat.st_size and h5path.stat().st_mtime_ns==hstat.st_mtime_ns)
    # One smallest raw library is reread, not all four matrices.
    manifest=pd.read_csv(NEXT/'organoid_sample_manifest.csv')
    lib=manifest.loc[manifest.gsm.eq('GSM6322577')].iloc[0]
    matrix_path=ORIG/'sources'/lib.matrix_url.rsplit('/',1)[-1]
    feature_path=NEXT/'organoid_inputs'/lib.features_url.rsplit('/',1)[-1]
    barcode_path=NEXT/'organoid_inputs'/lib.barcode_url.rsplit('/',1)[-1]
    with gzip.open(feature_path,'rt',encoding='utf-8') as f:
        features=[line.rstrip('\n').split('\t') for line in f]
    with gzip.open(barcode_path,'rt',encoding='utf-8') as f:
        barcodes=[line.strip() for line in f]
    with gzip.open(matrix_path,'rb') as f:
        matrix=mmread(f).tocsr()
    matrix.sum_duplicates(); matrix.eliminate_zeros()
    selected=cells.loc[cells.gsm.eq('GSM6322577')]
    positions=pd.Index(barcodes).get_indexer(selected.barcode)
    check('small_library_barcode_join',(positions>=0).all() and len(np.unique(positions))==3093)
    check('stored_raw_barcode_column',np.array_equal(positions,selected.raw_barcode_column_0based))
    raw=matrix[:,positions]
    check('small_library_raw_total_QC',np.array_equal(np.asarray(raw.sum(axis=0)).ravel(),selected.raw_total_counts))
    check('small_library_raw_detected_genes_QC',np.array_equal(raw.getnnz(axis=0),selected.raw_n_genes_by_counts))
    mitochondrial=[i for i,row in enumerate(features) if row[1].startswith('MT-')]
    check('small_library_MT_counts_QC',len(mitochondrial)==13 and np.array_equal(np.asarray(raw[mitochondrial,:].sum(axis=0)).ravel(),selected.raw_mito_counts))
    gene_ids=pd.Index([row[0] for row in features])
    for item in coverage.itertuples():
        pos=gene_ids.get_indexer([item.gene_id])[0]
        if pos<0:
            check('raw_AOC4P_absent',item.symbol=='AOC4P' and selected[item.symbol+'__counts'].isna().all())
        else:
            check('raw_GSM6322577_'+item.symbol,np.array_equal(raw[pos,:].toarray().ravel(),selected[item.symbol+'__counts']))
    aggregate_specs=[
       ('candidate25_library_subtype_all10.csv',['gsm','condition','celltype_category','celltype_scTypedb'],None,None),
       ('candidate25_subtype_all10_pooled.csv',['celltype_category','celltype_scTypedb'],None,None),
       ('candidate25_condition_subtype_all10.csv',['condition','celltype_category','celltype_scTypedb'],None,None),
       ('candidate25_hepatocyte_subtype_pooled.csv',['celltype_scTypedb'],hep,None),
       ('candidate25_library_hepatocyte_subtype.csv',['gsm','condition','celltype_scTypedb'],hep,None),
       ('spotlight4_category4_pooled.csv',['celltype_category'],None,['VNN1','AREG','JAG2','CXCL8']),
       ('spotlight4_subtype10_pooled.csv',['celltype_category','celltype_scTypedb'],None,['VNN1','AREG','JAG2','CXCL8'])]
    for file,keys,subset,genes in aggregate_specs:
        compare_groups(cells,coverage,file,keys,subset,genes)
    for file,category,subset in [
        ('cell_composition_library_category4.csv','celltype_category',cells),
        ('cell_composition_library_subtype10.csv','celltype_scTypedb',cells),
        ('cell_composition_library_hepatocyte_subtype5.csv','celltype_scTypedb',cells.loc[hep])]:
        data=pd.read_csv(ANALYSIS/file)
        for row in data.to_dict('records'):
            libdata=subset.loc[subset.gsm.eq(row['gsm'])]
            numerator=int(libdata[category].eq(row[category]).sum())
            check(file+'_'+row['gsm']+'_'+row[category],numerator==row['n_cells'] and len(libdata)==row['denominator_n_cells'] and math.isclose(row['cell_fraction_pct'],100*numerator/len(libdata),rel_tol=6e-12,abs_tol=1e-10))
        check(file+'_fractions_sum100',np.allclose(data.groupby('gsm').cell_fraction_pct.sum(),100,atol=1e-9,rtol=0))
    qc=pd.read_csv(ANALYSIS/'qc_library_summary.csv')
    for _,row in qc.iterrows():
        subset=cells.loc[cells.gsm.eq(row.gsm)]
        for metric in ['raw_total_counts','raw_n_genes_by_counts','raw_pct_mito_MT_prefix']:
            for pct in [0,10,25,50,75,90,100]:
                val=float(subset[metric].quantile(pct/100))
                check(f'QC_{row.gsm}_{metric}_q{pct}',math.isclose(val,row[f'{metric}_q{pct:02d}'],rel_tol=6e-12,abs_tol=1e-9))
    inputs=[{'file':str(file),'sha256':sha(file)} for file in sorted(ANALYSIS.glob('*')) if file.is_file()]
    passed={'status':'PASS_BOUNDED_TECHNICAL_CHECKS','date':'2026-10-09','checks':CHECKS,
            'interpreter':sys.executable,'source_scope':'All small derived tables, all source h5ad observations and coordinates, and one smallest raw library only.',
            'raw_library_rechecked':'GSM6322577','raw_matrix_sha256':sha(matrix_path),
            'all_cell_normalization_max_abs_serialization_error':max_normalization,
            'analysis_artifacts_at_review':inputs,
            'limits':['No reread of the other three raw matrices','No new embedding or clustering','Technical arithmetic review, not independent biological validation','Figure and manuscript review recorded separately when frozen']}
    (HERE/'numeric_review.json').write_text(json.dumps(passed,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':passed['status'],'checks':len(CHECKS),'raw_library':'GSM6322577','normalization_max_error':max_normalization}))

if __name__=='__main__':
    try:
        main()
    except Exception as error:
        (HERE/'numeric_review_failure.json').write_text(json.dumps({'error':repr(error),'checks':CHECKS},indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        raise
