"""Bounded reanalysis of deposited GSE115469 expression; see frozen protocol.

Run from this folder with Python >=3.10, NumPy and pandas:
  python analyze_external.py --input-dir inputs --output-dir results
No random algorithms, additional downloads, dependency installation, or raw-read processing.
"""
import argparse
import datetime as dt
import gzip
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def dump_json(path, obj):
    # JSON null rather than nonstandard NaN preserves missingness portably.
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')


def records(df):
    return json.loads(df.to_json(orient='records', double_precision=15))


def save_csv(df, path):
    df.to_csv(path, index=False, na_rep='NA', float_format='%.12g')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=Path('inputs'))
    parser.add_argument('--output-dir', type=Path, default=Path('results'))
    args = parser.parse_args()
    inp, out = args.input_dir, args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    here = Path(__file__).resolve().parent
    start = dt.datetime.now(dt.timezone.utc)
    timer = time.perf_counter()
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        if not condition:
            dump_json(out / 'calculation_checks.json', checks)
            raise AssertionError(name)

    frozen = json.loads((here / 'protocol_freeze.json').read_text('utf-8'))
    check('protocol_matches_freeze', sha(here/'protocol.md') == frozen['files']['protocol.md'])
    scale_freeze = json.loads((here/'source_scale_decision_freeze.json').read_text('utf-8'))
    check('source_scale_decision_matches_freeze', sha(here/'source_scale_decision.md') == scale_freeze['sha256'])
    candidates = pd.read_csv(inp/'human_candidates_all25.csv')[['gene_id', 'symbol']]
    xwalk = pd.read_csv(inp/'original_candidate_gene_coverage.csv', keep_default_na=False)
    check('fixed_25_unique_candidate_ids', len(candidates) == 25 and candidates.gene_id.is_unique)
    check('candidate_id_order_matches_original_crosswalk', candidates.gene_id.tolist() == xwalk.gene_id.tolist())
    for local in ['human_candidates_all25.csv', 'original_candidate_gene_coverage.csv']:
        expected = frozen['files'].get('inputs\\'+local, frozen['files'].get('inputs/'+local))
        check('frozen_'+local, sha(inp/local) == expected)
    annot = pd.read_csv(inp/'GSE115469_CellClusterType.txt.gz', sep='\t')
    annot.columns = ['cell_id', 'donor', 'barcode', 'cluster', 'author_cell_type']
    check('8444_author_cells_unique', len(annot) == 8444 and annot.cell_id.is_unique)
    check('five_expected_donors', set(annot.donor) == {f'P{i}TLH' for i in range(1, 6)})
    check('twenty_author_clusters_and_labels', annot.cluster.nunique() == 20 and annot.author_cell_type.nunique() == 20)
    check('cell_prefix_matches_donor', all(a.startswith(b+'_') for a, b in zip(annot.cell_id, annot.donor)))
    annot['broad_cell_type'] = annot.author_cell_type.str.replace(r'^Hepatocyte_\d+$', 'Hepatocytes', regex=True)
    donor_order = [f'P{i}TLH' for i in range(1,6)]
    annot['geo_sample'] = annot.donor.map({f'P{i}TLH': f'GSM{3178781+i}' for i in range(1,6)})
    geo_text = gzip.open(inp/'GSE115469_family.soft.gz', 'rt', encoding='utf-8').read()
    check('source_scale_documented_log2CPM', geo_text.count('log2CPM gene expression values') == 5)
    check('GEO_donor_titles_present', all(f'Patient {i} Total Liver Homogenate' in geo_text for i in range(1,6)))
    matrix_path = inp/'GSE115469_Data.csv.gz'
    header = pd.read_csv(matrix_path, nrows=0).columns.tolist()
    cells = header[1:]
    check('matrix_cells_unique_and_match_annotation', len(cells) == len(set(cells)) == 8444 and set(cells) == set(annot.cell_id))
    annot = annot.set_index('cell_id').loc[cells].reset_index()
    lookup = {}
    for row, alias in zip(candidates.itertuples(index=False), xwalk.itertuples(index=False)):
        lookup[row.gene_id] = set([row.gene_id, row.symbol] + ([alias.raw_feature_symbol] if alias.raw_feature_symbol not in ['', 'NA'] else []))
    union_alias = set().union(*lookup.values())
    feature_names, selected = [], {}
    matrix_min, matrix_max = math.inf, -math.inf
    exact_zeros = nonintegral_entries = matrix_entries = 0
    nonzero_genes = np.zeros(len(cells), dtype=np.int64)
    for chunk in pd.read_csv(matrix_path, chunksize=256):
        names = chunk.iloc[:,0].astype(str).tolist()
        arr = chunk.iloc[:,1:].to_numpy(dtype=np.float64)
        check('all_matrix_values_finite_so_far', bool(np.isfinite(arr).all()))
        matrix_min, matrix_max = min(matrix_min,float(arr.min())), max(matrix_max,float(arr.max()))
        exact_zeros += int(np.count_nonzero(arr == 0))
        nonintegral_entries += int(np.count_nonzero(arr != np.floor(arr)))
        matrix_entries += arr.size
        nonzero_genes += (arr > 0).sum(axis=0)
        for i, name in enumerate(names):
            if name in union_alias:
                selected.setdefault(name, []).append(arr[i].copy())
        feature_names.extend(names)
    check('source_matrix_nonnegative', matrix_min >= 0)
    check('source_matrix_not_raw_integer_counts', nonintegral_entries > 0)
    check('source_matrix_has_exact_zeros', exact_zeros > 0)
    check('source_feature_names_unique', len(feature_names) == len(set(feature_names)))
    print(f'Read {len(feature_names)} features x {len(cells)} cells; source range {matrix_min:g} to {matrix_max:g}', flush=True)
    scale_audit = {'n_features':len(feature_names), 'n_cells':len(cells), 'matrix_entries':matrix_entries,
                   'minimum':matrix_min, 'maximum':matrix_max, 'exact_zero_entries':exact_zeros,
                   'nonintegral_entries':nonintegral_entries,
                   'scale':'author-normalized log2CPM (GEO label)',
                   'normalization_applied_here':False,
                   'nonzero_metric':'percentage of deposited source values >0; raw-UMI equivalence not verified',
                   'raw_UMI_counts_available_in_downloaded_inputs':False,
                   'raw_UMI_depth_available':False}
    dump_json(out/'matrix_scale_audit.json', scale_audit)
    pd.DataFrame({'feature':feature_names}).to_csv(out/'feature_universe.csv',index=False)
    annot['n_nonzero_source_features'] = nonzero_genes
    mapping, expr = [], {}
    for candidate in candidates.itertuples(index=False):
        found = sorted(lookup[candidate.gene_id] & set(selected))
        status = 'matched' if len(found)==1 and len(selected[found[0]])==1 else 'absent_from_deposited_features' if not found else 'ambiguous'
        mapping.append({'gene_id':candidate.gene_id,'symbol':candidate.symbol,
                        'mapping_status':status,'matched_source_feature':';'.join(found),
                        'allowed_exact_aliases':';'.join(sorted(lookup[candidate.gene_id]))})
        expr[candidate.gene_id] = selected[found[0]][0] if status=='matched' else np.full(len(cells), np.nan)
    mapping = pd.DataFrame(mapping)
    save_csv(mapping,out/'candidate_feature_coverage_all25.csv')
    check('AREG_unambiguous_exact_match', mapping.loc[mapping.symbol.eq('AREG'),'matched_source_feature'].tolist()==['AREG'])
    full_cells = annot.copy()
    for gene_id in candidates.gene_id:
        full_cells[gene_id] = expr[gene_id]
    # Store complete frozen candidate observations to enable independent bounded reruns.
    with gzip.GzipFile(filename=str(out/'cell_source_expression_all25.csv.gz'), mode='wb', mtime=0) as gz:
        gz.write(full_cells.to_csv(index=False,na_rep='NA',float_format='%.17g').encode('utf-8'))
    save_csv(annot,out/'cell_metadata.csv')
    results = {}
    for label, type_col in [('author','author_cell_type'),('broad','broad_cell_type')]:
        types = sorted(annot[type_col].unique())
        coverage_rows, value_rows = [], []
        for donor in donor_order:
            for celltype in types:
                mask = (annot.donor.eq(donor) & annot[type_col].eq(celltype)).to_numpy()
                n = int(mask.sum())
                coverage_rows.append({'donor':donor,'geo_sample':f'GSM{3178781+int(donor[1])}',
                                      'cell_type':celltype,'n_cells':n,'sparse_lt20':n<20,
                                      'cell_fraction_within_donor':n/int(annot.donor.eq(donor).sum()),
                                      'mean_n_nonzero_source_features':float(nonzero_genes[mask].mean()) if n else np.nan,
                                      'median_n_nonzero_source_features':float(np.median(nonzero_genes[mask])) if n else np.nan})
                for candidate in candidates.itertuples(index=False):
                    values = expr[candidate.gene_id][mask]
                    available = n>0 and np.isfinite(values).all()
                    nonzero = int(np.count_nonzero(values >0)) if available else np.nan
                    value_rows.append({'gene_id':candidate.gene_id,'symbol':candidate.symbol,
                                       'donor':donor,'cell_type':celltype,'n_cells':n,
                                       'sparse_lt20':n<20,'available':available,
                                       'n_nonzero_source':nonzero,
                                       'pct_nonzero_source':100*nonzero/n if available else np.nan,
                                       'mean_source_log2CPM':float(values.mean()) if available else np.nan})
        coverage, table = pd.DataFrame(coverage_rows),pd.DataFrame(value_rows)
        save_csv(coverage,out/f'cell_coverage_{label}_by_donor.csv')
        save_csv(table,out/f'donor_{label}_results_all25.csv')
        check(f'{label}_coverage_sums_8444', coverage.n_cells.sum()==8444)
        check(f'{label}_complete_grid',len(table)==25*5*len(types))
        check(f'{label}_percent_in_bounds',table.pct_nonzero_source.dropna().between(0,100).all())
        check(f'{label}_donor_fractions_sum_one',np.allclose(coverage.groupby('donor').cell_fraction_within_donor.sum(),1,atol=1e-12))
        results[label]=table
    broad = results['broad']
    pair_rows=[]
    for candidate in candidates.itertuples(index=False):
        g = broad.loc[broad.gene_id.eq(candidate.gene_id)]
        for donor in donor_order:
            c = g.loc[g.donor.eq(donor)&g.cell_type.eq('Cholangiocytes')].iloc[0]
            h = g.loc[g.donor.eq(donor)&g.cell_type.eq('Hepatocytes')].iloc[0]
            pair_rows.append({'gene_id':candidate.gene_id,'symbol':candidate.symbol,'donor':donor,
                              'n_cholangiocytes':int(c.n_cells),'n_hepatocytes':int(h.n_cells),
                              'both_types_ge20':c.n_cells>=20 and h.n_cells>=20,
                              'paired_available':bool(c.available and h.available),
                              'chol_nonzero_pct':c.pct_nonzero_source,'hep_nonzero_pct':h.pct_nonzero_source,
                              'delta_nonzero_percentage_points':c.pct_nonzero_source-h.pct_nonzero_source,
                              'chol_mean_source_log2CPM':c.mean_source_log2CPM,'hep_mean_source_log2CPM':h.mean_source_log2CPM,
                              'delta_mean_source_log2CPM':c.mean_source_log2CPM-h.mean_source_log2CPM})
    pairs=pd.DataFrame(pair_rows)
    save_csv(pairs,out/'paired_cholangiocyte_minus_hepatocyte_all25.csv')
    check('125_pairs_preserved',len(pairs)==125)
    summary_rows=[]
    for candidate in candidates.itertuples(index=False):
        g = pairs.loc[pairs.gene_id.eq(candidate.gene_id)]
        for scope in ['all_five_donors','both_types_ge20']:
            inc = g.loc[g.paired_available & (True if scope=='all_five_donors' else g.both_types_ge20)]
            for metric in ['delta_nonzero_percentage_points','delta_mean_source_log2CPM']:
                v=inc[metric]
                summary_rows.append({'gene_id':candidate.gene_id,'symbol':candidate.symbol,'scope':scope,'metric':metric,
                                     'n_complete_donor_pairs':len(inc),'included_donors':';'.join(inc.donor),
                                     'excluded_or_missing_donors':';'.join(g.loc[~g.donor.isin(inc.donor),'donor']),
                                     'mean_difference':v.mean(),'median_difference':v.median(),'minimum_difference':v.min(),
                                     'maximum_difference':v.max(),'n_positive':int((v>0).sum()),'n_negative':int((v<0).sum()),'n_zero':int((v==0).sum())})
    summaries=pd.DataFrame(summary_rows)
    save_csv(summaries,out/'paired_descriptive_summary_all25.csv')
    context_rows=[]
    max_pooled_error=0.
    for (gene_id,symbol,celltype), t in broad.groupby(['gene_id','symbol','cell_type'],sort=False):
        valid=t.loc[t.available]
        mask=annot.broad_cell_type.eq(celltype).to_numpy()
        values=expr[gene_id][mask]
        valid_expr=np.isfinite(values).all()
        pooled_pct=float(100*np.mean(values>0)) if valid_expr else np.nan
        pooled_mean=float(values.mean()) if valid_expr else np.nan
        if valid_expr:
            max_pooled_error=max(max_pooled_error,abs(float(np.average(valid.mean_source_log2CPM,weights=valid.n_cells))-pooled_mean))
            check('pooled_nonzero_reconstructed',np.isclose(valid.n_nonzero_source.sum()/valid.n_cells.sum()*100,pooled_pct,atol=1e-12))
        context_rows.append({'gene_id':gene_id,'symbol':symbol,'cell_type':celltype,'n_donors_with_cells':int(t.n_cells.gt(0).sum()),
                             'n_donors_with_available_expression':len(valid),'n_cells':int(t.n_cells.sum()),
                             'pooled_nonzero_pct':pooled_pct,'pooled_mean_source_log2CPM':pooled_mean,
                             'equal_donor_nonzero_pct':valid.pct_nonzero_source.mean(),
                             'equal_donor_mean_source_log2CPM':valid.mean_source_log2CPM.mean()})
    context=pd.DataFrame(context_rows)
    save_csv(context,out/'all_celltype_context_all25.csv')
    check('pooled_mean_reconstructed',max_pooled_error<1e-12)
    check('all25_retained_context',context.gene_id.nunique()==25)
    absent=set(mapping.loc[mapping.mapping_status.ne('matched'),'gene_id'])
    check('missing_features_never_zero_imputed',broad.loc[broad.gene_id.isin(absent),['pct_nonzero_source','mean_source_log2CPM']].isna().all().all())
    # Independent scalar aggregation of every AREG donor/category using math.fsum.
    areg_id=candidates.loc[candidates.symbol.eq('AREG'),'gene_id'].iloc[0]
    scalar_rows=[]
    max_scalar_error=0.
    for donor in donor_order:
        for celltype in ['Cholangiocytes','Hepatocytes']:
            indices=[i for i,row in enumerate(annot.itertuples(index=False)) if row.donor==donor and row.broad_cell_type==celltype]
            vv=[float(expr[areg_id][i]) for i in indices]
            sc_mean=math.fsum(vv)/len(vv)
            sc_pct=sum(1 for v in vv if v>0)*100/len(vv)
            observed=broad.loc[broad.gene_id.eq(areg_id)&broad.donor.eq(donor)&broad.cell_type.eq(celltype)].iloc[0]
            max_scalar_error=max(max_scalar_error,abs(sc_mean-observed.mean_source_log2CPM),abs(sc_pct-observed.pct_nonzero_source))
            scalar_rows.append({'donor':donor,'cell_type':celltype,'n_cells':len(vv),'n_nonzero_source':sum(v>0 for v in vv),
                                'mean_scalar':sc_mean,'mean_vectorized':observed.mean_source_log2CPM,'pct_scalar':sc_pct})
    save_csv(pd.DataFrame(scalar_rows),out/'AREG_independent_scalar_checks.csv')
    check('AREG_all10_scalar_checks',max_scalar_error<1e-12)
    check('paired_direction_identity',np.allclose(pairs.delta_mean_source_log2CPM,pairs.chol_mean_source_log2CPM-pairs.hep_mean_source_log2CPM,equal_nan=True))
    areg_pairs=pairs.loc[pairs.symbol.eq('AREG')]
    check('five_complete_AREG_pairs',areg_pairs.paired_available.all() and len(areg_pairs)==5)
    areg_context=context.loc[context.symbol.eq('AREG')]
    save_csv(areg_pairs,out/'figure_AREG_paired_donors.csv')
    save_csv(areg_context,out/'figure_AREG_all_cell_context.csv')
    dump_json(out/'calculation_checks.json',{'all_passed':all(checks.values()),'checks':checks,
                'max_abs_scalar_difference':max_scalar_error,'max_abs_pooled_mean_reconstruction_error':max_pooled_error,
                'scope':'Bounded source-matrix arithmetic and identity checks; not raw FASTQ reproduction, biological validation, or normalization re-estimation.'})
    summary={'dataset':'GSE115469','doi':'10.1038/s41467-018-06318-7','status':'completed',
             'eligible_for':'descriptive external adult-liver cellular context',
             'not_validation_of':['sesamin response','lipotoxicity response','hepatoprotection','clinical benefit'],
             'n_donors':5,'n_cells':8444,'n_author_cell_types':20,'n_broad_types_after_hepatocyte_collapse':15,
             'n_candidates_fixed':25,'n_candidates_matched':int(mapping.mapping_status.eq('matched').sum()),
             'missing_candidates':records(mapping.loc[mapping.mapping_status.ne('matched')]),
             'source_scale_audit':scale_audit,'AREG_pairs':records(areg_pairs),
             'AREG_paired_summaries':records(summaries.loc[summaries.symbol.eq('AREG')]),
             'AREG_full_context':records(areg_context),'inferential_tests_performed':False,
             'independence_basis':'Primary adult donor tissues in a separate study versus W01 hPSC-derived GSE207889 HLOs; genotype identity not tested.',
             'limitation':'Only author-normalized expression, no raw UMI or verified UMI-detection equivalence; sparse cholangiocytes P1/P4, unequal donor/cell coverage, no treatment, possible brain-death inflammation.'}
    dump_json(out/'analysis_summary.json',summary)
    input_manifest=[{'path':str(p.relative_to(inp)), 'bytes':p.stat().st_size, 'sha256':sha(p)} for p in sorted(inp.iterdir()) if p.is_file()]
    dump_json(out/'execution_record.json',{'started_utc':start.isoformat(),'finished_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
                'elapsed_seconds':time.perf_counter()-timer,'command':[sys.executable,*sys.argv],
                'environment':{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__},
                'script_sha256':sha(Path(__file__)),'protocol_sha256':sha(here/'protocol.md'),
                'source_scale_decision_sha256':sha(here/'source_scale_decision.md'),'inputs':input_manifest,
                'reproducibility_scope':'Deposited processed matrix to donor-level and complete fixed-candidate descriptive tables; no raw-read or author normalization/clustering rerun.'})
    print(areg_pairs.to_string(index=False),flush=True)
    print('Missing candidates: '+', '.join(mapping.loc[mapping.mapping_status.ne('matched'),'symbol']),flush=True)
    print('All bounded checks passed.',flush=True)


if __name__=='__main__':
    main()
