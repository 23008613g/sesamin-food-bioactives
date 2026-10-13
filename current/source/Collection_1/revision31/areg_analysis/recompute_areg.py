"""Recompute AREG from deposited organoid raw matrices and adult normalized matrix.

Only writes beside this script. Original project inputs are read-only.
See protocol.md for the fixed post-hoc scope and limitations.
"""
from pathlib import Path
import csv
import datetime as dt
import gc
import gzip
import hashlib
import json
import math
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd
import scipy
from scipy.io import mmread

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
SRC = BASE / 'organoid_expression' / 'sources'
SMALL = BASE / 'organoid_inputs'
ADULT = BASE / 'manuscript_revision_24' / 'external_single_cell' / 'inputs'
OLD = BASE / 'pfhn_submission_28'
AREG_ID = 'ENSG00000109321'
CHECKS = {}
LOG = []


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def js(name, data):
    (HERE / name).write_text(json.dumps(data, indent=2, ensure_ascii=False,
                                      allow_nan=False), encoding='utf-8')


def save(name, data):
    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)
    data.to_csv(HERE / name, index=False, na_rep='NA', float_format='%.17g')


def log(message):
    line = dt.datetime.now(dt.timezone.utc).isoformat() + ' ' + message
    LOG.append(line)
    (HERE / 'run.log').write_text('\n'.join(LOG) + '\n', encoding='utf-8')
    print(line, flush=True)


def check(name, condition):
    CHECKS[name] = bool(condition)
    if not condition:
        js('checks.json', {'all_passed': False, 'checks': CHECKS})
        raise AssertionError(name)


def decode(node):
    if isinstance(node, h5py.Dataset):
        return node.asstr()[...].tolist() if h5py.check_string_dtype(node.dtype) else node[...].tolist()
    cats = decode(node['categories'])
    return [cats[int(x)] if x >= 0 else None for x in node['codes'][...]]


def metrics(df):
    n = len(df)
    if not n:
        return {'n_cells': 0, 'n_raw_nonzero': np.nan, 'raw_detection_pct': np.nan,
                'AREG_raw_UMI_sum': np.nan, 'all_feature_UMI_sum': np.nan,
                'mean_log1p_CP10k': np.nan, 'pooled_UMI_CPM': np.nan}
    raw = df.AREG_raw_UMI.to_numpy()
    total = df.all_feature_raw_UMI.to_numpy()
    values = df.AREG_log1p_CP10k.to_numpy()
    return {'n_cells': n, 'n_raw_nonzero': int((raw > 0).sum()),
            'raw_detection_pct': float(100 * np.count_nonzero(raw) / n),
            'AREG_raw_UMI_sum': int(raw.sum()), 'all_feature_UMI_sum': int(total.sum()),
            'mean_log1p_CP10k': float(values.mean()),
            'pooled_UMI_CPM': float(1e6 * raw.sum() / total.sum())}


def main():
    started = time.perf_counter()
    log('Post-hoc source reanalysis begins; old manuscript results were already known.')
    manifest = pd.read_csv(BASE / 'organoid_sample_manifest.csv', dtype=str)
    files = [HERE / 'protocol.md', HERE / 'recompute_areg.py',
             BASE / 'organoid_sample_manifest.csv', SRC / 'GSE207889_os.h5ad',
             BASE / 'organoid_expression/link_expression.py',
             ADULT / 'GSE115469_Data.csv.gz', ADULT / 'GSE115469_CellClusterType.txt.gz',
             ADULT / 'GSE115469_family.soft.gz']
    for row in manifest.itertuples(index=False):
        files.extend([SRC / row.matrix_url.rsplit('/', 1)[-1],
                      SMALL / row.barcode_url.rsplit('/', 1)[-1],
                      SMALL / row.features_url.rsplit('/', 1)[-1]])
    old_summary = OLD / 'analysis/matched_scale_areg/outputs/AREG_summary.csv'
    old_adult = OLD / 'external_single_cell/results/figure_AREG_paired_donors.csv'
    files.extend([old_summary, old_adult])
    frozen = [{'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in files]
    js('input_hashes_before_run.json', {'frozen_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                                      'post_hoc': True, 'files': frozen})
    log('Input/protocol/script hashes frozen before expression calculations.')

    # Direct read of original author obs. No imported local annotation extraction.
    with h5py.File(SRC / 'GSE207889_os.h5ad', 'r') as f:
        obs = f['obs']
        fields = ['sample', 'condition', 'detailed_condition', 'celltype_category', 'celltype_scTypedb']
        labels = pd.DataFrame({field: decode(obs[field]) for field in fields})
        labels['source_obs_index'] = decode(obs[obs.attrs.get('_index', '_index')])
        original_obs_keys = list(obs.keys())
        original_obs_n = len(labels)
    aliases = {}
    for row in manifest.itertuples(index=False):
        for alias in [row.author_sample_id, row.author_short_id_alias]:
            aliases[str(alias)] = row.gsm
    labels['library_accession'] = labels['sample'].astype(str).map(aliases)
    labels = labels[labels.library_accession.notna()].copy()
    labels = labels.rename(columns={'sample': 'source_sample_label', 'condition': 'source_condition',
                                    'detailed_condition': 'source_detailed_condition',
                                    'celltype_category': 'source_broad_label',
                                    'celltype_scTypedb': 'source_subtype_label'})
    barcode_matches = labels.source_obs_index.str.findall(r'(?<![ACGTN])[ACGTN]{16}-[0-9]+')
    check('one_barcode_per_author_obs_index', barcode_matches.map(len).eq(1).all())
    labels['barcode'] = barcode_matches.str[0]
    labels['analysis_broad_label'] = labels.source_broad_label
    labels['analysis_subtype_label'] = labels.source_subtype_label
    labels['analysis_chol_ductal_excluding_SMC'] = labels.source_subtype_label.isin(['Cholangiocytes', 'Ductal cells'])
    labels['analysis_condition'] = labels.library_accession.map(manifest.set_index('gsm')['group']).replace({'CTRL-PA': 'Control'})
    check('retained_accession_barcode_unique', not labels.duplicated(['library_accession', 'barcode']).any())
    check('all_four_libraries_represented', set(labels.library_accession) == set(manifest.gsm))
    check('condition_labels_match', ((labels.source_condition.eq('CTRL') & labels.analysis_condition.eq('Control')) |
                                   (labels.source_condition.eq('PA') & labels.analysis_condition.eq('PA'))).all())
    crosswalk = labels.groupby(['source_broad_label', 'source_subtype_label',
                               'analysis_broad_label', 'analysis_subtype_label'], dropna=False).size().reset_index(name='n_cells')
    save('label_crosswalk.csv', crosswalk)
    save('cell_label_mapping.csv', labels)
    js('label_origin_audit.json', {
        'source': str(SRC / 'GSE207889_os.h5ad'), 'obs_keys': original_obs_keys,
        'source_broad_field': 'obs/celltype_category', 'source_subtype_field': 'obs/celltype_scTypedb',
        'broad_Cholangiocytes_origin': 'Directly deposited author obs field; not newly merged by our script.',
        'earlier_local_code_evidence': 'organoid_expression/link_expression.py lines 135-137 reads both obs fields directly; lines 161-167 preserves them.',
        'new_analysis_defined_group': 'Cholangiocytes + Ductal cells excluding Smooth muscle cells; post-hoc sensitivity only.',
        'retained_n': len(labels), 'original_obs_total_n': original_obs_n,
        'clinical_or_cell_identity_validation': False})
    log('Direct author-label crosswalk retained; broad Cholangiocytes is original author category.')

    cells = []
    joins = []
    for row in manifest.itertuples(index=False):
        lib = labels[labels.library_accession.eq(row.gsm)].copy()
        feature = pd.read_csv(SMALL / row.features_url.rsplit('/', 1)[-1], sep='\t', header=None)
        with gzip.open(SMALL / row.barcode_url.rsplit('/', 1)[-1], 'rt') as handle:
            barcodes = [x.strip() for x in handle]
        check(row.gsm + '_raw_barcodes_unique', len(barcodes) == len(set(barcodes)))
        lookup = {b: i for i, b in enumerate(barcodes)}
        absent = set(lib.barcode) - set(lookup)
        check(row.gsm + '_all_retained_barcodes_matched', not absent)
        indices = np.asarray([lookup[b] for b in lib.barcode], dtype=np.int64)
        areg_indices = np.flatnonzero(feature.iloc[:, 0].eq(AREG_ID))
        check(row.gsm + '_AREG_unique_exact_id', len(areg_indices) == 1)
        check(row.gsm + '_AREG_symbol', feature.iloc[areg_indices[0], 1] == 'AREG')
        matrix = mmread(SRC / row.matrix_url.rsplit('/', 1)[-1]).tocsr()
        check(row.gsm + '_shape_matches_annotations', matrix.shape == (len(feature), len(barcodes)))
        check(row.gsm + '_all_33538_measured_features', len(feature) == 33538)
        check(row.gsm + '_counts_nonnegative_integer', (matrix.data >= 0).all() and np.equal(matrix.data, np.floor(matrix.data)).all())
        totals = np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64)
        raw = matrix[areg_indices[0]].toarray().ravel().astype(np.int64)
        check(row.gsm + '_all_retained_cells_have_positive_totals', (totals[indices] > 0).all())
        lib['raw_matrix_column_0based'] = indices
        lib['feature_count_in_T'] = len(feature)
        lib['all_feature_raw_UMI'] = totals[indices]
        lib['AREG_raw_UMI'] = raw[indices]
        lib['AREG_log1p_CP10k'] = np.log1p(10000.0 * raw[indices] / totals[indices])
        scalar = math.fsum(math.log1p(10000 * int(x) / int(t)) for x, t in zip(raw[indices], totals[indices])) / len(indices)
        check(row.gsm + '_scalar_mean_agrees', abs(scalar - lib.AREG_log1p_CP10k.mean()) < 1e-12)
        cells.append(lib)
        joins.append({'library_accession': row.gsm, 'raw_matrix_n_cells': len(barcodes),
                      'retained_cells': len(lib), 'unmatched_retained': len(absent),
                      'raw_not_retained': len(barcodes) - len(lib), 'n_features': len(feature),
                      'raw_barcode_duplicates': 0, 'joined_cell_key_duplicates': 0})
        log(f'Read raw {row.gsm}: {matrix.shape}, retained {len(lib)} cells, all-feature denominator used.')
        del matrix
        gc.collect()
    cells = pd.concat(cells, ignore_index=True)
    save('cell_label_mapping.csv', cells.drop(columns=['AREG_raw_UMI', 'AREG_log1p_CP10k']))
    save('organoid_areg_cells.csv.gz', cells)
    save('raw_annotation_join_audit.csv', joins)
    check('retained_counts_reconcile_crosswalk', crosswalk.n_cells.sum() == len(cells))

    # Complete grid including absent subtype/library combinations (NA expressions).
    subtypes = sorted(cells.source_subtype_label.unique())
    rows = []
    for row in manifest.itertuples(index=False):
        lib = cells[cells.library_accession.eq(row.gsm)]
        for subtype in subtypes:
            s = lib[lib.source_subtype_label.eq(subtype)]
            rows.append({'library_accession': row.gsm, 'condition': lib.analysis_condition.iloc[0],
                         'source_subtype_label': subtype, 'analysis_subtype_label': subtype,
                         'source_broad_label': '|'.join(sorted(cells.loc[cells.source_subtype_label.eq(subtype), 'source_broad_label'].unique())),
                         **metrics(s)})
    subtype = pd.DataFrame(rows)
    save('areg_subtype_by_library.csv', subtype)
    pooled = [{'source_subtype_label': st, **metrics(cells[cells.source_subtype_label.eq(st)])} for st in subtypes]
    save('areg_subtype_pooled.csv', pooled)
    sensitivity = []
    population_rows = []
    filters = {'Whole organoid': lambda d: d,
               'Hepatocytes': lambda d: d[d.source_broad_label.eq('Hepatocytes')],
               'Deposited broad Cholangiocytes': lambda d: d[d.source_broad_label.eq('Cholangiocytes')],
               'Cholangiocytes + Ductal cells excluding SMC': lambda d: d[d.analysis_chol_ductal_excluding_SMC]}
    for accession in manifest.gsm:
        lib = cells[cells.library_accession.eq(accession)]
        for population, fn in filters.items():
            result = {'library_accession': accession, 'condition': lib.analysis_condition.iloc[0],
                      'population': population, **metrics(fn(lib))}
            population_rows.append(result)
            if population in list(filters)[2:]:
                sensitivity.append(result)
    lib_results = pd.DataFrame(population_rows)
    save('areg_population_by_library.csv', lib_results)
    save('areg_excluding_smc.csv', sensitivity)
    check('all_40_subtype_library_rows_kept', len(subtype) == len(subtypes) * 4)
    check('raw_detection_bounded', subtype.raw_detection_pct.dropna().between(0, 100).all())
    check('subtype_counts_sum_to_retained', int(subtype.n_cells.sum()) == len(cells))
    check('subtype_raw_UMIs_reconcile', int(subtype.AREG_raw_UMI_sum.sum()) == int(cells.AREG_raw_UMI.sum()))
    check('subtype_total_UMIs_reconcile', int(subtype.all_feature_UMI_sum.sum()) == int(cells.all_feature_raw_UMI.sum()))
    summary = []
    pooled_umi = []
    for population, fn in filters.items():
        for weighting in ['equal_library', 'pooled_cell']:
            r = {'gene_id': AREG_ID, 'symbol': 'AREG', 'population': population, 'weighting': weighting}
            for cond, prefix in [('Control', 'control'), ('PA', 'PA')]:
                values = lib_results[lib_results.population.eq(population) & lib_results.condition.eq(cond)]
                check(population + '_' + cond + '_two_libraries', len(values) == 2)
                r[prefix + '_n_libraries'] = len(values)
                r[prefix + '_n_cells'] = int(values.n_cells.sum())
                if weighting == 'equal_library':
                    mean = values.mean_log1p_CP10k.mean(skipna=False)
                else:
                    mean = np.average(values.mean_log1p_CP10k, weights=values.n_cells)
                    direct = fn(cells[cells.analysis_condition.eq(cond)]).AREG_log1p_CP10k.mean()
                    check(population + '_' + cond + '_pooled_reconstruction', abs(mean - direct) < 1e-12)
                r[prefix + '_mean'] = float(mean)
            r['PA_minus_control'] = r['PA_mean'] - r['control_mean']
            summary.append(r)
        r = {'population': population, 'estimand': 'log2 ratio of condition mean library pooled-UMI CPM plus 1'}
        for cond, prefix in [('Control', 'control'), ('PA', 'PA')]:
            r[prefix + '_mean_library_CPM'] = float(lib_results.loc[lib_results.population.eq(population) & lib_results.condition.eq(cond), 'pooled_UMI_CPM'].mean())
        r['log2_PA_control_plus1'] = math.log2((r['PA_mean_library_CPM'] + 1) / (r['control_mean_library_CPM'] + 1))
        pooled_umi.append(r)
    summary = pd.DataFrame(summary)
    save('areg_matched_scale.csv', summary)
    save('areg_pooled_umi_different_estimand.csv', pooled_umi)
    # A posteriori reconciliation after all new organoid expression was calculated.
    old = pd.read_csv(old_summary)
    comparison = summary.merge(old, on=['gene_id', 'symbol', 'population', 'weighting'], suffixes=('_new', '_old'))
    for field in ['control_mean', 'PA_mean', 'PA_minus_control']:
        comparison[field + '_difference'] = comparison[field + '_new'] - comparison[field + '_old']
    save('organoid_old_new_comparison.csv', comparison)
    check('four_old_matched_scale_rows_reconciled', len(comparison) == 4)
    check('organoid_matched_scale_agrees_with_old_within_1e-12', comparison.filter(like='_difference').abs().to_numpy().max() < 1e-12)
    log('Organoid raw-source matched-scale and subtype results complete.')

    # Adult normalized source matrix: read exact AREG row; no rescaling or raw-UMI claim.
    adult_labels = pd.read_csv(ADULT / 'GSE115469_CellClusterType.txt.gz', sep='\t')
    adult_labels.columns = ['cell_id', 'donor', 'barcode', 'cluster', 'source_cell_type']
    adult_labels['analysis_broad_label'] = adult_labels.source_cell_type.str.replace(r'^Hepatocyte_\d+$', 'Hepatocytes', regex=True)
    check('adult_cell_ids_unique', adult_labels.cell_id.is_unique)
    adult_source = gzip.open(ADULT / 'GSE115469_family.soft.gz', 'rt', encoding='utf-8').read()
    check('adult_source_scale_log2CPM_documented', 'log2CPM gene expression values' in adult_source)
    matched = []
    feature_count = 0
    with gzip.open(ADULT / 'GSE115469_Data.csv.gz', 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.reader(handle)
        header = next(reader)
        matrix_cells = header[1:]
        for row in reader:
            feature_count += 1
            if row[0] == 'AREG':
                matched.append(np.asarray(row[1:], dtype=np.float64))
    check('adult_AREG_one_source_row', len(matched) == 1)
    check('adult_matrix_unique_cells', len(matrix_cells) == len(set(matrix_cells)))
    check('adult_matrix_annotation_one_to_one', set(matrix_cells) == set(adult_labels.cell_id))
    adult_cells = adult_labels.set_index('cell_id').loc[matrix_cells].reset_index()
    adult_cells['AREG_source_log2CPM'] = matched[0]
    check('adult_source_values_finite_nonnegative', np.isfinite(matched[0]).all() and (matched[0] >= 0).all())
    check('adult_five_donors_retained', set(adult_cells.donor) == {f'P{i}TLH' for i in range(1, 6)})
    save('adult_areg_cells.csv.gz', adult_cells)
    donor_context = []
    for donor in sorted(adult_cells.donor.unique()):
        for celltype in sorted(adult_cells.analysis_broad_label.unique()):
            d = adult_cells[adult_cells.donor.eq(donor) & adult_cells.analysis_broad_label.eq(celltype)]
            v = d.AREG_source_log2CPM.to_numpy()
            donor_context.append({'donor': donor, 'cell_type': celltype, 'n_cells': len(d),
                                  'low_coverage_lt20': len(d) < 20,
                                  'n_nonzero_source': int((v > 0).sum()) if len(v) else np.nan,
                                  'pct_nonzero_source': float(100 * (v > 0).mean()) if len(v) else np.nan,
                                  'mean_source_log2CPM': float(v.mean()) if len(v) else np.nan,
                                  'scale': 'author-normalized source log2CPM; nonzero is not verified raw detection'})
    dc = pd.DataFrame(donor_context)
    save('areg_adult_all_celltypes_by_donor.csv', dc)
    pairs = []
    for donor in sorted(adult_cells.donor.unique()):
        c = dc[(dc.donor == donor) & (dc.cell_type == 'Cholangiocytes')].iloc[0]
        h = dc[(dc.donor == donor) & (dc.cell_type == 'Hepatocytes')].iloc[0]
        pairs.append({'donor': donor, 'n_cholangiocytes': c.n_cells, 'n_hepatocytes': h.n_cells,
                      'both_types_ge20': c.n_cells >= 20 and h.n_cells >= 20,
                      'chol_nonzero_pct': c.pct_nonzero_source, 'hep_nonzero_pct': h.pct_nonzero_source,
                      'delta_nonzero_percentage_points': c.pct_nonzero_source - h.pct_nonzero_source,
                      'chol_mean_source_log2CPM': c.mean_source_log2CPM, 'hep_mean_source_log2CPM': h.mean_source_log2CPM,
                      'delta_mean_source_log2CPM': c.mean_source_log2CPM - h.mean_source_log2CPM})
    pairs = pd.DataFrame(pairs)
    save('areg_adult_donor_context.csv', pairs)
    adult_summary = []
    for scope, d in [('all_five_donors_primary', pairs), ('both_types_ge20_posthoc_sensitivity', pairs[pairs.both_types_ge20])]:
        adult_summary.append({'scope': scope, 'n_donors': len(d),
                              'mean_delta_source_log2CPM': d.delta_mean_source_log2CPM.mean(),
                              'mean_delta_nonzero_percentage_points': d.delta_nonzero_percentage_points.mean(),
                              'n_positive': int(d.delta_mean_source_log2CPM.gt(0).sum()),
                              'n_negative': int(d.delta_mean_source_log2CPM.lt(0).sum()),
                              'n_zero': int(d.delta_mean_source_log2CPM.eq(0).sum())})
    save('areg_adult_summary.csv', adult_summary)
    all_context = []
    for celltype, d in dc.groupby('cell_type'):
        selected = adult_cells[adult_cells.analysis_broad_label.eq(celltype)]
        v = selected.AREG_source_log2CPM.to_numpy()
        all_context.append({'cell_type': celltype, 'n_cells_total': len(v),
                            'n_donors_present': int(d.n_cells.gt(0).sum()),
                            'equal_donor_mean_source_log2CPM': d.mean_source_log2CPM.mean(),
                            'equal_donor_pct_nonzero_source': d.pct_nonzero_source.mean(),
                            'pooled_cells_mean_source_log2CPM': v.mean(),
                            'pooled_cells_pct_nonzero_source': 100 * (v > 0).mean()})
    save('areg_adult_immune_context.csv', all_context)
    old = pd.read_csv(old_adult)
    comp = pairs.merge(old, on='donor', suffixes=('_new', '_old'))
    for field in ['delta_nonzero_percentage_points', 'delta_mean_source_log2CPM']:
        comp[field + '_difference'] = comp[field + '_new'] - comp[field + '_old']
    save('adult_old_new_comparison.csv', comp)
    check('adult_old_new_agree_archived_precision_1e-9', comp.filter(like='_difference').abs().to_numpy().max() < 1e-9)
    check('adult_context_counts_reconcile', int(dc.n_cells.sum()) == len(adult_cells))
    js('adult_scale_audit.json', {'source_matrix': str(ADULT / 'GSE115469_Data.csv.gz'),
                                'source_scale': 'GEO-deposited author-normalized log2CPM',
                                'n_features_scanned': feature_count, 'n_cells': len(adult_cells),
                                'new_normalization_applied': False, 'raw_UMI_detection_equivalence': 'not established',
                                'only_AREG_values_loaded': True})
    log('Adult matrix reread complete: all five donors, contrary and low-coverage results retained.')
    # Confirm original source bytes were not modified by the run.
    for item in frozen:
        check('source_unchanged_' + str(Path(item['path']).name), sha(Path(item['path'])) == item['sha256'])
    js('checks.json', {'all_passed': all(CHECKS.values()), 'n_checks': len(CHECKS), 'checks': CHECKS,
                       'validation_scope': 'Raw-source joins and descriptive arithmetic; not independent biological validation.'})
    js('execution.json', {'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                          'elapsed_seconds': time.perf_counter() - started,
                          'python_executable': sys.executable, 'python': sys.version,
                          'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__, 'h5py': h5py.__version__,
                          'command': 'project .venv/Scripts/python.exe recompute_areg.py',
                          'cells_organoid': len(cells), 'cells_adult': len(adult_cells),
                          'random_algorithm': False, 'new_UMAP': False,
                          'P_values_or_biological_CIs_computed': False})
    log(f'Completed; {len(CHECKS)} checks passed.')
    output_files = [p for p in HERE.iterdir() if p.is_file() and p.name != 'output_hashes.json']
    js('output_hashes.json', [{'name': p.name, 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in sorted(output_files)])


if __name__ == '__main__':
    main()
