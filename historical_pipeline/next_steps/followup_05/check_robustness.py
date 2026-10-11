"""Read-only follow-up to stage04; all outputs are descriptive."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / 'sesamin_exploration_04'
ORG = HERE.parent / 'organoid_expression'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def effect(a, b, pseudocount):
    if a + pseudocount <= 0 or b + pseudocount <= 0:
        return np.nan
    return np.log2((a + pseudocount) / (b + pseudocount))


def main():
    started = datetime.now(timezone.utc).isoformat()
    freeze = json.loads((HERE / 'robustness_freeze.json').read_text())
    assert sha(HERE / 'robustness_protocol.json') == freeze['protocol_sha256']
    candidate_path = PREVIOUS / 'priority_all.csv'
    candidates = pd.read_csv(candidate_path, keep_default_na=False)
    assert len(candidates) == 25 and sha(candidate_path) == freeze['candidate_sha256']
    original_protocol = json.loads((PREVIOUS / 'protocol.json').read_text(encoding='utf-8'))
    raw_path = (PREVIOUS / original_protocol['input']).resolve()
    assert sha(raw_path) == original_protocol['input_sha256']
    columns = original_protocol['sample_order']
    raw = pd.read_csv(raw_path, sep='\t', usecols=['name'] + columns).set_index('name')
    qc_path = PREVIOUS / 'sample_qc.csv'
    qc = pd.read_csv(qc_path).set_index('sample').loc[columns]
    assert candidates['gene_id'].isin(raw.index).all()
    normalized = raw.loc[candidates['gene_id'], columns].div(qc['median_ratio_factor'], axis=1)
    sample_rows, stats, sensitivity = [], [], []
    group_columns = {g: list(qc.index[qc['group'] == g]) for g in ['Control', 'OA', 'OA + sesamin']}
    for gene in candidates.itertuples():
        group_means = {}
        for group, samples in group_columns.items():
            values = normalized.loc[gene.gene_id, samples].to_numpy()
            raw_values = raw.loc[gene.gene_id, samples].to_numpy()
            group_means[group] = float(values.mean())
            stats.append({'gene_id': gene.gene_id, 'symbol': gene.symbol, 'group': group,
                          'mean_normalized': values.mean(), 'median_normalized': np.median(values),
                          'min_normalized': values.min(), 'max_normalized': values.max(),
                          'raw_nonzero_libraries': int((raw_values > 0).sum()),
                          'max_library_share_normalized': values.max() / values.sum() if values.sum() else np.nan})
            for sample in samples:
                sample_rows.append({'gene_id': gene.gene_id, 'symbol': gene.symbol, 'sample': sample,
                                    'group': group, 'raw_count': int(raw.at[gene.gene_id, sample]),
                                    'normalized_count': normalized.at[gene.gene_id, sample]})
        for pc in [0, 0.5, 1, 5]:
            d = effect(group_means['OA'], group_means['Control'], pc)
            a = effect(group_means['OA + sesamin'], group_means['OA'], pc)
            r = effect(group_means['OA + sesamin'], group_means['Control'], pc)
            sensitivity.append({'gene_id': gene.gene_id, 'symbol': gene.symbol, 'pseudocount': pc,
                                'OA_vs_control': d, 'sesamin_vs_OA': a, 'sesamin_vs_control': r})
    pd.DataFrame(sample_rows).to_csv(HERE / 'bulk_candidate_samples.csv', index=False)
    pd.DataFrame(stats).to_csv(HERE / 'bulk_group_ranges.csv', index=False)
    pd.DataFrame(sensitivity).to_csv(HERE / 'bulk_pseudocount_sensitivity.csv', index=False)
    separations = []
    for gene in candidates.itertuples():
        ctrl, oa, tr = [normalized.loc[gene.gene_id, group_columns[g]].to_numpy() for g in group_columns]
        disease_direction = np.sign(gene.log2_OA_vs_control)
        treatment_direction = np.sign(gene.log2_sesamin_vs_OA)
        separations.append({'gene_id': gene.gene_id, 'symbol': gene.symbol,
                            'OA_control_pairs_expected_direction_of9': int((np.sign(oa[:, None] - ctrl) == disease_direction).sum()),
                            'sesamin_OA_pairs_expected_direction_of9': int((np.sign(tr[:, None] - oa) == treatment_direction).sum()),
                            'comparison_unit': '9 dependent library-pair comparisons, not independent n=9'})
    pd.DataFrame(separations).to_csv(HERE / 'bulk_group_separation.csv', index=False)

    metadata_path = ORG / 'library_subtype_summary.csv'
    context_path = PREVIOUS / 'organoid_candidate_context.csv'
    meta = pd.read_csv(metadata_path)
    context = pd.read_csv(context_path)
    subtypes = meta.loc[meta['subtype'] != 'ALL_HEPATOCYTE_LIKE', 'subtype'].drop_duplicates().tolist()
    libraries = meta['gsm'].drop_duplicates().tolist()
    assert len(subtypes) == 5 and len(libraries) == 4
    totals = meta.pivot(index='gsm', columns='subtype', values='total_counts').loc[libraries]
    assert np.array_equal(totals[subtypes].sum(axis=1), totals['ALL_HEPATOCYTE_LIKE'])
    fractions = totals[subtypes].div(totals['ALL_HEPATOCYTE_LIKE'], axis=0)
    weights = fractions.mean(axis=0)
    assert np.isclose(weights.sum(), 1)
    fractions.to_csv(HERE / 'organoid_library_UMI_fractions.csv')
    pd.DataFrame({'subtype': weights.index, 'reference_UMI_weight': weights.values,
                  'minimum_cells_across_four_libraries': [int(meta.loc[meta['subtype'] == s, 'cells_aggregated'].min()) for s in weights.index]}).to_csv(HERE / 'organoid_reference_weights.csv', index=False)
    standardized, by_subtype = [], []
    reconstruction_error = 0.0
    for gene in candidates.itertuples():
        rows = context.loc[context['gene_id'] == gene.gene_id]
        assert len(rows) == 24
        covered = rows['status'].ne('feature_not_covered').all()
        wide = rows.pivot(index='gsm', columns='subtype', values='CPM').loc[libraries]
        for subtype in subtypes:
            subset = rows.loc[rows['subtype'] == subtype]
            ctr = subset.loc[subset['condition'] == 'CTRL-PA', 'CPM']
            pa = subset.loc[subset['condition'] == 'PA', 'CPM']
            item = {'gene_id': gene.gene_id, 'symbol': gene.symbol, 'subtype': subtype,
                    'feature_covered': covered, 'minimum_cells': int(subset['cells_aggregated'].min()),
                    'any_library_below50_cells': bool((subset['cells_aggregated'] < 50).any()),
                    'control_mean_CPM': ctr.mean(), 'PA_mean_CPM': pa.mean(),
                    'control_min_CPM': ctr.min(), 'control_max_CPM': ctr.max(),
                    'PA_min_CPM': pa.min(), 'PA_max_CPM': pa.max(),
                    'log2_CPM_ratio_plus1': effect(pa.mean(), ctr.mean(), 1) if covered else np.nan}
            for row in subset.itertuples():
                item[row.gsm + '_CPM'] = row.CPM
                item[row.gsm + '_UMI'] = row.candidate_UMI
            by_subtype.append(item)
        for gsm in libraries:
            if covered:
                reconstructed = float((wide.loc[gsm, subtypes] * fractions.loc[gsm]).sum())
                observed = wide.at[gsm, 'ALL_HEPATOCYTE_LIKE']
                standard = float((wide.loc[gsm, subtypes] * weights).sum())
                reconstruction_error = max(reconstruction_error, abs(reconstructed - observed))
                assert np.isclose(reconstructed, observed, rtol=1e-12, atol=1e-12)
            else:
                observed = standard = np.nan
            standardized.append({'gene_id': gene.gene_id, 'symbol': gene.symbol, 'gsm': gsm,
                                 'condition': meta.loc[meta['gsm'] == gsm, 'group'].iloc[0],
                                 'feature_covered': covered, 'observed_CPM': observed, 'fixed_weight_CPM': standard})
    standardized = pd.DataFrame(standardized)
    standardized.to_csv(HERE / 'organoid_reweighted_samples.csv', index=False)
    pd.DataFrame(by_subtype).to_csv(HERE / 'organoid_within_subtype.csv', index=False)
    comparisons = []
    for gene in candidates.itertuples():
        values = standardized.loc[standardized['gene_id'] == gene.gene_id]
        row = {'gene_id': gene.gene_id, 'symbol': gene.symbol, 'feature_covered': bool(values['feature_covered'].all())}
        for column in ['observed_CPM', 'fixed_weight_CPM']:
            ctrl = values.loc[values['condition'] == 'CTRL-PA', column].mean()
            pa = values.loc[values['condition'] == 'PA', column].mean()
            row[column + '_CTRL_mean'] = ctrl
            row[column + '_PA_mean'] = pa
            row[column + '_log2_ratio_plus1'] = effect(pa, ctrl, 1) if row['feature_covered'] else np.nan
        comparisons.append(row)
    comparisons = pd.DataFrame(comparisons)
    comparisons.to_csv(HERE / 'organoid_composition_sensitivity.csv', index=False)
    artifacts = [HERE / name for name in ['bulk_candidate_samples.csv', 'bulk_group_ranges.csv',
                 'bulk_pseudocount_sensitivity.csv', 'bulk_group_separation.csv',
                 'organoid_library_UMI_fractions.csv', 'organoid_reference_weights.csv',
                 'organoid_reweighted_samples.csv', 'organoid_within_subtype.csv',
                 'organoid_composition_sensitivity.csv']]
    record = {'started_utc': started, 'completed_utc': datetime.now(timezone.utc).isoformat(),
              'script_sha256': sha(Path(__file__)), 'protocol_sha256': sha(HERE / 'robustness_protocol.json'),
              'source_hashes': {str(f.relative_to(HERE.parent)): sha(f) for f in [candidate_path, qc_path, metadata_path, context_path]},
              'bulk_sample_rows': len(sample_rows), 'candidate_n': 25, 'subtype_rows': len(by_subtype),
              'reweighted_rows': len(standardized), 'reconstruction_max_abs_error_CPM': reconstruction_error,
              'spotlight': comparisons.loc[comparisons['symbol'].isin(['CXCL8', 'VNN1'])].to_dict('records'),
              'limitations': json.loads((HERE / 'robustness_protocol.json').read_text())['limitations'],
              'artifact_hashes': {f.name: sha(f) for f in artifacts}}
    (HERE / 'robustness_summary.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'candidate_n': 25, 'bulk_sample_rows': len(sample_rows), 'spotlight': record['spotlight']}))


if __name__ == '__main__':
    main()
