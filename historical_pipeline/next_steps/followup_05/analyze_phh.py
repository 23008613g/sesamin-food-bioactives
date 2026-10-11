"""All-candidate, label-paired PHH reference analysis; no inferential tests."""
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def median_factors(matrix):
    positive = (matrix > 0).all(axis=1)
    ref = matrix[positive]
    assert len(ref) >= 100
    geometric = np.exp(np.log(ref).mean(axis=1))
    scale = np.median(ref / geometric[:, None], axis=0)
    scale /= np.exp(np.log(scale).mean())
    return scale, int(positive.sum())


def main():
    started = datetime.now(timezone.utc).isoformat()
    spec = json.loads((HERE / 'phh_protocol.json').read_text())
    freeze = json.loads((HERE / 'phh_freeze.json').read_text())
    assert sha(HERE / 'phh_protocol.json') == freeze['protocol_sha256']
    raw_path = (HERE / spec['input']).resolve()
    assert sha(raw_path) == freeze['raw_count_sha256']
    manifest_path = HERE / spec['manifest']
    meta = pd.read_csv(manifest_path, keep_default_na=False)
    columns = meta['matrix_column'].tolist()
    assert len(columns) == len(set(columns)) == len(set(meta['gsm'])) == 30
    assert set(meta['condition']) == {'Control', 'FFA'}
    assert meta.groupby(['donor_id', 'condition']).size().eq(3).all()
    donors = meta['donor_id'].drop_duplicates().tolist()
    assert len(donors) == 5
    table = pd.read_csv(raw_path, sep='\t', usecols=['gene_id', 'gene_name'] + columns, keep_default_na=False)
    assert table['gene_id'].is_unique
    counts = table[columns].to_numpy(dtype=float)
    assert counts.shape == (58735, 30) and np.isfinite(counts).all()
    assert (counts >= 0).all() and (counts % 1 == 0).all()
    totals = counts.sum(axis=0)
    factors, reference_n = median_factors(counts)
    norm = counts / factors
    cpm = counts / totals * 1e6
    identifiers = {gene: i for i, gene in enumerate(table['gene_id'])}
    candidate_path = HERE / spec['candidate_table']
    assert sha(candidate_path) == freeze['candidate_sha256']
    candidates = pd.read_csv(candidate_path, keep_default_na=False)
    assert len(candidates) == 25
    cluster_indices = {}
    pooled_keys = []
    for donor in donors:
        for condition in ['Control', 'FFA']:
            key = (donor, condition)
            cluster_indices[key] = np.flatnonzero((meta['donor_id'] == donor) & (meta['condition'] == condition))
            pooled_keys.append(key)
    pooled_counts = np.column_stack([counts[:, cluster_indices[key]].sum(axis=1) for key in pooled_keys])
    pool_factors, pool_ref_n = median_factors(pooled_counts)
    pooled_norm = pooled_counts / pool_factors
    pool_index = {key: i for i, key in enumerate(pooled_keys)}
    samples, effects = [], []
    for candidate in candidates.itertuples():
        row_id = identifiers.get(candidate.gene_id)
        covered = row_id is not None
        for j, label in enumerate(meta.itertuples()):
            samples.append({'gene_id': candidate.gene_id, 'candidate_symbol': candidate.symbol,
                            'reference_symbol': table.at[row_id, 'gene_name'] if covered else None,
                            'sample': label.matrix_column, 'gsm': label.gsm, 'donor_label': label.donor_id,
                            'condition': label.condition, 'feature_covered': covered,
                            'raw_count': int(counts[row_id, j]) if covered else None,
                            'normalized_count': norm[row_id, j] if covered else None,
                            'CPM': cpm[row_id, j] if covered else None})
        for donor in donors:
            row = {'gene_id': candidate.gene_id, 'symbol': candidate.symbol, 'donor_label': donor,
                   'feature_covered': covered, 'HepG2_OA_log2_effect': candidate.log2_OA_vs_control}
            if covered:
                means = {}
                for condition in ['Control', 'FFA']:
                    idx = cluster_indices[(donor, condition)]
                    nvalues, cvalues = norm[row_id, idx], cpm[row_id, idx]
                    row[condition + '_mean_norm'] = nvalues.mean()
                    row[condition + '_min_norm'] = nvalues.min()
                    row[condition + '_max_norm'] = nvalues.max()
                    row[condition + '_mean_CPM'] = cvalues.mean()
                    row[condition + '_CPM_ge1_libraries'] = int((cvalues >= 1).sum())
                    row[condition + '_nonzero_libraries'] = int((counts[row_id, idx] > 0).sum())
                    means[condition] = nvalues.mean()
                main = float(np.log2((means['FFA'] + 1) / (means['Control'] + 1)))
                simple = float(np.log2((row['FFA_mean_CPM'] + 1) / (row['Control_mean_CPM'] + 1)))
                pooled = float(np.log2((pooled_norm[row_id, pool_index[(donor, 'FFA')]] + 1) /
                                       (pooled_norm[row_id, pool_index[(donor, 'Control')]] + 1)))
                row.update({'log2_FFA_Control_mean_norm_plus1': main,
                            'log2_FFA_Control_mean_CPM_plus1': simple,
                            'log2_FFA_Control_pooled_norm_plus1': pooled,
                            'sign_matches_CPM': bool(np.sign(main) == np.sign(simple)),
                            'sign_matches_pooled': bool(np.sign(main) == np.sign(pooled)),
                            'same_disease_direction_as_HepG2': bool(np.sign(main) == np.sign(candidate.log2_OA_vs_control))})
            effects.append(row)
    pd.DataFrame(samples).to_csv(HERE / 'phh_candidate_samples.csv', index=False)
    effects = pd.DataFrame(effects)
    effects.to_csv(HERE / 'phh_donor_effects.csv', index=False)
    summaries, omit_rows = [], []
    for candidate in candidates.itertuples():
        df = effects.loc[effects['gene_id'] == candidate.gene_id]
        covered = bool(df['feature_covered'].all())
        row = {'gene_id': candidate.gene_id, 'symbol': candidate.symbol, 'feature_covered': covered}
        if covered:
            values = df['log2_FFA_Control_mean_norm_plus1']
            row.update({'label_n': len(df), 'median_log2_effect': float(values.median()),
                        'minimum_log2_effect': float(values.min()), 'maximum_log2_effect': float(values.max()),
                        'positive_labels': int((values > 0).sum()), 'negative_labels': int((values < 0).sum()),
                        'CPM_direction_agree_labels': int(df['sign_matches_CPM'].sum()),
                        'pool_direction_agree_labels': int(df['sign_matches_pooled'].sum())})
            for donor in donors:
                omit_rows.append({'gene_id': candidate.gene_id, 'symbol': candidate.symbol,
                                  'omitted_label': donor,
                                  'median_effect_other4': df.loc[df['donor_label'] != donor, 'log2_FFA_Control_mean_norm_plus1'].median()})
        summaries.append(row)
    summary = pd.DataFrame(summaries)
    summary.to_csv(HERE / 'phh_candidate_summary.csv', index=False)
    pd.DataFrame(omit_rows).to_csv(HERE / 'phh_label_omission.csv', index=False)
    qc = meta[['matrix_column', 'gsm', 'donor_id', 'condition']].copy()
    qc['assigned_counts'] = totals.astype(int)
    qc['median_ratio_factor'] = factors
    qc['detected_genes'] = (counts > 0).sum(axis=0)
    qc.to_csv(HERE / 'phh_library_qc.csv', index=False)
    pd.DataFrame({'donor_label': [x[0] for x in pooled_keys], 'condition': [x[1] for x in pooled_keys],
                  'pooled_factor': pool_factors, 'pooled_total': pooled_counts.sum(axis=0).astype(int)}).to_csv(HERE / 'phh_pooled_factors.csv', index=False)
    files = ['phh_candidate_samples.csv', 'phh_donor_effects.csv', 'phh_candidate_summary.csv',
             'phh_label_omission.csv', 'phh_library_qc.csv', 'phh_pooled_factors.csv']
    record = {'started_utc': started, 'completed_utc': datetime.now(timezone.utc).isoformat(),
              'python': sys.version, 'interpreter': sys.executable,
              'script_sha256': sha(Path(__file__)), 'protocol_sha256': sha(HERE / 'phh_protocol.json'),
              'input_sha256': sha(raw_path), 'manifest_sha256': sha(manifest_path), 'candidate_sha256': sha(candidate_path),
              'source_shape': list(counts.shape), 'donor_labels': donors, 'reference_genes_30_libraries': reference_n,
              'reference_genes_10_pooled_profiles': pool_ref_n,
              'matched_gene_ids': int(summary['feature_covered'].sum()),
              'missing_gene_ids': summary.loc[~summary['feature_covered'], 'gene_id'].tolist(),
              'sample_rows': len(samples), 'donor_effect_rows': len(effects),
              'spotlight': summary.loc[summary['symbol'].isin(['CXCL8', 'VNN1'])].to_dict('records'),
              'no_inference': True, 'no_sesamin_validation': True,
              'limitations': spec['claim_limits'],
              'artifact_hashes': {name: sha(HERE / name) for name in files}}
    (HERE / 'phh_summary.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'matched': record['matched_gene_ids'], 'missing': record['missing_gene_ids'], 'spotlight': record['spotlight']}))


if __name__ == '__main__':
    main()
