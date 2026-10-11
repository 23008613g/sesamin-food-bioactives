"""Descriptive exploration of deposited GSE248903 libraries, never inference."""
import hashlib
import itertools
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def factors(counts):
    positive = (counts > 0).all(axis=1)
    if positive.sum() < 100:
        raise ValueError('Too few positive reference genes for median-ratio normalization')
    refs = counts[positive]
    geometric = np.exp(np.log(refs).mean(axis=1))
    sf = np.median(refs / geometric[:, None], axis=0)
    sf /= np.exp(np.log(sf).mean())
    return sf, int(positive.sum())


def effects(norm, groups):
    means = np.column_stack([norm[:, idx].mean(axis=1) for idx in groups])
    c, m, t = means.T
    d = np.log2((m + 1) / (c + 1))
    a = np.log2((t + 1) / (m + 1))
    r = np.log2((t + 1) / (c + 1))
    between = (t > np.minimum(c, m)) & (t < np.maximum(c, m))
    return means, d, a, r, between


def main():
    started = datetime.now(timezone.utc).isoformat()
    protocol_path = HERE / 'protocol.json'
    spec = json.loads(protocol_path.read_text(encoding='utf-8'))
    freeze = json.loads((HERE / 'protocol_freeze.json').read_text())
    assert sha(protocol_path) == freeze['protocol_sha256']
    source = (HERE / spec['input']).resolve()
    assert sha(source) == spec['input_sha256']
    order = spec['sample_order']
    table = pd.read_csv(source, sep='\t', usecols=['name', 'GeneName'] + order,
                        dtype={'name': str, 'GeneName': str}, keep_default_na=False)
    counts = table[order].to_numpy(dtype=float)
    assert counts.shape == (60605, 9) and table['name'].is_unique
    assert np.isfinite(counts).all() and (counts >= 0).all() and (counts % 1 == 0).all()
    totals = counts.sum(axis=0)
    cpm = counts / totals[None, :] * 1e6
    keep = (cpm >= 1).sum(axis=1) >= 3
    sf, nref = factors(counts)
    norm = counts / sf
    groups = [np.arange(0, 3), np.arange(3, 6), np.arange(6, 9)]
    means, d, a, residual, between = effects(norm, groups)
    _, dcpm, acpm, rcpm, bcpm = effects(cpm, groups)
    threshold = spec['descriptive_priority_rule']['minimum_abs_log2_effect_each']
    base = keep & (abs(d) >= threshold) & (abs(a) >= threshold) & (d * a < 0) & between
    norm_agree = (np.sign(d) == np.sign(dcpm)) & (np.sign(a) == np.sign(acpm)) & bcpm

    loo_pass = np.zeros(len(table), dtype=int)
    loo_d, loo_a = [], []
    perturbations = []
    for dropped in itertools.product(*groups):
        retained = [j for j in range(9) if j not in dropped]
        active = counts[:, retained]
        new_sf, ref_n = factors(active)
        _, di, ai, _, bi = effects(active / new_sf, [np.arange(0, 2), np.arange(2, 4), np.arange(4, 6)])
        success = (np.sign(di) == np.sign(d)) & (np.sign(ai) == np.sign(a)) & bi
        success &= (abs(di) >= threshold) & (abs(ai) >= threshold)
        loo_pass += success
        loo_d.append(di)
        loo_a.append(ai)
        perturbations.append({'dropped': [order[i] for i in dropped],
                              'retained': [order[i] for i in retained],
                              'normalization_reference_genes': ref_n,
                              'size_factors': new_sf.tolist()})
    fraction = loo_pass / len(perturbations)
    priority = base & norm_agree & (fraction >= 0.8)
    output = pd.DataFrame({'gene_id': table['name'], 'symbol': table['GeneName'],
                           'retained_cpm_filter': keep, 'cpm_ge1_libraries': (cpm >= 1).sum(axis=1),
                           'mean_norm_control': means[:, 0], 'mean_norm_OA': means[:, 1],
                           'mean_norm_sesamin_OA': means[:, 2], 'log2_OA_vs_control': d,
                           'log2_sesamin_vs_OA': a, 'log2_sesamin_vs_control': residual,
                           'cpm_log2_OA_vs_control': dcpm, 'cpm_log2_sesamin_vs_OA': acpm,
                           'cpm_log2_sesamin_vs_control': rcpm,
                           'between_control_and_OA': between, 'base_descriptive_rule': base,
                           'cpm_direction_and_ordering_agree': norm_agree,
                           'loo_pass_n_of_27': loo_pass, 'loo_fraction': fraction,
                           'loo_d_min': np.min(loo_d, axis=0), 'loo_d_max': np.max(loo_d, axis=0),
                           'loo_a_min': np.min(loo_a, axis=0), 'loo_a_max': np.max(loo_a, axis=0),
                           'descriptive_priority': priority,
                           'minimum_abs_primary_effect': np.minimum(abs(d), abs(a)),
                           'prespecified_known_panel': table['GeneName'].isin(spec['known_claim_panel'])})
    for name, idx in zip(['control', 'OA', 'sesamin_OA'], groups):
        output['sd_norm_' + name] = norm[:, idx].std(axis=1, ddof=1)
        output['median_norm_' + name] = np.median(norm[:, idx], axis=1)
    output.to_csv(HERE / 'all_gene_descriptive.csv.gz', index=False, compression='gzip')
    candidates = output.loc[priority].sort_values(['minimum_abs_primary_effect', 'gene_id'], ascending=[False, True])
    candidates.to_csv(HERE / 'priority_all.csv', index=False)
    candidates.head(20).to_csv(HERE / 'priority_top20.csv', index=False)
    output.loc[output['prespecified_known_panel']].to_csv(HERE / 'known_panel_summary.csv', index=False)

    rows = []
    for i, rec in table.loc[table['GeneName'].isin(spec['known_claim_panel'])].iterrows():
        for j, sample in enumerate(order):
            rows.append({'gene_id': rec['name'], 'symbol': rec['GeneName'], 'sample': sample,
                         'group': ['Control', 'OA', 'OA + sesamin'][j // 3],
                         'raw_count': int(counts[i, j]), 'normalized_count': norm[i, j],
                         'log2_normalized_plus1': np.log2(norm[i, j] + 1)})
    pd.DataFrame(rows).to_csv(HERE / 'known_panel_samples.csv', index=False)
    qc = pd.DataFrame({'sample': order, 'group': np.repeat(['Control', 'OA', 'OA + sesamin'], 3),
                       'assigned_counts': totals.astype(int), 'detected_genes': (counts > 0).sum(axis=0),
                       'median_ratio_factor': sf,
                       'top10_gene_fraction': np.sort(counts, axis=0)[-10:].sum(axis=0) / totals})
    qc.to_csv(HERE / 'sample_qc.csv', index=False)
    log_norm = np.log2(norm[keep] + 1)
    corr = np.corrcoef(log_norm.T)
    pd.DataFrame(corr, index=order, columns=order).to_csv(HERE / 'sample_correlations.csv')
    variance = log_norm.var(axis=1, ddof=1)
    top = np.argsort(-variance, kind='stable')[:min(500, len(variance))]
    selected = log_norm[top].T
    centered = selected - selected.mean(axis=0)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    scores = u[:, :2] * s[:2]
    explained = s ** 2 / (s ** 2).sum()
    pca = qc[['sample', 'group']].copy()
    pca['PC1'], pca['PC2'] = scores[:, 0], scores[:, 1]
    pca.to_csv(HERE / 'pca_scores.csv', index=False)
    pd.DataFrame({'gene_id': table.loc[keep, 'name'].to_numpy()[top],
                  'log_normalized_variance': variance[top]}).to_csv(HERE / 'pca_genes.csv', index=False)
    pairs = []
    for i, j in itertools.combinations(range(9), 2):
        pairs.append({'sample1': order[i], 'sample2': order[j],
                      'exact_duplicate': bool(np.array_equal(counts[:, i], counts[:, j])),
                      'log_normalized_pearson': float(corr[i, j])})
    pd.DataFrame(pairs).to_csv(HERE / 'sample_pair_checks.csv', index=False)
    (HERE / 'loo_normalization.json').write_text(json.dumps(perturbations, indent=2) + '\n')
    artifacts = ['all_gene_descriptive.csv.gz', 'priority_all.csv', 'priority_top20.csv',
                 'known_panel_summary.csv', 'known_panel_samples.csv', 'sample_qc.csv',
                 'sample_correlations.csv', 'pca_scores.csv', 'pca_genes.csv',
                 'sample_pair_checks.csv', 'loo_normalization.json']
    summary = {'started_utc': started, 'completed_utc': datetime.now(timezone.utc).isoformat(),
               'interpreter': sys.executable, 'script_sha256': sha(Path(__file__)),
               'protocol_sha256': sha(protocol_path), 'input_sha256': sha(source),
               'gene_rows': len(table), 'filtered_gene_rows': int(keep.sum()),
               'excluded_low_expression_rows': int((~keep).sum()),
               'normalization_reference_genes': nref,
               'sample_n': 9, 'sample_exclusions': 0, 'duplicate_sample_pairs': sum(p['exact_duplicate'] for p in pairs),
               'correlation_range_off_diagonal': [float(corr[np.triu_indices(9, 1)].min()), float(corr[np.triu_indices(9, 1)].max())],
               'pca_explained_variance': explained.tolist(),
               'base_descriptive_rule_n': int(base.sum()),
               'base_and_CPM_agreement_n': int((base & norm_agree).sum()),
               'priority_rule_n': int(priority.sum()),
               'prespecified_panel_found': sorted(set(output.loc[output['prespecified_known_panel'], 'symbol'])),
               'prespecified_panel_missing': sorted(set(spec['known_claim_panel']) - set(table['GeneName'])),
               'panel_gene_id_counts': output.loc[output['prespecified_known_panel']].groupby('symbol').size().to_dict(),
               'priority_top20': candidates.head(20)[['gene_id', 'symbol', 'log2_OA_vs_control', 'log2_sesamin_vs_OA', 'loo_fraction']].to_dict('records'),
               'inference': 'No p-values or FDR. Priority is descriptive and provisional. LOO is a stability perturbation, not an independent sample or CI.',
               'actual_AI_nomination_calls': 0,
               'artifacts': [{'file': f, 'sha256': sha(HERE / f)} for f in artifacts]}
    (HERE / 'analysis_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({k: summary[k] for k in ['filtered_gene_rows', 'base_descriptive_rule_n',
                                             'base_and_CPM_agreement_n', 'priority_rule_n', 'duplicate_sample_pairs']}, ensure_ascii=True))


if __name__ == '__main__':
    main()
