"""Stage 15 bounded post-selection sensitivity; no inferential statistics."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import itertools
import json
import platform
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HERE = Path(sys.argv[1]).resolve()
INPUT = ROOT / 'inputs'
MOUSE = INPUT / 'mouse'
GROUPS = [list(range(0, 3)), list(range(3, 6)), list(range(6, 9))]
CONTRASTS = [('model_vs_control', 1, 0), ('treatment_vs_model', 2, 1), ('treatment_vs_control', 2, 0)]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def write_json(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')


def factors(a):
    good = (a > 0).all(axis=1)
    assert good.sum() > 100
    ref = a[good]
    geometric = np.exp(np.log(ref).mean(axis=1))
    sf = np.median(ref / geometric[:, None], axis=0)
    sf /= np.exp(np.log(sf).mean())
    return sf, int(good.sum())


def subset_specs():
    choices = [[group] + [[j for j in group if j != dropped] for dropped in group] for group in GROUPS]
    return list(itertools.product(*choices))


def qc(species, ids, raw, norm, keep, names, labels, normalization, unit):
    transformed = np.log2(norm[keep] + 1)
    centered = transformed.T - transformed.mean(axis=1)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    for k in range(len(s)):
        anchor = np.argmax(np.abs(vt[k]))
        if vt[k, anchor] < 0:
            vt[k] *= -1
            u[:, k] *= -1
    scores = u * s
    explained = s**2 / np.sum(s**2)
    scoretable = pd.DataFrame({'library': names, 'group': labels})
    for k in range(8):
        scoretable[f'PC{k+1}'] = scores[:, k]
    scoretable.to_csv(HERE / f'{species}_pca_scores.csv', index=False)
    pd.DataFrame({'component': [f'PC{k+1}' for k in range(9)], 'singular_value': s,
                  'explained_variance_ratio': explained}).to_csv(HERE / f'{species}_pca_variance.csv', index=False)
    loads = pd.DataFrame(vt[:8].T, columns=[f'PC{k+1}_loading' for k in range(8)])
    loads.insert(0, 'gene_id', np.asarray(ids)[keep])
    loads.to_csv(HERE / f'{species}_pca_loadings.csv.gz', index=False)
    pd.DataFrame({'gene_id': ids, 'expression_filter_pass': keep}).to_csv(HERE / f'{species}_qc_filter.csv.gz', index=False)
    corr = np.corrcoef(transformed.T)
    spear = pd.DataFrame(transformed, columns=names).corr(method='spearman')
    pd.DataFrame(corr, index=names, columns=names).to_csv(HERE / f'{species}_sample_pearson.csv')
    spear.to_csv(HERE / f'{species}_sample_spearman.csv')
    samples = pd.DataFrame({'library': names, 'group': labels, 'assigned_count_sum': raw.sum(axis=0),
                            'detected_genes_raw_positive': (raw > 0).sum(axis=0),
                            'qc_retained_genes': int(keep.sum()),
                            'top10_gene_count_fraction': np.sort(raw, axis=0)[-10:].sum(axis=0) / raw.sum(axis=0)})
    samples.to_csv(HERE / f'{species}_sample_qc.csv', index=False)
    pairs = []
    for i, j in itertools.combinations(range(9), 2):
        pairs.append({'library1': names[i], 'library2': names[j], 'same_submitted_group': labels[i] == labels[j],
                      'exact_duplicate_profile': bool(np.array_equal(raw[:, i], raw[:, j])),
                      'pearson': corr[i, j], 'spearman': spear.iloc[i, j]})
    pd.DataFrame(pairs).to_csv(HERE / f'{species}_sample_pair_qc.csv', index=False)
    # Reconstruction and centering are numerical checks, not evidence of biology.
    assert np.allclose((u * s) @ vt, centered, atol=1e-10)
    assert np.allclose(scores.mean(axis=0), 0, atol=1e-10)
    return {'gene_rows': len(ids), 'retained_genes': int(keep.sum()), 'normalization': normalization,
            'unit': unit, 'PC1_PC2_explained': explained[:2].tolist(),
            'correlation_range': [float(corr[np.triu_indices(9, 1)].min()), float(corr[np.triu_indices(9, 1)].max())],
            'duplicate_profiles': sum(p['exact_duplicate_profile'] for p in pairs), 'excluded_libraries': 0}


def main():
    started = datetime.now(timezone.utc).isoformat()
    freeze_path = HERE / 'protocol_freeze.json'
    if freeze_path.exists():
        freeze = json.loads(freeze_path.read_text(encoding='utf-8'))
        assert freeze['sha256'] == sha(HERE / 'protocol.json'), 'Protocol changed after freeze'
    else:
        write_json('protocol_freeze.json', {'utc': started, 'sha256': sha(HERE / 'protocol.json'),
                                          'stage': 'before input expression matrices were loaded by this script'})
    sources = {
        'human_counts': INPUT / 'human/counts.txt.gz',
        'human_candidates': INPUT / 'human/priority_all.csv',
        'human_prior_qc': INPUT / 'human/sample_qc.csv',
        'human_prior_loo': INPUT / 'human/loo_normalization.json',
        'mouse_tpm': MOUSE / 'gene_tpm.tsv.gz',
        'mouse_estimated_counts': MOUSE / 'gene_est_counts.tsv.gz',
        'mouse_candidates': MOUSE / 'candidate_mouse_descriptive.tsv',
        'mouse_libraries': MOUSE / 'library_manifest.tsv',
    }
    manifest = [{'key': key, 'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size} for key, p in sources.items()]
    oldhuman = pd.read_csv(sources['human_candidates'])
    oldmouse = pd.read_csv(sources['mouse_candidates'], sep='\t')
    primary = oldmouse.loc[oldmouse['primary_comparison_eligible'].astype(str).str.lower().eq('true')].copy()
    assert len(primary) == 15 and primary.original_human_id.is_unique and primary.mouse_gene_id.is_unique
    primary = primary[['original_human_id', 'human_symbol', 'mouse_gene_id', 'mouse_symbol', 'all_libraries_TPM_lt1']]
    primary.to_csv(HERE / 'frozen_primary_candidates.csv', index=False)
    hn = [f'{g}_{i}' for g in ['Control', 'Model', 'Treat'] for i in range(1, 4)]
    ml = pd.read_csv(sources['mouse_libraries'], sep='\t')
    mn = ml.run.tolist()
    assert mn == [f'DRR{x}' for x in range(1010532, 1010541)]
    ht = pd.read_csv(sources['human_counts'], sep='\t', usecols=['name', 'GeneName'] + hn)
    mt = pd.read_csv(sources['mouse_tpm'], sep='\t')
    mc = pd.read_csv(sources['mouse_estimated_counts'], sep='\t')
    assert ht.name.is_unique and mt.gene_id.is_unique and mc.gene_id.is_unique
    assert mt.gene_id.equals(mc.gene_id)
    h = ht[hn].to_numpy(float)
    tpm = mt[mn].to_numpy(float)
    counts = mc[mn].to_numpy(float)
    for a in [h, tpm, counts]:
        assert np.isfinite(a).all() and (a >= 0).all()
    assert (h % 1 == 0).all() and h.shape == (60605, 9)
    hsf, hnref = factors(h)
    msf, mnref = factors(counts)
    normh = h / hsf
    normm = counts / msf
    cpm = h / h.sum(axis=0) * 1e6
    hi = ht.set_index('name').index.get_indexer(primary.original_human_id)
    mi = mt.set_index('gene_id').index.get_indexer(primary.mouse_gene_id)
    assert (hi >= 0).all() and (mi >= 0).all()
    assert np.allclose(hsf, pd.read_csv(sources['human_prior_qc']).median_ratio_factor, rtol=1e-12)
    assert np.allclose(tpm[mi], oldmouse.set_index('original_human_id').loc[primary.original_human_id,
                       [x + '_TPM' for x in mn]].to_numpy(float), rtol=1e-12)
    candidate_samples = []
    for k, r in enumerate(primary.itertuples(index=False)):
        for j in range(9):
            candidate_samples.extend([
                {'species': 'human', 'human_symbol': r.human_symbol, 'gene_id': r.original_human_id,
                 'library': hn[j], 'group': ['Control', 'Model', 'Treat'][j // 3],
                 'raw_or_estimated_count': h[hi[k], j], 'median_ratio_count': normh[hi[k], j],
                 'CPM': cpm[hi[k], j], 'TPM': np.nan},
                {'species': 'mouse', 'human_symbol': r.human_symbol, 'gene_id': r.mouse_gene_id,
                 'library': mn[j], 'group': ml.provisional_group.iloc[j],
                 'raw_or_estimated_count': counts[mi[k], j], 'median_ratio_count': normm[mi[k], j],
                 'CPM': np.nan, 'TPM': tpm[mi[k], j]}])
    pd.DataFrame(candidate_samples).to_csv(HERE / 'candidate_library_values.csv', index=False)
    configs = subset_specs()
    assert len(configs) == 64 and len({tuple(tuple(g) for g in config) for config in configs}) == 64
    rows, cfgrows, sfrows = [], [], []
    for cid, config in enumerate(configs):
        active = sum(config, [])
        hs, hr = factors(h[:, active])
        ms, mr = factors(counts[:, active])
        local = [[active.index(j) for j in group] for group in config]
        cfgrows.append({'configuration': cid, 'group_sizes': '/'.join(str(len(g)) for g in config),
                        'retained_indices': '|'.join(','.join(map(str, g)) for g in config),
                        'human_retained_libraries': ';'.join(hn[j] for j in active),
                        'mouse_retained_libraries': ';'.join(mn[j] for j in active),
                        'human_factor_reference_genes': hr, 'mouse_factor_reference_genes': mr})
        for q, j in enumerate(active):
            sfrows += [{'configuration': cid, 'species': 'human', 'library': hn[j], 'factor': hs[q]},
                       {'configuration': cid, 'species': 'mouse', 'library': mn[j], 'factor': ms[q]}]
        metrics = [
            ('human', 'median_ratio_reestimated', h[hi][:, active] / hs, 1.),
            ('human', 'median_ratio_frozen', normh[hi][:, active], 1.),
            ('human', 'CPM', cpm[hi][:, active], 1.),
            ('mouse', 'TPM', tpm[mi][:, active], .1),
            ('mouse', 'estimated_count_median_ratio', counts[mi][:, active] / ms, 1.),
        ]
        for species, metric, array, pc in metrics:
            for aggregation in ['mean', 'median']:
                centers = np.column_stack([getattr(np, aggregation)(array[:, g], axis=1) for g in local])
                for k, gene in enumerate(primary.human_symbol):
                    c, m, t = centers[k]
                    for contrast, numerator, denominator in CONTRASTS:
                        a, b = centers[k, numerator], centers[k, denominator]
                        logratio = np.log2((a + pc) / (b + pc))
                        rows.append({'configuration': cid, 'species': species, 'metric': metric,
                                     'aggregation': aggregation, 'human_symbol': gene, 'contrast': contrast,
                                     'control_summary': c, 'model_summary': m, 'treatment_summary': t,
                                     'log2_ratio': logratio, 'direction': int(np.sign(a-b)), 'pseudocount': pc,
                                     'percent_difference_no_pseudocount': (a-b) / b * 100 if b > 0 else np.nan})
    pd.DataFrame(cfgrows).to_csv(HERE / 'enumerated_configurations.csv', index=False)
    pd.DataFrame(sfrows).to_csv(HERE / 'configuration_size_factors.csv', index=False)
    results = pd.DataFrame(rows)
    results.to_csv(HERE / 'all_configuration_effects.csv.gz', index=False)
    baseline = results.loc[results.configuration == 0].copy()
    baseline.to_csv(HERE / 'all_library_mean_median_effects.csv', index=False)
    summaryrows = []
    keys = ['species', 'metric', 'aggregation', 'human_symbol', 'contrast']
    for key, d in results.groupby(keys, sort=False):
        base = d.loc[d.configuration == 0].iloc[0]
        summaryrows.append(dict(zip(keys, key)) | {
            'full_data_log2_ratio': base.log2_ratio, 'full_data_direction': int(base.direction),
            'full_data_percent_difference': base.percent_difference_no_pseudocount,
            'configurations': len(d), 'baseline_direction_retained_n': int(d.direction.eq(base.direction).sum()),
            'positive_n': int(d.direction.gt(0).sum()), 'negative_n': int(d.direction.lt(0).sum()),
            'unchanged_n': int(d.direction.eq(0).sum()), 'min_log2_ratio': d.log2_ratio.min(),
            'max_log2_ratio': d.log2_ratio.max(),
            'min_percent_difference': d.percent_difference_no_pseudocount.min(),
            'max_percent_difference': d.percent_difference_no_pseudocount.max()})
    sums = pd.DataFrame(summaryrows)
    sums.to_csv(HERE / 'configuration_effect_summary.csv', index=False)
    cross = []
    for aggregation in ['mean', 'median']:
        for gene in primary.human_symbol:
            d = results.loc[(results.human_symbol == gene) & (results.aggregation == aggregation)]
            hd = d.loc[(d.species == 'human') & (d.metric == 'median_ratio_reestimated')].pivot(index='configuration', columns='contrast', values='direction')
            md = d.loc[(d.species == 'mouse') & (d.metric == 'TPM')].pivot(index='configuration', columns='contrast', values='direction')
            joint = []
            rec = {'human_symbol': gene, 'aggregation': aggregation, 'configuration_pairs': 4096}
            for contrast, _, _ in CONTRASTS:
                ha, ma = hd[contrast].to_numpy(), md[contrast].to_numpy()
                agree = (ha[:, None] == ma[None, :]) & (ha[:, None] != 0)
                baseline_agree = int(ha[0] == ma[0] and ha[0] != 0)
                baseline_sign_preserved = agree & (ha[:, None] == ha[0]) & bool(baseline_agree)
                rec.update({contrast + '_baseline_same_nonzero': bool(baseline_agree),
                            contrast + '_concordant_n': int(agree.sum()),
                            contrast + '_baseline_sign_preserved_n': int(baseline_sign_preserved.sum())})
                if contrast != 'treatment_vs_control':
                    joint.append(agree)
            both = joint[0] & joint[1]
            rec['both_model_and_treatment_concordant_n'] = int(both.sum())
            rec['both_concordant_configuration_share'] = float(both.mean())
            rec['both_baseline_same_nonzero'] = rec['model_vs_control_baseline_same_nonzero'] and rec['treatment_vs_model_baseline_same_nonzero']
            cross.append(rec)
    crossdf = pd.DataFrame(cross)
    crossdf.to_csv(HERE / 'cross_species_configuration_consensus.csv', index=False)
    # Independent direct baseline reproductions, and old 27 factor sets.
    hp = baseline.loc[(baseline.species == 'human') & (baseline.metric == 'median_ratio_reestimated') & (baseline.aggregation == 'mean')]
    oldmap = oldhuman.set_index('symbol')
    for r in hp.itertuples(index=False):
        oldcol = {'model_vs_control': 'log2_OA_vs_control', 'treatment_vs_model': 'log2_sesamin_vs_OA',
                  'treatment_vs_control': 'log2_sesamin_vs_control'}[r.contrast]
        assert np.isclose(r.log2_ratio, oldmap.loc[r.human_symbol, oldcol], rtol=1e-10, atol=1e-10)
    sfs = pd.DataFrame(sfrows)
    oldloo = json.loads(sources['human_prior_loo'].read_text(encoding='utf-8'))
    for old in oldloo:
        cid = next(r['configuration'] for r in cfgrows if r['human_retained_libraries'].split(';') == old['retained'])
        cur = sfs.loc[(sfs.configuration == cid) & (sfs.species == 'human'), 'factor'].to_numpy()
        assert np.allclose(cur, old['size_factors'], rtol=1e-12, atol=1e-12)
    hqc = qc('human', ht.name.to_numpy(), h, normh, (cpm >= 1).sum(axis=1) >= 3, hn,
             ['Control']*3 + ['Model']*3 + ['Treat']*3, 'Full-data median-ratio count', 'Raw counts')
    mqc = qc('mouse', mt.gene_id.to_numpy(), counts, tpm, (tpm >= 1).sum(axis=1) >= 3, mn,
             ml.provisional_group.tolist(), 'TPM', 'Summed transcript estimated counts; fractional')
    assert hqc['retained_genes'] == 13014
    means = crossdf.loc[crossdf.aggregation == 'mean']
    assert set(means.loc[means.both_baseline_same_nonzero, 'human_symbol']) == {'VNN1', 'AREG', 'JAG2'}
    assert len(means.loc[means.treatment_vs_model_baseline_same_nonzero]) == 8
    for rec in manifest:
        assert sha(rec['path']) == rec['sha256']
    outputs = [{'file': p.name, 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in HERE.iterdir()
               if p.is_file() and p.suffix in ['.csv', '.gz']]
    conclusion = {
        'full_data_both_concordant': means.loc[means.both_baseline_same_nonzero, 'human_symbol'].tolist(),
        'all4096_both_concordant': means.loc[means.both_model_and_treatment_concordant_n.eq(4096), 'human_symbol'].tolist(),
        'all4096_treatment_concordant': means.loc[means.treatment_vs_model_concordant_n.eq(4096), 'human_symbol'].tolist(),
        'all4096_treatment_baseline_preserved': means.loc[means.treatment_vs_model_baseline_sign_preserved_n.eq(4096), 'human_symbol'].tolist(),
        'spotlight': means.loc[means.human_symbol.isin(['VNN1', 'AREG', 'JAG2'])].to_dict('records')}
    write_json('analysis_summary.json', {'status': 'PASS_DESCRIPTIVE_EXECUTION', 'started_utc': started,
        'completed_utc': datetime.now(timezone.utc).isoformat(), 'python': sys.version,
        'executable': sys.executable, 'platform': platform.platform(), 'numpy': np.__version__, 'pandas': pd.__version__,
        'script_sha256': sha(__file__), 'protocol_sha256': sha(HERE / 'protocol.json'),
        'invocation': [sys.executable, str(Path(__file__).resolve()), str(HERE)],
        'seed': 'Not applicable; complete deterministic enumeration and exact SVD',
        'inputs': manifest, 'outputs': outputs, 'candidate_count': 15,
        'configuration_count_per_species': 64, 'cross_species_configuration_pairs': 4096,
        'self_checks': {'human_baseline_all45_contrasts_reproduced': True,
                        'old27_human_factor_sets_reproduced': True, 'mouse135_TPM_values_reproduced': True,
                        'pca_reconstruction': True, 'input_hashes_unchanged': True},
        'human_qc': hqc, 'mouse_qc': mqc, 'conclusions': conclusion,
        'limits': ['Configuration shares are enumeration fractions, not probabilities or CIs.',
                   'No independent held-out validation; human candidates selected on these data.',
                   'Independent biological units unresolved; no inferential testing.',
                   'No sample removal. Mouse groups provisional; no matched-treatment causal comparison.',
                   'Mouse normalized estimated counts are method sensitivity, not tximport length-adjusted DE.']})
    print(json.dumps({'status': 'PASS', 'conclusions': conclusion, 'human_qc': hqc, 'mouse_qc': mqc}, indent=2))


if __name__ == '__main__':
    main()
