"""Independent checks from the deposited count file; does not import analysis code."""
import csv
import gzip
import hashlib
import itertools
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
SOURCE_STAGE = STAGE.parent / 'candidate_screen_03'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(x):
    eligible = np.min(x, axis=1) > 0
    reference = np.exp(np.sum(np.log(x[eligible]), axis=1) / x.shape[1])
    ratios = x[eligible] / reference[:, None]
    size = np.array([np.median(ratios[:, j]) for j in range(x.shape[1])])
    size = size / np.prod(size) ** (1 / x.shape[1])
    return x / size, size, int(eligible.sum())


def summarize(x, n):
    c, m, t = [x[:, i * n:(i + 1) * n].sum(axis=1) / n for i in range(3)]
    # Difference of logs is independently expressed versus the production ratio.
    dm = np.log2(m + 1) - np.log2(c + 1)
    dt = np.log2(t + 1) - np.log2(m + 1)
    residual = np.log2(t + 1) - np.log2(c + 1)
    interior = ((c < t) & (t < m)) | ((m < t) & (t < c))
    return np.array([c, m, t]).T, dm, dt, residual, interior


def main():
    spec = json.loads((STAGE / 'protocol.json').read_text(encoding='utf-8'))
    source = (STAGE / spec['input']).resolve()
    samples = spec['sample_order']
    summary = json.loads((STAGE / 'analysis_summary.json').read_text(encoding='utf-8'))
    freeze = json.loads((STAGE / 'protocol_freeze.json').read_text(encoding='utf-8'))
    checks, differences = {}, {}
    checks['input_hash'] = digest(source) == spec['input_sha256']
    checks['protocol_frozen_hash'] = digest(STAGE / 'protocol.json') == freeze['protocol_sha256']
    checks['analysis_code_hash'] = digest(STAGE / 'run_analysis.py') == summary['script_sha256']
    checks['analysis_output_hashes'] = all(digest(STAGE / x['file']) == x['sha256'] for x in summary['artifacts'])

    ids, symbols, values = [], [], []
    with gzip.open(source, 'rt', encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            ids.append(row['name'])
            symbols.append(row['GeneName'])
            values.append([int(row[sample]) for sample in samples])
    x = np.array(values, dtype=float)
    checks['raw_shape_and_keys'] = x.shape == (60605, 9) and len(set(ids)) == len(ids)
    manifest = pd.read_csv(SOURCE_STAGE / 'sesamin_sample_manifest.csv')
    checks['source_mapping'] = True
    for row in manifest.itertuples():
        soft = (SOURCE_STAGE / 'target_metadata' / (row.gsm + '.soft.txt')).read_text(encoding='utf-8')
        checks['source_mapping'] &= ('!Sample_description = ' + row.matrix_column) in soft
    checks['sample_order_mapping'] = set(samples) == set(manifest.matrix_column) and len(samples) == 9
    output = pd.read_csv(STAGE / 'all_gene_descriptive.csv.gz', keep_default_na=False)
    checks['output_gene_mapping'] = output.gene_id.tolist() == ids and output.symbol.tolist() == symbols

    def close(name, actual, expected, atol=1e-10, rtol=1e-10):
        actual, expected = np.array(actual), np.array(expected)
        differences[name] = float(np.max(np.abs(actual - expected)))
        checks[name] = bool(np.allclose(actual, expected, atol=atol, rtol=rtol))

    def exact(name, actual, expected):
        checks[name] = bool(np.array_equal(np.array(actual), np.array(expected)))

    totals = x.sum(axis=0)
    cpm = x * (1e6 / totals)
    keep = np.count_nonzero(cpm >= 1, axis=1) >= 3
    z, factors, nref = normalize(x)
    means, d, a, r, between = summarize(z, 3)
    _, dc, ac, rc, bc = summarize(cpm, 3)
    qc = pd.read_csv(STAGE / 'sample_qc.csv')
    exact('qc_sample_order', qc['sample'], samples)
    exact('library_totals', qc.assigned_counts, totals)
    exact('detected_genes', qc.detected_genes, np.count_nonzero(x, axis=0))
    close('size_factors', qc.median_ratio_factor, factors)
    exact('filter', output.retained_cpm_filter, keep)
    for name, expected in [('mean_norm_control', means[:, 0]), ('mean_norm_OA', means[:, 1]),
                           ('mean_norm_sesamin_OA', means[:, 2]), ('log2_OA_vs_control', d),
                           ('log2_sesamin_vs_OA', a), ('log2_sesamin_vs_control', r),
                           ('cpm_log2_OA_vs_control', dc), ('cpm_log2_sesamin_vs_OA', ac),
                           ('cpm_log2_sesamin_vs_control', rc)]:
        close(name, output[name], expected)
    exact('between_order', output.between_control_and_OA, between)

    threshold = np.log2(1.5)
    close('threshold_log2_1_5', spec['descriptive_priority_rule']['minimum_abs_log2_effect_each'], threshold)
    base = keep & (np.abs(d) >= threshold) & (np.abs(a) >= threshold) & (d * a < 0) & between
    agreement = (np.sign(d) == np.sign(dc)) & (np.sign(a) == np.sign(ac)) & bc
    exact('base_rule', output.base_descriptive_rule, base)
    exact('cpm_agreement', output.cpm_direction_and_ordering_agree, agreement)
    passes = np.zeros(len(x), dtype=int)
    ds, ats = [], []
    recorded = json.loads((STAGE / 'loo_normalization.json').read_text(encoding='utf-8'))
    for i, dropped in enumerate(itertools.product(range(3), range(3, 6), range(6, 9))):
        retained = [j for j in range(9) if j not in dropped]
        newz, sf, reference_n = normalize(x[:, retained])
        _, di, ai, _, bi = summarize(newz, 2)
        success = ((np.sign(di) == np.sign(d)) & (np.sign(ai) == np.sign(a)) & bi &
                   (np.abs(di) >= threshold) & (np.abs(ai) >= threshold))
        passes += success
        ds.append(di)
        ats.append(ai)
        checks['loo_' + str(i) + '_metadata'] = (recorded[i]['dropped'] == [samples[j] for j in dropped] and
                                               recorded[i]['normalization_reference_genes'] == reference_n)
        close('loo_' + str(i) + '_factors', recorded[i]['size_factors'], sf)
    exact('loo_pass_counts', output.loo_pass_n_of_27, passes)
    close('loo_fractions', output.loo_fraction, passes / 27)
    for name, expected in [('loo_d_min', np.min(ds, axis=0)), ('loo_d_max', np.max(ds, axis=0)),
                           ('loo_a_min', np.min(ats, axis=0)), ('loo_a_max', np.max(ats, axis=0))]:
        close(name, output[name], expected)
    selected = base & agreement & (passes >= 22)
    exact('priority_membership', output.descriptive_priority, selected)
    rankings = sorted(np.flatnonzero(selected), key=lambda i: (-min(abs(d[i]), abs(a[i])), ids[i]))
    priorities = pd.read_csv(STAGE / 'priority_all.csv', keep_default_na=False)
    exact('priority_order', priorities.gene_id, [ids[i] for i in rankings])
    top20 = pd.read_csv(STAGE / 'priority_top20.csv', keep_default_na=False)
    exact('top20_order', top20.gene_id, [ids[i] for i in rankings[:20]])

    log = np.log2(1 + z[keep])
    corr = np.corrcoef(log.T)
    close('sample_correlations', pd.read_csv(STAGE / 'sample_correlations.csv', index_col=0), corr)
    variance = ((log - log.mean(axis=1)[:, None]) ** 2).sum(axis=1) / 8
    indices = sorted(range(len(variance)), key=lambda i: (-variance[i], i))[:500]
    feature_ids = np.array(ids)[keep][indices]
    exact('pca_feature_ids', pd.read_csv(STAGE / 'pca_genes.csv').gene_id, feature_ids)
    centered = log[indices].T - log[indices].mean(axis=1)
    eigenvalues, eigenvectors = np.linalg.eigh(centered @ centered.T)
    ei = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.maximum(eigenvalues[ei], 0)
    close('pca_explained_variance', summary['pca_explained_variance'], eigenvalues / eigenvalues.sum())
    independently_projected = eigenvectors[:, ei[:2]] * np.sqrt(eigenvalues[:2])
    pca = pd.read_csv(STAGE / 'pca_scores.csv')[['PC1', 'PC2']].to_numpy()
    # Projection Gram matrices remove the arbitrary signs of eigenvectors.
    close('pca_two_component_gram', pca @ pca.T, independently_projected @ independently_projected.T, atol=1e-8)
    exact_duplicates = sum(np.array_equal(x[:, i], x[:, j]) for i, j in itertools.combinations(range(9), 2))
    checks['duplicate_pair_summary'] = exact_duplicates == summary['duplicate_sample_pairs']
    checks['summary_counts'] = (int(keep.sum()) == summary['filtered_gene_rows'] and
                               int(base.sum()) == summary['base_descriptive_rule_n'] and
                               int((base & agreement).sum()) == summary['base_and_CPM_agreement_n'] and
                               int(selected.sum()) == summary['priority_rule_n'] and
                               nref == summary['normalization_reference_genes'])

    result = {'executed_utc': datetime.now(timezone.utc).isoformat(), 'executable': sys.executable,
              'python': sys.version, 'numpy_version': np.__version__, 'pandas_version': pd.__version__,
              'review_script_sha256': digest(Path(__file__)), 'reviewed_analysis_sha256': digest(STAGE / 'run_analysis.py'),
              'reviewed_summary_sha256': digest(STAGE / 'analysis_summary.json'),
              'method': 'Independent raw TSV parsing; no import or execution of production functions. Full gene effect/filter/LOO recomputation; alternate PCA eigendecomposition.',
              'checks': checks, 'failed_checks': [k for k, value in checks.items() if not value],
              'max_absolute_differences': differences,
              'verified_counts': {'genes': len(x), 'libraries': 9, 'filtered': int(keep.sum()),
                                  'normalization_reference': nref, 'base_rule': int(base.sum()),
                                  'cpm_agree_among_base': int((base & agreement).sum()), 'priority': int(selected.sum()),
                                  'loo_combinations': 27, 'minimum_successes_for_80_percent': 22},
              'scope_limits': ['Computational checks do not establish independent cultures or exposure-time matching.',
                               'Shared Model group induces negative dependence of successive contrasts; LOO does not remove it.',
                               'Priority is a heuristic expression-pattern list, not significant genes, an effect estimate with CI, recovery evidence, novelty or independent validation.',
                               'No biological annotation audit, novel mechanism test, raw-read alignment audit or AI nomination experiment performed.',
                               'Figures and final narrative not yet inspected by this numerical reviewer.'],
              'design_notes': ['Main factors use all input genes positive in all nine libraries; each LOO refits on all input genes positive in six retained libraries. The CPM feature filter alone stays fixed.',
                               'Primary pseudocount is one normalized count; CPM output pseudocount is one CPM. CPM direction/order check does not require effect magnitude agreement.',
                               'No p or q values are computed. Full output retains low-expression rows with eligibility flags.'],
              'overall_assessment': 'Share with caveats' if all(checks.values()) else 'Needs revision'}
    (HERE / 'review.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({'checks': len(checks), 'failed': result['failed_checks'], 'counts': result['verified_counts']}))


if __name__ == '__main__':
    main()
