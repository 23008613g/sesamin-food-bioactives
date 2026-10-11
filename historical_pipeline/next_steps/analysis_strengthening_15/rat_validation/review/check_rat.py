"""Bounded independent recomputation without importing production functions."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import io
import itertools
import json
import math
import statistics
import sys

import numpy as np
import pandas as pd
import scipy
from scipy import stats

HERE = Path(__file__).resolve().parent
R = HERE.parent
N = R.parents[1]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def get_table(path, table_begin, table_end):
    lines = gzip.open(path, 'rt', encoding='utf-8').read().splitlines()
    start = lines.index(table_begin) + 1
    end = lines.index(table_end)
    return pd.read_csv(io.StringIO('\n'.join(lines[start:end])), sep='\t', dtype=str, keep_default_na=False), lines[:start-1]


def main():
    checklog = []
    def check(name, result, detail=None):
        checklog.append({'name': name, 'pass': bool(result), 'detail': detail})
        if not result:
            raise AssertionError(name)
    def close(name, actual, expected):
        aa, bb = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
        check(name, np.allclose(aa, bb, atol=2e-9, rtol=2e-10, equal_nan=True),
              {'values': int(aa.size), 'max_abs_difference': float(np.nanmax(np.abs(aa-bb))) if aa.size else None})
    summ = json.loads((R / 'analysis_summary.json').read_text(encoding='utf-8'))
    all_input_hashes = []
    for record in summ['inputs']:
        check('initial input hash ' + Path(record['path']).name, sha(record['path']) == record['sha256'])
        all_input_hashes.append(record)
    original_outputs = {f['file']: sha(R / f['file']) for f in summ['output_files']}
    for record in summ['output_files']:
        check('initial output hash ' + record['file'], original_outputs[record['file']] == record['sha256'])
    matrixpath = Path(summ['inputs'][0]['path'])
    matrix, metadata = get_table(matrixpath, '!series_matrix_table_begin', '!series_matrix_table_end')
    annotation, ametadata = get_table(R / 'sources/GPL1355.annot.gz', '!platform_table_begin', '!platform_table_end')
    annotation = annotation.set_index('ID')
    matrix = matrix.set_index('ID_REF').astype(float)
    assert annotation.index.is_unique and matrix.index.is_unique
    animals = [f'GSM{i}' for i in range(1048036, 1048056)]
    check('matrix exact 31099 by20 and expected order', matrix.shape == (31099, 20) and matrix.columns.tolist() == animals)
    check('unique probes and platform coverage', matrix.index.is_unique and set(matrix.index) == set(annotation.index))
    metadata_rows = {}
    for line in metadata:
        if not line.strip():
            continue
        row = next(csv.reader([line], delimiter='\t'))
        metadata_rows.setdefault(row[0], []).append(row[1:])
    titles = metadata_rows['!Sample_title'][0]
    check('primary5+5 exact exposures, othercompound only remaining10',
          all(x.startswith('liver_control diet_') for x in titles[:5]) and
          all(x.startswith('liver_sesamin_') for x in titles[5:10]) and
          all(x.startswith('liver_episesamin_') for x in titles[10:15]) and
          all(x.startswith('liver_sesamolin_') for x in titles[15:]))
    design = metadata_rows['!Series_overall_design'][0][0]
    check('source specifies original7 selected5 and weight exclusions',
          all(x in design for x in ['7 animals', '5 rats', 'highest and lowest body weights']))
    growth = metadata_rows['!Sample_growth_protocol_ch1'][0][0]
    processing = metadata_rows['!Sample_data_processing'][0][0]
    check('distinct individually housed animals', 'housed individually' in growth)
    check('deposited processed global scaling not RMA', 'global scaling' in processing and 'GeneChip Operating Software' in processing)
    candidates = pd.read_csv(N / 'publication_package_14/data/human_candidates_all25.csv', dtype=str)
    ortho = pd.read_csv(R / 'sources/human_rat_orthology_ensembl116.tsv', sep='\t', dtype=str, keep_default_na=False)
    xr = pd.read_csv(R / 'sources/rat_entrez_ensembl116.tsv', sep='\t', dtype=str, keep_default_na=False)
    outputmap = pd.read_csv(R / 'outputs/candidate_mapping_all25.csv', dtype=str, keep_default_na=False).set_index('human_id')
    # Build independent dictionaries keyed by original annotation fields.
    human_targets = {}
    for row in ortho.itertuples(index=False, name=None):
        if row[2]:
            human_targets.setdefault(row[0], []).append(row)
    ratxref = {}
    for ratid, symbol, entrez in xr.itertuples(index=False, name=None):
        if entrez:
            ratxref.setdefault(ratid, set()).add(entrez)
    unique_probes, ambiguous_probes = {}, {}
    global_ambiguous = 0
    for probe, row in annotation.iterrows():
        identifiers = {v.strip() for v in row['Gene ID'].split('///') if v.strip() not in ['', '---']}
        if len(identifiers) == 1:
            unique_probes.setdefault(next(iter(identifiers)), []).append(probe)
        elif len(identifiers) > 1:
            global_ambiguous += 1
            for identifier in identifiers:
                ambiguous_probes.setdefault(identifier, []).append((probe, identifiers))
    mapping = {}
    exclusions = []
    for candidate in candidates.itertuples(index=False):
        target = human_targets.get(candidate.gene_id, [])
        if len(target) != 1 or target[0][4] != 'ortholog_one2one':
            status = 'no_rat_orthologue' if len(target) == 0 else 'non_one_to_one_orthology'
            check('unavailable mapping ' + candidate.symbol, outputmap.loc[candidate.gene_id, 'status'] == status)
            continue
        t = target[0]
        identifiers = ratxref.get(t[2], set())
        check('unique rat Entrez ' + candidate.symbol, len(identifiers) == 1)
        eid = next(iter(identifiers))
        probes = sorted(unique_probes.get(eid, []))
        row = outputmap.loc[candidate.gene_id]
        check('mapping fields ' + candidate.symbol,
              row.rat_id == t[2] and row.rat_symbol == t[3] and row.rat_entrez == eid and
              row.orthology_type == t[4] and row.orthology_confidence == t[5] and row.probe_ids.split(';') == probes)
        mapping[candidate.gene_id] = {'symbol': candidate.symbol, 'probes': probes, 'rat_id': t[2], 'entrez': eid}
        for probe, ids in ambiguous_probes.get(eid, []):
            exclusions.append((candidate.gene_id, probe, '///'.join(sorted(ids))))
    check('15 mapped candidates eachoneprobe', len(mapping) == 15 and all(len(x['probes']) == 1 for x in mapping.values()))
    excluded = pd.read_csv(R / 'outputs/excluded_ambiguous_probes.csv')
    check('candidate multiEntrez exclusion matches annotation', len(exclusions) == len(excluded) == 0,
          {'global_platform_multiEntrez_probes': global_ambiguous, 'candidate_multiEntrez_probes': len(exclusions)})
    # Entire matrix rank-based normalization, independent of production sorted tie blocks.
    raw = matrix.to_numpy()
    logged = np.log1p(raw) / math.log(2)
    ordered = np.sort(logged, axis=0)
    reference = np.mean(ordered, axis=1)
    prefix = np.r_[0.0, np.cumsum(reference)]
    quantile = np.zeros_like(logged)
    max_tie = 0
    for j in range(20):
        lo = stats.rankdata(logged[:, j], method='min').astype(int)-1
        hi = stats.rankdata(logged[:, j], method='max').astype(int)
        quantile[:, j] = (prefix[hi]-prefix[lo])/(hi-lo)
        max_tie = max(max_tie, int(np.max(hi-lo)))
    check('tied values receive equal rankblock means', max_tie > 1, {'max_tie_size': max_tie})
    production_results = pd.read_csv(R / 'outputs/candidate_results_all25_two_normalizations.csv')
    production_values = pd.read_csv(R / 'outputs/candidate_sample_values.csv')
    production_probes = pd.read_csv(R / 'outputs/candidate_probe_values.csv')
    rebuilt, vnn = [], []
    for method, array in [('deposited_log2', logged), ('quantile_sensitivity', quantile)]:
        pvals = []
        indices = []
        for candidate in candidates.itertuples(index=False):
            if candidate.gene_id not in mapping:
                pvals.append(1.0)
                indices.append(None)
                prow = production_results[(production_results.human_id == candidate.gene_id) & (production_results.method == method)].iloc[0]
                check('unavailable statistics remain NA ' + candidate.symbol + method,
                      all(pd.isna(prow[f]) for f in ['effect_log2_signal', 'p_value', 'q_bh_fixed25']))
                continue
            probeids = mapping[candidate.gene_id]['probes']
            locs = matrix.index.get_indexer(probeids)
            x = np.median(array[locs], axis=0)
            control, treatment = x[:5].tolist(), x[5:10].tolist()
            difference = statistics.mean(treatment)-statistics.mean(control)
            vcontrol = statistics.variance(control)/5
            vtreatment = statistics.variance(treatment)/5
            standard_error = math.sqrt(vcontrol+vtreatment)
            degrees = (vcontrol+vtreatment)**2/(vcontrol*vcontrol/4+vtreatment*vtreatment/4)
            pvalue = 2*stats.t.sf(abs(difference)/standard_error, degrees)
            radius = stats.t.isf(.025, degrees)*standard_error
            pr = production_results[(production_results.human_id == candidate.gene_id) & (production_results.method == method)].iloc[0]
            metrics = {'effect_log2_signal': difference, 'p_value': pvalue, 'ci95_low': difference-radius,
                       'ci95_high': difference+radius, 'welch_df': degrees,
                       'control_mean_log2': statistics.mean(control), 'sesamin_mean_log2': statistics.mean(treatment),
                       'control_sd_log2': statistics.stdev(control), 'sesamin_sd_log2': statistics.stdev(treatment),
                       'median_group_difference': statistics.median(treatment)-statistics.median(control)}
            close('all effect/Welch/summary values ' + candidate.symbol + method, list(metrics.values()),
                  [pr[k] for k in metrics])
            pv = production_values[(production_values.human_id == candidate.gene_id) & (production_values.method == method)].set_index('sample').loc[animals]
            close('all20 sample values ' + candidate.symbol + method, x, pv.aggregate_log2_signal.to_numpy())
            check('primary comparison flag excludes othercompound ' + candidate.symbol + method,
                  pv.primary_comparison.tolist() == [True]*10+[False]*10)
            for probe, i in zip(probeids, locs):
                pp = production_probes[(production_probes.human_id == candidate.gene_id) & (production_probes.method == method) &
                                       (production_probes.probe_id == probe)].iloc[0]
                close('all20 probe values ' + candidate.symbol + method, array[i], pp[animals].to_numpy(float))
            if candidate.symbol == 'VNN1':
                for j in range(10):
                    vnn.append({'method': method, 'sample': animals[j], 'group': 'control' if j < 5 else 'sesamin',
                                'raw_deposited_scaled_signal': raw[locs[0], j], 'recomputed_log2_signal': x[j],
                                'reported_log2_signal': pv.aggregate_log2_signal.iloc[j]})
            indices.append(len(rebuilt))
            rebuilt.append({'human_id': candidate.gene_id, 'human_symbol': candidate.symbol, 'method': method, **metrics})
            pvals.append(pvalue)
        # Independent maintained SciPy implementation, full25 placeholders included.
        adjusted = stats.false_discovery_control(pvals, method='bh')
        for idx, q in zip(indices, adjusted):
            if idx is None:
                continue
            row = rebuilt[idx]
            row['q_bh_fixed25'] = float(q)
            pr = production_results[(production_results.human_id == row['human_id']) & (production_results.method == method)].iloc[0]
            close('BH25 ' + row['human_symbol'] + method, [q], [pr.q_bh_fixed25])
    pd.DataFrame(rebuilt).to_csv(HERE / 'independently_recomputed_all15.csv', index=False)
    pd.DataFrame(vnn).to_csv(HERE / 'VNN1_all10_per_normalization.csv', index=False)
    rframe = pd.DataFrame(rebuilt)
    significant = {method: rframe.loc[(rframe.method == method) & (rframe.q_bh_fixed25 < .05), 'human_symbol'].tolist()
                   for method in ['deposited_log2', 'quantile_sensitivity']}
    check('VNN1 sole q<.05 in both fixed analyses', all(v == ['VNN1'] for v in significant.values()))
    check('VNN1 UP both analyses', (rframe.loc[rframe.human_symbol == 'VNN1', 'effect_log2_signal'] > 0).all())
    separation = {}
    for method in ['deposited_log2', 'quantile_sensitivity']:
        values = [v for v in vnn if v['method'] == method]
        control = [v['recomputed_log2_signal'] for v in values if v['group'] == 'control']
        sesamin = [v['recomputed_log2_signal'] for v in values if v['group'] == 'sesamin']
        separation[method] = {'all_sesamin_above_all_control': bool(min(sesamin) > max(control)),
                              'max_control_log2': max(control), 'min_sesamin_log2': min(sesamin),
                              'minimum_between_group_gap_log2': min(sesamin)-max(control)}
    check('VNN1 complete separation in both normalizations', all(x['all_sesamin_above_all_control'] for x in separation.values()))
    check('AREG JAG2 retained despite low signals', set(rframe.loc[rframe.human_symbol.isin(['AREG', 'JAG2']), 'method']) == {'deposited_log2', 'quantile_sensitivity'})
    exact_summary = json.loads((R / 'exact_sensitivity_summary.json').read_text(encoding='utf-8'))
    exact_path = R / 'outputs/posthoc_exact_sensitivity_all25.csv'
    check('posthoc output matches executed summary', sha(exact_path) == exact_summary['output_sha256'])
    for file, digest in exact_summary['input_sha256'].items():
        check('posthoc source matches executed summary ' + Path(file).name, sha(file) == digest)
    exact_inputs = [R / 'exact_sensitivity.py', R / 'exact_sensitivity_summary.json', R / 'posthoc_sensitivity_plan.md', exact_path]
    all_input_hashes.extend([{'path': str(file), 'sha256': sha(file), 'bytes': file.stat().st_size} for file in exact_inputs])
    exact = pd.read_csv(exact_path)
    allocations = list(itertools.combinations(range(10), 5))
    check('all252 unique label allocations includes observedandcomplement',
          len(allocations) == len(set(allocations)) == 252 and tuple(range(5)) in allocations and tuple(range(5,10)) in allocations)
    permutation_rows = []
    for method, array in [('deposited_log2', logged), ('quantile_sensitivity', quantile)]:
        probabilities, record_indices = [], []
        for candidate in candidates.itertuples(index=False):
            reported = exact[(exact.human_id == candidate.gene_id) & (exact.method == method)].iloc[0]
            if candidate.gene_id not in mapping:
                probabilities.append(1.)
                record_indices.append(None)
                check('posthoc unavailable remainsNA '+candidate.symbol+method, pd.isna(reported.permutation_p) and pd.isna(reported.permutation_q_bh25))
                continue
            locs = matrix.index.get_indexer(mapping[candidate.gene_id]['probes'])
            x = np.median(array[locs], axis=0)[:10].tolist()
            observed = statistics.mean(x[5:])-statistics.mean(x[:5])
            effects = []
            for labels in allocations:
                selected = set(labels)
                aa = [x[i] for i in range(10) if i in selected]
                bb = [x[i] for i in range(10) if i not in selected]
                effects.append(statistics.mean(aa)-statistics.mean(bb))
            exceed = sum(abs(value) >= abs(observed)-1e-12 for value in effects)
            p = exceed/252
            close('independent exactpermutation '+candidate.symbol+method,
                  [observed, p, exceed], [reported.observed_mean_difference, reported.permutation_p, reported.exceedances])
            probabilities.append(p)
            record_indices.append(len(permutation_rows))
            permutation_rows.append({'human_id': candidate.gene_id, 'human_symbol': candidate.symbol, 'method': method,
                                     'observed_mean_difference': observed, 'exceedances': exceed, 'allocations':252,'permutation_p':p})
        qs = stats.false_discovery_control(probabilities, method='bh')
        for index, q in zip(record_indices, qs):
            if index is None:
                continue
            rr = permutation_rows[index]
            rr['permutation_q_bh25'] = float(q)
            reported = exact[(exact.human_id == rr['human_id']) & (exact.method == method)].iloc[0]
            close('independent exactBH25 '+rr['human_symbol']+method, [q], [reported.permutation_q_bh25])
    pframe = pd.DataFrame(permutation_rows)
    pframe.to_csv(HERE / 'independently_recomputed_exact_sensitivity.csv', index=False)
    check('posthoc no FDR05 and VNN1 minimal2of252',
          not pframe.permutation_q_bh25.lt(.05).any() and pframe.loc[pframe.human_symbol.eq('VNN1'),'exceedances'].eq(2).all())
    for rec in all_input_hashes:
        check('final input unchanged ' + Path(rec['path']).name, sha(rec['path']) == rec['sha256'])
    for file, digest in original_outputs.items():
        check('final production output unchanged ' + file, sha(R / file) == digest)
    review = {
        'status': 'PASS_BOUNDED_REVIEW_WITH_CONDITIONAL_INFERENCE_LIMITS',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'reviewer': 'Separate Codex internal agent; not Claude or a human statistical expert',
        'scope': 'Independent mapping, all 15 mapped candidate statistics in both normalizations, all 20 sample and probe values, VNN1 primary10 values, BH25, all252 exact label allocations for every mappedcandidate and both normalizations, metadata, code review. No new network verification or raw CEL reprocessing.',
        'production_script_sha256': sha(R / 'analyze_rat.py'), 'check_script_sha256': sha(__file__),
        'protocol_sha256': sha(R / 'protocol.md'), 'production_summary_sha256': sha(R / 'analysis_summary.json'),
        'environment': {'python': sys.executable, 'version': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
        'input_sources': all_input_hashes, 'source_metadata_evidence': {'overall_design': design, 'processing': processing, 'growth': growth},
        'checks': checklog, 'checked_candidates': 15, 'candidate_family': 25, 'methods': 2,
        'blocking_findings': [],
        'required_reporting_limits': [
            'Selected distinct animals are biological units, but original outcome-time body-weight selection is unadjusted and can affect inference.',
            'Welch CIs are nominal small-sample approximations for deposited animals, not selection-adjusted intervals or guaranteed coverage for original seven-animal cohorts.',
            'Processed GCOS globally scaled array signals, not raw sequencing counts, RMA or CEL reprocessing.',
            'Healthy whole-liver dietary exposure differs from fatty-acid challenged cultured hepatocytes; opposite VNN1 direction is not matched-context replication.',
            'Quantile sensitivity uses all20 arrays as specified; other compounds affect normalization target but are excluded from the primary5v5 effect/test.',
            'Annotations are used as deposited; no probe sequence remapping or detection-call filtering. Low signal is not treated as biological absence.',
            'Nominal intervals are per-candidate, not simultaneous FDR-adjusted intervals; adjusted q values control the fixed25 test family under BH assumptions.'
        ],
        'findings': {'mapped': 15, 'unmapped_or_nonone2one': 10, 'candidate_multiEntrez_exclusions': 0,
                     'platform_multiEntrez_probes': global_ambiguous, 'all_mapped_single_probe': True, 'significant_symbols': significant,
                     'VNN1_direction': 'up', 'VNN1_complete_separation': separation, 'selection_by_effect_detected': False},
        'posthoc_sensitivity_review': {'status': 'executed_and_independently_recomputed',
            'all_candidate_tests_per_method':25,'mapped_candidate_tests_per_method':15,'allocations_per_test':252,
            'VNN1_exact_p':2/252,'VNN1_exact_BH25':float(pframe.loc[pframe.human_symbol.eq('VNN1'),'permutation_q_bh25'].iloc[0]),
            'FDR05_symbols':[],
            'interpretation':'This explicitly post-hoc sensitivity is informative because very small Welch P values rely on a parametric log-scale model at n=5 per group. The same fixed25 family is retained. Direction and magnitude are unchanged, while FDR05 status depends on the inferential method. Exchangeability is an assumption; original group-wise body-weight selection and undocumented randomization prevent a causal randomization-test interpretation. With a two-sided absolute-mean-difference statistic, complementary allocations give a theoretical minimum nonzero exact P of2/252, which constrains resolution under multiplicity.'},
        'review_outputs': [{'file': p.name, 'sha256': sha(p)} for p in HERE.glob('*.csv')]
    }
    (HERE / 'review.json').write_text(json.dumps(review, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'status': review['status'], 'checks': len(checklog), 'significant': significant,
                      'platform_ambiguous_probes': global_ambiguous}))


if __name__ == '__main__':
    main()
