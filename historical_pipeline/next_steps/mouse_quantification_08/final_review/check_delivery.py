"""Independent final audit, run only after all nine quantifications and aggregation.

No parent execution/aggregation helper is imported. Source abundance values are
summed with Decimal, then every gene-matrix cell and every candidate quantity is
checked. The 31GB FASTQ payload is not reread: existing integrity/QC provenance is
cross-checked against source metadata and file stats, with that limit disclosed.
All writes are confined to this final_review directory. Stdlib only.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import traceback

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
STEPS = STAGE.parent
PREV = STEPS / 'mouse_readiness_07'
DATA = Path('D:/ResearchData/sesamin_mouse_2026-10-08')
REL_TOL, ABS_TOL = 2e-12, 1e-9
CHECKS, FINGERPRINTS, NOTES = [], {}, []
KEY = {}


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    path = Path(path).resolve()
    key = str(path)
    if key not in FINGERPRINTS:
        with path.open('rb') as f:
            FINGERPRINTS[key] = hashlib.file_digest(f, 'sha256').hexdigest()
    return FINGERPRINTS[key]


def read_json(path):
    sha(path)
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def read_table(path, delimiter='\t'):
    sha(path)
    if Path(path).suffix == '.gz':
        stream = gzip.open(path, 'rt', encoding='utf-8-sig', newline='')
    else:
        stream = Path(path).open(encoding='utf-8-sig', newline='')
    with stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def check(name, passed, detail=None, critical=False):
    item = {'name': name, 'passed': bool(passed)}
    if detail is not None:
        # A corrupt NaN/Inf matrix cell must remain reportable in valid JSON.
        item['detail'] = json_safe(detail)
    CHECKS.append(item)
    if critical and not passed:
        raise ValueError(name + ': ' + str(detail))


def json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def close(a, b, rel=REL_TOL, absolute=ABS_TOL):
    return math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=rel, abs_tol=absolute)


def norm(path):
    return str(Path(path).resolve()).casefold()


def keyed_hashes(record, field):
    return {norm(p): h for p, h in record.get(field, {}).items()}


def numeric_raw_run_info(path):
    """Parse numeric/version members independently of malformed Windows call JSON."""
    sha(path)
    text = Path(path).read_text(encoding='utf-8')
    result = {}
    for name in ['n_targets', 'n_bootstraps', 'n_processed', 'n_pseudoaligned', 'n_unique', 'index_version', 'k-mer length']:
        matches = re.findall(r'^\s*"' + re.escape(name) + r'"\s*:\s*(\d+)\s*,?\s*$', text, re.M)
        if len(matches) != 1:
            raise ValueError(f'Expected one integer {name} in {path}')
        result[name] = int(matches[0])
    version = re.findall(r'^\s*"kallisto_version"\s*:\s*"([^"]+)"\s*,?\s*$', text, re.M)
    if len(version) != 1:
        raise ValueError('Cannot recover kallisto version without touching raw JSON')
    result['kallisto_version'] = version[0]
    return result


def direction(number):
    return 'NA' if number is None else 'up' if number > 0 else 'down' if number < 0 else 'unchanged'


def agreement(a, b):
    if a is None or b is None:
        return 'NA'
    if a == 0 or b == 0:
        return 'unchanged_in_at_least_one'
    return 'same' if (a > 0) == (b > 0) else 'opposite'


def field_matches(cell, wanted):
    if wanted is None:
        return cell == 'NA'
    if isinstance(wanted, bool):
        return cell == ('true' if wanted else 'false')
    if isinstance(wanted, float):
        try:
            return close(float(cell), wanted)
        except (ValueError, TypeError):
            return False
    return cell == str(wanted)


def ready_paths(data, matrices):
    runs = [f'DRR{n}' for n in range(1010532, 1010541)]
    required = [STAGE / 'tool_validation.json', STAGE / 'index_validation.json',
                STAGE / 'pilot_gate.json', data / 'reference/tx2gene.tsv', matrices / 'aggregation_summary.json',
                matrices / 'gene_est_counts.tsv.gz', matrices / 'gene_tpm.tsv.gz',
                matrices / 'candidate_coverage25.tsv', matrices / 'candidate_mouse_descriptive.tsv']
    for run in runs:
        run_result_path = STAGE / 'run_results' / (run + '.json')
        execution_key = 'quant_' + run
        if run_result_path.is_file():
            execution_key = json.loads(run_result_path.read_text(encoding='utf-8-sig')).get('quant_execution_key', execution_key)
        if re.fullmatch('quant_' + run + r'(?:_retry\d+)?', execution_key) is None:
            raise ValueError('Unexpected quantification execution key: ' + execution_key)
        required.extend([run_result_path, STAGE / 'qc' / (run + '.json'),
                         STAGE / 'execution' / (execution_key + '.json'),
                         STAGE / 'execution' / ('qc_' + run + '.json'),
                         data / 'quant' / run / 'abundance.tsv', data / 'quant' / run / 'run_info.json'])
    missing = [str(p) for p in required if not p.is_file()]
    check('All required completed-study artifacts exist', not missing, missing, critical=True)


def audit_sources(data, matrices):
    plan_path = PREV / 'processing_plan.json'
    plan = read_json(plan_path)
    manifest = read_json(PREV / 'fastq_manifest.json')
    orthology = read_json(PREV / 'orthology/orthology.json')
    design = read_json(PREV / 'design/sample_design.json')
    freeze = read_json(STAGE / 'execution_rules_before_pilot.json')
    summary = read_json(matrices / 'aggregation_summary.json')
    check('Real completed aggregation, not synthetic', summary['status'] == 'passed' and summary['synthetic'] is False, critical=True)
    check('Scientific plan unchanged since pre-pilot freeze', sha(plan_path) == freeze['scientific_plan_sha256'])
    for rel, expected in plan['input_snapshot_sha256'].items():
        check('Frozen stage07 input ' + rel, sha(PREV / rel) == expected)
    check('Mouse source does not establish independent animal count', design['study']['known_independent_animal_count'] is None)
    entries = manifest['entries']
    run_entries = defaultdict(list)
    for e in entries:
        run_entries[e['run']].append(e)
    runs = sorted(run_entries)
    check('Exactly nine frozen runs', runs == [f'DRR{n}' for n in range(1010532, 1010541)], runs, critical=True)
    check('Eighteen unique FASTQ names/URLs/checksums', len(entries) == 18 and all(len({e[k] for e in entries}) == 18 for k in ['filename', 'url', 'expected_md5']))
    check('Frozen raw size total', sum(e['bytes'] for e in entries) == 31230644073)
    check('Exactly two distinct mates per run', all(sorted(e['mate'] for e in pair) == [1, 2] for pair in run_entries.values()), critical=True)
    metadata = {r: dict(sorted(es, key=lambda e: e['mate'])[0]) for r, es in run_entries.items()}
    check('Three submitted libraries per inferred label group', Counter(m['provisional_group'] for m in metadata.values()) == {'ctrl': 3, 'HFD': 3, 'sess': 3})
    design_by_run = {s['fields']['mirror_run']['value']: s['fields'] for s in design['samples']}
    for run, m in metadata.items():
        check('Library label matches source design ' + run, m['original_label'] == design_by_run[run]['sample_label']['value'])
        check('Inferred group matches source prefix ' + run, m['provisional_group'] == design_by_run[run]['label_prefix']['value'] and 'inferred' in m['condition_mapping'])
        check('Biological independence remains unknown ' + run, design_by_run[run]['independent_biological_replicate']['status'] == 'unknown')
    human_path = STEPS / 'sesamin_exploration_04/priority_all.csv'
    human_rows = read_table(human_path, delimiter=',')
    human = {r['gene_id']: r for r in human_rows}
    check('Original25 human input hash retained', sha(human_path) == orthology['input_sha256'])
    check('Original25 human/candidate IDs match', len(human_rows) == len(human) == len(orthology['candidates']) == 25 and set(human) == {c['original_human_id'] for c in orthology['candidates']}, critical=True)
    for path, expected in summary['input_sha256'].items():
        check('Aggregation input hash ' + Path(path).name, sha(path) == expected)
    for rel, expected in summary['output_sha256'].items():
        check('Aggregation output hash ' + rel, sha(matrices / rel) == expected)
    check('Aggregation script matches executed hash', sha(STAGE / 'aggregation/aggregate_quantification.py') == summary['script_sha256'])
    return plan, manifest, orthology, metadata, human, summary, freeze


def audit_integrity_records(data, manifest):
    paths = [STAGE / 'pilot_download_integrity.json', STAGE / 'remaining_download_integrity.json']
    paths += sorted((STAGE / 'per_run_download_integrity').glob('*.json'))
    records = defaultdict(list)
    for path in paths:
        if not path.exists():
            continue
        record = read_json(path)
        check('Integrity record passed ' + path.name, record.get('passed') is True)
        for f in record['files']:
            records[Path(f['file']).name].append(f)
    needed = {e['filename'] for e in manifest['entries']}
    check('Integrity evidence covers exactly18 FASTQ', set(records) == needed, sorted(set(records) ^ needed), critical=True)
    evidence = {}
    for e in manifest['entries']:
        name = e['filename']
        for f in records[name]:
            check('Saved complete-file checksum agrees with ENA ' + name,
                  f['passed'] is True and f['expected_md5'] == f['actual_md5'] == e['expected_md5']
                  and f['expected_bytes'] == f['actual_bytes'] == e['bytes'] and f['url'] == e['url'])
        check('Repeated integrity records agree ' + name, len({r['sha256'] for r in records[name]}) == 1)
        path = data / 'fastq' / name
        check('FASTQ current size/complete transfer ' + name,
              path.is_file() and path.stat().st_size == e['bytes'] and not Path(str(path) + '.aria2').exists())
        evidence[norm(path)] = records[name][0]
    NOTES.append('The18 original FASTQ payloads were not reread or rehashed by this final checker; saved full-file MD5/SHA256, per-run input fingerprints, exact sizes and QC timestamps are audited instead.')
    return evidence


def audit_tool_and_index(data, plan):
    tool = read_json(STAGE / 'tool_validation.json')
    index = read_json(STAGE / 'index_validation.json')
    ref = read_json(STAGE / 'reference_build/reference_validation.json')
    check('Pinned portable kallisto passed version/toy checks', tool['passed'] is True and tool['version_output'] == 'kallisto, version 0.51.1' and tool['executable_run'] is True)
    check('Tool source archive pin', tool['archive_sha256'] == plan['tool']['archive_sha256'])
    check('Current executable matches actual tested file', sha(tool['executable']) == tool['executable_sha256'])
    check('Index completed with k31 and pinned version', index['passed'] is True and index['k'] == 31 and index['version'] == tool['version_output'])
    check('Fixed actual index hash/size', sha(index['index']) == index['index_sha256'] and Path(index['index']).stat().st_size == index['index_bytes'])
    check('Index input FASTA matches unchanged reference', sha(index['fasta']) == index['fasta_sha256'])
    for output in ref['outputs']:
        check('Reference output remains unchanged ' + Path(output['path']).name, sha(output['path']) == output['sha256'])
    check('All15 primary genes covered in reference', ref['all_15_primary_genes_covered'] is True)
    actual_index = read_json(STAGE / 'execution/mouse_index.json')
    wanted = [tool['executable'], 'index', '-i', index['index'], '-k', '31', index['fasta']]
    check('Recorded index argv/exit match fixed design', actual_index['argv'] == wanted and actual_index['status'] == 'passed' and actual_index['exit_code'] == 0)
    return tool, index, ref


def audit_runs(data, metadata, evidence, tool, index, freeze):
    results = {}
    for run in sorted(metadata):
        result_path = STAGE / 'run_results' / (run + '.json')
        result = read_json(result_path)
        qc_path = STAGE / 'qc' / (run + '.json')
        qc = read_json(qc_path)
        quant_key = result.get('quant_execution_key', 'quant_' + run)
        check('Explicit successful quant execution key belongs to run ' + run,
              re.fullmatch('quant_' + run + r'(?:_retry\d+)?', quant_key) is not None, quant_key, critical=True)
        quant_exec = read_json(STAGE / 'execution' / (quant_key + '.json'))
        prior_attempts = []
        for attempt_path in sorted((STAGE / 'execution').glob('quant_' + run + '*.json')):
            if attempt_path.stem != quant_key:
                attempt = read_json(attempt_path)
                prior_attempts.append({'execution_key': attempt_path.stem, 'status': attempt['status'], 'exit_code': attempt.get('exit_code')})
        if prior_attempts:
            NOTES.append({'run': run, 'retained_other_quant_attempts': prior_attempts, 'successful_execution_key_checked': quant_key})
        qc_exec = read_json(STAGE / 'execution' / ('qc_' + run + '.json'))
        abundance = data / 'quant' / run / 'abundance.tsv'
        info_path = data / 'quant' / run / 'run_info.json'
        info = numeric_raw_run_info(info_path)
        check('Quantification technical gate passed ' + run, result['passed'] is True and result['numerical_and_structure_checks_passed'] is True and result['flags_requiring_review'] == [], critical=True)
        check('Whole-pair gzip QC passed ' + run, qc['passed'] is True and qc['full_scan_complete'] is True and qc['gzip_crc_and_isize_validation'] == 'passed_end_of_both_streams')
        check('All processed pairs equal full paired QC ' + run, info['n_processed'] == qc['total_pairs'] == qc['mates']['R1']['records_validated'] == qc['mates']['R2']['records_validated'] and info['n_processed'] > 0)
        check('Quantification numeric order/count/version ' + run, 0 <= info['n_unique'] <= info['n_pseudoaligned'] <= info['n_processed'] and info['n_pseudoaligned'] > 0 and info['n_bootstraps'] == 0 and info['kallisto_version'] == '0.51.1' and info['k-mer length'] == 31)
        check('Raw numeric info matches reported result ' + run, all(result['run_info'][k] == v for k, v in info.items()))
        mates = [data / 'fastq' / f'{run}_{m}.fastq.gz' for m in (1, 2)]
        wanted = [tool['executable'], 'quant', '-i', index['index'], '-o', str(data / 'quant' / run), '-t', '8', '--plaintext', *map(str, mates)]
        check('Quantification actual argv/exit ' + run, quant_exec['argv'] == wanted and quant_exec['status'] == 'passed' and quant_exec['exit_code'] == 0)
        check('QC actual argv/exit ' + run, qc_exec['status'] == 'passed' and qc_exec['exit_code'] == 0 and qc_exec['argv'][1:] == ['-X', 'utf8', str(STAGE / 'fastq_qc/fastq_qc.py'), *map(str, mates), '--output', str(qc_path)])
        check('Per-run output hash records ' + run,
              result['abundance_sha256'] == sha(abundance) and result['raw_run_info_sha256'] == sha(info_path) and result['qc_sha256'] == sha(qc_path))
        for execrec, expected_outputs in [(quant_exec, [abundance, info_path]), (qc_exec, [qc_path])]:
            for path in expected_outputs:
                check('Execution binds actual output ' + run + ' ' + path.name,
                      keyed_hashes(execrec, 'output_sha256').get(norm(path)) == sha(path))
        quant_inputs = keyed_hashes(quant_exec, 'input_fingerprints')
        qc_inputs = keyed_hashes(qc_exec, 'input_fingerprints')
        qc_script = STAGE / 'fastq_qc/fastq_qc.py'
        check('QC implementation fingerprint retained ' + run,
              quant_inputs.get(norm(qc_script)) == qc_inputs.get(norm(qc_script)) == sha(qc_script))
        for path in mates:
            check('QC/quant input fingerprint agrees with full-file integrity ' + path.name,
                  quant_inputs.get(norm(path)) == qc_inputs.get(norm(path)) == evidence[norm(path)]['sha256'])
        check('Quantification bound to same fixed executable/index/QC ' + run,
              quant_inputs.get(norm(tool['executable'])) == tool['executable_sha256']
              and quant_inputs.get(norm(index['index'])) == result['index_sha256'] == index['index_sha256']
              and quant_inputs.get(norm(qc_path)) == sha(qc_path))
        for recorded in qc['input_files']:
            path = Path(recorded['path'])
            check('FASTQ stats unchanged since full QC ' + path.name,
                  path.stat().st_size == recorded['compressed_bytes'] and path.stat().st_mtime_ns == recorded['mtime_ns'])
        rate = info['n_pseudoaligned'] / info['n_processed']
        check('Pseudoalignment rate recomputed ' + run, close(rate, result['pseudoalignment_fraction']))
        thresholds = freeze['diagnostic_thresholds']
        flags = []
        if rate < thresholds['pseudoalignment_fraction_below']:
            flags.append('pseudoalignment_below_50_percent_requires_diagnostic_review')
        for mate, stat in qc['mates'].items():
            sample = stat['sample_first_pairs_only']
            check('Only first100k quality/adapter observations ' + run + ' ' + mate,
                  sample['reads'] == min(100000, qc['total_pairs']) and sample['phred_offset_assumption'] == 33)
            if sample['q30_base_fraction_assuming_33'] < thresholds['sampled_Q30_fraction_below']:
                flags.append(mate + '_first100k_q30_below_80pct')
            if sample['fraction_containing_motif'] > thresholds['sampled_exact_adapter_fraction_above']:
                flags.append(mate + '_first100k_exact_adapter_motif_over_5pct')
        check('Predeclared diagnostic flags recomputed ' + run, flags == result['flags_requiring_review'], flags)
        results[run] = {'info': info, 'reported': result}
    gate = read_json(STAGE / 'pilot_gate.json')
    check('Pilot gate bound to deterministic first run', gate['passed'] is True and gate['pilot_run'] == sorted(metadata)[0] and gate['run_result_sha256'] == sha(STAGE / 'run_results' / (gate['pilot_run'] + '.json')) and gate['biological_validation'] is False)
    KEY['total_processed_read_pairs'] = sum(r['info']['n_processed'] for r in results.values())
    KEY['pseudoalignment_fraction_minmax'] = [min(r['reported']['pseudoalignment_fraction'] for r in results.values()), max(r['reported']['pseudoalignment_fraction'] for r in results.values())]
    return results


def recompute_genes(data, matrices, metadata, results, ref):
    # Independent decimal accumulation avoids copying the aggregation function.
    mapping_rows = read_table(data / 'reference/tx2gene.tsv')
    mapping = {r['transcript_id']: (r['gene_id'], int(r['sequence_length'])) for r in mapping_rows}
    genes = {g for g, length in mapping.values()}
    check('Unique versioned transcript mapping', len(mapping) == len(mapping_rows) == ref['combined_transcripts'], critical=True)
    check('Quantifiable gene denominator matches FASTA reference', len(genes) == ref['combined_genes'], critical=True)
    computed = {}
    for run in sorted(metadata):
        acc = {g: [Decimal(0), Decimal(0)] for g in genes}
        seen = set()
        total_tpm, total_counts = Decimal(0), Decimal(0)
        input_path = data / 'quant' / run / 'abundance.tsv'
        sha(input_path)
        with localcontext() as ctx:
            ctx.prec = 40
            with input_path.open(encoding='utf-8-sig', newline='') as f:
                for row in csv.DictReader(f, delimiter='\t'):
                    tid = row['target_id']
                    if tid in seen or tid not in mapping:
                        raise ValueError('Duplicate/unmapped abundance transcript ' + run + '/' + tid)
                    seen.add(tid)
                    gid, length = mapping[tid]
                    count, tpm = Decimal(row['est_counts']), Decimal(row['tpm'])
                    eff = Decimal(row['eff_length'])
                    if not all(v.is_finite() and v >= 0 for v in (count, tpm, eff)):
                        raise ValueError('Nonfinite/negative abundance ' + run + '/' + tid)
                    if Decimal(row['length']) != length or (eff == 0 and tpm > 0):
                        raise ValueError('Invalid target/effective length ' + run + '/' + tid)
                    acc[gid][0] += count
                    acc[gid][1] += tpm
                    total_counts += count
                    total_tpm += tpm
        check('Every reference transcript retained ' + run, seen == set(mapping) and len(seen) == results[run]['info']['n_targets'])
        check('Original abundance TPM normalization ' + run, abs(total_tpm - Decimal(1000000)) <= Decimal('5.000001'), str(total_tpm))
        check('Original abundance fragment mass vs pseudoaligned ' + run,
              abs(float(total_counts) - results[run]['info']['n_pseudoaligned']) < max(10, results[run]['info']['n_pseudoaligned'] * 0.001))
        check('Recomputed totals match per-run report ' + run,
              close(float(total_tpm), results[run]['reported']['total_tpm']) and close(float(total_counts), results[run]['reported']['total_estimated_counts']))
        computed[run] = {g: (float(v[0]), float(v[1])) for g, v in acc.items()}
        print(json.dumps({'event': 'recomputed_gene_sums', 'run': run, 'genes': len(acc)}, ensure_ascii=False), flush=True)
    runs = sorted(metadata)
    for metric, index in [('est_counts', 0), ('tpm', 1)]:
        path = matrices / f'gene_{metric}.tsv.gz'
        sha(path)
        failures, seen = [], set()
        cells, largest_error = 0, 0.0
        with gzip.open(path, 'rt', encoding='utf-8-sig', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            check('Exact gene-matrix library columns ' + metric, reader.fieldnames == ['gene_id', *runs], reader.fieldnames, critical=True)
            for row in reader:
                gid = row['gene_id']
                if gid in seen or gid not in genes:
                    raise ValueError('Duplicate/unexpected emitted gene ' + gid)
                seen.add(gid)
                for run in runs:
                    observed, expected = float(row[run]), computed[run][gid][index]
                    cells += 1
                    error = abs(observed - expected)
                    largest_error = max(largest_error, error)
                    if not close(observed, expected) and len(failures) < 20:
                        failures.append({'gene_id': gid, 'run': run, 'observed': observed, 'independent_expected': expected})
        check('Exact gene-matrix universe ' + metric, seen == genes)
        check('All gene-matrix cells independently reproduced ' + metric, not failures,
              {'cells': cells, 'maximum_absolute_error': largest_error, 'relative_tolerance': REL_TOL, 'absolute_tolerance': ABS_TOL, 'first_discrepancies': failures})
    KEY.update(libraries=len(runs), transcripts=len(mapping), quantifiable_genes=len(genes), gene_matrix_cells_checked=2 * len(runs) * len(genes))
    return genes, computed


def audit_candidates(matrices, orthology, plan, metadata, human, genes, computed, summary):
    candidates = orthology['candidates']
    coverage = read_table(matrices / 'candidate_coverage25.tsv')
    rows = read_table(matrices / 'candidate_mouse_descriptive.tsv')
    check('Exactly25 coverage rows retained', len(coverage) == len({r['original_human_id'] for r in coverage}) == 25)
    check('Exactly26 descriptive rows including both ACAT2 homologues and NA rows', len(rows) == 26)
    expected_keys = {(c['original_human_id'], h['mouse_id']) for c in candidates for h in c['mouse_orthologues']}
    expected_keys |= {(c['original_human_id'], 'NA') for c in candidates if not c['mouse_orthologues']}
    actual_keys = {(r['original_human_id'], r['mouse_gene_id']) for r in rows}
    check('Exact candidate/homologue key set without duplicate rows', actual_keys == expected_keys and len(actual_keys) == len(rows), critical=True)
    coverage_by_id = {r['original_human_id']: r for r in coverage}
    descriptive = {(r['original_human_id'], r['mouse_gene_id']): r for r in rows}
    primary = {c['original_human_id'] for c in candidates if c['primary_comparison_eligible']}
    check('Frozen primary denominator is15', primary == set(plan['primary_human_gene_ids']) and len(primary) == 15)
    check('Exact frozen primary mouse genes', {h['mouse_id'] for c in candidates if c['primary_comparison_eligible'] for h in c['mouse_orthologues']} == set(plan['primary_mouse_gene_ids']))
    runs = sorted(metadata)
    findings = []
    na_count, both_same, low_count = 0, 0, 0
    for candidate in candidates:
        hid = candidate['original_human_id']
        homologues = candidate['mouse_orthologues']
        available = [h['mouse_id'] in genes for h in homologues]
        wanted_coverage = {'original_human_id': hid, 'human_symbol': candidate['input_symbol'],
                           'orthology_status': candidate['status'], 'primary_comparison_eligible': hid in primary,
                           'mouse_orthologue_count': len(homologues), 'mouse_gene_ids': ';'.join(h['mouse_id'] for h in homologues) or None,
                           'mouse_symbols': ';'.join(h['mouse_symbol'] for h in homologues) or None,
                           'reference_gene_count': sum(available), 'coverage_status': 'no_mouse_orthologue' if not homologues else 'all_targets_in_reference' if all(available) else 'target_missing_from_reference'}
        errors = {k: {'observed': coverage_by_id[hid].get(k), 'expected': v} for k, v in wanted_coverage.items() if not field_matches(coverage_by_id[hid].get(k), v)}
        check('Candidate coverage metadata ' + candidate['input_symbol'], not errors, errors)
        for homologue in homologues or [None]:
            gid = homologue['mouse_id'] if homologue else None
            present = gid is not None and gid in genes
            row = descriptive[(hid, gid or 'NA')]
            wanted = {'human_symbol': candidate['input_symbol'], 'mouse_symbol': homologue['mouse_symbol'] if homologue else None,
                      'primary_comparison_eligible': hid in primary, 'orthology_status': candidate['status'],
                      'comparison_role': 'primary_one_to_one' if hid in primary else 'secondary_one_to_many' if homologue else 'not_comparable_no_orthologue',
                      'coverage_status': 'quantified' if present else 'no_mouse_orthologue' if not homologue else 'target_missing_from_reference'}
            vals = {r: computed[r][gid][1] if present else None for r in runs}
            wanted.update({run + '_TPM': value for run, value in vals.items()})
            means = {group: math.fsum(vals[r] for r in runs if metadata[r]['provisional_group'] == group) / 3 if present else None for group in ['ctrl', 'HFD', 'sess']}
            wanted.update({'mean_' + g + '_TPM_provisional': v for g, v in means.items()})
            wanted['all_libraries_TPM_lt1'] = all(v < 1 for v in vals.values()) if present else None
            for suffix, pc in [('0p1', 0.1), ('0', 0.0), ('0p5', 0.5), ('1', 1.0)]:
                for a, b in [('HFD', 'ctrl'), ('sess', 'HFD')]:
                    value = None
                    if present and means[a] + pc > 0 and means[b] + pc > 0:
                        value = math.log2((means[a] + pc) / (means[b] + pc))
                    wanted[f'log2_{a}_vs_{b}_pc{suffix}'] = value
            md, ma = wanted['log2_HFD_vs_ctrl_pc0p1'], wanted['log2_sess_vs_HFD_pc0p1']
            hd, ha = float(human[hid]['log2_OA_vs_control']), float(human[hid]['log2_sesamin_vs_OA'])
            wanted.update(mouse_HFD_vs_ctrl_direction=direction(md), mouse_sess_vs_HFD_direction=direction(ma),
                          human_OA_vs_ctrl_direction=direction(hd), human_sesamin_vs_OA_direction=direction(ha),
                          cross_species_disease_direction=agreement(md, hd), cross_species_sesamin_direction=agreement(ma, ha))
            wanted['both_cross_species_directions_same'] = (agreement(md, hd) == agreement(ma, ha) == 'same') if present else None
            wanted['mouse_two_contrasts_opposite_direction'] = (md * ma < 0) if present else None
            errors = {k: {'observed': row.get(k), 'expected': v} for k, v in wanted.items() if not field_matches(row.get(k), v)}
            check('All retained TPM, means, ratios and directions ' + candidate['input_symbol'] + '/' + str(gid), not errors, errors)
            if not present:
                na_count += 1
            if hid in primary:
                both_same += wanted['both_cross_species_directions_same'] is True
                low_count += wanted['all_libraries_TPM_lt1'] is True
                findings.append({'human_symbol': candidate['input_symbol'], 'mouse_gene_id': gid,
                                 'means_provisional_TPM': means, 'mouse_log2_HFD_vs_ctrl_pc0p1': md,
                                 'mouse_log2_sess_vs_HFD_pc0p1': ma,
                                 'both_cross_species_directions_same': wanted['both_cross_species_directions_same'],
                                 'all_libraries_TPM_lt1': wanted['all_libraries_TPM_lt1']})
    check('All nine absent-orthologue rows remain NA', na_count == 9)
    acat_rows = [r for r in rows if r['human_symbol'] == 'ACAT2']
    check('ACAT2 retains both predeclared secondary homologues', {r['mouse_gene_id'] for r in acat_rows} == {'ENSMUSG00000023832', 'ENSMUSG00000062480'} and all(r['primary_comparison_eligible'] == 'false' for r in acat_rows))
    check('CXCL8 has no substituted mouse gene', len([r for r in rows if r['human_symbol'] == 'CXCL8' and r['mouse_gene_id'] == 'NA']) == 1)
    for key, expected in [('library_count', 9), ('human_candidate_count', 25), ('candidate_homologue_rows_including_NA', 26), ('primary_comparison_gene_count', 15), ('primary_both_cross_species_directions_same_count', both_same), ('primary_all_libraries_TPM_lt1_count', low_count), ('transcript_count', KEY['transcripts']), ('gene_count', len(genes))]:
        check('Published aggregation headline ' + key, summary[key] == expected, expected)
    check('Published low-expression candidate ID list', set(summary['primary_low_expression_gene_ids']) == {r['mouse_gene_id'] for r in findings if r['all_libraries_TPM_lt1']})
    check('Published complete coverage counts', summary['coverage_counts'] == dict(Counter(r['coverage_status'] for r in coverage)))
    KEY.update(human_candidates=25, candidate_rows_with_NA=26, primary_denominator=15, absent_orthologue_NA_rows=na_count,
               primary_both_cross_species_directions_same_count=both_same, primary_all_libraries_TPM_lt1_count=low_count,
               primary_findings=findings)


def write_report(result):
    HERE.mkdir(parents=True, exist_ok=True)
    (HERE / 'review.json').write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    failed = result['failed_checks']
    text = [
        '# 第08步獨立交付核查', '',
        '**' + result['assessment'] + '**；檢查時間 ' + result['completed_utc'] + '。', '',
        f"已執行 {len(result['checks'])} 項分組檢查，失敗 {len(failed)} 項。", '',
        '本檢查直接從九份原始 abundance.tsv 與 tx2gene，以 Decimal 重新累計基因估計計數與 TPM，逐格核對基因矩陣，重算完整候選表。', '',
        '## 關鍵數值', '',
    ]
    for k, v in result['key_numbers'].items():
        if k != 'primary_findings':
            text.append(f'- `{k}`：{v}')
    text += ['', '## 必須保留的限制', ''] + ['- ' + v for v in result['caveats']]
    text += ['', '## 未通過項目', '']
    text += ['- ' + item['name'] + ': ' + json.dumps(item.get('detail'), ensure_ascii=False) for item in failed] or ['無。']
    if result.get('exception'):
        text += ['', '程序中止資訊：', '', '```text', result['exception'], '```']
    text += ['', '完整檢查、输入 SHA256、候選重算值及容差見 [review.json](review.json)。']
    (HERE / 'report.md').write_text('\n'.join(text) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=DATA)
    parser.add_argument('--matrices', type=Path)
    args = parser.parse_args()
    started = utc()
    exception = None
    ready = False
    try:
        matrices = args.matrices or args.data_root / 'matrices'
        ready_paths(args.data_root, matrices)
        ready = True
        plan, manifest, ortho, metadata, human, summary, freeze = audit_sources(args.data_root, matrices)
        evidence = audit_integrity_records(args.data_root, manifest)
        tool, index, ref = audit_tool_and_index(args.data_root, plan)
        runs = audit_runs(args.data_root, metadata, evidence, tool, index, freeze)
        genes, computed = recompute_genes(args.data_root, matrices, metadata, runs, ref)
        audit_candidates(matrices, ortho, plan, metadata, human, genes, computed, summary)
    except Exception:
        exception = traceback.format_exc()
        if not any(not c['passed'] for c in CHECKS):
            check('Checker completed without execution error', False, exception)
    failed = [c for c in CHECKS if not c['passed']]
    passed = ready and not failed and exception is None
    result = {'passed': passed, 'status': 'passed' if passed else 'failed' if ready else 'not_ready',
              'assessment': 'Share with caveats' if passed else 'Needs revision / incomplete artifacts',
              'started_utc': started, 'completed_utc': utc(), 'python_executable': sys.executable,
              'script_sha256': sha(Path(__file__)), 'input_sha256': FINGERPRINTS,
              'checks': CHECKS, 'failed_checks': failed, 'key_numbers': KEY, 'exception': exception,
              'caveats': [
                  '分組 ctrl/HFD/sess 仍依樣本標籤及研究摘要推定；九份文庫不代表九隻獨立動物或九次獨立培養。',
                  '這是描述性探索；方向一致不等於獨立生物學驗證、新機制、因果作用或臨床效益。',
                  '基因 estimated counts 是可含小數的轉錄本估计之和，不是原始整數讀數；沒有計算生物學 P 值、FDR 或信賴區間。',
                  '跨物種僅比較方向，沒有合併人源 normalized counts 與鼠 TPM 的大小。',
                  '品質及接頭指標僅來自固定前100,000對，假設 Phred+33；完整gzip與結構檢查則涵蓋全部讀對。',
                  '本次覆核不重讀31GB FASTQ；核對已執行的完整檔案雜湊、QC與實際命令／輸出鏈。',
                  '基因矩陣的每一格均與原始 abundance.tsv 重算結果比較；轉錄本矩陣本步核對產物雜湊，未另外逐格重算。',
                  '3,262個GTF-only TEC轉錄本／3,219基因不在提供的FASTA中；定量分母限於參考中75,129基因。',
              ], 'notes': NOTES, 'relative_numeric_tolerance': REL_TOL, 'absolute_numeric_tolerance': ABS_TOL,
              'scientific_independence_established': False, 'external_model_APIs_or_Claude_called_by_checker': False,
              'authorship_note': 'Internal Codex agents authored and reviewed scripts; this numerical checker itself calls no model or external service.'}
    write_report(result)
    print(json.dumps({'passed': passed, 'status': result['status'], 'checks': len(CHECKS), 'failed_checks': len(failed), 'report': str(HERE / 'report.md')}, ensure_ascii=False))
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
