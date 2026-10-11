"""Download, validate and combine pinned Ensembl116/GRCm39 reference; stdlib only.

Run using the project .venv. Raw source files are never rewritten after completion.
The generated tx2gene table is restricted to exactly the combined FASTA records.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
STEPS = HERE.parents[1]
PLAN_PATH = STEPS / 'mouse_readiness_07' / 'processing_plan.json'
REF = Path(r'D:\ResearchData\sesamin_mouse_2026-10-08\reference')
ATTR = re.compile(r'([A-Za-z0-9_]+) "([^"]*)";')


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def events(kind, **fields):
    with (HERE / 'events.jsonl').open('a', encoding='utf-8') as f:
        f.write(json.dumps({'utc': utc(), 'event': kind, **fields}, ensure_ascii=False) + '\n')


def download(plan):
    pending = []
    records = []
    for spec in plan['reference']['files']:
        path = REF / Path(spec['destination']).name
        if path.exists():
            if path.stat().st_size != int(spec['bytes']):
                raise ValueError(f'Existing raw input has unexpected size; not modifying: {path}')
            events('reuse_existing_raw_input', path=str(path), size=path.stat().st_size)
        else:
            pending.append(spec['url'] + '\n  out=' + path.name + '\n')
        records.append({**spec, 'destination': str(path)})
    input_path = HERE / 'reference_downloads.aria2'
    input_path.write_text('\n'.join(pending), encoding='utf-8')
    if pending:
        exe = shutil.which('aria2c')
        if exe is None:
            raise RuntimeError('Existing aria2c was not found; no dependencies installed')
        args = [exe, '--input-file=' + str(input_path), '--dir=' + str(REF),
                '--max-concurrent-downloads=3', '--max-connection-per-server=8', '--split=8',
                '--continue=true', '--allow-overwrite=false', '--auto-file-renaming=false',
                '--file-allocation=none', '--check-certificate=true', '--max-tries=5',
                '--retry-wait=5', '--summary-interval=10', '--console-log-level=notice']
        started = utc()
        with (HERE / 'aria2_download.log').open('a', encoding='utf-8') as log:
            log.write('\nSTART ' + started + '\n')
            log.flush()
            result = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, check=False)
        events('download_command', args=args, started_utc=started, completed_utc=utc(), exit_code=result.returncode)
        if result.returncode:
            raise RuntimeError(f'aria2 download failed ({result.returncode}); retained inputs/log for investigation')
    for record in records:
        path = Path(record['destination'])
        actual = path.stat().st_size
        if actual != int(record['bytes']):
            raise ValueError(f'Reference size mismatch: {path} {actual} != {record["bytes"]}')
        record.update(actual_bytes=actual, sha256=sha(path), upstream_sha256_available=False,
                      local_hash_interpretation='Local provenance fingerprint; not an upstream checksum/signature verification')
    write_json(HERE / 'download_records.json', records)
    return {r['name']: Path(r['destination']) for r in records}, records


def parse_gtf(path):
    mapping, gene_meta = {}, {}
    headers = []
    line_count = 0
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        for line in f:
            line_count += 1
            if line.startswith('#'):
                headers.append(line.rstrip())
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'GTF malformed line {line_count}')
            if fields[2] != 'transcript':
                continue
            a = dict(ATTR.findall(fields[8]))
            for required in ['gene_id', 'transcript_id', 'transcript_version', 'gene_version']:
                if required not in a:
                    raise ValueError(f'GTF missing {required} line {line_count}')
            tid = a['transcript_id'] + '.' + a['transcript_version']
            row = {'transcript_id': tid, 'gene_id': a['gene_id'],
                   'gene_version': a['gene_version'], 'gene_name': a.get('gene_name', ''),
                   'gene_biotype': a.get('gene_biotype', ''),
                   'transcript_biotype': a.get('transcript_biotype', '')}
            if tid in mapping and mapping[tid] != row:
                raise ValueError(f'Conflicting GTF mapping for {tid}')
            mapping[tid] = row
            g = tuple(row[k] for k in ['gene_version', 'gene_name', 'gene_biotype'])
            if row['gene_id'] in gene_meta and gene_meta[row['gene_id']] != g:
                raise ValueError('Conflicting GTF gene metadata ' + row['gene_id'])
            gene_meta[row['gene_id']] = g
    if not any('genome-build GRCm39' in h for h in headers):
        raise ValueError('GTF genome build does not confirm GRCm39')
    events('gtf_parsed_full_gzip', lines=line_count, transcripts=len(mapping), genes=len(gene_meta))
    return mapping, {'lines': line_count, 'transcripts': len(mapping), 'genes': len(gene_meta), 'headers': headers, 'gzip_integrity_passed': True}


def fasta_records(path):
    head, chunks = None, []
    with gzip.open(path, 'rt', encoding='ascii') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if head is not None:
                    yield head, ''.join(chunks)
                head, chunks = line[1:], []
            elif line:
                if head is None:
                    raise ValueError('FASTA sequence before header')
                chunks.append(line)
    if head is not None:
        yield head, ''.join(chunks)


def combine(paths, mapping, plan):
    combined = REF / 'combined.fa'
    tx2gene = REF / 'tx2gene.tsv'
    # Refuse to overwrite completed products: rerunning validates the same completed build.
    if combined.exists() or tx2gene.exists():
        raise FileExistsError('Completed generated reference exists; inspect reference_validation.json, do not overwrite it')
    partial = REF / 'combined.fa.building'
    if partial.exists():
        raise FileExistsError('Partial reference exists from a prior attempt; retain and investigate before rerun')
    seen = {}
    duplicates, source_counts, biotypes = [], {}, {}
    gene_counts = {}
    total_bases = 0
    with partial.open('w', encoding='ascii', newline='\n') as out:
        for source in ['cdna', 'ncrna']:
            total, unique = 0, 0
            for header, sequence in fasta_records(paths[source]):
                total += 1
                tid = header.split()[0]
                if not sequence or re.search('[^ACGTNacgtn]', sequence):
                    raise ValueError(f'Unexpected FASTA sequence alphabet/empty sequence for {tid}')
                if tid not in mapping:
                    raise ValueError(f'FASTA transcript absent in same-release GTF: {tid}')
                row = mapping[tid]
                header_fields = dict(field.split(':', 1) for field in header.split()[1:] if ':' in field)
                fasta_gene = header_fields.get('gene')
                expected_gene = row['gene_id'] + '.' + row['gene_version']
                if fasta_gene != expected_gene:
                    raise ValueError(f'FASTA/GTF gene discrepancy for {tid}: {fasta_gene} != {expected_gene}')
                for key in ['gene_biotype', 'transcript_biotype']:
                    if header_fields.get(key) != row[key]:
                        raise ValueError(f'FASTA/GTF {key} discrepancy for {tid}')
                if tid in seen:
                    # Store full sequence: assert actual sequence equality, not only hash equality.
                    if seen[tid][1] != sequence:
                        raise ValueError(f'Conflicting duplicate transcript sequence {tid}')
                    duplicates.append({'transcript_id': tid, 'retained_source': seen[tid][0],
                                       'duplicate_source': source, 'length': len(sequence),
                                       'sequence_sha256': hashlib.sha256(sequence.encode('ascii')).hexdigest()})
                    continue
                seen[tid] = (source, sequence)
                out.write('>' + header + '\n')
                for start in range(0, len(sequence), 60):
                    out.write(sequence[start:start+60] + '\n')
                unique += 1
                total_bases += len(sequence)
                gene_counts[row['gene_id']] = gene_counts.get(row['gene_id'], 0) + 1
                biotypes[row['transcript_biotype']] = biotypes.get(row['transcript_biotype'], 0) + 1
            source_counts[source] = {'total_records': total, 'unique_records_added': unique, 'gzip_integrity_passed': True}
    missing_primary = sorted(set(plan['primary_mouse_gene_ids']) - set(gene_counts))
    if missing_primary:
        raise ValueError('Primary gene(s) missing from FASTA: ' + str(missing_primary))
    tsv_partial = REF / 'tx2gene.tsv.building'
    fields = ['transcript_id', 'gene_id', 'gene_version', 'gene_name', 'gene_biotype', 'transcript_biotype', 'sequence_length', 'source_fasta']
    with tsv_partial.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for tid, (source, sequence) in seen.items():
            writer.writerow({**mapping[tid], 'sequence_length': len(sequence), 'source_fasta': source})
    os.rename(partial, combined)
    os.rename(tsv_partial, tx2gene)
    write_json(HERE / 'identical_duplicate_transcripts.json', duplicates)
    return {'combined_fasta': str(combined), 'tx2gene': str(tx2gene),
            'combined_transcripts': len(seen), 'combined_genes': len(gene_counts),
            'combined_sequence_bases': total_bases, 'source_counts': source_counts,
            'identical_duplicates_removed': len(duplicates), 'conflicting_duplicates': 0,
            'fasta_gtf_join_complete': True, 'transcript_ids_include_versions': True,
            'transcript_biotype_counts': biotypes,
            'gtf_transcripts_not_in_combined_fasta': len(set(mapping) - set(seen)),
            'primary_gene_coverage': [{'gene_id': g, 'n_transcripts': gene_counts.get(g, 0)} for g in plan['primary_mouse_gene_ids']],
            'all_15_primary_genes_covered': len(plan['primary_mouse_gene_ids']) == 15 and not missing_primary,
            'outputs': [{'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)} for p in [combined, tx2gene]]}


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    REF.mkdir(parents=True, exist_ok=True)
    plan = json.loads(PLAN_PATH.read_text(encoding='utf-8-sig'))
    result = {'started_utc': utc(), 'python_executable': sys.executable, 'python_version': sys.version,
              'plan_sha256': sha(PLAN_PATH), 'script_sha256': sha(Path(__file__)),
              'release': 116, 'assembly': 'GRCm39', 'species': 'Mus musculus', 'passed': False,
              'kallisto_index_built_here': False, 'skill_used': 'validate-data'}
    try:
        paths, records = download(plan)
        result['raw_sources'] = records
        mapping, gtf = parse_gtf(paths['gtf'])
        result['gtf'] = gtf
        result.update(combine(paths, mapping, plan))
        result['passed'] = True
    except Exception:
        result['error'] = traceback.format_exc()
        events('failure', error=result['error'])
        raise
    finally:
        result['completed_utc'] = utc()
        write_json(HERE / 'reference_validation.json', result)
        if result['passed']:
            write_json(REF / 'reference_validation.json', result)
    print(json.dumps({k: result[k] for k in ['passed', 'combined_transcripts', 'combined_genes', 'identical_duplicates_removed', 'all_15_primary_genes_covered']}, ensure_ascii=False))


if __name__ == '__main__':
    main()
