"""Independently re-read final FASTA/TSV and inventory annotation-only transcripts."""
import collections
import csv
import gzip
import json
from pathlib import Path
import re
import sys
import build_reference as build

HERE, REF = build.HERE, build.REF
validation_path = HERE / 'reference_validation.json'
original = json.loads(validation_path.read_text(encoding='utf-8'))
checks = []


def check(name, condition, actual=None):
    checks.append({'name': name, 'passed': bool(condition), 'actual': actual})
    if not condition:
        raise AssertionError(f'{name}: {actual}')


rows = {}
with (REF / 'tx2gene.tsv').open(encoding='utf-8', newline='') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        tid = row['transcript_id']
        if tid in rows:
            raise AssertionError('duplicate output transcript ' + tid)
        rows[tid] = row
seen = set()
lengths = {}
tid = None
with (REF / 'combined.fa').open(encoding='ascii') as f:
    for line in f:
        if line.startswith('>'):
            tid = line[1:].split()[0]
            if tid in seen or tid not in rows:
                raise AssertionError('FASTA duplicate or unjoined transcript ' + tid)
            if not re.fullmatch(r'ENSMUST\d+\.\d+', tid):
                raise AssertionError('unexpected transcript ID ' + tid)
            seen.add(tid)
            lengths[tid] = 0
        else:
            lengths[tid] += len(line.strip())
check('FASTA/TSV exact identifier set', seen == set(rows), len(seen))
check('FASTA/TSV exact sequence lengths', all(int(rows[t]['sequence_length']) == n for t, n in lengths.items()))
check('transcript count matches build', len(seen) == original['combined_transcripts'], len(seen))
genes = set(r['gene_id'] for r in rows.values())
check('gene count matches build', len(genes) == original['combined_genes'], len(genes))
check('base count matches build', sum(lengths.values()) == original['combined_sequence_bases'], sum(lengths.values()))
for output in original['outputs']:
    p = Path(output['path'])
    check('output bytes ' + p.name, p.stat().st_size == output['bytes'], p.stat().st_size)
    check('output sha256 ' + p.name, build.sha(p) == output['sha256'])
plan = json.loads(build.PLAN_PATH.read_text(encoding='utf-8-sig'))
gene_counts = collections.Counter(r['gene_id'] for r in rows.values())
for g in plan['primary_mouse_gene_ids']:
    check('primary gene covered ' + g, gene_counts[g] > 0, gene_counts[g])

gtf_path = REF / 'Mus_musculus.GRCm39.116.gtf.gz'
missing = []
joined = 0
with gzip.open(gtf_path, 'rt', encoding='utf-8') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if fields[2] != 'transcript':
            continue
        a = dict(build.ATTR.findall(fields[8]))
        t = a['transcript_id'] + '.' + a['transcript_version']
        if t not in rows:
            missing.append({'transcript_id': t, 'gene_id': a['gene_id'], 'gene_name': a.get('gene_name', ''),
                            'transcript_biotype': a.get('transcript_biotype', ''), 'seqname': fields[0],
                            'start': fields[3], 'end': fields[4]})
        else:
            r = rows[t]
            if any(r[k] != a[k] for k in ['gene_id', 'gene_version', 'gene_biotype', 'transcript_biotype']):
                raise AssertionError('GTF output row discrepancy ' + t)
            joined += 1
check('GTF join count equals all FASTA transcripts', joined == len(rows), joined)
check('GTF-only transcript count matches build', len(missing) == original['gtf_transcripts_not_in_combined_fasta'], len(missing))
with (HERE / 'gtf_transcripts_not_in_fasta.tsv').open('w', encoding='utf-8', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['transcript_id', 'gene_id', 'gene_name', 'transcript_biotype', 'seqname', 'start', 'end'], delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(missing)

result = {'passed': all(c['passed'] for c in checks), 'utc': build.utc(), 'python_executable': sys.executable,
          'checks': checks, 'check_count': len(checks), 'script_sha256': build.sha(Path(__file__)),
          'annotation_only_transcripts': {'count': len(missing),
                                        'genes_with_any_annotation_only_transcript': len(set(r['gene_id'] for r in missing)),
                                        'genes_absent_from_quantifiable_reference': len(set(r['gene_id'] for r in missing) - genes),
                                        'genes_partially_represented_in_reference': len(set(r['gene_id'] for r in missing) & genes),
                                        'biotypes': dict(collections.Counter(r['transcript_biotype'] for r in missing)),
                                        'seqnames': dict(collections.Counter(r['seqname'] for r in missing)),
                                        'interpretation': 'Annotation-only transcripts were not present in the supplied cDNA/ncRNA FASTA. No sequences were inferred or silently created; these are outside the quantifiable reference.'}}
build.write_json(HERE / 'output_validation.json', result)

# Separate inherited stage07 HEAD metadata from actual stage08 payload validation.
# Save the original record before making this metadata-only clarification.
before_path = HERE / 'reference_validation.before_metadata_clarification.json'
if not before_path.exists():
    build.write_json(before_path, original)
for source in original['raw_sources']:
    if 'prior_stage07_head_probe' not in source:
        source['prior_stage07_head_probe'] = {k: source.pop(k) for k in ['method', 'checked_utc', 'payload_bytes_downloaded', 'status', 'content_validated'] if k in source}
    source['payload_bytes_downloaded'] = source['actual_bytes']
    source['content_validated'] = True
    source['gzip_integrity_passed'] = True
    source['content_validation_method'] = 'Exact pinned size, local SHA256, full gzip stream parse, and release-matched structural checks; no upstream checksum supplied'
initial_download_records = HERE / 'download_records.before_metadata_clarification.json'
if not initial_download_records.exists():
    initial_download_records.write_bytes((HERE / 'download_records.json').read_bytes())
build.write_json(HERE / 'download_records.json', original['raw_sources'])
original['post_build_validation'] = {'path': str(HERE / 'output_validation.json'), 'sha256': build.sha(HERE / 'output_validation.json'),
                                   'passed': result['passed'], 'checks': len(checks)}
original['annotation_only_transcript_inventory'] = result['annotation_only_transcripts']
original['metadata_clarified_utc'] = build.utc()
original['metadata_clarification_script_sha256'] = result['script_sha256']
original['scope_note'] = 'Reference coverage is not measured expression or cross-species biological validation. No kallisto index or quantification was run by this reference task.'
build.write_json(validation_path, original)
build.write_json(REF / 'reference_validation.json', original)
print(json.dumps({'passed': result['passed'], 'checks': len(checks), 'annotation_only_transcripts': result['annotation_only_transcripts']}, ensure_ascii=False))
