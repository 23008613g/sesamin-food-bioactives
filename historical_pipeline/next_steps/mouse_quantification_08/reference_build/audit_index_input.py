"""Locate ambiguity bases and long terminal poly-A runs without modifying input."""
import collections
import csv
import json
from pathlib import Path
import re
import sys
import urllib.request
import build_reference as b

plan = json.loads(b.PLAN_PATH.read_text(encoding='utf-8-sig'))
primary = set(plan['primary_mouse_gene_ids'])
with (b.REF / 'tx2gene.tsv').open(encoding='utf-8', newline='') as f:
    rows = {r['transcript_id']: r for r in csv.DictReader(f, delimiter='\t')}
ambiguous = []
poly_a = []
head, seq = None, []


def inspect():
    if head is None:
        return
    tid = head.split()[0]
    sequence = ''.join(seq)
    for match in re.finditer('[^ACGUTacgut]', sequence):
        row = rows[tid]
        ambiguous.append({'transcript_id': tid, 'gene_id': row['gene_id'], 'gene_name': row['gene_name'],
                          'primary_candidate': row['gene_id'] in primary, 'one_based_sequence_position': match.start() + 1,
                          'base': match.group(), 'sequence_length': len(sequence)})
    tail = len(sequence) - len(sequence.rstrip('Aa'))
    if tail >= 10:
        poly_a.append({'transcript_id': tid, 'gene_id': rows[tid]['gene_id'],
                       'primary_candidate': rows[tid]['gene_id'] in primary, 'tail_length': tail})


with (b.REF / 'combined.fa').open(encoding='ascii') as f:
    for line in f:
        if line.startswith('>'):
            inspect()
            head, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
inspect()
source_url = 'https://raw.githubusercontent.com/pachterlab/kallisto/v0.51.1/src/KmerIndex.cpp'
source_path = b.HERE / 'kallisto_v0.51.1_KmerIndex.cpp'
if not source_path.exists():
    with urllib.request.urlopen(source_url, timeout=60) as response:
        source_path.write_bytes(response.read())
source_code = source_path.read_text(encoding='utf-8')
assert 'std::mt19937 gen(42);' in source_code
assert 'str.size() >= 10 && str.substr(str.size()-10,10) == "AAAAAAAAAA"' in source_code
b.write_json(b.HERE / 'index_input_sequence_audit.json', {
    'utc': b.utc(), 'python_executable': sys.executable, 'script_sha256': b.sha(Path(__file__)),
    'source_fasta_sha256': b.sha(b.REF / 'combined.fa'),
    'non_acgut_characters': ambiguous, 'non_acgut_character_count': len(ambiguous),
    'terminal_poly_a_gt10_transcripts': sum(r['tail_length'] > 10 for r in poly_a),
    'terminal_poly_a_exactly10_transcripts': sum(r['tail_length'] == 10 for r in poly_a),
    'terminal_poly_a_ge10_transcripts': len(poly_a), 'poly_a_transcripts': poly_a,
    'pinned_source': {'url': source_url, 'local_path': str(source_path), 'sha256': b.sha(source_path),
                      'index_sequence_rng_seed': 42, 'index_sequence_rng': 'std::mt19937',
                      'poly_a_boundary_in_source': 'at least 10 trailing A bases, despite warning text longer than 10',
                      'verification_type': 'Static read of official pinned source, not compiled binary introspection'},
    'interpretation': 'This scans the original reference. The kallisto index process is executed by root separately. The official pinned source uses fixed seed 42 for ambiguous-base replacement and clips terminal runs of at least 10 A bases. Retain the fixed resulting index and SHA256 for all nine samples; do not claim original FASTA is byte-equivalent to indexed sequences or prove whole-index byte reproducibility across environments.',
    'reference_modified': False,
})
print(json.dumps({'non_acgut_characters': ambiguous, 'terminal_poly_a_ge10_transcripts': len(poly_a), 'terminal_poly_a_exactly10_transcripts': sum(r['tail_length'] == 10 for r in poly_a)}, ensure_ascii=False))
