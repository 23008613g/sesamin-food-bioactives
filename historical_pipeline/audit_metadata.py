"""Bounded public GEO metadata audit; no differential-expression analysis.

Reuses the project's standard-library header reader without modifying it.
All retrieved studies are development/exploration candidates, not blind tests.
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--accessions', nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent.parent / 'audit_candidate_headers.py'
    spec = importlib.util.spec_from_file_location('existing_header_audit', source)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    accessions = list(dict.fromkeys(args.accessions))
    if not all(a.startswith('GSE') and a[3:].isdigit() for a in accessions):
        parser.error('Each accession must be GSE followed by digits')
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(helper.audit, accessions))
    result = {
        'retrieved_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'decision_date_taipei': '2026-10-07',
        'command': sys.argv,
        'interpreter': sys.executable,
        'python_version': platform.python_version(),
        'environment_scope': 'Existing stdlib tool environment, not an analysis environment',
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'helper_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'scope': 'Series-matrix metadata only, stops before expression rows. Not QC or validation.',
        'split_status': 'All records are exploration; abstracts and sample metadata accessed.',
        'records': records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps([{k:r.get(k) for k in ['accession','status','sample_count','error']} for r in records], ensure_ascii=True, indent=2))


if __name__ == '__main__':
    main()
