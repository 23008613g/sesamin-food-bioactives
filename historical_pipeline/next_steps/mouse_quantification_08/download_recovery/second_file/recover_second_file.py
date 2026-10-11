"""Recover only DRR1010540 mate2, preserving the original failed download.

Uses existing project Python and aria2. Writes only this script directory and the
dedicated DATA/download_recovery/DRR1010540_2.fastq.gz file/control/log paths.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback
import urllib.request

HERE = Path(__file__).resolve().parent
STAGE = HERE.parents[1]
DATA = Path('D:/ResearchData/sesamin_mouse_2026-10-08')
NAME = 'DRR1010540_2.fastq.gz'
OLD = DATA / 'fastq' / NAME
NEW = DATA / 'download_recovery' / NAME
ENA = ('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=DRR1010540&result=read_run'
       '&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv')


def now():
    return datetime.now(timezone.utc).isoformat()


def write(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def hashes(path):
    before = path.stat()
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            md5.update(block)
            sha.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError('File changed during independent hash: ' + str(path))
    return {'actual_bytes': before.st_size, 'actual_md5': md5.hexdigest(),
            'actual_sha256': sha.hexdigest(), 'mtime_ns': before.st_mtime_ns}


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    NEW.parent.mkdir(parents=True, exist_ok=True)
    existing = list(HERE.glob('attempt_*.json'))
    attempt_id = max([int(p.stem.split('_')[-1]) for p in existing] or [0]) + 1
    key = f'recovery_DRR1010540_2_attempt{attempt_id}'
    rec_path = HERE / f'attempt_{attempt_id:02d}.json'
    log_path = NEW.parent / f'{NAME}.attempt{attempt_id:02d}.log'
    rec = {'execution_key': key, 'started_utc': now(), 'status': 'preparing',
           'python_executable': sys.executable, 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'raw_original_preserved': str(OLD), 'recovery_target': str(NEW),
           'log': str(log_path), 'actual_pid': None}
    write(rec_path, rec)
    try:
        manifest = json.loads((STAGE.parent / 'mouse_readiness_07/fastq_manifest.json').read_text(encoding='utf-8'))
        candidates = [e for e in manifest['entries'] if e['filename'] == NAME]
        assert len(candidates) == 1
        expected = candidates[0]
        request = urllib.request.Request(ENA, headers={'User-Agent': 'research-data-integrity/1.0'})
        with urllib.request.urlopen(request, timeout=90) as response:
            body = response.read()
            fresh = {'url': ENA, 'fetched_utc': now(), 'status': response.status,
                     'final_url': response.url, 'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest()}
        (HERE / f'ena_DRR1010540_attempt{attempt_id:02d}.tsv').write_bytes(body)
        rows = list(csv.DictReader(io.StringIO(body.decode('utf-8-sig')), delimiter='\t'))
        assert len(rows) == 1 and rows[0]['run_accession'] == 'DRR1010540'
        data = rows[0]
        triples = list(zip(data['fastq_ftp'].split(';'), data['fastq_md5'].split(';'), data['fastq_bytes'].split(';'), strict=True))
        matches = [(url, md5, int(size)) for url, md5, size in triples if url.rsplit('/', 1)[-1] == NAME]
        assert len(matches) == 1
        url, official_md5, official_bytes = matches[0]
        url = url if url.startswith('https://') else 'https://' + url.removeprefix('ftp://')
        fresh.update(file=NAME, official_url=url, official_md5=official_md5, official_bytes=official_bytes,
                     matches_frozen_manifest=(official_md5 == expected['expected_md5'] and official_bytes == expected['bytes'] and url == expected['url']))
        write(HERE / f'fresh_ena_metadata_attempt{attempt_id:02d}.json', fresh)
        assert fresh['matches_frozen_manifest'], 'Fresh official ENA metadata disagrees with frozen source'
        failure_path = HERE / 'failure_record.json'
        if not failure_path.exists():
            failed = hashes(OLD)
            failed.update(file=str(OLD), expected_md5=official_md5, expected_bytes=official_bytes,
                          recorded_utc=now(), checksum_matches_official=failed['actual_md5'] == official_md5,
                          cause='Unknown cause of checksum mismatch. No attribution to server, network, software or storage established.',
                          raw_file_modified=False, entered_QC=False)
            write(failure_path, failed)
        else:
            failed = json.loads(failure_path.read_text(encoding='utf-8'))
            assert failed['file'] == str(OLD)
            assert OLD.stat().st_size == failed['actual_bytes'] and OLD.stat().st_mtime_ns == failed['mtime_ns'], 'Preserved original changed'
        rec['failed_original_sha256'] = failed['actual_sha256']
        rec['fresh_official_md5'] = official_md5
        rec['fresh_official_bytes'] = official_bytes
        exe = shutil.which('aria2c')
        assert exe, 'Existing aria2 unavailable; no installation permitted here'
        args = [exe, '--dir=' + str(NEW.parent), '--out=' + NAME,
                '--max-concurrent-downloads=1', '--max-connection-per-server=4', '--split=4',
                '--continue=true', '--allow-overwrite=false', '--auto-file-renaming=false',
                '--file-allocation=none', '--check-certificate=true', '--check-integrity=true',
                '--checksum=md5=' + official_md5, '--max-tries=5', '--retry-wait=5',
                '--connect-timeout=30', '--timeout=60', '--summary-interval=10',
                '--console-log-level=notice', url]
        rec.update(status='running', argv=args, download_started_utc=now())
        write(rec_path, rec)
        started = time.monotonic()
        with log_path.open('wb') as log:
            process = subprocess.Popen(args, stdout=log, stderr=subprocess.STDOUT)
            rec['actual_pid'] = process.pid
            write(rec_path, rec)
            print(json.dumps({'execution_key': key, 'pid': process.pid, 'log': str(log_path)}), flush=True)
            code = process.wait()
        rec.update(exit_code=code, download_completed_utc=now(), wall_seconds=time.monotonic() - started)
        write(rec_path, rec)
        if code != 0:
            raise RuntimeError(f'aria2 returned {code}; partial recovery and log retained')
        current = hashes(NEW)
        assert current['actual_bytes'] == official_bytes, 'Recovered size mismatch'
        assert current['actual_md5'] == official_md5, 'Independent recovered MD5 mismatch'
        assert not Path(str(NEW) + '.aria2').exists(), 'aria2 control file still present'
        assert OLD.stat().st_size == failed['actual_bytes'] and OLD.stat().st_mtime_ns == failed['mtime_ns'], 'Original changed during recovery'
        replacement = {'passed': True, 'new_file': str(NEW), 'actual_bytes': current['actual_bytes'],
                       'actual_md5': current['actual_md5'], 'sha256': current['actual_sha256'],
                       'official_md5': official_md5, 'execution_key': key,
                       'verified_utc': now(), 'official_url': url, 'official_bytes': official_bytes,
                       'independent_full_file_md5_sha256_completed': True,
                       'raw_original_file': str(OLD), 'raw_original_unchanged_by_size_mtime': True,
                       'no_qc_or_quantification_executed': True}
        write(HERE / 'replacement_integrity.json', replacement)
        rec.update(status='passed', completed_utc=now(), replacement_integrity=str(HERE / 'replacement_integrity.json'))
        write(rec_path, rec)
        print(json.dumps(replacement, ensure_ascii=False), flush=True)
    except Exception:
        rec.update(status='failed', completed_utc=now(), exception=traceback.format_exc())
        write(rec_path, rec)
        raise


if __name__ == '__main__':
    main()
