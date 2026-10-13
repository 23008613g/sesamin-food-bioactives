"""Rebuild all three bounded analyses from package-local inputs, then compare tables."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results'),
                        help='New output folder, relative to the package unless absolute.')
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    output = output.resolve()
    if output.exists():
        raise SystemExit('Output already exists; supply --output with a new folder name. No overwrite was performed.')
    manifest = json.loads((ROOT/'provenance/source_manifest.json').read_text(encoding='utf-8'))
    for record in manifest:
        path = ROOT / record['file']
        if sha(path) != record['sha256']:
            raise AssertionError('Input/protocol/reference hash mismatch: ' + record['file'])
    adaptations = json.loads((ROOT/'provenance/script_adaptations.json').read_text(encoding='utf-8'))
    for record in adaptations:
        assert sha(ROOT/'code'/record['portable']) == record['portable_sha256'], record['portable']
    output.mkdir(parents=True)
    (output/'logs').mkdir()
    started = datetime.now(timezone.utc).isoformat()
    runs = []
    for name, module, script in [('human_mouse','robustness','robustness.py'),
                                 ('human_mouse_formula_check','robustness','robustness_check.py'),
                                 ('rat_welch','rat','rat.py'),
                                 ('rat_exact','rat','rat_exact.py'),
                                 ('organoid','organoid','organoid.py')]:
        folder = output/module
        folder.mkdir(exist_ok=True)
        for protocol in (ROOT/'protocols'/module).iterdir():
            if not (folder/protocol.name).exists():
                shutil.copy2(protocol, folder/protocol.name)
        command = [sys.executable, str(ROOT/'code'/script), str(folder)]
        begin = time.perf_counter()
        with (output/'logs'/(name+'.stdout.txt')).open('w', encoding='utf-8') as stdout, \
             (output/'logs'/(name+'.stderr.txt')).open('w', encoding='utf-8') as stderr:
            execution = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr)
        runs.append({'module': name, 'command': command, 'return_code': execution.returncode,
                     'seconds': time.perf_counter()-begin, 'script_sha256': sha(ROOT/'code'/script)})
        print(f'{name}: return code {execution.returncode}', flush=True)
        if execution.returncode:
            raise RuntimeError(f'{name} failed; see {output}/logs/{name}.stderr.txt')
    command = [sys.executable, str(ROOT/'code/external/analyze_external.py'),
               '--input-dir', str(ROOT/'inputs/external'), '--output-dir', str(output/'external')]
    begin = time.perf_counter()
    with (output/'logs/external.stdout.txt').open('w', encoding='utf-8') as stdout, (output/'logs/external.stderr.txt').open('w', encoding='utf-8') as stderr:
        execution = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr)
    runs.append({'module':'external_atlas','command':command,'return_code':execution.returncode,
                 'seconds':time.perf_counter()-begin,'script_sha256':sha(ROOT/'code/external/analyze_external.py')})
    print(f'external_atlas: return code {execution.returncode}', flush=True)
    if execution.returncode:
        raise RuntimeError('External atlas failed; see output logs.')
    command = [sys.executable, str(ROOT/'code/external/verify_results.py'),
               '--input-dir', str(ROOT/'inputs/external'), '--output-dir', str(output/'external')]
    begin = time.perf_counter()
    with (output/'logs/external_focused.stdout.txt').open('w', encoding='utf-8') as stdout, (output/'logs/external_focused.stderr.txt').open('w', encoding='utf-8') as stderr:
        execution = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr)
    runs.append({'module':'external_focused_check','command':command,'return_code':execution.returncode,
                 'seconds':time.perf_counter()-begin,'script_sha256':sha(ROOT/'code/external/verify_results.py')})
    print(f'external_focused_check: return code {execution.returncode}', flush=True)
    if execution.returncode:
        raise RuntimeError('External focused check failed; see output logs.')
    command = [sys.executable, str(ROOT/'verify.py'), '--results', str(output)]
    subprocess.run(command, cwd=ROOT, check=True)
    assert all(sha(ROOT/r['file']) == r['sha256'] for r in manifest)
    execution_record = {
        'status':'PASS_BOUNDED_PROCESSED_DATA_REPRODUCTION',
        'started_utc': started, 'completed_utc': datetime.now(timezone.utc).isoformat(),
        'package_root': str(ROOT), 'output_root': str(output), 'runs': runs,
        'environment': {'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
                        **{name:importlib.metadata.version(name) for name in ['numpy','pandas','scipy']}},
        'input_protocol_reference_hashes_unchanged':True,
        'source_project_required_at_runtime':False,'internet_used':False,
        'installations_performed':False,'raw_read_reprocessing':False,
        'single_cell_all_gene_UMI_reprocessing':False,
        'scope':'23 packaged inputs; quantitative rerun from processed matrices, organoid cell-level candidate counts, and the external author-normalized adult liver matrix. See README.md.',
    }
    (output/'execution.json').write_text(json.dumps(execution_record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: all reference tables compared; see verification.json and execution.json.',flush=True)


if __name__ == '__main__':
    main()
