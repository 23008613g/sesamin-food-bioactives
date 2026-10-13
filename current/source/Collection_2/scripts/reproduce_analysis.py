"""Offline reproduction of frozen rat, limma and cross-species analyses.

Original analysis scripts are loaded unchanged; only their input/output globals
are redirected. An already-existing output directory is always rejected.
"""
from __future__ import annotations

import argparse
import contextlib
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import traceback

sys.dont_write_bytecode = True


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f'Cannot load analysis module: {path}')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def compare_tables(output, reference):
    import numpy as np
    import pandas as pd

    rows = []
    tables = sorted(p for p in (output / 'results').rglob('*')
                    if p.is_file() and (p.name.endswith('.csv') or p.name.endswith('.csv.gz')))
    if not tables:
        raise AssertionError('No reproduced result tables found')
    for actual_path in tables:
        rel = actual_path.relative_to(output)
        expected_path = reference / rel
        if not expected_path.is_file():
            raise FileNotFoundError(f'Missing reference table: {rel}')
        # Read identifiers, categories, booleans and missing markers as strings.
        # Exactly matching columns need no numeric coercion; any differing
        # column must be an ordinary numeric measure, never an identifier.
        actual = pd.read_csv(actual_path, dtype=str, keep_default_na=False)
        expected = pd.read_csv(expected_path, dtype=str, keep_default_na=False)
        if actual.shape != expected.shape or actual.columns.tolist() != expected.columns.tolist():
            raise AssertionError(f'Row/column structure differs: {rel}')
        numeric_errors = {}
        for col in actual.columns:
            a, b = actual[col], expected[col]
            if a.equals(b):
                continue
            key = col.lower()
            if key in {'rat_entrez', 'id', 'id_ref', 'sample', 'group'} or key.endswith('_id') or 'sha256' in key:
                raise AssertionError(f'Identifier/order differs: {rel} / {col}')
            try:
                av = pd.to_numeric(a.replace({'NA': np.nan, '': np.nan, 'NaN': np.nan, 'nan': np.nan}), errors='raise').to_numpy(dtype=float)
                bv = pd.to_numeric(b.replace({'NA': np.nan, '': np.nan, 'NaN': np.nan, 'nan': np.nan}), errors='raise').to_numpy(dtype=float)
            except (ValueError, TypeError) as exc:
                raise AssertionError(f'Non-numeric contents differ: {rel} / {col}') from exc
            if not np.allclose(av, bv, rtol=1e-10, atol=1e-10, equal_nan=True):
                raise AssertionError(f'Numeric result differs: {rel} / {col}')
            finite = np.isfinite(av) & np.isfinite(bv)
            numeric_errors[col] = float(np.max(np.abs(av[finite] - bv[finite]))) if finite.any() else 0.0
        rows.append({'path': rel.as_posix(), 'rows': len(actual), 'columns': len(actual.columns),
                     'status': 'passed', 'identifier_and_row_order': 'exact',
                     'max_abs_numeric_error': max(numeric_errors.values(), default=0.0),
                     'numeric_columns_with_text_rounding_differences': numeric_errors})
    return {'status': 'passed', 'reference': str(reference), 'n_tables': len(rows),
            'comparison': 'All CSV/CSV.GZ tables produced by rat, limma and cross-species main functions; gzip timestamps are ignored.',
            'rtol': 1e-10, 'atol': 1e-10,
            'max_abs_numeric_error': max(row['max_abs_numeric_error'] for row in rows),
            'tables': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True, help='A directory that does not yet exist')
    parser.add_argument('--rscript', type=Path, default=Path('C:/Program Files/R/R-4.6.1/bin/Rscript.exe'),
                        help='Existing Rscript executable, R 4.6.1 with limma 3.68.4; no installation is performed')
    parser.add_argument('--reference-dir', type=Path, help='Reference bundle for numeric comparison; defaults to this bundle')
    args = parser.parse_args()
    bundle = Path(__file__).resolve().parents[1]
    output = args.output_dir.expanduser().resolve()
    reference = (args.reference_dir or bundle).expanduser().resolve()
    rscript = args.rscript.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite an existing directory: {output}')
    if output == bundle or bundle in output.parents:
        raise ValueError('Choose a new output directory outside the frozen bundle')
    if not rscript.is_file():
        raise FileNotFoundError('Pass --rscript with an existing Rscript executable; no software will be installed')
    # Fail before creating output if Python dependencies or source inputs are missing.
    import numpy as np
    import pandas as pd
    import scipy

    project_source = bundle / 'sources/project_inputs'
    project_manifest = json.loads((project_source / 'source_manifest.json').read_text(encoding='utf-8'))
    project_roles = {}
    for item in project_manifest['files']:
        path = project_source / item['relative_path']
        if not path.is_file() or sha(path) != item['sha256']:
            raise ValueError(f'Bundled original/archived input hash mismatch: {path}')
        project_roles[item['role']] = item['relative_path']
    needed_roles = {'rat_matrix', 'platform_annotation', 'human_count', 'archived_all_gene_descriptive', 'archived_sample_qc'}
    if set(project_roles) != needed_roles:
        raise ValueError('Bundled source manifest must contain exactly the five fixed original/archived inputs')
    freeze = json.loads((bundle / 'protocol_freeze.json').read_text(encoding='utf-8'))
    addendum = json.loads((bundle / 'PROTOCOL_ADDENDUM_1.freeze.json').read_text(encoding='utf-8'))
    if sha(bundle / 'PROTOCOL.txt') != freeze['protocol_sha256'] or sha(bundle / 'PROTOCOL_ADDENDUM_1.txt') != addendum['sha256']:
        raise ValueError('Protocol freeze hashes do not match')
    go_manifest = json.loads((bundle / 'sources/go/source_manifest.json').read_text(encoding='utf-8'))
    if sha(bundle / 'sources/go/go-basic.obo') != go_manifest['sha256']:
        raise ValueError('Frozen GO ontology hash mismatch')
    orthology = json.loads((bundle / 'sources/orthology/orthology_freeze_report.json').read_text(encoding='utf-8'))
    for item in orthology['outputs']:
        if sha(bundle / 'sources/orthology' / item['file']) != item['sha256']:
            raise ValueError(f'Frozen orthology hash mismatch: {item["file"]}')

    output.mkdir(parents=True, exist_ok=False)
    (output / 'scripts').mkdir()
    (output / 'review').mkdir()
    (output / 'logs').mkdir()
    record_path = output / 'review/portable_reproduction.json'
    record = {'status': 'started', 'started_utc': utc(), 'command': sys.argv,
              'bundle': str(bundle), 'output_dir': str(output), 'reference_dir': str(reference),
              'python': sys.executable, 'python_version': platform.python_version(),
              'versions': {'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
              'network_access_in_wrapper': False, 'dependency_installations': False,
              'scripts_modified': False, 'wrapper_sha256': sha(Path(__file__)), 'phases': []}
    write_json(record_path, record)
    try:
        for name in ('PROTOCOL.txt', 'protocol_freeze.json', 'PROTOCOL_ADDENDUM_1.txt', 'PROTOCOL_ADDENDUM_1.freeze.json'):
            shutil.copy2(bundle / name, output / name)
        for name in ('go', 'orthology', 'project_inputs'):
            shutil.copytree(bundle / 'sources' / name, output / 'sources' / name)
        script_names = ('analyze_rat_genome.py', 'limma_sensitivity.R', 'compare_cross_species.py')
        for name in script_names:
            shutil.copy2(bundle / 'scripts' / name, output / 'scripts' / name)
            if sha(bundle / 'scripts' / name) != sha(output / 'scripts' / name):
                raise AssertionError(f'Copied script hash mismatch: {name}')
        copied = [p for p in (output / 'sources').rglob('*') if p.is_file()]
        record['copied_inputs'] = [{'path': p.relative_to(output).as_posix(), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(copied)]
        record['scripts'] = [{'path': name, 'sha256': sha(output / 'scripts' / name)} for name in script_names]
        for item in record['copied_inputs']:
            if sha(bundle / item['path']) != item['sha256']:
                raise AssertionError(f'Copied source differs: {item["path"]}')
        project = output / 'sources/project_inputs'
        rat = load_module('stage39_reproduction_rat', output / 'scripts/analyze_rat_genome.py')
        rat.ROOT, rat.OUT, rat.PROJECT, rat.NEXT = output, output / 'results', project, project / 'next_steps'
        rat.INPUTS = {'rat_matrix': project / project_roles['rat_matrix'],
                      'platform_annotation': project / project_roles['platform_annotation'],
                      'go_ontology': output / 'sources/go/go-basic.obo',
                      'protocol': output / 'PROTOCOL.txt', 'protocol_addendum_1': output / 'PROTOCOL_ADDENDUM_1.txt'}
        with (output / 'logs/rat_stdout.log').open('w', encoding='utf-8') as stdout, (output / 'logs/rat_stderr.log').open('w', encoding='utf-8') as stderr:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                rat.main()
        record['phases'].append({'name': 'rat', 'status': 'completed', 'entrypoint': 'unchanged analyze_rat_genome.main(); redirected ROOT/OUT/PROJECT/NEXT/INPUTS'})
        write_json(record_path, record)

        command = [str(rscript), str(output / 'scripts/limma_sensitivity.R'), str(output)]
        r = subprocess.run(command, capture_output=True, check=False)
        (output / 'logs/limma_stdout.log').write_bytes(r.stdout)
        (output / 'logs/limma_stderr.log').write_bytes(r.stderr)
        record['phases'].append({'name': 'limma', 'status': 'completed' if r.returncode == 0 else 'failed',
                                 'command': command, 'returncode': r.returncode})
        write_json(record_path, record)
        if r.returncode:
            raise RuntimeError(f'Rscript failed with return code {r.returncode}; see logs/limma_stderr.log')
        record['R_sessionInfo'] = (output / 'results/limma/sessionInfo.txt').read_text(encoding='utf-8')

        cross = load_module('stage39_reproduction_cross', output / 'scripts/compare_cross_species.py')
        cross.ROOT, cross.OUT, cross.REVIEW = output, output / 'results/cross_species', output / 'review/cross_species_execution.json'
        cross.PROJECT = project
        with (output / 'logs/cross_stdout.log').open('w', encoding='utf-8') as stdout, (output / 'logs/cross_stderr.log').open('w', encoding='utf-8') as stderr:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                cross.main(project)
        record['phases'].append({'name': 'cross_species', 'status': 'completed', 'entrypoint': 'unchanged compare_cross_species.main(project); redirected ROOT/OUT/REVIEW/PROJECT'})
        write_json(record_path, record)
        record['numeric_comparison'] = compare_tables(output, reference)
        record['source_bundle_inputs_unchanged'] = all(sha(bundle / item['path']) == item['sha256'] for item in record['copied_inputs'])
        if not record['source_bundle_inputs_unchanged']:
            raise AssertionError('An input in the source bundle changed during reproduction')
        record['status'], record['completed_utc'] = 'completed', utc()
        record['scope'] = 'Complete rat/limma/cross-species computation and all generated table comparisons. Figure packaging and separate original independent-audit scripts are outside this wrapper.'
        write_json(record_path, record)
        print(json.dumps({'status': record['status'], 'output_dir': str(output),
                          'n_result_tables_verified': record['numeric_comparison']['n_tables'],
                          'max_abs_numeric_error': record['numeric_comparison']['max_abs_numeric_error'],
                          'record': str(record_path)}, indent=2))
    except Exception as exc:
        record.update(status='failed', failed_utc=utc(), error=repr(exc), traceback=traceback.format_exc())
        write_json(record_path, record)
        raise


if __name__ == '__main__':
    main()
