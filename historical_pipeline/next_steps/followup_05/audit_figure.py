"""Check the current single-panel plot using existing local runtimes."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
QA = Path('<LOCAL_PATH>')
PDF_PYTHON = '<LOCAL_PATH>'


def main():
    commands = [
        ('source', [sys.executable, '-X', 'utf8', str(QA / 'validate_figure.py'), str(HERE / 'plot_phh.py'), '--json']),
        ('glyphs', [sys.executable, '-X', 'utf8', str(QA / 'audit_pdf_text.py'), str(HERE / 'phh_spotlight.pdf'), '--min-pt', '5', '--json']),
        ('collision', [PDF_PYTHON, '-X', 'utf8', str(QA / 'audit_figure_collisions.py'), str(HERE / 'phh_spotlight.pdf'), '--json-out', str(HERE / 'phh_spotlight.collision.json')])]
    checks = []
    for label, cmd in commands:
        result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8')
        (HERE / ('figure_' + label + '.log')).write_text(result.stdout + '\n' + result.stderr, encoding='utf-8')
        checks.append({'check': label, 'command': cmd, 'returncode': result.returncode})
    record = {'executed_utc': datetime.now(timezone.utc).isoformat(), 'checks': checks,
              'commands_passed': all(x['returncode'] == 0 for x in checks),
              'runtime_scope': 'Analysis/plotting in project venv; already installed PyMuPDF tool Python only inspects PDF geometry; no installation or system mutation',
              'warnings_reviewed': ['PNG 300 dpi is a report preview; editable vector formats also supplied, no submission TIFF needed.',
                                    'Log ratio is guarded by nonnegative checks and +1 in both numerator and denominator; static text-label warnings do not indicate an unguarded log.',
                                    'Single-panel alignment is not applicable; raw label effects shown without inferential intervals.'],
              'visual_review': 'Recorded after actual inspection in visual_review.json',
              'file_hashes': {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest() for name in
                             ['plot_phh.py', 'phh_donor_effects.csv', 'phh_spotlight.png', 'phh_spotlight.pdf', 'phh_spotlight.svg']}}
    (HERE / 'figure_qa.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'commands_passed': record['commands_passed'], 'codes': [x['returncode'] for x in checks]}))
    if not record['commands_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
