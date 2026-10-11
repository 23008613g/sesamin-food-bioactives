"""Inspect final figure files; PDF tool runtime is separate from analysis."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
QA = Path('<LOCAL_PATH>')
PDF_PYTHON = '<LOCAL_PATH>'


def main():
    commands = [('source', [sys.executable, '-X', 'utf8', str(QA / 'validate_figure.py'),
                            str(HERE / 'plot_results.py'), '--json'])]
    for name in ['sample_pca', 'known_gene_panel']:
        commands.extend([
            (name + '.glyphs', [sys.executable, '-X', 'utf8', str(QA / 'audit_pdf_text.py'),
                               str(HERE / (name + '.pdf')), '--min-pt', '5', '--json']),
            (name + '.collision', [PDF_PYTHON, '-X', 'utf8', str(QA / 'audit_figure_collisions.py'),
                                  str(HERE / (name + '.pdf')), '--json-out', str(HERE / (name + '.collision.json'))])])
    checks = []
    for label, cmd in commands:
        result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8')
        (HERE / (label + '.audit.log')).write_text(result.stdout + '\n' + result.stderr, encoding='utf-8')
        checks.append({'label': label, 'command': cmd, 'exit_code': result.returncode})
    record = {'executed_utc': datetime.now(timezone.utc).isoformat(), 'checks': checks,
              'all_commands_exit0': all(x['exit_code'] == 0 for x in checks),
              'plotting_runtime': sys.executable, 'pdf_geometry_only_runtime': PDF_PYTHON,
              'scope': 'Existing PyMuPDF 1.27.2.2 used only for PDF geometry. Analysis/rendering remain in project venv; no installs.',
              'previous_failure': 'Static source validator could not recognize loop-built vector suffixes and mm conversions; explicit exports/inch dimensions resolved flags.',
              'manual_visual_review': 'Tracked separately in visual_review.json after inspecting the final PNG files.',
              'warnings_reviewed': ['PNG 300 dpi is a report preview; editable PDF/SVG also supplied. No TIFF needed for this exploratory report.',
                                    'LOG-GUARD is a static false positive on the axis label; code checks nonnegative normalized values and reconstructs log2(count+1).',
                                    'Single-panel figures: multi-panel alignment is not applicable, not a pass for an untested multi-panel layout.'],
              'figure_hashes': {f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                                for f in HERE.iterdir() if f.suffix in ['.pdf', '.svg', '.png']}}
    (HERE / 'figure_qa.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'all_commands_exit0': record['all_commands_exit0'], 'checks': len(checks)}))
    if not record['all_commands_exit0']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
