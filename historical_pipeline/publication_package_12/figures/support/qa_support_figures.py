"""Audit final support figures without replacing the project drawing environment."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
from PIL import Image

HERE = Path(__file__).resolve().parent
SKILL = Path.home() / '.agents/skills/nature-figure/scripts'
PDF_PYTHON = Path.home() / 'AppData/Local/Programs/Python/Python313/python.exe'
STEMS = ['fig1_design', 'figS1_quality', 'figS2_sensitivity']

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(command, log, json_output=False):
    result = subprocess.run([str(x) for x in command], capture_output=True, encoding='utf-8',
                            errors='replace', env={**__import__('os').environ,'PYTHONUTF8':'1'})
    (HERE/log).write_text(result.stdout+result.stderr,encoding='utf-8')
    return {'command':[str(x) for x in command], 'exit_code':result.returncode,
            'log':log, 'result':json.loads(result.stdout) if json_output and result.stdout.startswith('{') else None}

def main():
    checks = []
    checks.append(run([sys.executable,SKILL/'validate_figure.py',HERE/'make_support_figures.py','--json'],
                      'source_preflight.json',True))
    geometry_code = "import fitz,json,sys; from pathlib import Path; p=Path(sys.argv[1]); print(json.dumps({f.stem:{'page_size_mm':[fitz.open(f)[0].rect.width*25.4/72,fitz.open(f)[0].rect.height*25.4/72],'selectable_characters':len(fitz.open(f)[0].get_text())} for f in p.glob('*.pdf') if not f.name.endswith('overlay.pdf')}))"
    geometry = run([PDF_PYTHON,'-c',geometry_code,HERE], 'pdf_geometry.json', True)
    checks.append(geometry)
    figure_checks = {}
    for stem in STEMS:
        checks.append(run([sys.executable,SKILL/'audit_pdf_text.py',HERE/f'{stem}.pdf','--min-pt','5','--json'],
                          f'{stem}.text-audit.json',True))
        checks.append(run([PDF_PYTHON,SKILL/'audit_figure_collisions.py',HERE/f'{stem}.pdf',
                          '--json-out',HERE/f'{stem}.collision.json',
                          '--overlay-pdf',HERE/f'{stem}.collision-overlay.pdf'],f'{stem}.collision.log'))
        with Image.open(HERE/f'{stem}.tiff') as im:
            tiff = {'pixels':list(im.size),'dpi':[float(x) for x in im.info['dpi']],
                    'compression':im.info.get('compression')}
        assert tiff['compression'] == 'tiff_lzw'
        assert all(abs(d-600)<.01 for d in tiff['dpi'])
        assert abs(tiff['pixels'][0]/600*25.4-180)<.05
        assert abs(geometry['result'][stem]['page_size_mm'][0]-180)<.01
        assert geometry['result'][stem]['selectable_characters'] > 50
        figure_checks[stem]={'raster':tiff,
            'alignment':json.loads((HERE/f'{stem}.alignment.json').read_text(encoding='utf-8')),
            'export_hashes':{suffix:sha(HERE/f'{stem}.{suffix}') for suffix in ['png','pdf','svg','tiff']}}
    record = json.loads((HERE/'execution_record.json').read_text(encoding='utf-8'))
    project = HERE.parents[3]
    input_checks = [{**item, 'unchanged':sha(project/item['path']) == item['sha256']}
                    for item in record['input_fingerprints']]
    assert all(item['unchanged'] for item in input_checks)
    output = {'date':'2026-10-09','passed':all(c['exit_code']==0 for c in checks),
              'plotting_interpreter':sys.executable,
              'pdf_collision_auditor_interpreter':str(PDF_PYTHON),
              'independent_auditor_scope':'Existing PyMuPDF environment used only for PDF geometry checks; no drawing or dependency installation.',
              'commands':checks,'figure_checks':figure_checks,'input_checks':input_checks,
              'source_preflight_attempt_note':'Initial preflight rejected the explicit DejaVu Sans family because it requires a generic sans-serif family declaration. rcParams now declares sans-serif with DejaVu Sans fallback; rendered typeface is unchanged. Final width is verified from raster metadata and PDF page geometry. First QA summary serialization failed on the TIFF rational DPI object after the audits completed; casting it to float repaired the report writer without changing figures.',
              'visual_review_file':'visual_review.md'}
    (HERE/'validation.json').write_text(json.dumps(output,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'passed':output['passed'],'audits':len(checks),'input_hashes':len(input_checks)}))
    if not output['passed']:
        raise SystemExit(1)

if __name__ == '__main__':
    main()
