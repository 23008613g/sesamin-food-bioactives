"""Check rendered figures, unchanged input hashes and exported numerical tables."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
from PIL import Image
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
SKILL=Path.home()/'.agents/skills/nature-figure/scripts'
PDF_PYTHON=Path.home()/'AppData/Local/Programs/Python/Python313/python.exe'
STEMS=['figS3_phh_all25','figS4_mouse_all15']

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def run(command,log,json_output=False):
    result=subprocess.run([str(x) for x in command],capture_output=True,encoding='utf-8',
                          errors='replace',env={**os.environ,'PYTHONUTF8':'1'})
    (HERE/log).write_text(result.stdout+result.stderr,encoding='utf-8')
    return {'command':[str(x) for x in command],'exit_code':result.returncode,'log':log,
            'result':json.loads(result.stdout) if json_output and result.stdout.startswith('{') else None}

def main():
    checks=[run([sys.executable,SKILL/'validate_figure.py',HERE/'make_new_support.py','--json'],'source_preflight.json',True)]
    code="import fitz,json,sys; from pathlib import Path; p=Path(sys.argv[1]); print(json.dumps({f.stem:{'page_size_mm':[fitz.open(f)[0].rect.width*25.4/72,fitz.open(f)[0].rect.height*25.4/72],'selectable_characters':len(fitz.open(f)[0].get_text())} for f in p.glob('*.pdf') if not f.name.endswith('overlay.pdf')}))"
    geometry=run([PDF_PYTHON,'-c',code,HERE],'pdf_geometry.json',True)
    checks.append(geometry)
    rendered={}
    for stem in STEMS:
        checks.append(run([sys.executable,SKILL/'audit_pdf_text.py',HERE/f'{stem}.pdf','--min-pt','5','--json'],f'{stem}.text-audit.json',True))
        checks.append(run([PDF_PYTHON,SKILL/'audit_figure_collisions.py',HERE/f'{stem}.pdf','--json-out',HERE/f'{stem}.collision.json','--overlay-pdf',HERE/f'{stem}.collision-overlay.pdf'],f'{stem}.collision.log'))
        with Image.open(HERE/f'{stem}.tiff') as im:
            info={'pixels':list(im.size),'dpi':[float(x) for x in im.info['dpi']],'compression':im.info.get('compression')}
        assert info['compression']=='tiff_lzw' and info['dpi']==[600,600]
        assert abs(info['pixels'][0]/600*25.4-180)<.05
        assert abs(geometry['result'][stem]['page_size_mm'][0]-180)<.01
        assert geometry['result'][stem]['selectable_characters']>500
        rendered[stem]={'raster':info,'PDF_geometry':geometry['result'][stem],
            'alignment':json.loads((HERE/f'{stem}.alignment.json').read_text(encoding='utf-8')),
            'exports':{ext:sha(HERE/f'{stem}.{ext}') for ext in ['png','pdf','svg','tiff']}}
    record=json.loads((HERE/'execution_record.json').read_text(encoding='utf-8'))
    hashes=[{**r,'unchanged':sha(Path(r['path']))==r['sha256']} for r in record['input_fingerprints']]
    assert all(r['unchanged'] for r in hashes)
    for entry in record['source_data_fingerprints']:
        assert sha(HERE/entry['file'])==entry['sha256']
    h=pd.read_csv(HERE/'source_data/figS3_phh_125_values.csv')
    m=pd.read_csv(HERE/'source_data/figS4_mouse_135_values.csv')
    assert len(h)==125 and len(m)==135
    assert h.gene_id.nunique()==25 and h.donor_label.nunique()==5
    assert m.original_human_id.nunique()==15 and m.run.nunique()==9
    assert (m.TPM==0).sum()==32 and (m.log2_TPM_plus1==0).sum()==32
    assert np.allclose(m.log2_TPM_plus1,np.log2(m.TPM+1),atol=2e-15,rtol=0)
    status={'date':'2026-10-09','passed':all(r['exit_code']==0 for r in checks),
            'commands':checks,'rendered':rendered,'input_hash_checks':hashes,
            'numeric_source_export_checks':'125 PHH and 135 mouse entries preserved; 32 mouse zeros; exported log2 arithmetic passes; original frozen order verified in plotting script.',
            'auditor_environment':'Existing separate Python/PyMuPDF is used only for PDF geometry inspection; all scientific figures are drawn in the project environment.',
            'source_warnings_reviewed':'Named width variable is checked against the actual 180 mm PDF and raster geometry. Nonnegative inputs and additive +1 guards are explicit before the log transforms. If dropna is flagged, it is used only in an assertion, never to remove any plotted data.',
            'render_iteration_note':'Initial S4 x-axis label overlapped the explanatory footer; retained first collision report, increased bottom margin and re-exported. Final figures must have no collision FAIL.',
            'visual_review':'visual_review.md'}
    (HERE/'validation.json').write_text(json.dumps(status,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'passed':status['passed'],'audits':len(checks),'source_hashes':len(hashes),'entries':[len(h),len(m)]}))
    if not status['passed']:
        raise SystemExit(1)

if __name__=='__main__':
    main()
