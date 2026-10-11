"""Source, rendered PDF, alignment and export gates for stage 15 figures."""
from pathlib import Path
import hashlib,json,os,subprocess,sys
from PIL import Image
HERE=Path(__file__).resolve().parent
SKILLS=Path.home()/'.agents/skills/nature-figure/scripts'
PDF_PY=Path.home()/'AppData/Local/Programs/Python/Python313/python.exe'
STEMS=['fig15_rat_external','fig15_robustness','fig15_cell_context','fig15_global_pca']

def run(cmd,filename):
    env=dict(os.environ);env['PYTHONIOENCODING']='utf-8'
    result=subprocess.run([str(x) for x in cmd],capture_output=True,encoding='utf-8',env=env)
    report={'command':[str(x) for x in cmd],'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr}
    try:report['parsed']=json.loads(result.stdout)
    except json.JSONDecodeError:pass
    (HERE/filename).write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return report

results={'source':run([sys.executable,SKILLS/'validate_figure.py',HERE/'build_figures.py','--json'],'source_audit.json'),'figures':{}}
for stem in STEMS:
    item={}
    item['text']=run([sys.executable,SKILLS/'audit_pdf_text.py',HERE/(stem+'.pdf'),'--min-pt','7','--json'],stem+'.text_audit.json')
    item['collision']=run([PDF_PY,SKILLS/'audit_figure_collisions.py',HERE/(stem+'.pdf'),'--json-out',HERE/(stem+'.collision.json'),'--overlay-pdf',HERE/(stem+'.collision_overlay.pdf')],stem+'.collision_execution.json')
    item['alignment']=json.loads((HERE/(stem+'.alignment.json')).read_text(encoding='utf-8'))
    code="import fitz,json,sys; d=fitz.open(sys.argv[1]); p=d[0]; print(json.dumps({'width_mm':p.rect.width/72*25.4,'height_mm':p.rect.height/72*25.4,'text_characters':len(p.get_text()),'pages':len(d),'pymupdf':fitz.VersionBind}))"
    item['geometry']=run([PDF_PY,'-c',code,HERE/(stem+'.pdf')],stem+'.pdf_geometry.json')
    item['raster']={}
    for ext in ['png','tiff']:
        with Image.open(HERE/(stem+'.'+ext)) as im:
            item['raster'][ext]={'pixels':list(im.size),'dpi':[float(x) for x in im.info.get('dpi',[])],'compression':im.info.get('compression')}
    results['figures'][stem]=item
sources=json.loads((HERE/'source_data.json').read_text(encoding='utf-8'))
results['sources_unchanged']=all(hashlib.sha256(Path(x['source']).read_bytes()).hexdigest()==x['sha256'] for x in sources['inputs'])
results['plot_environment']=sys.executable;results['existing_PDF_QA_environment']=str(PDF_PY)
(HERE/'qa_results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
brief={stem:{'text':r['text']['returncode'],'collision':r['collision']['returncode'],'alignment':r['alignment']['verdict']} for stem,r in results['figures'].items()}
print(json.dumps(brief))
assert results['sources_unchanged']
assert results['source']['returncode']==0,'Source audit failed'
for stem,r in results['figures'].items():
    assert r['text']['returncode']==r['collision']['returncode']==0,stem
    assert r['alignment']['verdict']=='PASS',stem
    assert abs(r['geometry']['parsed']['width_mm']-180)<.001,stem
    assert r['raster']['tiff']['dpi']==[600.,600.] and r['raster']['tiff']['compression']=='tiff_lzw',stem
