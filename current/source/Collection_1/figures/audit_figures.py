"""Inspect final PDFs with the existing document Python; no scientific analysis."""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import fitz

HERE=Path(__file__).resolve().parent
QA=HERE/'qa';SCRIPTS=HERE/'scripts'
parser=argparse.ArgumentParser();parser.add_argument('--only',choices=['all','1','2','3','4','5'],default='all');args=parser.parse_args()
source=subprocess.run([sys.executable,'-X','utf8',str(SCRIPTS/'validate_figure.py'),str(HERE/'build_figures.py'),'--json'],capture_output=True,text=True,encoding='utf-8')
(QA/'source_validation.json').write_text(source.stdout,encoding='utf-8')
record=[{'kind':'source','returncode':source.returncode}]
for n in range(1,6):
    if args.only not in ['all',str(n)]:continue
    stem=f'Fig{n}';pdf=HERE/'main'/f'{stem}.pdf'
    glyph=subprocess.run([sys.executable,'-X','utf8',str(SCRIPTS/'audit_pdf_text.py'),str(pdf),'--min-pt','8','--json'],capture_output=True,text=True,encoding='utf-8')
    (QA/f'{stem}.glyphs.json').write_text(glyph.stdout,encoding='utf-8')
    collision=subprocess.run([sys.executable,'-X','utf8',str(SCRIPTS/'audit_figure_collisions.py'),str(pdf),'--json-out',str(QA/f'{stem}.collision.json'),'--overlay-pdf',str(QA/f'{stem}.collision.pdf')],capture_output=True,text=True,encoding='utf-8')
    (QA/f'{stem}.collision.txt').write_text(collision.stdout+'\n'+collision.stderr,encoding='utf-8')
    with fitz.open(pdf) as doc:
        page=doc[0];w=page.rect.width/72*25.4;h=page.rect.height/72*25.4;text=page.get_text()
        assert abs(w-174)<.01 and h<234 and len(text)>30
        page.get_pixmap(matrix=fitz.Matrix(150/72,150/72),alpha=False).save(QA/f'{stem}.pdf-review.png')
    record.append({'figure':stem,'glyph_returncode':glyph.returncode,'collision_returncode':collision.returncode,'width_mm':w,'height_mm':h,'selectable_characters':len(text)})
    print(stem,'glyph',glyph.returncode,'collision',collision.returncode)
(QA/f'pdf_audit_summary_{args.only}.json').write_text(json.dumps({'qa_environment':sys.executable,'pymupdf':fitz.VersionBind,'results':record},indent=2)+'\n',encoding='utf-8')
raise SystemExit(1 if any(r.get('returncode',0) or r.get('glyph_returncode',0) or r.get('collision_returncode',0) for r in record) else 0)
