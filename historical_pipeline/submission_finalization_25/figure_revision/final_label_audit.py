"""Read-only final label, numeric-text, export-size and source audits."""
from pathlib import Path
import hashlib,json,re,sys
import xml.etree.ElementTree as ET
import fitz
from PIL import Image

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
report={"audit_interpreter":sys.executable,"fitz_version":fitz.VersionBind,"figures":{},"legacy_svg_labels":[]}
for name in ('Figure_S01','Figure_S02','Figure_S13'):
    old=HERE/'before'/f'{name}.pdf'
    new=ROOT/'figures/supplement'/f'{name}.pdf'
    with fitz.open(old) as a,fitz.open(new) as b:
        before='\n'.join(page.get_text() for page in a)
        after='\n'.join(page.get_text() for page in b)
        assert re.findall(r'-?\d+(?:\.\d+)?',before)==re.findall(r'-?\d+(?:\.\d+)?',after)
        expected=before.replace('OA + Sesamin / Control','Treat / Control').replace('OA + Sesamin / OA','Treat / Model').replace('OA / Control','Model / Control')
        expected=expected.replace('OA + sesamin','Treat')
        expected=re.sub(r'^OA$','Model',expected,flags=re.M)
        assert expected==after,(name,'unexpected PDF text edit')
        width,height=b[0].rect.width/72*25.4,b[0].rect.height/72*25.4
        assert abs(width-180)<0.001
    raster={}
    for ext in ('png','tiff'):
        with Image.open(ROOT/'figures/supplement'/f'{name}.{ext}') as im:
            raster[ext]={"pixels":list(im.size),"dpi":[float(x) for x in im.info['dpi']]}
            assert min(im.info['dpi'])>=599
    stage24=ROOT.parent/'manuscript_revision_24/figures/supplement'
    original_unchanged=all(sha(HERE/'before'/f'{name}.{ext}')==sha(stage24/f'{name}.{ext}') for ext in ('pdf','svg','png','tiff'))
    assert original_unchanged
    report['figures'][name]={"pdf_width_mm":width,"pdf_height_mm":height,"raster":raster,
        "only_authorized_pdf_text_changed":True,"all_rendered_numeric_tokens_unchanged":True,
        "before_copies_equal_read_only_stage24":True}
for path in (ROOT/'figures').rglob('*.svg'):
    for node in ET.parse(path).getroot().iter():
        if node.tag.rsplit('}',1)[-1] in ('text','title','desc'):
            value=''.join(node.itertext())
            if re.search(r'OA\s*(?:\+|plus)\s*Sesamin|OA / Control|共同处理|油酸.{0,5}芝麻素',value,re.I):
                report['legacy_svg_labels'].append({'file':str(path.relative_to(ROOT)),'text':value})
assert not report['legacy_svg_labels']
(HERE/'final_label_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
