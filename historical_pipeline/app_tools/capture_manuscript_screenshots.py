"""Capture unmodified public UI for manuscript illustration, excluding footer/contact.

Existing isolated Playwright runtime; screenshots are browser captures, not AI images.
No pixel editing, injected styles, hidden caveats, data mutations or profile access.
"""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import hashlib,json,math,struct
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
deployment=json.loads((ROOT/'audit/deployment.json').read_text(encoding='utf-8'))
VERSION=deployment['app_version']
OUT=ROOT.parents[1]/'APP交付_2026-10-09'/('论文截图_v'+VERSION)
AUDIT=ROOT/'audit'/('screenshots_v'+VERSION)
OUT.mkdir(exist_ok=True)
AUDIT.mkdir(exist_ok=True)
URL='https://fc28-sesamin-atlas.static.hf.space/'
manifest=json.loads((ROOT/'audit/publication_manifest.json').read_text(encoding='utf-8'))
checks=[]; captures=[]; errors=[]
def check(name,condition,detail=None):
    checks.append({'name':name,'passed':bool(condition),'detail':detail})
    if not condition:raise RuntimeError(name)
def record(path,kind,link,clip=None):
    data=path.read_bytes(); dimensions=struct.unpack('>II',data[16:24])
    captures.append({'file':path.name,'kind':kind,'url':link,'pixel_dimensions':dimensions,
        'sha256':hashlib.sha256(data).hexdigest(),'css_clip':clip})

with sync_playwright() as pw:
    browser=pw.chromium.launch(executable_path=r'<LOCAL_PATH>',headless=True)
    ctx=browser.new_context(viewport={'width':1440,'height':1100},device_scale_factor=2)
    info=ctx.request.get('https://huggingface.co/api/spaces/fc28/sesamin-atlas').json()
    check('current reviewed version',info['sha']==deployment['commit_oid'],info['sha'])
    source_hashes={}
    for file in ['app.js','charts.js','data/evidence.json','data/cells.json']:
        response=ctx.request.get(URL+file)
        sha=hashlib.sha256(response.body()).hexdigest()
        expected=next(row['sha256'] for row in manifest['files'] if row['path']==file)
        check('public source hash '+file,response.ok and sha==expected)
        source_hashes[file]=sha
    page=ctx.new_page(); page.on('pageerror',lambda e:errors.append(str(e)))
    for stem,view,gene,end_selector in [
        ('01_Cell_atlas_interface_EN','atlas','AREG','.atlas-panel'),
        ('02_Cross_model_interface_EN','evidence','VNN1','.evidence-charts')]:
        page.set_viewport_size({'width':1440,'height':1100})
        page.goto(URL+f'#view={view}&gene={gene}&lang=en',wait_until='networkidle')
        page.locator('#app-content').wait_for(state='visible')
        if view=='atlas':
            page.wait_for_function("document.querySelector('#atlas-loading').classList.contains('hidden') && !document.querySelector('#download-atlas').disabled")
            page.locator('#reset-zoom').click()
            check('all organoid cells displayed','24,034' in page.locator('#cell-count').inner_text())
            check('no sesamin qualification visible','no sesamin treatment' in page.locator('.atlas-intro').inner_text())
            check('technical replication note present','Four technical libraries come from one differentiation' in page.locator('.atlas-footnote').inner_text())
        else:
            check('both rat q values retained',all(s in page.locator('#rat-stats').inner_text() for s in ['Welch','exact']))
            check('nine mouse libraries plotted',page.locator('#mouse-chart circle').count()==9)
        page.evaluate('window.scrollTo(0,0)')
        page.evaluate('document.fonts.ready')
        page.wait_for_timeout(200)
        bounds=page.locator(end_selector).bounding_box()
        height=math.ceil(bounds['y']+bounds['height']+8)
        page.set_viewport_size({'width':1440,'height':height})
        page.wait_for_timeout(150)
        bounds=page.locator(end_selector).bounding_box()
        height=math.ceil(bounds['y']+bounds['height']+8)
        left=math.floor(page.locator('.workspace').bounding_box()['x'])
        clip={'x':left,'y':0,'width':1440-left,'height':height}
        brand=page.locator('.sidebar .brand').bounding_box()
        check(stem+' excludes company header logo',brand['x']+brand['width']<=left)
        footer=page.locator('.page-footer').bounding_box()
        check(stem+' excludes institutional footer',footer['y']>height)
        check(stem+' full scientific panel retained',bounds['y']+bounds['height']<=height)
        check(stem+' institutional contact hidden',page.locator('.institution-attribution').is_hidden())
        check(stem+' English gene state',page.locator('html').get_attribute('lang')=='en' and page.locator('#gene-select').input_value()==gene)
        path=OUT/(stem+'.png'); page.screenshot(path=str(path),clip=clip,animations='disabled')
        record(path,'unaltered live UI capture at 2x device scale',page.url,clip)
        if view=='atlas':
            with page.expect_download() as download:page.locator('#download-atlas').click()
            path=OUT/'03_AREG_UMAP_native_export_EN.png';download.value.save_as(str(path))
            record(path,'native in-app PNG export',page.url)
    check('no JavaScript errors',not errors,errors)
    ctx.close();browser.close()

report={'status':'PASS','captured_at':datetime.now(ZoneInfo('Asia/Taipei')).isoformat(),
    'app_version':deployment['app_version'],'commit_oid':info['sha'],'app_url':URL,
    'device_scale_factor':2,'viewport_css_width':1440,'source_hashes':source_hashes,
    'captures':captures,'checks':checks,'visual_review':'Pending root inspection',
    'not_changed':['DOM','CSS','data','legends','scientific text','source documents'],
    'crop_scope':'UI screenshots show the main workspace, excluding the company sidebar. They stop after the complete principal scientific panel; subsequent panels and footer are outside frame.'}
(AUDIT/'capture_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'PASS','captures':len(captures),'checks':len(checks),'dimensions':[c['pixel_dimensions'] for c in captures]}))
