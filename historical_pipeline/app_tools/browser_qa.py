"""Functional/visual capture with existing Playwright and a fresh isolated Chrome.

No personal profile, cookies, credentials or external analytics are used.
Reports failures honestly; screenshots still require a human/agent visual review.
"""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse
from datetime import datetime, timezone
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base-url', default='http://127.0.0.1:8765/')
    parser.add_argument('--label', default='local')
    args = parser.parse_args()
    out = ROOT / 'audit' / ('browser_' + args.label)
    out.mkdir(exist_ok=True)
    checks, errors, console_errors, requests, failed_http = [], [], [], [], []
    def check(name, value, detail=None):
        checks.append({'name': name, 'passed': bool(value), 'detail': detail})
    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=r'<LOCAL_PATH>', headless=True)
        context = browser.new_context(viewport={'width': 1440, 'height': 1100}, device_scale_factor=1)
        page = context.new_page()
        page.on('pageerror', lambda e: errors.append(str(e)))
        page.on('console', lambda e: console_errors.append(e.text) if e.type == 'error' else None)
        hosting_asset_redirects=[]
        logo=ROOT/'space/assets/infinia-logo.png'
        logo_hash=hashlib.sha256(logo.read_bytes()).hexdigest() if logo.exists() else None
        def record_request(r):
            first=r
            while first.redirected_from:first=first.redirected_from
            target=urlparse(r.url)
            # The Space redirects its local PNG path to HF-managed content storage.
            # Permit only this exact image's redirect chain and content hash, not a CDN wildcard.
            if (first.url==args.base_url.rstrip('/')+'/assets/infinia-logo.png'
                and r.resource_type=='image' and target.scheme=='https'
                and target.hostname=='cdn-lfs-us-1.hf.co' and target.path.endswith('/'+str(logo_hash))):
                hosting_asset_redirects.append(target._replace(query='',fragment='').geturl())
            else:requests.append(r.url)
        page.on('request', record_request)
        page.on('response', lambda r: failed_http.append({'url':r.url,'status':r.status}) if r.status >= 400 else None)
        try:
            response = page.goto(args.base_url, wait_until='networkidle', timeout=60000)
            check('initial HTTP 200', response.status == 200)
            page.locator('#app-content').wait_for(state='visible', timeout=20000)
            check('app loaded', page.locator('#load-error').is_hidden())
            opts = page.locator('#gene-select option').evaluate_all('(xs)=>xs.map(x=>({value:x.value,label:x.textContent}))')
            check('all 25 candidates selectable', len(opts) == 25, len(opts))
            def select_gene(symbol):
                option = next(x for x in opts if x['label'] == symbol or x['label'].split(' ')[0] == symbol)
                page.locator('#gene-select').select_option(option['value'])
                page.wait_for_timeout(180)
            select_gene('VNN1')
            check('initial VNN1 narrative', 'VNN1' in page.locator('#gene-insight-title').inner_text())
            stats = page.locator('#rat-stats').inner_text()
            check('both inference methods visible', 'Welch' in stats and ('exact' in stats.lower() or '精確' in stats), stats)
            check('exact q not hidden', '0.198' in stats or '0.20' in stats, stats)
            check('all 9 mouse libraries plotted', page.locator('#mouse-chart circle').count() == 9)
            check('all 10 rat animals plotted', page.locator('#rat-chart circle').count() == 10)
            page.screenshot(path=str(out / '00_first_screen.png'))
            page.screenshot(path=str(out / '01_evidence_zh.png'), full_page=True)
            before = page.locator('#rat-stats').inner_text()
            page.locator('#normalization').select_option('quantile')
            page.wait_for_timeout(150)
            after = page.locator('#rat-stats').inner_text()
            check('rat normalization changes statistics', before != after, {'before':before,'after':after})
            page.locator('#normalization').select_option('deposited')
            page.locator('[data-lang="en"]').click()
            check('English language applies', page.locator('html').get_attribute('lang').startswith('en'))
            check('English content translated', 'mouse' in page.locator('#study-cards').inner_text().lower())
            page.screenshot(path=str(out / '02_evidence_en.png'), full_page=True)
            page.locator('#candidate-search').fill('AREG')
            rows = page.locator('#candidate-table tbody tr')
            check('table search filters to AREG', rows.count() == 1 and 'AREG' in rows.first.inner_text())
            rows.first.click()
            check('table click selects gene', 'AREG' in page.locator('#gene-insight-title').inner_text())
            page.locator('#candidate-search').fill('')
            check('all candidates retained after clear', page.locator('#candidate-table tbody tr').count() == 25)
            select_gene('JAG2')
            jag_text = page.locator('#view-evidence').inner_text()
            check('low abundance caveat visible for JAG2', any(x in jag_text.lower() for x in ['low abundance', 'low-abundance', 'below 1', '< 1', '<1']))
            select_gene('CXCL8')
            cx = page.locator('#study-cards').inner_text()
            check('missing mouse contrast is NA', 'NA' in cx)
            select_gene('NEGR1')
            zero_class = page.locator('#study-cards .mouse .effect-value').get_attribute('class')
            check('true zero uses neutral color', 'neutral' in zero_class or 'zero' in zero_class, zero_class)
            select_gene('AREG')
            page.locator('[data-view="atlas"]').click()
            check('view navigation returns to page top', page.evaluate('window.scrollY') == 0)
            page.wait_for_function("document.querySelector('#atlas-loading').classList.contains('hidden')", timeout=30000)
            page.wait_for_timeout(250)
            check('all 24034 cells displayed', '24,034' in page.locator('#cell-count').inner_text())
            canvas = page.locator('#umap-canvas')
            check('canvas initialized', canvas.evaluate('(c)=>c.width>500&&c.height>200'))
            check('canvas has plotted pixels', canvas.evaluate("(c)=>{const d=c.getContext('2d').getImageData(0,0,c.width,c.height).data;let n=0;for(let i=0;i<d.length;i+=4)if(d[i+3]>0 && Math.abs(d[i]-d[i+1])>12)n++;return n>1000;}"))
            page.screenshot(path=str(out / '03_atlas_areg.png'), full_page=True)
            h1 = hashlib.sha256(canvas.screenshot()).hexdigest()
            page.locator('#zoom-in').click()
            h2 = hashlib.sha256(canvas.screenshot()).hexdigest()
            check('zoom changes rendering', h1 != h2)
            page.locator('#reset-zoom').click()
            h3 = hashlib.sha256(canvas.screenshot()).hexdigest()
            check('reset restores rendering', h1 == h3)
            libs = page.locator('#library-filter option').evaluate_all('(xs)=>xs.map(x=>x.value)')
            check('four library filters', len(libs) == 5)
            page.locator('#library-filter').select_option(libs[1])
            page.wait_for_function("!document.querySelector('#cell-count').textContent.includes('24,034')", timeout=5000)
            one_count = page.locator('#cell-count').inner_text()
            check('library filter changes count', '24,034' not in one_count, one_count)
            page.locator('#library-filter').select_option('all')
            page.wait_for_function("document.querySelector('#cell-count').textContent.includes('24,034')", timeout=5000)
            page.locator('#color-by').select_option('subtype')
            check('categorical legend visible', len(page.locator('#atlas-legend').inner_text()) > 30)
            page.screenshot(path=str(out / '04_atlas_subtype.png'), full_page=True)
            page.locator('#color-by').select_option('expression')
            # Hover an actual saved coordinate; use the same documented axes/fit,
            # not hidden browser state, to locate a real cell on the visible canvas.
            raw = context.request.get(args.base_url.rstrip('/') + '/data/cells.json').json()
            box = canvas.bounding_box()
            xs, ys = raw['coordinates']['x'], raw['coordinates']['y']
            scale = min((box['width'] - 82) / (max(xs) - min(xs)), (box['height'] - 76) / (max(ys) - min(ys)))
            px = box['x'] + box['width']/2 + (xs[0] - (max(xs)+min(xs))/2)*scale
            py = box['y'] + box['height']/2 - (ys[0] - (max(ys)+min(ys))/2)*scale
            page.mouse.move(px, py)
            page.locator('#cell-tooltip').wait_for(state='visible', timeout=5000)
            check('real-cell tooltip contains expression and library', 'AREG' in page.locator('#cell-tooltip').inner_text() and 'log1p' in page.locator('#cell-tooltip').inner_text())
            canvas.focus()
            pan_before = hashlib.sha256(canvas.screenshot()).hexdigest()
            canvas.press('ArrowRight')
            pan_after = hashlib.sha256(canvas.screenshot()).hexdigest()
            check('keyboard pan changes plot', pan_before != pan_after)
            page.locator('#reset-zoom').click()
            select_gene('AOC4P')
            na_text = page.locator('#atlas-legend').inner_text()
            check('AOC4P marked unavailable, not zero', any(x in na_text.lower() for x in ['not measured', 'not covered', 'not present', 'missing', 'unavailable', 'not in', '未測', '未收錄', '缺失', 'na']), na_text)
            page.screenshot(path=str(out / '05_atlas_missing.png'), full_page=True)
            select_gene('AREG')
            page.locator('#weighting').select_option('pooled')
            pooled = page.locator('#composition-chart').inner_text()
            page.locator('#weighting').select_option('equal_library')
            equal = page.locator('#composition-chart').inner_text()
            check('weighting changes decomposition', pooled != equal)
            with page.expect_download(timeout=10000) as event:
                page.locator('#download-atlas').click()
            download = event.value
            download.save_as(str(out / 'exported_atlas.png'))
            check('PNG export succeeds', (out / 'exported_atlas.png').stat().st_size > 5000)
            page.locator('[data-view="methods"]').click()
            links = page.locator('#download-list a').evaluate_all('(xs)=>xs.map(x=>({href:x.href,label:x.textContent}))')
            check('12 source downloads present', len(links) == 12, len(links))
            for link in links:
                result = context.request.get(link['href'])
                check('download ' + link['href'].split('/')[-1], result.status == 200 and len(result.body()) > 50)
            check('dataset source links present', page.locator('#source-cards a').count() >= 4)
            page.screenshot(path=str(out / '06_methods.png'), full_page=True)
            page.locator('[data-lang="zh"]').click()
            page.locator('[data-view="evidence"]').click()
            select_gene('VNN1')
            for width in [390, 768, 1440]:
                page.set_viewport_size({'width':width,'height':900})
                page.wait_for_timeout(150)
                dims = page.evaluate('({scroll:document.documentElement.scrollWidth, client:document.documentElement.clientWidth})')
                check(f'no page overflow at {width}px', dims['scroll'] <= dims['client'] + 1, dims)
                if width == 390:
                    page.screenshot(path=str(out / '07_mobile_evidence.png'), full_page=True)
                    page.locator('[data-view="atlas"]').click()
                    page.wait_for_timeout(200)
                    page.screenshot(path=str(out / '08_mobile_atlas.png'), full_page=True)
                    dims = page.evaluate('({scroll:document.documentElement.scrollWidth, client:document.documentElement.clientWidth})')
                    check('mobile atlas no page overflow', dims['scroll'] <= dims['client'] + 1, dims)
                    page.locator('[data-view="evidence"]').click()
            check('gene included in share URL', 'VNN1' in page.url, page.url)
            page.reload(wait_until='networkidle')
            page.locator('#app-content').wait_for(state='visible')
            check('URL state survives refresh', 'VNN1' in page.locator('#gene-insight-title').inner_text())
        except Exception as exc:
            check('browser workflow completed', False, str(exc))
            page.screenshot(path=str(out / 'failure.png'), full_page=True)
        check('no JavaScript exceptions', not errors, errors)
        check('no failed app HTTP requests', not failed_http, failed_http)
        check('no console errors', not console_errors, console_errors)
        origin = urlparse(args.base_url).netloc
        external = sorted({urlparse(u)._replace(query='',fragment='').geturl() for u in requests if urlparse(u).scheme in {'http','https'} and urlparse(u).netloc != origin})
        # HF may inject its own standard iframe host script; record all such requests explicitly.
        check('no app-owned third-party requests', not external, external)
        report = {'status':'PASS' if all(c['passed'] for c in checks) else 'FAIL',
                  'base_url':args.base_url, 'timestamp':datetime.now(timezone.utc).isoformat(),
                  'browser':browser.version, 'isolated_context':True, 'checks':checks,
                  'hosting_asset_redirects':sorted(set(hosting_asset_redirects)),
                  'visual_review':'Screenshots captured; root must inspect before release.'}
        (out / 'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
        print(json.dumps({'status':report['status'], 'passed':sum(c['passed'] for c in checks), 'total':len(checks),
                          'failures':[c for c in checks if not c['passed']]},ensure_ascii=False))
        browser.close()
    return report['status'] == 'PASS'

if __name__ == '__main__':
    raise SystemExit(0 if run() else 1)
