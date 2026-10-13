"""Bounded retrieval of public, primary metadata; no expression data altered."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib, json, urllib.request

ROOT = Path(__file__).resolve().parent
SOURCES = ROOT / 'sources'
SOURCES.mkdir(parents=True, exist_ok=True)
targets = []
for a in ['GSE248903', 'GSE42693', 'GSE247407', 'GSE207889', 'GSE115469', *[f'GSM{i}' for i in range(7922593,7922602)], *[f'GSM{i}' for i in range(6322577,6322581)]]:
    targets.append((a, f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={a}&targ=self&form=text&view=full', f'{a}.soft.txt'))
targets += [
 ('PRJCA021826','https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA021826','PRJCA021826.html'),
 ('CRA013797','https://ngdc.cncb.ac.cn/gsa/browse/CRA013797','CRA013797.html'),
 ('DRP018055','https://www.ebi.ac.uk/ena/browser/api/xml/DRP018055','DRP018055.xml'),
 ('PRJDB40955_DDBJ','https://ddbj.nig.ac.jp/resource/bioproject/PRJDB40955','PRJDB40955_DDBJ.html'),
 ('PRJDB40955_runs','https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB40955&result=read_run&fields=run_accession,experiment_accession,sample_accession,secondary_sample_accession,sample_title,experiment_title,library_name,library_source,library_strategy&format=json','PRJDB40955_runs.json'),
 ('rat_paper_PubMed','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=19352061&retmode=xml','rat_paper_PubMed.xml'),
 ('Hess_correction','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11649770/fullTextXML','PMC11649770.xml'),
 ('Hess_original','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10711666/fullTextXML','PMC10711666.xml'),
]
for i in range(3214402,3214411):
    targets.append((f'SAMC{i}',f'https://ngdc.cncb.ac.cn/biosample/browse/SAMC{i}',f'SAMC{i}.html'))
for i in range(1872804,1872813):
    a=f'SAMD0{i}'
    targets.append((a,f'https://www.ebi.ac.uk/ena/browser/api/xml/{a}',f'{a}.xml'))

def fetch(target):
    ident,url,filename=target
    record=dict(id=ident,url=url,request_started_utc=datetime.now(timezone.utc).isoformat())
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 (bounded public research metadata review)'})
        with urllib.request.urlopen(req,timeout=25) as response:
            body=response.read()
            record.update(http_status=response.status,final_url=response.url,bytes=len(body))
        path=SOURCES/filename
        path.write_bytes(body)
        txt=body.decode('utf-8',errors='replace')
        challenge=('GOAPPLE' in txt or 'security verification' in txt.lower() or 'captcha' in txt.lower() or '正在进行安全检查' in txt) and len(body)<20000
        record.update(status='retrieved_security_challenge' if challenge else 'retrieved',path=path.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(body).hexdigest())
    except Exception as e:
        record.update(status='failed',error=f'{type(e).__name__}: {e}')
    record['completed_utc']=datetime.now(timezone.utc).isoformat()
    print(json.dumps(record,ensure_ascii=True),flush=True)
    return record

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=3) as pool:
        result=list(pool.map(fetch,targets))
    (ROOT/'source_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
