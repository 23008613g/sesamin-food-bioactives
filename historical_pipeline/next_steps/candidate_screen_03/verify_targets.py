"""Verify selected public GEO metadata and linked small-file availability."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
import urllib.request
import urllib.parse

HERE = Path(__file__).resolve().parent
OUT = HERE / 'target_metadata'
OUT.mkdir(exist_ok=True)
from search_geo import get, calls


def fetch(acc, target='self'):
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?' + urllib.parse.urlencode({'acc': acc, 'targ': target, 'form': 'text', 'view': 'full'})
    path = OUT / (acc + '.soft.txt')
    if not path.exists():
        with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent':'PublicResearchMetadataAudit/1.0'}), timeout=35) as response:
            data = response.read(5_000_001)
        if len(data) > 5_000_000 or not data.startswith(b'^'):
            raise ValueError('Unexpected SOFT response')
        path.write_bytes(data)
        time.sleep(.4)
    data = path.read_bytes()
    text = data.decode('utf-8')
    fields = {}
    for line in text.splitlines():
        if line.startswith('!') and ' = ' in line:
            key, value = line.split(' = ', 1)
            fields.setdefault(key, []).append(value)
    return {'accession':acc, 'url':url, 'local_file':str(path.relative_to(HERE)),
            'sha256':hashlib.sha256(data).hexdigest(), 'fields':fields}


def main():
    results=[]
    for acc in ['GSE248903','GSM7922593','GSM7922596','GSM7922599','GSE45961','GSE93819','GSE303090']:
        try:
            result=fetch(acc)
            results.append(result)
            print(json.dumps(result),flush=True)
        except Exception as exc:
            results.append({'accession':acc,'error':str(exc)})
    for label,query in [('sesamin','sesamin AND gse[ETYP]'),('linalool','linalool AND gse[ETYP]')]:
        data=get('esearch.fcgi',{'db':'gds','term':query,'retmode':'json','retmax':100},label+'_search.json')['esearchresult']
        ids=data['idlist']
        if ids:
            summary=get('esummary.fcgi',{'db':'gds','id':','.join(ids),'retmode':'json'},label+'_summaries.json')['result']
            results.append({'search':query,'count':int(data['count']), 'studies':[summary[x] for x in summary['uids']]})
    (HERE/'target_metadata.json').write_text(json.dumps({'executed_utc':datetime.now(timezone.utc).isoformat(),
        'records':results,'extra_api_calls':calls},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
