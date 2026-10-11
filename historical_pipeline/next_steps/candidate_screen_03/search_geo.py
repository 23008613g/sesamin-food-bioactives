"""Bounded official GEO DataSets retrieval, public metadata only, stdlib only."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request

HERE = Path(__file__).resolve().parent
CACHE = HERE / "geo_api"
CACHE.mkdir(exist_ok=True)
QUERIES = [
    '(palmitate OR palmitic OR lipotoxicity) AND (hepatocyte OR HepG2 OR HepaRG OR liver) AND gse[ETYP]',
    '("oleic acid" OR "free fatty acids") AND (HepG2 OR HepaRG OR hepatocytes) AND gse[ETYP]',
    '(ergothioneine OR urolithin OR spermidine OR astaxanthin) AND (liver OR hepatic OR hepatocyte OR HepG2) AND gse[ETYP]',
]
calls = []


def get(endpoint, params, name):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/' + endpoint + '?' + urllib.parse.urlencode(params)
    path = CACHE / name
    started = datetime.now(timezone.utc).isoformat()
    if path.exists():
        body = path.read_bytes()
        status = 'cached_read_only'
    else:
        request = urllib.request.Request(url, headers={'User-Agent': 'PublicResearchMetadataAudit/1.0'})
        for attempt in range(2):
            try:
                with urllib.request.urlopen(request, timeout=35) as response:
                    body = response.read(10_000_001)
                if len(body) > 10_000_000:
                    raise ValueError('Metadata response over 10 MB budget')
                json.loads(body)
                path.write_bytes(body)
                status = 'downloaded'
                break
            except Exception as exc:
                calls.append({'url': url, 'attempt': attempt+1, 'error': str(exc), 'started_utc': started})
                if attempt == 1:
                    raise
                time.sleep(1)
        time.sleep(.4)
    calls.append({'url': url, 'file': str(path.relative_to(HERE)), 'bytes': len(body),
                  'sha256': hashlib.sha256(body).hexdigest(), 'status': status, 'started_utc': started})
    return json.loads(body)


def main():
    results = []
    ids = []
    try:
        for i, query in enumerate(QUERIES):
            obj = get('esearch.fcgi', {'db': 'gds', 'term': query, 'retmode': 'json', 'retmax': 300}, f'search_{i+1}.json')
            record = obj['esearchresult']
            results.append({'query': query, 'total_count': int(record['count']), 'retrieved_count': len(record['idlist']),
                            'query_translation': record.get('querytranslation'), 'warnings': record.get('warninglist')})
            ids.extend(record['idlist'])
            print(json.dumps(results[-1]), flush=True)
        unique = list(dict.fromkeys(ids))
        studies = []
        for i in range(0, len(unique), 50):
            obj = get('esummary.fcgi', {'db': 'gds', 'id': ','.join(unique[i:i+50]), 'retmode': 'json'}, f'summary_{i//50+1}.json')
            result = obj['result']
            studies.extend(result[uid] for uid in result['uids'])
        (HERE / 'geo_study_summaries.json').write_text(json.dumps(studies, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
        print(json.dumps({'unique_studies': len(studies)}), flush=True)
    finally:
        (HERE / 'geo_query_log.json').write_text(json.dumps({'executed_utc': datetime.now(timezone.utc).isoformat(),
            'queries': results, 'calls': calls, 'scope': 'Bounded metadata retrieval, no expression analysis; first successful esearch is endpoint preflight'}, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
