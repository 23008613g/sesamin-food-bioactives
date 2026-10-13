"""Retrieve and freeze the official GO basic ontology; no expression analysis."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, urllib.request

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'sources/go'
OUT.mkdir(parents=True, exist_ok=True)
URL = 'https://current.geneontology.org/ontology/go-basic.obo'
target = OUT/'go-basic.obo'
if target.exists():
    raise FileExistsError('Do not silently replace a frozen ontology')
started = datetime.now(timezone.utc).isoformat()
req = urllib.request.Request(URL, headers={'User-Agent':'PFHN-reproducible-public-data-reanalysis/1.0'})
try:
    with urllib.request.urlopen(req, timeout=45) as response:
        data = response.read()
        meta = {'requested_url':URL, 'resolved_url':response.url,
                'last_modified':response.headers.get('Last-Modified'),
                'retrieved_utc':datetime.now(timezone.utc).isoformat(),
                'started_utc':started, 'bytes':len(data),
                'sha256':hashlib.sha256(data).hexdigest()}
    assert data.startswith(b'format-version:') and b'[Term]' in data
    meta['header'] = data.decode('utf-8').split('[Term]',1)[0]
    target.write_bytes(data)
    (OUT/'source_manifest.json').write_text(json.dumps(meta,indent=2),encoding='utf8')
    print(json.dumps({k:v for k,v in meta.items() if k!='header'}))
except Exception as err:
    (OUT/'retrieval_failure.json').write_text(json.dumps({'url':URL,'started_utc':started,'error':repr(err)},indent=2),encoding='utf8')
    raise
