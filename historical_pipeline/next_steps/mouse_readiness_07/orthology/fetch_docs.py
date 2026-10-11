"""Save the official documentation consulted; no packages or keys required."""
import hashlib
import json
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / 'docs'
RAW.mkdir(exist_ok=True)
urls = {
    'homology_api': 'https://rest.ensembl.org/documentation/info/homology_species_gene_id',
    'orthology_quality': 'https://jun2026.archive.ensembl.org/info/genome/compara/Ortholog_qc_manual.html',
    'biomart_api': 'https://jun2026.archive.ensembl.org/info/data/biomart/biomart_restful.html',
    'archive_api': 'https://rest.ensembl.org/documentation/info/archive_id_get',
}
records = []
for key, url in urls.items():
    dest, meta = RAW / (key + '.html'), RAW / (key + '.meta.json')
    if dest.exists() and meta.exists():
        records.append(json.loads(meta.read_text(encoding='utf-8')))
        continue
    started = datetime.now(timezone.utc).isoformat()
    try:
        with urllib.request.urlopen(url, timeout=25) as r:
            body = r.read(1024 * 1024)
            code, final = r.status, r.url
    except urllib.error.HTTPError as e:
        body, code, final = e.read(1024 * 1024), e.code, e.url
    except Exception as e:
        body, code, final = repr(e).encode(), None, url
    record = {'url': url, 'final_url': final, 'started_utc': started, 'http_status': code,
              'file': dest.relative_to(HERE).as_posix(), 'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest()}
    dest.write_bytes(body)
    meta.write_text(json.dumps(record, indent=2), encoding='utf-8')
    records.append(record)
    print(key, code, len(body), flush=True)
(HERE / 'documentation_records.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
