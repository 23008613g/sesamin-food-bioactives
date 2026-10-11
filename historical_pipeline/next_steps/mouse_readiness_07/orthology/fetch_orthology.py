"""Fetch a small, frozen 25-candidate orthology audit using stdlib only.

No expression is queried. Successful and failed response bodies are immutable
cache entries. Delete nothing: use --refresh-suffix for a future collection.
"""
from __future__ import annotations
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import platform
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
INPUT = PROJECT / 'next_steps/sesamin_exploration_04/priority_all.csv'
REST = 'https://rest.ensembl.org'
ap = argparse.ArgumentParser()
ap.add_argument('--refresh-suffix', default='')
args = ap.parse_args()
CACHE = HERE / ('raw' + args.refresh_suffix)
CACHE.mkdir(parents=True, exist_ok=True)
records = []

def sha(data):
    return hashlib.sha256(data).hexdigest()

def fetch(key, url, json_response=True):
    meta = CACHE / (key + '.meta.json')
    bodypath = CACHE / (key + ('.json' if json_response else '.txt'))
    if meta.exists() and bodypath.exists():
        rec = json.loads(meta.read_text(encoding='utf-8'))
        body = bodypath.read_bytes()
        assert sha(body) == rec['sha256'], key
    else:
        start = datetime.now(timezone.utc).isoformat()
        req = urllib.request.Request(url, headers={'User-Agent': 'Public-data-orthology-audit/1.0', 'Accept': 'application/json' if json_response else '*/*'})
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                body = r.read(8 * 1024 * 1024 + 1)
                if len(body) > 8 * 1024 * 1024:
                    raise ValueError('Metadata response exceeded 8 MiB cap')
                status, final_url, headers = r.status, r.url, dict(r.headers)
        except urllib.error.HTTPError as e:
            body = e.read(1024 * 1024)
            status, final_url, headers = e.code, e.url, dict(e.headers)
        except Exception as e:
            body = json.dumps({'transport_error': repr(e)}).encode()
            status, final_url, headers = None, url, {}
        rec = {'key': key, 'requested_url': url, 'final_url': final_url, 'started_utc': start,
               'completed_utc': datetime.now(timezone.utc).isoformat(), 'http_status': status,
               'bytes': len(body), 'sha256': sha(body), 'response_headers': headers,
               'file': bodypath.relative_to(HERE).as_posix()}
        bodypath.write_bytes(body)
        meta.write_text(json.dumps(rec, indent=2, ensure_ascii=False), encoding='utf-8')
        time.sleep(.15)
    records.append(rec)
    parsed = None
    if json_response:
        try: parsed = json.loads(body)
        except (ValueError, UnicodeDecodeError): pass
    else:
        parsed = body.decode('utf-8', errors='replace')
    print(key, rec['http_status'], len(body), flush=True)
    return parsed, rec

def rest(key, path):
    return fetch(key, REST + path + ('&' if '?' in path else '?') + 'content-type=application/json')

candidate_bytes = INPUT.read_bytes()
rows = list(csv.DictReader(candidate_bytes.decode('utf-8-sig').splitlines()))
assert len(rows) == 25 and len({r['gene_id'] for r in rows}) == 25
freeze = HERE / 'mapping_rule.json'
rule = {
    'created_utc': datetime.now(timezone.utc).isoformat(),
    'input': str(INPUT.relative_to(PROJECT)), 'input_sha256': sha(candidate_bytes),
    'candidate_count': 25,
    'primary_rule': 'Original human ENSG resolves as human gene; exactly one mouse target with official ortholog_one2one, matching mouse gene lookup and official high-confidence flag 1. No symbol-based replacement or functional analog substitution.',
    'unavailable_confidence_rule': 'Leave primary eligibility false/pending if orthology confidence unavailable; report identity percentages without invented confidence threshold.',
    'expression_usage': 'none; this audit does not inspect mouse expression or select favorable effects',
    'retired_ids': 'retain original ID and archive result; do not silently remap',
}
if freeze.exists():
    frozen = json.loads(freeze.read_text(encoding='utf-8'))
    assert frozen['input_sha256'] == rule['input_sha256']
else:
    freeze.write_text(json.dumps(rule, indent=2, ensure_ascii=False), encoding='utf-8')

metadata = {}
for key, path in [('info_data', '/info/data'), ('info_rest', '/info/rest'), ('info_software', '/info/software'), ('assembly_human', '/info/assembly/homo_sapiens'), ('assembly_mouse', '/info/assembly/mus_musculus')]:
    metadata[key] = rest(key, path)[0]

def fetch_candidate(row):
    gid = row['gene_id']
    lookup, lookup_rec = rest('lookup_' + gid, '/lookup/id/' + gid)
    hom, hom_rec = rest('homology_' + gid, '/homology/id/homo_sapiens/' + gid + '?target_species=mus_musculus&type=orthologues&sequence=none&format=full')
    result = {'original_human_id': gid, 'input_symbol': row['symbol'], 'human_lookup': lookup,
              'human_lookup_status': lookup_rec['http_status'], 'homology_response': hom,
              'homology_http_status': hom_rec['http_status'], 'mouse_lookups': {}}
    if not isinstance(lookup, dict) or 'id' not in lookup:
        result['archive_response'] = rest('archive_' + gid, '/archive/id/' + gid)[0]
    if isinstance(hom, dict):
        for block in hom.get('data', []):
            for h in block.get('homologies', []):
                mid = h.get('target', {}).get('id')
                if mid:
                    result['mouse_lookups'][mid] = rest('lookup_' + mid, '/lookup/id/' + mid)[0]
    return result

with ThreadPoolExecutor(max_workers=3) as pool:
    results = list(pool.map(fetch_candidate, rows))

# A small BioMart query supplies confidence/GOC/WGA fields absent from many REST responses.
attrs = ['ensembl_gene_id', 'external_gene_name', 'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name', 'mmusculus_homolog_orthology_type', 'mmusculus_homolog_orthology_confidence', 'mmusculus_homolog_perc_id', 'mmusculus_homolog_perc_id_r1', 'mmusculus_homolog_goc_score', 'mmusculus_homolog_wga_coverage']
xml = '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Filter name="ensembl_gene_id" value="' + ','.join(r['gene_id'] for r in rows) + '"/>' + ''.join('<Attribute name="' + a + '"/>' for a in attrs) + '</Dataset></Query>'
(HERE / ('biomart_query' + args.refresh_suffix + '.xml')).write_text(xml, encoding='utf-8')
mart, mart_rec = fetch('biomart_25', 'https://www.ensembl.org/biomart/martservice?query=' + urllib.parse.quote(xml), False)
registry, registry_rec = fetch('biomart_registry', 'https://www.ensembl.org/biomart/martservice?type=registry', False)
dataset, dataset_rec = fetch('biomart_datasets', 'https://www.ensembl.org/biomart/martservice?type=datasets&mart=ENSEMBL_MART_ENSEMBL', False)
out = {'metadata': metadata, 'runtime': {'python': sys.version, 'platform': platform.platform()},
       'candidate_sha256': sha(candidate_bytes), 'mapping_rule_sha256': sha(freeze.read_bytes()),
       'candidates': results, 'biomart_attributes': attrs, 'biomart_status': mart_rec['http_status']}
(HERE / ('collected' + args.refresh_suffix + '.json')).write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding='utf-8')
(HERE / ('source_records' + args.refresh_suffix + '.json')).write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps({'candidates': len(results), 'requests': len(records), 'bytes': sum(r['bytes'] for r in records)}))
