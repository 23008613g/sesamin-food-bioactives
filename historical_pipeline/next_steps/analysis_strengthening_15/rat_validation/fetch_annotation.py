"""Public annotation acquisition only; no candidate expression effects inspected."""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import csv, hashlib, json, urllib.request, urllib.parse, urllib.error, sys

S=Path(__file__).resolve().parent
N=S.parents[1]
C=S/'sources';C.mkdir(exist_ok=True)
def write(p,j): p.write_text(json.dumps(j,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def fetch(key,url,limit=50_000_000):
    dest=C/key; meta=dest.with_name(dest.name+'.meta.json')
    if dest.exists() and meta.exists():
        rec=json.loads(meta.read_text(encoding='utf-8'));body=dest.read_bytes()
        assert hashlib.sha256(body).hexdigest()==rec['sha256']
        return body,rec
    t=datetime.now(timezone.utc).isoformat()
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Public-research-reanalysis/1.0','Accept':'application/json, text/plain, */*'}),timeout=40) as r:
            body=r.read(limit+1);status=r.status; final=r.url
            if len(body)>limit: raise ValueError('Response exceeds metadata cap')
    except urllib.error.HTTPError as e: body=e.read(1_000_000);status=e.code;final=e.url
    except Exception as e: body=json.dumps({'error':repr(e)}).encode();status=None;final=url
    rec={'url':url,'final_url':final,'http_status':status,'started_utc':t,'completed_utc':datetime.now(timezone.utc).isoformat(),'sha256':hashlib.sha256(body).hexdigest(),'bytes':len(body),'file':str(dest)}
    dest.write_bytes(body);write(meta,rec)
    print(key,status,len(body),flush=True)
    return body,rec

def biomart(key,dataset,filtername,ids,attrs):
    xml='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6"><Dataset name="'+dataset+'" interface="default"><Filter name="'+filtername+'" value="'+','.join(ids)+'"/>'+''.join('<Attribute name="'+a+'"/>' for a in attrs)+'</Dataset></Query>'
    (S/(key+'.query.xml')).write_text(xml,encoding='utf-8')
    return fetch(key+'_ensembl116.tsv','https://jun2026.archive.ensembl.org/biomart/martservice?query='+urllib.parse.quote(xml))

rows=list(csv.DictReader((N/'publication_package_14/data/human_candidates_all25.csv').open(encoding='utf-8-sig')))
ids=[r['gene_id'] for r in rows];assert len(ids)==25
items=[('GPL1355.annot.gz','https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1355/annot/GPL1355.annot.gz'),
       ('ensembl_release.json','https://rest.ensembl.org/info/data?content-type=application/json')]
with ThreadPoolExecutor(max_workers=2) as pool: list(pool.map(lambda x:fetch(*x),items))
attrs=['ensembl_gene_id','external_gene_name','rnorvegicus_homolog_ensembl_gene','rnorvegicus_homolog_associated_gene_name','rnorvegicus_homolog_orthology_type','rnorvegicus_homolog_orthology_confidence','rnorvegicus_homolog_perc_id','rnorvegicus_homolog_perc_id_r1']
body,rec=biomart('human_rat_orthology','hsapiens_gene_ensembl','ensembl_gene_id',ids,attrs)
assert rec['http_status']==200 and not body.startswith(b'Query ERROR'), body[:300]
parsed=list(csv.reader(body.decode().splitlines(),delimiter='\t'))
ratids=sorted({r[2] for r in parsed[1:] if len(r)>2 and r[2].startswith('ENSRNOG')})
assert ratids,'No official rat IDs'
body,rec=biomart('rat_entrez','rnorvegicus_gene_ensembl','ensembl_gene_id',ratids,['ensembl_gene_id','external_gene_name','entrezgene_id'])
assert rec['http_status']==200 and not body.startswith(b'Query ERROR'),body[:300]
write(S/'fetch_summary.json',{'status':'ANNOTATION_FETCH_COMPLETE','rat_ids':ratids,'candidate_ids':ids,'protocol_sha256':hashlib.sha256((S/'protocol.md').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.executable,'expression_effects_inspected':False})
