"""Bounded retrieval of official public annotation and paper, without expression analysis."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,urllib.request,urllib.error

HERE=Path(__file__).resolve().parent
OUT=HERE/'sources'
OUT.mkdir(exist_ok=True)
TARGETS=[
 ('GPL1355.annot.gz','https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1355/annot/GPL1355.annot.gz'),
 ('Ide2009_jnsv55_31.pdf','https://www.jstage.jst.go.jp/article/jnsv/55/1/55_1_31/_pdf')]
rows=[]
for name,url in TARGETS:
 p=OUT/name
 row={'name':name,'url':url,'retrieved_utc':datetime.now(timezone.utc).isoformat()}
 try:
  if p.exists():
   row['retrieval']='reused_current_stage_cache'
  else:
   req=urllib.request.Request(url,headers={'User-Agent':'ResearchDataAudit/1.0'})
   with urllib.request.urlopen(req,timeout=30) as response:
    row.update(status=response.status,final_url=response.url,content_type=response.headers.get('Content-Type'))
    data=response.read()
   if name.endswith('.gz'): assert data[:2]==b'\x1f\x8b','not gzip'
   if name.endswith('.pdf'): assert data.startswith(b'%PDF'),'not PDF'
   p.write_bytes(data)
  row.update(ok=True,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),local_file=str(p))
 except Exception as error:
  row.update(ok=False,error=str(error))
 rows.append(row)
(HERE/'fetch_log.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(rows,ensure_ascii=True))
