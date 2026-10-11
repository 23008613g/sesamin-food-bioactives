"""Freeze provenance and query log for the bounded stage-15 source review."""
from pathlib import Path
import hashlib,json,sys
from datetime import datetime,timezone
D=Path(__file__).resolve().parent
P=D.parents[2]
now=datetime.now(timezone.utc).isoformat()
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
records=[]
for f in sorted((D/"sources").glob("web_*.json")):
    j=json.loads(f.read_text(encoding="utf-8"))
    records.append({"cache":str(f.relative_to(D)),"request":j.get("query"),"sha256":sha(f)})
log={"generated_utc":now,"scope":"16 initial bounded source/data queries; 8 additional VNN1 novelty queries explicitly requested by root after the rat result","search_query_count":sum(len(x["request"].get("search_query",[])) for x in records),"requests":records,"actual_direct_downloads":json.loads((D/"fetch_log.json").read_text(encoding="utf8")),"failure_boundary":["ScienceDirect original human/mouse article direct opens failed with 403 and were not retried repeatedly.","Sage article direct open timed out; no usable public matrix/accession established.","Some GEO HTML requests returned recaptcha/access errors; existing official cached records and successful official annotation FTP were used.","No access challenge bypass, author contact, installation or paid service."],"limitations":["Not a systematic search or proof of exhaustive absence.","Web-index excerpts for human/mouse paper identity do not replace complete Methods.","Non-primary hits were navigational only; cited conclusions use original sources."]}
assert log["search_query_count"]==24
(D/"search_log.json").write_text(json.dumps(log,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
e=json.loads((D/"evidence.json").read_text(encoding="utf8"))
e["generated_utc"]=now
e["additional_prior_reporting_check"]="vnn1_prior_reporting.json"
(D/"evidence.json").write_text(json.dumps(e,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
inputs=[
P/"next_steps/design_closure_10/design_comparison.json",
P/"next_steps/design_closure_10/human/evidence.json",
P/"next_steps/design_closure_10/mouse/evidence.json",
P/"next_steps/evidence_closure_06/intervention_access/findings.json",
P/"next_steps/evidence_closure_06/intervention_access/sources/GSE42693_series_matrix.txt.gz",
P/"next_steps/evidence_closure_06/intervention_access/sources/GSE42693_series.txt",
]
files=inputs+sorted(x for x in D.rglob("*") if x.is_file() and x.name!="source_manifest.json")
manifest={"generated_utc":now,"project_python_executable":sys.executable,"document_python_for_PDF_read_only":"<LOCAL_PATH>","sha_method":"streaming SHA256","files":[{"path":str(f),"bytes":f.stat().st_size,"sha256":sha(f),"role":"prior_read_only_input" if f in inputs else "stage15_source_or_report"} for f in files]}
(D/"source_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
audit=json.loads((D/"GSE42693_metadata_audit.json").read_text(encoding="utf8"))
assert audit["matrix_ids_missing_annotation"]==[]
assert audit["matrix_probe_rows"]==audit["annotation_rows"]==31099
assert audit["human_GSE248903_GSM_intersection"]==[]
for j in D.glob("*.json"):json.loads(j.read_text(encoding="utf8"))
print(json.dumps({"status":"PASS","search_queries":log["search_query_count"],"manifest_files":len(files),"evidence_sha256":sha(D/"evidence.json"),"prior_reporting_sha256":sha(D/"vnn1_prior_reporting.json"),"meaning":"File integrity and metadata checks only; no new biological validation"},ensure_ascii=True))

