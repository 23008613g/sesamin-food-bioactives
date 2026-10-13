"""Run finite Figure 1 export and source QA with the existing document runtime."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import fitz
from PIL import Image

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
p=argparse.ArgumentParser()
p.add_argument("--audit-scripts",type=Path,default=HERE/"qa_tools")
p.add_argument("--output-dir",type=Path,default=ROOT/"figures"/"main")
a=p.parse_args()
out=a.output_dir.resolve()
reports={}
for name,script,extra in [
    ("source", "validate_figure.py", [str(HERE/"draw_figure1.py"),"--json"]),
    ("pdf_text", "audit_pdf_text.py", [str(out/"Fig1.pdf"),"--min-pt","8","--json"]),
    ("collision", "audit_figure_collisions.py", [str(out/"Fig1.pdf"),"--json",
       "--json-out",str(HERE/"Fig1.collision-audit.json"),"--overlay-pdf",str(HERE/"Fig1.collision-audit.pdf")]),
]:
    result=subprocess.run([sys.executable,str(a.audit_scripts/script)]+extra,capture_output=True,text=True,encoding="utf-8")
    (HERE/f"Fig1.{name}-audit.json").write_text(result.stdout,encoding="utf-8")
    assert result.returncode==0,(name,result.stdout,result.stderr)
    reports[name]=json.loads(result.stdout)

pdf=fitz.open(out/"Fig1.pdf")
text=" ".join(pdf[0].get_text().split())
required=["Human HepG2 cells","GSE248903","9 libraries","25 human candidates",
          "Timing and independence unresolved", "Dietary rat liver", "GSE42693", "5 control + 5 sesamin",
          "Healthy animals; 15-day exposure", "Original body-weight selection retained",
          "Mouse source-label comparison","CRA013797","ctrl / HFD / sess: exploration only",
          "Exposure mapping remains unconfirmed", "Primary human hepatocytes", "GSE247407",
          "30 libraries; 5 labels", "Fatty-acid expression background", "Donor-label relationships unresolved",
          "Human liver organoids", "GSE207889", "24,034 retained cells", "AREG subtypes and expression summaries",
          "One differentiation; 4 technical libraries", "Adult human liver atlas", "GSE115469",
          "8,444 cells; 5 donors", "AREG donor and immune-cell context", "Author-normalized source expression",
          "no sesamin intervention", "Library and cell counts do not establish independent biological replication."]
for item in required:
    assert item in text,item
assert len(pdf)==1
assert abs(pdf[0].rect.width-493.2)<.02
assert abs(pdf[0].rect.height-338.4)<.02
svg=(out/"Fig1.svg").read_text(encoding="utf-8")
assert "<text" in svg
rasters={}
for ext in ["png","tiff"]:
    with Image.open(out/f"Fig1.{ext}") as im:
        assert im.size==(4110,2820),im.size
        dpi=im.info["dpi"]
        assert all(abs(d-600)<.1 for d in dpi),dpi
        rasters[ext]={"size_px":im.size,"dpi":[float(d) for d in dpi],"mode":im.mode}

report={"status":"PASS","source_preflight":reports["source"]["summary"],
        "source_warning_resolution":"174 mm continues the manuscript double-column figure width; common 89/183 mm defaults are not the target requirement.",
        "pdf_glyph_minimum_pt":reports["pdf_text"]["minimum_found_pt"],
        "collision_verdict":reports["collision"]["verdict"],
        "required_label_count":len(required),"all_labels_present":True,
        "editable_svg_text":True,"pdf_text_extractable":True,
        "rasters":rasters,"pdf_dimensions_pt":[pdf[0].rect.width,pdf[0].rect.height],
        "sha256":{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(out.glob("Fig1.*"))}}
(HERE/"figure1_export_checks.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print(json.dumps(report,indent=2))
