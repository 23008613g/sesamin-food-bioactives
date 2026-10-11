"""Audit final figure exports and their data provenance; no scientific analysis."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from PIL import Image

HERE = Path(__file__).resolve().parent
SKILLS = Path.home() / ".agents/skills/nature-figure/scripts"
PDF_QA_PYTHON = Path.home() / "AppData/Local/Programs/Python/Python313/python.exe"
STEMS = ["fig2_human_all25", "fig6_organoid_subtypes"]


def run_json(command, filename):
    environment = dict(os.environ)
    environment["PYTHONIOENCODING"] = "utf-8"
    result = subprocess.run([str(x) for x in command], capture_output=True, encoding="utf-8", errors="strict", env=environment)
    payload = {"command": [str(x) for x in command], "returncode": result.returncode, "stdout": result.stdout, "stderr": result.stderr}
    try:
        payload["parsed"] = json.loads(result.stdout)
    except json.JSONDecodeError:
        pass
    (HERE / filename).write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return payload


results = {"source": run_json([sys.executable, SKILLS / "validate_figure.py", HERE / "build_figures.py", "--json"], "source_audit.json")}
for stem in STEMS:
    results[stem] = {
        "text": run_json([sys.executable, SKILLS / "audit_pdf_text.py", HERE / f"{stem}.pdf", "--min-pt", "5", "--json"], f"{stem}.text_audit.json"),
        "collision": run_json([PDF_QA_PYTHON, SKILLS / "audit_figure_collisions.py", HERE / f"{stem}.pdf", "--json-out", HERE / f"{stem}.collision.json", "--overlay-pdf", HERE / f"{stem}.collision_overlay.pdf"], f"{stem}.collision_execution.json"),
    }
    code = "import fitz,json,sys; d=fitz.open(sys.argv[1]); p=d[0]; print(json.dumps({'width_mm':p.rect.width/72*25.4,'height_mm':p.rect.height/72*25.4,'text_characters':len(p.get_text()),'pages':len(d),'pymupdf':fitz.VersionBind}))"
    results[stem]["geometry"] = run_json([PDF_QA_PYTHON, "-c", code, HERE / f"{stem}.pdf"], f"{stem}.pdf_geometry.json")
    raster = {}
    for extension in ["png", "tiff"]:
        with Image.open(HERE / f"{stem}.{extension}") as image:
            raster[extension] = {"pixels": list(image.size), "dpi": [float(x) for x in image.info.get("dpi", [])], "compression": image.info.get("compression")}
    results[stem]["raster"] = raster
    results[stem]["alignment"] = json.loads((HERE / f"{stem}.alignment.json").read_text(encoding="utf-8"))

data = json.loads((HERE / "source_data.json").read_text(encoding="utf-8"))
expected = {
    "human": "2d1ad658c782ee0558bbf4663f0cec1187604b10a4e7dbe8ce57f42391893796",
    "organoid_subtypes": "3356cbd6ebaae029c19f06fdd27a80137f19e1e24c8eda2260856b8af91b5a8e",
    "organoid_umi": "d33e5bb8429d87e4d30e9722938536b1706d583642ac06290fe6f215eaf3d15f",
    "organoid_library_metadata": "bf16aea91afd550c857014389bf25351cd9a1816e1f3ca254fff6c6a085f340f",
}
source_checks = []
for record in data["inputs"]:
    actual = hashlib.sha256(Path(record["path"]).read_bytes()).hexdigest()
    source_checks.append({"key": record["key"], "unchanged": actual == record["sha256"],
                          "historical_analysis_hash_matches": actual == expected[record["key"]] if record["key"] in expected else None})
assert all(r["unchanged"] and r["historical_analysis_hash_matches"] is not False for r in source_checks)
summary = {"date": "2026-10-09", "figures": {}, "source_checks": source_checks,
           "plotting_environment": sys.executable, "pdf_QA_environment": str(PDF_QA_PYTHON),
           "dependencies_installed": False, "results": results}
for stem in STEMS:
    result = results[stem]
    summary["figures"][stem] = {"text_returncode": result["text"]["returncode"], "collision_returncode": result["collision"]["returncode"],
                                "alignment": result["alignment"]["verdict"], "geometry": result["geometry"].get("parsed"), "raster": result["raster"]}
(HERE / "qa_results.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
manifest = [{"file": p.name, "bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in sorted(HERE.iterdir()) if p.is_file() and p.name != "artifact_manifest.json"]
(HERE / "artifact_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"figures": summary["figures"], "source_hashes_passed": True}))
for stem in STEMS:
    result = results[stem]
    assert result["text"]["returncode"] == result["collision"]["returncode"] == 0
    assert abs(result["geometry"]["parsed"]["width_mm"] - 180) < .001
    assert result["raster"]["tiff"]["dpi"] == [600., 600.]
    assert result["raster"]["tiff"]["compression"] == "tiff_lzw"
    assert result["alignment"]["verdict"] in {"PASS", "NOT APPLICABLE"}
