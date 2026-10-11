"""Audit already-rendered assets without changing scientific geometry or data."""
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
STEMS = ["fig2_expression", "fig3_concordance", "fig4_background"]


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


results = {}
results["source"] = run_json([sys.executable, SKILLS / "validate_figure.py", HERE / "build_figures.py", "--json"], "source_audit.json")
for stem in STEMS:
    results[stem] = {
        "text": run_json([sys.executable, SKILLS / "audit_pdf_text.py", HERE / f"{stem}.pdf", "--min-pt", "5", "--json"], f"{stem}.text_audit.json"),
        "collision": run_json([PDF_QA_PYTHON, SKILLS / "audit_figure_collisions.py", HERE / f"{stem}.pdf", "--json-out", HERE / f"{stem}.collision.json", "--overlay-pdf", HERE / f"{stem}.collision_overlay.pdf"], f"{stem}.collision_execution.json"),
    }
    # PDF inspection alone uses the existing independent PyMuPDF environment.
    code = "import fitz,json,sys; d=fitz.open(sys.argv[1]); p=d[0]; print(json.dumps({'width_mm':p.rect.width/72*25.4,'height_mm':p.rect.height/72*25.4,'text_characters':len(p.get_text()),'pages':len(d),'pymupdf':fitz.VersionBind}))"
    results[stem]["pdf_geometry"] = run_json([PDF_QA_PYTHON, "-c", code, HERE / f"{stem}.pdf"], f"{stem}.pdf_geometry.json")
    raster = {}
    for extension in ["png", "tiff"]:
        with Image.open(HERE / f"{stem}.{extension}") as image:
            raster[extension] = {"pixels": list(image.size), "dpi": [float(x) for x in image.info.get("dpi", [])], "compression": image.info.get("compression")}
    results[stem]["raster"] = raster
    results[stem]["alignment"] = json.loads((HERE / f"{stem}.alignment.json").read_text(encoding="utf-8"))

source_data = json.loads((HERE / "source_data.json").read_text(encoding="utf-8"))
hash_checks = []
for record in source_data["source_inputs"]:
    path = Path(record["path"])
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    hash_checks.append({"key": record["key"], "unchanged": actual == record["sha256"]})
historical_expected = {
    "human_library_counts": "2d1ad658c782ee0558bbf4663f0cec1187604b10a4e7dbe8ce57f42391893796",
    "mouse_library_tpm": "bad3fd0a665ae2690f18813976c6f19dcd5dddb2c61dbe04cfcb4caa98ad34ce",
    "phh_label_effects": "788b5cac44ca24b7e6ac3de8e8f1a340b4a8d250f65e0d27256717f10964492f",
    "organoid_sensitivity": "e7cd42b8aaf0fb1453fad7a76d78dd3a152492d1e9cb5e7b94f2f5ee072442f1",
}
historical_checks = [{"key": r["key"], "matches_historical_analysis_sha256": r["sha256"] == historical_expected[r["key"]]}
                     for r in source_data["source_inputs"] if r["key"] in historical_expected]
assert all(r["matches_historical_analysis_sha256"] for r in historical_checks)
all_outputs = [{"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size}
               for path in sorted(HERE.iterdir()) if path.is_file() and path.name not in {"qa_results.json", "artifact_manifest.json"}]
(HERE / "artifact_manifest.json").write_text(json.dumps(all_outputs, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
summary = {
    "date": "2026-10-09", "figures": {}, "source_static_returncode": results["source"]["returncode"],
    "source_input_hash_checks": hash_checks,
    "historical_analysis_hash_checks": historical_checks,
    "plotting_environment": sys.executable, "pdf_geometry_environment": str(PDF_QA_PYTHON),
    "dependencies_installed": False, "results": results,
}
for stem in STEMS:
    r = results[stem]
    summary["figures"][stem] = {"text_returncode": r["text"]["returncode"], "collision_returncode": r["collision"]["returncode"],
                                "alignment": r["alignment"]["verdict"], "geometry": r["pdf_geometry"].get("parsed"), "raster": r["raster"]}
    assert abs(r["pdf_geometry"]["parsed"]["width_mm"] - 180) < .001
    assert r["raster"]["tiff"]["dpi"] == [600.0, 600.0]
    assert r["raster"]["tiff"]["compression"] == "tiff_lzw"
    assert r["text"]["returncode"] == r["collision"]["returncode"] == 0
    assert r["alignment"]["verdict"] in {"PASS", "NOT APPLICABLE"}
(HERE / "qa_results.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"figures": summary["figures"], "all_sources_unchanged": all(r["unchanged"] for r in hash_checks)}))
