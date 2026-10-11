"""Rendered figure QA and export provenance; no reanalysis or model execution."""
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
STEMS = ["fig7_sc_umap", "fig8_sc_features", "fig9_sc_dotplot", "figS5_sc_library_umap", "figS6_sc_quality", "figS7_sc_composition"]


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
        "text": run_json([sys.executable, SKILLS / "audit_pdf_text.py", HERE / f"{stem}.pdf", "--min-pt", "7", "--json"], f"{stem}.text_audit.json"),
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
source_checks = []
for record in data["inputs"]:
    actual = hashlib.sha256(Path(record["source"]).read_bytes()).hexdigest()
    source_checks.append({"name": record["name"], "unchanged": actual == record["sha256"]})
assert all(r["unchanged"] for r in source_checks)
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
print(json.dumps({"figures": {stem:{"text":summary["figures"][stem]["text_returncode"],"collision":summary["figures"][stem]["collision_returncode"],"alignment":summary["figures"][stem]["alignment"]} for stem in STEMS}, "source_hashes_passed": True}))
for stem in STEMS:
    result = results[stem]
    assert result["text"]["returncode"] == result["collision"]["returncode"] == 0
    assert abs(result["geometry"]["parsed"]["width_mm"] - 180) < .001
    assert result["raster"]["tiff"]["dpi"] == [600., 600.]
    assert result["raster"]["tiff"]["compression"] == "tiff_lzw"
    assert result["alignment"]["verdict"] in {"PASS", "NOT APPLICABLE"}
