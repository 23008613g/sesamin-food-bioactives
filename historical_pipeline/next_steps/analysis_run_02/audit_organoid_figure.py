"""Check the actual final organoid figure without changing source data."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
QA = Path("<LOCAL_PATH>")
PDF_PYTHON = "<LOCAL_PATH>"
name = "organoid_reference"
pdf = HERE / (name + ".pdf")
checks = []
for label, command, destination in [
    ("source", [sys.executable, "-X", "utf8", str(QA / "validate_figure.py"), str(HERE / "plot_organoid.py"), "--json"], "organoid_figure_source_audit.json"),
    ("glyphs", [sys.executable, "-X", "utf8", str(QA / "audit_pdf_text.py"), str(pdf), "--min-pt", "5", "--json"], name + ".glyphs.json"),
    ("collision", [PDF_PYTHON, "-X", "utf8", str(QA / "audit_figure_collisions.py"), str(pdf), "--json-out", str(HERE / (name + ".collision.json"))], None),
]:
    result = subprocess.run(command, capture_output=True, text=True, encoding="utf-8")
    if destination:
        (HERE / destination).write_text(result.stdout, encoding="utf-8")
    checks.append({"check": label, "returncode": result.returncode, "command": command})
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
record = {"checked_at_utc": datetime.now(timezone.utc).isoformat(),
          "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(), "checks": checks,
          "alignment": name + ".alignment.json",
          "source_warnings_accepted": ["300 dpi PNG is a report preview; vector PDF/SVG supplied", "No TIFF production requirement for this pilot"],
          "visual_review": "Root viewed the final PNG: labels/numbers readable, no clipping, four libraries shown, low nonzero CPM remains visibly nonzero using two significant digits, no significance marks. Raw subtype naming is explained in caption.",
          "pdf_runtime": "Existing user-site PyMuPDF 1.27.2.2, read-only geometric QA; analysis/rendering use project venv"}
(HERE / "organoid_figure_qa.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
print(json.dumps(record))
