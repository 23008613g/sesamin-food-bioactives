"""Audit actual final exports; keep PDF-only QA separate from numerical runtime."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qa-scripts", type=Path, required=True)
    parser.add_argument("--pdf-python", required=True)
    args = parser.parse_args()
    source = subprocess.run([sys.executable, "-X", "utf8", str(args.qa_scripts/"validate_figure.py"), str(HERE/"run_bulk_analysis.py"), "--json"], capture_output=True, text=True, encoding="utf-8")
    source_result = json.loads(source.stdout)
    (HERE/"figure_source_audit.json").write_text(json.dumps(source_result, indent=2)+"\n", encoding="utf-8")
    assert source.returncode == 0
    results = []
    for name in ["sample_quality", "literature_panel"]:
        pdf = HERE/(name+".pdf")
        glyph = subprocess.run([sys.executable, "-X", "utf8", str(args.qa_scripts/"audit_pdf_text.py"), str(pdf), "--min-pt", "5", "--json"], capture_output=True, text=True, encoding="utf-8")
        (HERE/(name+".glyphs.json")).write_text(glyph.stdout, encoding="utf-8")
        assert glyph.returncode == 0, glyph.stdout+glyph.stderr
        collision = subprocess.run([args.pdf_python, "-X", "utf8", str(args.qa_scripts/"audit_figure_collisions.py"), str(pdf), "--json-out", str(HERE/(name+".collision.json"))], capture_output=True, text=True, encoding="utf-8")
        assert collision.returncode == 0, collision.stdout+collision.stderr
        results.append({"figure": name, "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
                        "glyph_audit_passed": True, "collision_audit_passed": True,
                        "alignment_report": name+".alignment.json"})
    record = {"checked_at_utc": datetime.now(timezone.utc).isoformat(), "figures": results,
        "numerical_and_plotting_runtime": sys.executable, "pdf_only_audit_runtime": args.pdf_python,
        "pdf_runtime_scope": "Existing Python with installed PyMuPDF 1.27.2.2 used only to inspect final PDF geometry; no package or system change. Numerical analysis and rendering remain entirely in project venv.",
        "source_warnings_reviewed": ["300 dpi PNG is an analysis preview; vector PDF/SVG are also provided", "No TIFF because this is a research pilot/report, not a journal production submission"],
        "visual_review": "Root inspected every panel of PNG previews: all sample points shown; panel and axis labels readable; no visible clipping; mean and sample SD stated; figures do not assert efficacy or original protein-result replication.",
        "previous_failures": ["First static source validator could not recognize loop-built export suffixes, converted millimetre expressions or a passed alignment callable; explicit equivalent calls resolved false flags", "PDF-only audit with Python -I hid existing user-site PyMuPDF; regular existing tool runtime was verified and used without installation", "Initial audit wrapper decoded Windows child output as UTF-8 without setting child UTF-8 mode; added -X utf8 to child commands and reran"]}
    (HERE/"figure_qa.json").write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(record))


if __name__ == "__main__":
    main()
