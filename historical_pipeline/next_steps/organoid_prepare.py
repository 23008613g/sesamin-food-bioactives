"""Fetch small public metadata for GSE207889; Python standard library only.

Files are capped at 10 MiB each, cached without overwriting, and hashed.
No expression matrices are downloaded or analysed by this script.
"""
from __future__ import annotations

import argparse
import ast
import csv
import gzip
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET
import zipfile

BASE = Path(__file__).resolve().parent
INPUTS = BASE / "organoid_inputs"
LIMIT = 10 * 1024 * 1024
SAMPLES = ["GSM6322577", "GSM6322578", "GSM6322579", "GSM6322580"]


def fetch(url: str, name: str, records: list) -> Path | None:
    path = INPUTS / name
    if path.parent != INPUTS or path.name != name:
        raise ValueError("Input filenames must remain directly within organoid_inputs")
    entry = {"url": url, "local_file": str(path.relative_to(BASE))}
    try:
        if path.exists():
            blob = path.read_bytes()
            if len(blob) > LIMIT:
                raise ValueError("Cached file exceeds per-file limit")
            entry["retrieval"] = "existing_file_read_without_overwrite"
        else:
            request = Request(url, headers={"User-Agent": "Research metadata audit; Python urllib"})
            with urlopen(request, timeout=40) as response:
                size = response.headers.get("Content-Length")
                if size and int(size) > LIMIT:
                    raise ValueError("Content-Length exceeds 10 MiB limit")
                blob = response.read(LIMIT + 1)
                if len(blob) > LIMIT:
                    raise ValueError("Response exceeds 10 MiB limit")
                entry["http_status"] = response.status
                entry["content_type"] = response.headers.get("Content-Type")
                entry["final_url"] = response.url
            path.write_bytes(blob)
            entry["retrieval"] = "downloaded"
        entry.update(status="read", bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest())
        records.append(entry)
        return path
    except Exception as exc:
        entry.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        records.append(entry)
        return None


def xlsx_rows(path: Path) -> dict:
    ns = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    result = {}
    with zipfile.ZipFile(path) as archive:
        if sum(info.file_size for info in archive.infolist()) > 50 * 1024 * 1024:
            raise ValueError("XLSX uncompressed content exceeds metadata cap")
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            strings = ["".join(node.itertext()) for node in root.findall("s:si", ns)]
        for name in archive.namelist():
            if not (name.startswith("xl/worksheets/sheet") and name.endswith(".xml")):
                continue
            root = ET.fromstring(archive.read(name))
            rows = []
            for row in root.findall("s:sheetData/s:row", ns):
                cells = {}
                for cell in row.findall("s:c", ns):
                    value = cell.find("s:v", ns)
                    text = value.text if value is not None else "".join(cell.itertext())
                    if cell.attrib.get("t") == "s" and value is not None:
                        text = strings[int(value.text)]
                    cells[cell.attrib["r"]] = text or ""
                if cells:
                    rows.append(cells)
            result[name] = rows
    return result


def parse_soft(path: Path) -> dict:
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("!") and " = " in line:
            key, value = line.split(" = ", 1)
            out.setdefault(key.removeprefix("!"), []).append(value)
    return out


def assignment_literal(path: Path, name: str):
    """Read author mapping as literals, never execute downloaded source."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise KeyError(name)


def export_manifest(soft: dict, records: list) -> list:
    rows = []
    for gsm in SAMPLES:
        meta = soft[gsm]
        characteristics = dict(x.split(": ", 1) for x in meta["Sample_characteristics_ch1"])
        urls = [u.replace("ftp://", "https://", 1) for k, values in meta.items() if k.startswith("Sample_supplementary_file_") for u in values]
        matrix = next(u for u in urls if u.endswith("matrix.mtx.gz"))
        sample = re.search(r"_(SM-[A-Z0-9]+)_matrix", matrix).group(1)
        counts = {}
        for suffix in ["barcodes.tsv.gz", "features.tsv.gz"]:
            url = next(u for u in urls if u.endswith(suffix))
            path = fetch(url, url.rsplit("/", 1)[1], records)
            if path is None:
                raise RuntimeError(f"Required small index file unavailable: {gsm} {suffix}")
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                lines = [line.rstrip("\n") for line in handle]
            counts[suffix] = {"rows": len(lines), "duplicates": len(lines)-len(set(lines))}
        mapping = assignment_literal(INPUTS / "github_src_sc_constants.py", "sample2condition")
        expected = "CTRL" if characteristics["treatment"] == "Control-PA" else "PA"
        if mapping[sample] != expected:
            raise ValueError(f"GEO/author condition mismatch: {gsm}")
        rows.append({
            "geo_series": "GSE207889", "gsm": gsm, "author_sample_id": sample,
            "author_short_id_alias": sample.removeprefix("SM-KBV"),
            "geo_title": meta["Sample_title"][0], "group": "CTRL-PA" if expected == "CTRL" else "PA",
            "replicate_label": int(re.search(r"replicate (\d+)", meta["Sample_title"][0]).group(1)),
            "cell_line": characteristics["cell line"], "culture": "orbital shaker HLO",
            "harvest_day": 25, "pa_uM": 0 if expected == "CTRL" else 500,
            "duration_days_geo": 5, "duration_days_corrected_paper": 4,
            "duration_resolution": "unresolved_source_conflict; corrected_paper_day21_to25",
            "differentiation_batch_id": "not_provided_in_accessed_metadata",
            "replicate_level": "technical_culture_well_library; independent_differentiation_not_established",
            "matching_control": "CTRL-PA group; same-replicate-number pairing not established",
            "barcodes_in_public_index": counts["barcodes.tsv.gz"]["rows"],
            "duplicate_barcodes_within_library": counts["barcodes.tsv.gz"]["duplicates"],
            "feature_rows_in_public_index": counts["features.tsv.gz"]["rows"],
            "matrix_url": matrix,
            "barcode_url": next(u for u in urls if u.endswith("barcodes.tsv.gz")),
            "features_url": next(u for u in urls if u.endswith("features.tsv.gz")),
            "annotation_object": "GSE207889_os.h5ad.gz; listed_not_downloaded",
            "role": "disease_stimulus_reference_only_not_ingredient_validation",
        })
    with (BASE / "organoid_sample_manifest.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def readiness(rows: list, records: list, workbook_rows: dict) -> dict:
    source_root = "https://raw.githubusercontent.com/anjahess/hlo-scrnaseq/"
    commit = json.loads((INPUTS / "github_commit.json").read_text(encoding="utf-8"))["sha"]
    notebook = json.loads((INPUTS / "github_hlo_masld.ipynb").read_text(encoding="utf-8"))
    os_cell = notebook["cells"][4]
    saved_output = "".join("".join(x.get("text", [])) for x in os_cell.get("outputs", []) if "AnnData object" in "".join(x.get("text", [])))
    return {
        "date_taipei": "2026-10-07", "accession": "GSE207889",
        "status": "manifest_and_small_indexes_ready; descriptive_reference_only; large_matrix_connection_pending",
        "sample_count": len(rows), "groups": {"CTRL-PA": 2, "PA": 2},
        "known_hpsc_lines": ["W01, WiCell"], "independent_donors": "one cell-line identity; not four independent donors",
        "differentiation_replication": {
            "finding": "Original Discussion calls the two scRNA-seq replicates technical; Statistical analysis defines replicate suspensions as different wells from the same differentiation experiment.",
            "contradictory_wording": "Library-preparation paragraph uses biological replicates; do not let this generic wording override explicit well-level definition.",
            "not_established": ["Unique differentiation batch ID of each library", "Whether control and PA groups share the same differentiation batch", "Independent differentiation replicates per condition", "Replicate-number pairing across groups"],
            "inference_allowed": "Per-library descriptive aggregation and within-study consistency only; no biological n=2/group or donor-generalizable P-values.",
            "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10711666/"
        },
        "duration": {"geo_days": 5, "geo_start_day": 20, "paper_and_2024_correction_days": 4, "paper_start_day": 21, "harvest_day_agrees": 25, "resolution": "Retain both sources. Use corrected-paper protocol description with explicit discrepancy; do not silently claim GEO was corrected."},
        "control": {"group": "matched PA carrier control", "corrected_preparation": "2024 correction describes PA stock in 100% ethanol at 500 mM, then conjugation with 10% BSA in HCM before final dilution; control prepared without active PA.", "final_ethanol_bsa_concentrations": "Not inferred from stock recipe.", "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11649770/"},
        "correction": {"doi": "10.1038/s44318-024-00308-w", "changes": "Figure 4B labels and omitted PA/OA preparation details; corrected protocol retains four days.", "impact_on_reuse": "Use corrected cell labels/protocol; verify on-object labels against corrected figure before biological interpretation."},
        "connection": {
            "manifest": "organoid_sample_manifest.csv", "sample_join": "GSM to SM-KBV62/63/64/65 is directly established by GEO file names and agrees with author constants; aliases 62/63/64/65 occur in author mappings.",
            "cell_join_key": "(library/sample identity, 10x barcode); never barcode alone across libraries",
            "small_indexes_downloaded": True, "barcode_counts_are_cells_not_biological_replicates": True,
            "large_matrices_downloaded": False, "expression_values_analyzed": False,
            "cell_annotations": "Author notebook statically documents sample, detailed_condition, celltype_category and celltype_scTypedb fields in processed os.h5ad. No exported per-cell annotation table was found among accessed small files; live H5AD contents not yet verified.",
            "raw_counts": "Individual MTX/feature/barcode files are listed by GEO; matrices not read. Processed object has only 5000 main variables in author's saved output. Inspect raw/layers or use original MTX before aggregation; scaled/Harmony values must not be summed as counts.",
            "objects": {
                "processed_os": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207889/suppl/GSE207889_os.h5ad.gz",
                "unprocessed_merged": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207889/suppl/GSE207889_masld-hlos_merged.h5ad.gz"
            },
            "object_scope_warning": "Downloaded GEO workbook maps merged object to 17 GSMs, while author's saved notebook output for merged object also shows Day34 cells. All-object contents need actual inspection; select the four manifest libraries explicitly.",
            "next_analysis": ["Read only required samples/cell types from processed annotation object after authorized large-data acquisition", "Verify sample IDs and barcode suffix conventions", "Check integer counts and feature identity before summing by library x hepatocyte-like subtype", "Report descriptive log-expression/direction consistency without biological replicate significance", "Keep disease-reference support separate from ingredient-effect replication"]
        },
        "author_code": {"commit": commit, "static_source": source_root+commit+"/src/sc_constants.py", "notebook_source": source_root+commit+"/hlo_masld.ipynb", "executed": False, "saved_notebook_output_not_our_result": saved_output},
        "spreadsheet_role": "File-to-GSM mapping and deposited checksums only; not differentiation-batch metadata.",
        "allowed_claim": "Candidate molecular response has within-study organoid PA-reference relevance, subject to cell annotation and count QC.",
        "disallowed_claims": ["Food ingredient protects organoids", "Independent replication in four donors", "Two independent differentiation experiments per arm", "Biological validation established from cell-level P-values"],
        "download_log": "organoid_metadata_audit.json", "files_recorded": len(records),
        "no_dependency_installation": True, "no_model_calls": True,
        "runtime_scope": "Existing standard-library tool environment; not a validated single-cell analysis environment"
    }


def main() -> None:
    INPUTS.mkdir(parents=True, exist_ok=True)
    records = []
    workbook = fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207889/suppl/GSE207889_GEO_supplementary_files.xlsx", "GSE207889_GEO_supplementary_files.xlsx", records)
    soft = {}
    for gsm in SAMPLES:
        path = fetch(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=full", f"{gsm}.soft.txt", records)
        if path:
            soft[gsm] = parse_soft(path)
    paper = fetch("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10711666/fullTextXML", "PMC10711666.xml", records)
    fetch("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11649770/fullTextXML", "PMC11649770_correction.xml", records)
    fetch("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207889&targ=self&form=text&view=full", "GSE207889.soft.txt", records)
    repo = fetch("https://api.github.com/repos/anjahess/hlo-scrnaseq", "github_repository.json", records)
    if repo:
        repo_data = json.loads(repo.read_text(encoding="utf-8"))
        branch = repo_data.get("default_branch", "main")
        fetch(f"https://api.github.com/repos/anjahess/hlo-scrnaseq/git/trees/{branch}?recursive=1", "github_tree.json", records)
        commit = fetch(f"https://api.github.com/repos/anjahess/hlo-scrnaseq/commits/{branch}", "github_commit.json", records)
        if commit:
            sha = json.loads(commit.read_text(encoding="utf-8"))["sha"]
            for remote in ["README.md", "src/sc_constants.py", "src/data_parsing.py", "hlo_masld.ipynb"]:
                fetch(f"https://raw.githubusercontent.com/anjahess/hlo-scrnaseq/{sha}/{remote}", "github_"+remote.replace("/", "_"), records)
    rows = export_manifest(soft, records)
    workbook_data = xlsx_rows(workbook) if workbook else {}
    report = {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "interpreter": sys.executable,
        "python_version": sys.version.split()[0],
        "downloads": records,
        "sample_metadata": soft,
        "spreadsheet_rows": workbook_data,
        "scope": "Public metadata only; no expression matrices or analysis",
    }
    (BASE / "organoid_metadata_audit.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    (BASE / "organoid_readiness.json").write_text(json.dumps(readiness(rows, records, workbook_data), ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"downloads": len(records), "failed": [x["url"] for x in records if x["status"] == "failed"], "samples": len(soft)}, ensure_ascii=True))


if __name__ == "__main__":
    main()
