"""Audit public candidate files with stdlib; source files are never overwritten."""
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

from run_descriptive_pilot import read_table

HERE = Path(__file__).resolve().parent
DOWNLOADS = HERE / "additional_inputs"
S = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def fetch(url, name, cap):
    path = DOWNLOADS / name
    if not path.exists():
        with urllib.request.urlopen(url, timeout=90) as response:
            data = response.read(cap + 1)
        if len(data) > cap:
            raise ValueError("Download exceeded cap")
        with path.open("xb") as stream:
            stream.write(data)
    return path, {"url": url, "local_file": str(path.relative_to(HERE)), "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def xlsx_rows(path):
    with zipfile.ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            with archive.open("xl/sharedStrings.xml") as stream:
                for _, elem in ET.iterparse(stream, events=("end",)):
                    if elem.tag == S + "si":
                        strings.append("".join(t.text or "" for t in elem.iter(S + "t")))
                        elem.clear()
        sheets = ET.fromstring(archive.read("xl/workbook.xml")).find(S + "sheets")
        first = list(sheets)[0]
        rel_id = first.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
        relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target = next(e.attrib["Target"] for e in relations if e.attrib["Id"] == rel_id)
        member = target.lstrip("/") if target.startswith("/") else "xl/" + target
        with archive.open(member) as stream:
            for _, elem in ET.iterparse(stream, events=("end",)):
                if elem.tag != S + "row":
                    continue
                result = {}
                for cell in elem.findall(S + "c"):
                    letters = "".join(c for c in cell.attrib["r"] if c.isalpha())
                    index = 0
                    for c in letters:
                        index = index * 26 + ord(c) - ord("A") + 1
                    value = cell.find(S + "v")
                    raw = value.text if value is not None else ""
                    if cell.attrib.get("t") == "s":
                        raw = strings[int(raw)]
                    elif cell.attrib.get("t") == "inlineStr":
                        raw = "".join(t.text or "" for t in cell.iter(S + "t"))
                    result[index - 1] = raw
                yield result
                elem.clear()


def main():
    DOWNLOADS.mkdir(exist_ok=True)
    vit_e, re = fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195619/matrix/GSE195619_series_matrix.txt.gz", "GSE195619_series_matrix.txt.gz", 5_000_000)
    butyrate, rb = fetch("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325422/suppl/GSE325422_Allgene_info.xlsx", "GSE325422_Allgene_info.xlsx", 15_000_000)
    ah, ar = read_table(vit_e, "!series_matrix_table_begin", "!series_matrix_table_end")
    bh, br = read_table(HERE.parent / "downloads/GSE200765_series_matrix.txt.gz", "!series_matrix_table_begin", "!series_matrix_table_end")
    a = {r[0]: list(map(float, r[1:])) for r in ar}
    b = {r[0]: list(map(float, r[1:])) for r in br}
    ids = sorted(set(a) & set(b))
    pairs, exact = [], []
    for i, ga in enumerate(ah[1:]):
        for j, gb in enumerate(bh[1:]):
            eq = sum(a[k][i] == b[k][j] for k in ids)
            rmse = math.sqrt(sum((a[k][i] - b[k][j])**2 for k in ids) / len(ids))
            row = {"GSE195619_gsm": ga, "GSE200765_gsm": gb, "matching_features": len(ids), "exactly_equal_features": eq, "rmse": rmse}
            pairs.append(row)
            if eq == len(ids):
                exact.append(row)
    with (HERE / "cross_study_sample_comparison.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(pairs[0]))
        writer.writeheader()
        writer.writerows(pairs)
    reader = xlsx_rows(butyrate)
    header = next(reader)
    counts = {i: s for i, s in header.items() if s.startswith("Count_")}
    fpkm = {i: s for i, s in header.items() if s.startswith("FPKM_")}
    gene_col = next(i for i, s in header.items() if s == "gene_id")
    assert len(counts) == 15 and len(fpkm) == 15
    stat = {s: {"sum": 0, "nonzero": 0, "missing": 0, "nonnumeric": 0, "nonfinite": 0, "negative": 0, "fractional": 0} for s in counts.values()}
    ids_seen, duplicate_genes, num_rows = set(), 0, 0
    for row in reader:
        if not row:
            continue
        num_rows += 1
        gene = row.get(gene_col, "")
        duplicate_genes += gene in ids_seen
        ids_seen.add(gene)
        for i, label in counts.items():
            rec = stat[label]
            raw = row.get(i, "")
            if raw == "":
                rec["missing"] += 1
                continue
            try:
                value = float(raw)
            except ValueError:
                rec["nonnumeric"] += 1
                continue
            if not math.isfinite(value):
                rec["nonfinite"] += 1
                continue
            rec["negative"] += value < 0
            rec["fractional"] += value != int(value)
            rec["sum"] += value
            rec["nonzero"] += value != 0
    out = {"executed_at_utc": datetime.now(timezone.utc).isoformat(), "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "downloads": [re, rb],
           "vitamin_E_vs_D": {"features_compared": len(ids), "column_pairs_compared": len(pairs), "exact_duplicate_pairs": exact, "inference": "Identical deposited control profiles require these series to be grouped as sharing source data; do not count as independent validation. Nonidentical values alone do not establish independence."},
           "butyrate": {"feature_rows": num_rows, "duplicate_gene_ids": duplicate_genes, "count_columns": list(counts.values()), "other_expression_columns": list(fpkm.values()), "count_value_checks": stat, "status": "Structural count check only; not biological QC or confirmed suitable ingredient", "metadata_issues": ["GEO dose says PA 200mM; do not silently substitute units", "Series mentions 24/48h but sample records say 48h", "File has FPKM columns while protocol says TPM", "No butyrate-only group", "No verified independent matched PA+butyrate partner", "Sodium butyrate is a nutrition-associated metabolite; food-origin scope must be defined before promoting it to main ingredient"]}}
    (HERE / "additional_input_audit.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"duplicate_pairs": exact, "butyrate_rows": num_rows, "count_columns": len(counts), "count_format_issues": sum(v[k] for v in stat.values() for k in ["missing", "nonnumeric", "nonfinite", "negative", "fractional"])}))


if __name__ == "__main__":
    main()
