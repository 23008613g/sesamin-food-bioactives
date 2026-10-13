"""Freeze whole-genome Ensembl 116 human-rat orthology and rat Entrez maps.

Uses only the Python standard library. Does not inspect expression data or
candidate genes. Raw downloads are cached with SHA256; no release fallback.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import platform
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "sources" / "orthology"
BASE = "https://jun2026.archive.ensembl.org/biomart/martservice"
RELEASE = 116
HUMAN_ATTRS = ["ensembl_gene_id", "external_gene_name",
    "rnorvegicus_homolog_ensembl_gene", "rnorvegicus_homolog_associated_gene_name",
    "rnorvegicus_homolog_orthology_type", "rnorvegicus_homolog_orthology_confidence",
    "rnorvegicus_homolog_perc_id", "rnorvegicus_homolog_perc_id_r1"]
RAT_ATTRS = ["ensembl_gene_id", "external_gene_name", "entrezgene_id"]
HUMAN_COLS = ["human_ensembl_gene_id", "human_symbol", "rat_ensembl_gene_id",
    "rat_symbol", "orthology_type", "orthology_confidence",
    "rat_perc_identity", "human_perc_identity"]
RAT_COLS = ["rat_ensembl_gene_id", "rat_symbol", "rat_entrez_id"]


def utc():
    return datetime.now(timezone.utc).isoformat()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def fetch(name, url):
    dest = OUT / name
    meta_path = dest.with_name(dest.name + ".meta.json")
    if dest.exists() and meta_path.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        data = dest.read_bytes()
        if meta.get("url") == url and meta.get("sha256") == digest(data):
            print("Verified cached", name, len(data), flush=True)
            return data
        raise RuntimeError(f"Existing cache mismatch: {name}; no overwrite performed")
    started = utc()
    attempts = []
    for attempt in range(1, 4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PFHN-reproducible-analysis/39"})
            with urllib.request.urlopen(req, timeout=240) as resp:
                data = resp.read()
                meta = {"url": url, "final_url": resp.url, "http_status": resp.status,
                        "headers": dict(resp.headers), "started_utc": started,
                        "completed_utc": utc(), "release": RELEASE,
                        "sha256": digest(data), "bytes": len(data), "attempts": attempts}
            if b"Query ERROR" in data[:1000]:
                raise RuntimeError(data[:2000].decode("utf-8", errors="replace"))
            dest.write_bytes(data)
            dump(meta_path, meta)
            print("Downloaded", name, len(data), flush=True)
            return data
        except Exception as exc:
            attempts.append({"attempt": attempt, "utc": utc(), "exception": repr(exc)})
            dump(OUT / (name + ".attempts.json"), {"url": url, "release": RELEASE, "attempts": attempts})
            if attempt == 3:
                raise
            time.sleep(2 * attempt)


def query(dataset, attrs, name):
    xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
           '<Query virtualSchemaName="default" formatter="TSV" header="1" '
           'uniqueRows="1" count="" datasetConfigVersion="0.6" completionStamp="1">'
           f'<Dataset name="{dataset}" interface="default">' +
           ''.join(f'<Attribute name="{escape(a)}"/>' for a in attrs) +
           '</Dataset></Query>')
    path = OUT / (name + ".query.xml")
    if path.exists() and path.read_text(encoding="utf-8") != xml:
        raise RuntimeError("Query XML conflict: " + str(path))
    path.write_text(xml, encoding="utf-8")
    data = fetch(name + ".raw.tsv", BASE + "?" + urllib.parse.urlencode({"query": xml}))
    text = data.decode("utf-8")
    lines = text.splitlines()
    if not lines or lines[-1] != "[success]":
        raise RuntimeError(f"Missing BioMart completion stamp: {name}")
    rows = list(csv.reader(lines[:-1], delimiter="\t"))
    if not rows or any(len(row) != len(attrs) for row in rows):
        raise RuntimeError(f"Unexpected column count: {name}")
    return rows[0], rows[1:]


def write_tsv(name, cols, rows):
    path = OUT / name
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(cols)
        writer.writerows(rows)
    return {"file": name, "rows": len(rows), "sha256": digest(path.read_bytes()), "columns": cols}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    started = utc()
    metadata_jobs = [
        ("registry.xml", BASE + "?type=registry"),
        ("datasets.tsv", BASE + "?type=datasets&mart=ENSEMBL_MART_ENSEMBL"),
        ("human_attributes.tsv", BASE + "?type=attributes&dataset=hsapiens_gene_ensembl"),
        ("rat_attributes.tsv", BASE + "?type=attributes&dataset=rnorvegicus_gene_ensembl"),
        ("archive_release.html", "https://jun2026.archive.ensembl.org/info/website/archives/index.html"),
    ]
    with ThreadPoolExecutor(max_workers=2) as pool:
        raw_meta = list(pool.map(lambda args: fetch(*args), metadata_jobs))
    registry = raw_meta[0].decode("utf-8")
    if 'ensembl_mart_116' not in registry and 'Ensembl Genes 116' not in registry:
        raise RuntimeError("Registry does not confirm release 116")
    attr_headers = {}
    for label, required, raw in [("human", HUMAN_ATTRS, raw_meta[2]), ("rat", RAT_ATTRS, raw_meta[3])]:
        table = {row[0]: row for row in csv.reader(io.StringIO(raw.decode("utf-8")), delimiter="\t") if row}
        absent = sorted(set(required) - set(table))
        if absent:
            raise RuntimeError(f"Missing {label} attributes {absent}")
        attr_headers[label] = {key: table[key] for key in required}
    dump(OUT / "selected_attribute_definitions.json", attr_headers)
    with ThreadPoolExecutor(max_workers=2) as pool:
        hfuture = pool.submit(query, "hsapiens_gene_ensembl", HUMAN_ATTRS, "human_rat_all_ensembl116")
        rfuture = pool.submit(query, "rnorvegicus_gene_ensembl", RAT_ATTRS, "rat_entrez_all_ensembl116")
        hheader, hrows = hfuture.result()
        rheader, rrows = rfuture.result()
    if len(hrows) < 15000 or len(rrows) < 15000:
        raise RuntimeError("Unexpectedly short full-genome response")
    # Every row is retained. Symbols are labels and never joining/filtering keys.
    hraw_count, rraw_count = len(hrows), len(rrows)
    hrows = sorted(set(tuple(row) for row in hrows))
    rrows = sorted(set(tuple(row) for row in rrows))
    outputs = [write_tsv("human_rat_all_ensembl116.tsv", HUMAN_COLS, hrows),
               write_tsv("rat_entrez_all_ensembl116.tsv", RAT_COLS, rrows)]
    human_to_rat, rat_to_human = defaultdict(set), defaultdict(set)
    for row in hrows:
        if row[2]:
            human_to_rat[row[0]].add(row[2])
            rat_to_human[row[2]].add(row[0])
    rat_to_entrez, entrez_to_rat = defaultdict(set), defaultdict(set)
    for rat, symbol, entrez in rrows:
        if entrez:
            rat_to_entrez[rat].add(entrez)
            entrez_to_rat[entrez].add(rat)
    rat_unique = []
    rat_audit = []
    for row in rrows:
        rat, symbol, entrez = row
        reasons = []
        if not entrez:
            reasons.append("no_entrez")
        if len(rat_to_entrez[rat]) > 1:
            reasons.append("rat_ensembl_multiple_entrez")
        if entrez and len(entrez_to_rat[entrez]) > 1:
            reasons.append("entrez_multiple_rat_ensembl")
        status = "excluded:" + ";".join(reasons) if reasons else "included_unique_bidirectional"
        rat_audit.append(row + (len(rat_to_entrez[rat]), len(entrez_to_rat[entrez]) if entrez else 0, status))
        if not reasons:
            rat_unique.append(row)
    # Collapse label differences only if stable-ID keys are uniquely identified;
    # conflicting labels remain visible in raw/full tables and audit tables.
    unique_rat_by_id = {}
    for rat, symbol, entrez in rat_unique:
        unique_rat_by_id.setdefault(rat, [rat, set(), entrez])[1].add(symbol)
    rat_unique = [(rat, "|".join(sorted(symbols)), entrez) for rat, symbols, entrez in unique_rat_by_id.values()]
    rat_unique.sort()
    outputs.append(write_tsv("rat_entrez_unique_bidirectional.tsv", RAT_COLS, rat_unique))
    outputs.append(write_tsv("rat_entrez_mapping_audit.tsv", RAT_COLS + ["n_entrez_for_rat", "n_rat_for_entrez", "mapping_status"], rat_audit))
    orth_audit, orth_one2one, orth_highconf, analysis_map = [], [], [], []
    unique_entrez = {row[0]: row[2] for row in rat_unique}
    for row in hrows:
        human, hsym, rat, rsym, otype, confidence, rid, hid = row
        reasons = []
        if not rat:
            reasons.append("no_rat_ortholog")
        elif otype != "ortholog_one2one":
            reasons.append("not_one2one_type")
        if rat and len(human_to_rat[human]) != 1:
            reasons.append("human_multiple_rat_ids")
        if rat and len(rat_to_human[rat]) != 1:
            reasons.append("rat_multiple_human_ids")
        status = "excluded:" + ";".join(reasons) if reasons else "included_one2one"
        orth_audit.append(row + (len(human_to_rat[human]), len(rat_to_human[rat]) if rat else 0, status,
                                "unique_bidirectional" if rat in unique_entrez else "no_unique_entrez"))
        if not reasons:
            orth_one2one.append(row)
            if confidence == "1":
                orth_highconf.append(row)
                if rat in unique_entrez:
                    analysis_map.append(row + (unique_entrez[rat],))
    # Per-gene output must not contain conflicting duplicate records.
    for name, rows in [("one2one", orth_one2one), ("highconfidence", orth_highconf), ("analysis", analysis_map)]:
        if len({r[0] for r in rows}) != len(rows) or len({r[2] for r in rows}) != len(rows):
            raise RuntimeError(f"Ambiguous duplicate records remain in {name}; raw/audit files retained")
    if len({r[-1] for r in analysis_map}) != len(analysis_map):
        raise RuntimeError("Nonunique rat Entrez ID in strict analysis mapping")
    outputs.append(write_tsv("human_rat_one2one.tsv", HUMAN_COLS, orth_one2one))
    outputs.append(write_tsv("human_rat_one2one_highconfidence.tsv", HUMAN_COLS, orth_highconf))
    outputs.append(write_tsv("human_rat_one2one_highconfidence_unique_entrez.tsv", HUMAN_COLS + ["rat_entrez_id"], analysis_map))
    outputs.append(write_tsv("human_rat_mapping_audit.tsv", HUMAN_COLS + ["n_rat_ids_for_human", "n_human_ids_for_rat", "orthology_mapping_status", "entrez_mapping_status"], orth_audit))
    summary = {
        "release": RELEASE, "archive": "jun2026.archive.ensembl.org", "started_utc": started, "completed_utc": utc(),
        "script": str(Path(__file__).resolve()), "script_sha256": digest(Path(__file__).read_bytes()),
        "python": sys.version, "python_executable": sys.executable, "platform": platform.platform(),
        "scope": "Unfiltered entire BioMart human/rat datasets; no expression or candidate gene filtering.",
        "symbol_policy": "Symbols are display labels only; never primary joining or disambiguation keys.",
        "strict_mapping_rule": "ortholog_one2one, unique human/rat IDs across full orthology response, confidence=1, and unique rat Ensembl/Entrez mapping in BOTH directions across complete rat dataset.",
        "percent_identity_columns": {"rat_perc_identity": hheader[6], "human_perc_identity": hheader[7]},
        "source_headers": {"human_rat": hheader, "rat_entrez": rheader},
        "counts": {
            "human_rat_raw_rows": hraw_count, "rat_entrez_raw_rows": rraw_count,
            "human_rat_rows": len(hrows), "human_ids": len({r[0] for r in hrows}),
            "human_ids_with_rat_ortholog": sum(bool(v) for v in human_to_rat.values()),
            "rat_ids_with_human_ortholog": sum(bool(v) for v in rat_to_human.values()),
            "orthology_types": dict(Counter(r[4] or "no_rat_ortholog" for r in hrows)),
            "human_ids_with_multiple_rat_ids": sum(len(v) > 1 for v in human_to_rat.values()),
            "rat_ids_with_multiple_human_ids": sum(len(v) > 1 for v in rat_to_human.values()),
            "rat_entrez_rows": len(rrows), "rat_ids_total": len({r[0] for r in rrows}),
            "rat_ids_with_entrez": sum(bool(v) for v in rat_to_entrez.values()),
            "entrez_ids_total": sum(bool(v) for v in entrez_to_rat.values()),
            "rat_ids_with_multiple_entrez": sum(len(v) > 1 for v in rat_to_entrez.values()),
            "entrez_ids_with_multiple_rat_ids": sum(len(v) > 1 for v in entrez_to_rat.values()),
            "rat_entrez_unique_bidirectional": len(rat_unique), "one2one": len(orth_one2one),
            "one2one_highconfidence": len(orth_highconf), "strict_analysis_map": len(analysis_map),
            "orthology_audit_statuses": dict(Counter(r[-2] for r in orth_audit)),
            "rat_entrez_audit_statuses": dict(Counter(r[-1] for r in rat_audit))},
        "outputs": outputs,
        "official_sources": ["https://jun2026.archive.ensembl.org/info/website/archives/index.html",
                             "https://training.ensembl.org/exercise/BioMart_homologues",
                             "https://jun2026.archive.ensembl.org/info/genome/compara/homology_types.html",
                             "https://jun2026.archive.ensembl.org/info/genome/compara/Ortholog_qc_manual.html"]}
    dump(OUT / "orthology_freeze_report.json", summary)
    lines = ["# Ensembl 116 full-genome orthology freeze", "", "All source data queried from the June 2026 archive without a gene-list filter.", "",
             "Strict analysis map: `human_rat_one2one_highconfidence_unique_entrez.tsv`.",
             "Use stable Ensembl and Entrez IDs as keys; symbols are labels only.",
             "All rejected/ambiguous/unmapped rows remain in raw data and mapping audit tables.", "",
             "## Counts", ""] + [f"- {k}: {v}" for k, v in summary["counts"].items() if not isinstance(v, dict)]
    lines += ["", "## Provenance", "", "Queries, raw TSV (including completion stamp), source attribute definitions, BioMart registry, dataset definitions, HTTP metadata, retrieval UTC times, and SHA256 hashes accompany this report.",
              "No expression statistics or candidate-gene selection were performed by this script.", ""]
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(summary["counts"], ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
