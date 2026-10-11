"""Metadata and probe-ID availability audit only; no expression-effect analysis."""
from pathlib import Path
import csv, gzip, hashlib, json, platform
from collections import Counter
from datetime import datetime, timezone

D = Path(__file__).resolve().parent
P = D.parents[2]
M = P / "next_steps/evidence_closure_06/intervention_access/sources/GSE42693_series_matrix.txt.gz"
A = D / "sources/GPL1355.annot.gz"

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for b in iter(lambda: fh.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()

metadata = {}
matrix_ids = []
with gzip.open(M, "rt", encoding="utf-8") as f:
    table = False
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            table = True
            columns = next(csv.reader([next(f)], delimiter="\t"))[1:]
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if table:
            matrix_ids.append(next(csv.reader([line], delimiter="\t"))[0])
        elif line.startswith("!"):
            k, *vals = next(csv.reader([line], delimiter="\t"))
            metadata[k] = vals
with gzip.open(A, "rt", encoding="utf-8") as f:
    annotation_metadata = []
    for line in f:
        if line.startswith("!platform_table_begin"):
            break
        if line.startswith("!"):
            annotation_metadata.append(line.rstrip())
    fields = next(f).rstrip("\n").split("\t")
    rows = []
    for line in f:
        if line.startswith("!platform_table_end"):
            break
        if line.strip():
            vals = line.rstrip("\n").split("\t")
            rows.append(dict(zip(fields, vals)))
annot_ids = [r["ID"] for r in rows]
assert len(matrix_ids) == len(set(matrix_ids))
assert len(annot_ids) == len(set(annot_ids))
assert columns == metadata["!Sample_geo_accession"]
samples = [
    {"accession": x, "title": metadata["!Sample_title"][i],
     "source": metadata["!Sample_source_name_ch1"][i],
     "organism": metadata["!Sample_organism_ch1"][i],
     "included_sesamin_vs_control": i < 10,
     "group": "Control" if i < 5 else "Sesamin" if i < 10 else "Episesamin" if i < 15 else "Sesamolin"}
    for i, x in enumerate(columns)
]
assert [r["accession"] for r in samples[:10]] == [f"GSM{i}" for i in range(1048036,1048046)]
human_soft = P / "next_steps/design_closure_10/human/sources/GSE248903_family.soft"
human_ids = [x.split("=", 1)[1].strip() for x in human_soft.read_text(encoding="utf-8-sig").splitlines() if x.startswith("^SAMPLE")]
assert len(human_ids) == 9
result = {
    "generated_utc": datetime.now(timezone.utc).isoformat(),
    "execution": {"python": platform.python_version(), "executable": __import__("sys").executable, "scope": "metadata, identifiers, annotation availability; no expression values summarized or tested"},
    "sources": [
        {"path": str(x), "sha256": sha(x), "bytes": x.stat().st_size}
        for x in [M,A,human_soft]
    ],
    "matrix_probe_rows": len(matrix_ids),
    "matrix_sample_columns": len(columns),
    "matrix_unique_probe_ids": len(set(matrix_ids)),
    "annotation_rows": len(rows),
    "annotation_unique_ids": len(set(annot_ids)),
    "matrix_ids_in_annotation": len(set(matrix_ids) & set(annot_ids)),
    "matrix_ids_missing_annotation": sorted(set(matrix_ids) - set(annot_ids)),
    "annotation_ids_absent_matrix": sorted(set(annot_ids)-set(matrix_ids)),
    "annotation_metadata": annotation_metadata,
    "annotation_columns": fields,
    "annotation_symbol_nonempty_rows": sum(bool(r.get("Gene symbol","").strip()) for r in rows),
    "annotation_symbol_multi_label_rows": sum("///" in r.get("Gene symbol","") for r in rows),
    "sample_groups": dict(Counter(s["group"] for s in samples)),
    "primary_comparison_samples": [s["accession"] for s in samples if s["included_sesamin_vs_control"]],
    "deposited_processing_unique": sorted(set(metadata["!Sample_data_processing"])),
    "series_identifiers": {k: metadata.get(k,[]) for k in ["!Series_geo_accession","!Series_pubmed_id","!Series_relation","!Series_platform_id"]},
    "human_samples": human_ids,
    "human_GSE248903_GSM_intersection": sorted(set(human_ids)&set(columns)),
    "overlap_check_scope": "Accession and species metadata only. No genotype or contamination analysis.",
    "candidate_mapping_status": "Not performed here; root freezes orthologue/probe mapping independently of expression outcomes."
}
(D/"GSE42693_metadata_audit.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
with (D/"GSE42693_sample_selection.csv").open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(samples[0]));w.writeheader();w.writerows(samples)
print(json.dumps({k:result[k] for k in ["matrix_probe_rows","annotation_rows","matrix_ids_in_annotation","matrix_ids_missing_annotation","annotation_symbol_nonempty_rows","annotation_symbol_multi_label_rows","sample_groups","human_GSE248903_GSM_intersection","series_identifiers"]},ensure_ascii=True))

