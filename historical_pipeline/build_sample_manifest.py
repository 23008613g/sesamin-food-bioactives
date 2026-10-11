"""Rebuild explicit sample mapping from cached GEO metadata and table headers.

No network, dependency installation, or biological analysis. Input files are read-only.
Run from any directory; output is this script's sibling sample_manifest.csv.
"""
import csv
import hashlib
import json
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent
FIELDS = ["accession", "matrix_column", "gsm", "title", "group",
          "donor_or_cell_line", "replicate_label", "role", "split",
          "independence_status", "local_file"]


def main():
    metadata = json.loads((BASE / "metadata_audit.json").read_text(encoding="utf-8-sig"))
    tables = json.loads((BASE / "processed_table_audit.json").read_text(encoding="utf-8-sig"))
    sources = {r["accession"]: r for r in metadata["records"]}
    output = []
    for table in tables["records"]:
        accession = table["accession"]
        source = sources[accession]
        assert len(source["sample_accessions"]) == len(source["titles"][0])
        mapped = {}
        for i, (gsm, title) in enumerate(zip(source["sample_accessions"], source["titles"][0], strict=True)):
            match = re.search(r"([1-3])$", title)
            assert match, title
            rep = match.group(1)
            if accession == "GSE200765":
                col, donor, role = gsm, "cell_line_HepaRG", "known_result_pilot"
                if "untreated" in title:
                    group = "vehicle"
                elif "pre-treated only" in title:
                    group = "FFA"
                elif "1,25" in title:
                    group = "active_vitamin_D_plus_FFA"
                elif "SM-VD2" in title:
                    group = "SM_VD2_plus_FFA"
                else:
                    assert "SV-VD3" in title, title
                    group = "SV_VD3_plus_FFA"
            elif accession == "GSE270357":
                group = source["sample_characteristics"][-1][i].split(": ", 1)[1]
                col = {"PA": "PA", "Cur + PA": "PA_Cur", "Cur-Cu + PA": "PA_Cur_CuCl2"}[group] + rep
                donor, role = "cell_line_HepG2", "limited_intervention_candidate"
            elif accession == "GSE247407":
                donor = source["sample_characteristics"][0][i].split(": ", 1)[1]
                group = source["sample_characteristics"][-1][i].split(": ", 1)[1]
                col = donor + {"Control": "C", "FFA": "F"}[group] + rep
                role = "fatty_acid_reference_only"
            else:
                raise ValueError(f"No explicit mapping for {accession}")
            assert col not in mapped, (accession, col)
            mapped[col] = dict(zip(FIELDS, [
                accession, col, gsm, title, group, donor, rep, role, "exploration",
                "biological_vs_culture_vs_technical_replication_requires_source_review",
                table["local_file"],
            ], strict=True))
        assert set(mapped) == set(table["sample_columns"]), accession
        output.extend(mapped[col] for col in table["sample_columns"])
    assert len(output) == 53
    destination = BASE / "sample_manifest.csv"
    with destination.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    print(json.dumps({"rows": len(output), "output": destination.name,
                      "sha256": hashlib.sha256(destination.read_bytes()).hexdigest()}))


if __name__ == "__main__":
    main()
