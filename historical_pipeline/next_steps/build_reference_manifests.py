"""Build explicit donor and count-column mappings from cached audit records."""
import csv
import json
from pathlib import Path
import re

BASE = Path(__file__).resolve().parent


def write(name, rows):
    with (BASE / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    human = json.loads((BASE.parent / "single_cell_datasets.json").read_text(encoding="utf-8"))["candidates"][0]
    assert human["accession"] == "GSE289173"
    hr = [{"accession": "GSE289173", **r, "unit": "human_donor", "role": "disease_reference_only", "split": "exploration", "expression_matrix_read": False, "cell_type_annotation_verified": False} for r in human["donor_library_mapping"]]
    assert len(hr) == len({r["donor_id"] for r in hr}) == 48
    write("human_snRNA_sample_manifest.csv", hr)
    metadata = json.loads((BASE / "additional_metadata.json").read_text(encoding="utf-8"))
    source = next(r for r in metadata["records"] if r["accession"] == "GSE325422")
    observed = json.loads((BASE / "additional_input_audit.json").read_text(encoding="utf-8"))["butyrate"]["count_columns"]
    prefixes = {"Control": "con", "palmitate": "IR", "palmitate + low sodium butyrate": "LNaB", "palmitate + medium sodium butyrate": "MNaB", "palmitate + high sodium butyrate": "HNaB"}
    mapped = {}
    for i, gsm in enumerate(source["sample_accessions"]):
        title = source["titles"][0][i]
        group = source["sample_characteristics"][-2][i].split(": ", 1)[1]
        batch = source["sample_characteristics"][-1][i].split(": ", 1)[1]
        rep = re.search(r"replicate ([1-3])$", title).group(1)
        col = "Count_ " + prefixes[group] + "-" + rep
        assert col not in mapped
        mapped[col] = {"accession": "GSE325422", "matrix_column_exact": col, "gsm": gsm, "title": title, "group": group, "batch": batch, "replicate_label": rep, "model": "HepG2", "role": "candidate_pending_metadata_resolution", "split": "exploration"}
    assert set(mapped) == set(observed)
    write("butyrate_sample_manifest.csv", [mapped[col] for col in observed])
    print(json.dumps({"human_donors_mapped": len(hr), "butyrate_columns_mapped": len(mapped), "butyrate_column_order_taken_from_file_not_GEO": True}))


if __name__ == "__main__":
    main()
