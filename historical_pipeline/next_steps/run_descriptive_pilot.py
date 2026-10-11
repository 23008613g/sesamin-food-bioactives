"""Stdlib-only descriptive pilot; does not perform hypothesis tests or AI calls."""
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics as st
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OUT = HERE / "pilot_descriptive"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_table(path, begin, end):
    with gzip.open(path, "rt", encoding="utf-8-sig") as stream:
        for line in stream:
            if line.strip() == begin:
                break
        else:
            raise ValueError(f"Missing {begin}")
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        rows = []
        for row in reader:
            if row and row[0] == end:
                return header, rows
            if len(row) != len(header):
                raise ValueError(f"Width mismatch in {path.name}")
            rows.append(row)
        raise ValueError(f"Missing {end}")


def write_csv(name, columns, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def quantile(values, fraction):
    ordered = sorted(values)
    pos = (len(ordered) - 1) * fraction
    lo = int(pos)
    hi = min(lo + 1, len(ordered) - 1)
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (pos - lo)


def main():
    OUT.mkdir(exist_ok=True)
    protocol = json.loads((HERE / "pilot_protocol.json").read_text(encoding="utf-8"))
    matrix = BASE / "downloads/GSE200765_series_matrix.txt.gz"
    expected = next(r for r in json.loads((BASE / "processed_table_audit.json").read_text(encoding="utf-8"))["records"] if r["accession"] == "GSE200765")
    assert sha(matrix) == expected["sha256"]
    header, all_rows = read_table(matrix, "!series_matrix_table_begin", "!series_matrix_table_end")
    columns = header[1:]
    with (BASE / "sample_manifest.csv").open(encoding="utf-8-sig", newline="") as stream:
        sample_records = [r for r in csv.DictReader(stream) if r["accession"] == "GSE200765"]
    by_gsm = {r["gsm"]: r for r in sample_records}
    assert set(columns) == set(by_gsm)
    groups = {}
    for i, gsm in enumerate(columns):
        groups.setdefault(by_gsm[gsm]["group"], []).append(i)
    assert sorted(map(len, groups.values())) == [2, 3, 3, 3, 3]
    selected = [r for r in all_rows if r[0].startswith("TC") and ".hg." in r[0]]
    ids = [r[0] for r in selected]
    values = [list(map(float, r[1:])) for r in selected]
    assert len(ids) == len(set(ids)) == 19911
    assert all(math.isfinite(v) for row in values for v in row)
    annotations = {}
    platform = HERE / "GPL31262_platform.soft.txt.gz"
    if platform.exists():
        ah, ar = read_table(platform, "!platform_table_begin", "!platform_table_end")
        for fields in ar:
            row = dict(zip(ah, fields, strict=True))
            symbols = set()
            for assignment in row.get("SPOT_ID", "").split(" /// "):
                parts = assignment.split(" // ")
                if len(parts) > 1 and parts[1].strip() not in ["", "---"]:
                    symbols.add(parts[1].strip())
            annotations[row["ID"]] = sorted(symbols)
    write_csv("feature_annotation.csv", ["feature_id", "gene_symbols", "mapping_status"], [
        [i, "|".join(annotations.get(i, [])), "unmapped" if not annotations.get(i) else "unique_symbol" if len(annotations[i]) == 1 else "ambiguous_symbols"] for i in ids])
    vectors = list(zip(*values))
    qc = []
    for i, gsm in enumerate(columns):
        v = vectors[i]
        qc.append([gsm, by_gsm[gsm]["group"], len(v), min(v), quantile(v, .25), st.median(v), quantile(v, .75), max(v), st.mean(v), st.stdev(v)])
    write_csv("sample_qc.csv", ["gsm", "group", "features", "min", "q25", "median", "q75", "max", "mean", "sd"], qc)
    corrs = []
    duplicates = []
    for i, a in enumerate(columns):
        for j in range(i + 1, len(columns)):
            b = columns[j]
            corr = st.correlation(vectors[i], vectors[j])
            exact = vectors[i] == vectors[j]
            if exact:
                duplicates.append([a, b])
            corrs.append([a, b, by_gsm[a]["group"], by_gsm[b]["group"], corr, exact])
    write_csv("sample_correlations.csv", ["sample_a", "sample_b", "group_a", "group_b", "pearson_r_across_features", "identical_values"], corrs)
    group_columns = list(groups)
    per_probe = []
    contrasts = protocol["contrasts"]
    panel = []
    for feature, row in zip(ids, values, strict=True):
        means = {g: st.mean([row[i] for i in ix]) for g, ix in groups.items()}
        sds = {g: st.stdev([row[i] for i in ix]) for g, ix in groups.items()}
        result = [feature, "|".join(annotations.get(feature, []))]
        for g in group_columns:
            result.extend([len(groups[g]), means[g], sds[g]])
        for treated, reference in contrasts:
            result.extend([means[treated] - means[reference], math.sqrt(sds[treated]**2/len(groups[treated]) + sds[reference]**2/len(groups[reference]))])
        per_probe.append(result)
        for symbol in annotations.get(feature, []):
            if symbol in protocol["fixed_literature_panel"]:
                for i, gsm in enumerate(columns):
                    panel.append([symbol, feature, "|".join(annotations[feature]), gsm, by_gsm[gsm]["group"], row[i]])
    names = ["feature_id", "gene_symbols"]
    for g in group_columns:
        names.extend([g + "_n", g + "_mean", g + "_sd"])
    for a, b in contrasts:
        names.extend([a + "_minus_" + b + "_signal_difference", a + "_minus_" + b + "_descriptive_SE"])
    write_csv("all_feature_descriptives.csv", names, per_probe)
    write_csv("literature_panel_samples.csv", ["gene_symbol", "feature_id", "all_mapped_symbols", "gsm", "group", "submitted_signal"], panel)
    found = sorted({r[0] for r in panel})
    summary = {
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv, "interpreter": sys.executable, "python_version": sys.version,
        "dependencies": "Python standard library only; existing skills-smoke tool environment, not a new formal analysis environment",
        "script_sha256": sha(Path(__file__)), "protocol_sha256": sha(HERE / "pilot_protocol.json"),
        "matrix_sha256": sha(matrix), "platform_sha256": sha(platform) if platform.exists() else None,
        "raw_rows": len(all_rows), "human_transcript_clusters": len(ids), "excluded_nonhuman_or_control_rows": len(all_rows)-len(ids),
        "groups": {g: [columns[i] for i in ix] for g, ix in groups.items()},
        "unique_symbol_features": sum(len(annotations.get(i, [])) == 1 for i in ids),
        "ambiguous_features": sum(len(annotations.get(i, [])) > 1 for i in ids),
        "unmapped_features": sum(not annotations.get(i) for i in ids),
        "fixed_panel_found": found, "fixed_panel_missing": sorted(set(protocol["fixed_literature_panel"])-set(found)),
        "minimum_pairwise_correlation": min(r[4] for r in corrs), "maximum_pairwise_correlation": max(r[4] for r in corrs),
        "within_study_exact_duplicate_pairs": duplicates,
        "limits": ["Descriptive outputs only; no statistical significance or effect replication claim", "No samples dropped; high expression correlation does not prove biological replicate quality", "Panel based on published protein assays does not assert matching RNA directions", "Shared FFA controls prevent interpreting inverse contrasts as independent rescue evidence", "No AI calls, no independent ingredient validation, no new discovery claimed"],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
