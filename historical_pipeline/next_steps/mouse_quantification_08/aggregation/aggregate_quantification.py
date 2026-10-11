"""Aggregate complete kallisto abundances; no quantification or inference is performed here.

Run this only after the external per-run FASTQ/quantification QC gate has passed.
Uses Python standard library. Synthetic tests live separately in test_aggregation.py.
"""
from __future__ import annotations

import argparse
from array import array
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import sys

HERE = Path(__file__).resolve().parent
STEPS = HERE.parent.parent
PREV = STEPS / "mouse_readiness_07"
DEFAULT_DATA = Path("D:/ResearchData/sesamin_mouse_2026-10-08")
GROUPS = ("ctrl", "HFD", "sess")
METRICS = ("est_counts", "tpm", "eff_length")
PSEUDOCOUNTS = (("0p1", 0.1), ("0", 0.0), ("0p5", 0.5), ("1", 1.0))
# kallisto plain text uses roughly six significant digits. An aggregate absolute
# deviation of five TPM (5 ppm of 1e6) covers worst-case six-digit rounding.
TPM_SUM_ABS_TOLERANCE = 5.000001


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                                   allow_nan=False) + "\n", encoding="utf-8")


def finite_nonnegative(value, context):
    result = float(value)
    require(math.isfinite(result) and result >= 0, f"Invalid {context}: {value!r}")
    return result


def read_run_manifest(path):
    manifest = load_json(path)
    entries = manifest["entries"]
    runs = {}
    for entry in entries:
        run = entry["run"]
        expected = {key: entry[key] for key in
                    ("run", "sample", "original_label", "provisional_group", "condition_mapping")}
        if run in runs:
            require(runs[run]["metadata"] == expected, f"Conflicting run metadata: {run}")
        else:
            runs[run] = {"metadata": expected, "mates": []}
        runs[run]["mates"].append(entry["mate"])
    require(len(runs) == 9 and len(entries) == 18, "Expected all nine runs and eighteen mates")
    require(all(sorted(rec["mates"]) == [1, 2] for rec in runs.values()), "Duplicated/missing mate")
    result = [runs[run]["metadata"] for run in sorted(runs)]
    require(Counter(r["provisional_group"] for r in result) == {g: 3 for g in GROUPS},
            "Expected three submitted libraries per provisional group")
    require(len({r["original_label"] for r in result}) == 9, "Duplicate library labels")
    return result


def read_mapping(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"transcript_id", "gene_id", "gene_version", "gene_name", "gene_biotype",
                    "transcript_biotype", "sequence_length", "source_fasta"}
        require(required.issubset(reader.fieldnames or []), "Unexpected tx2gene columns")
        rows = list(reader)
    require(bool(rows), "Empty tx2gene mapping")
    transcripts, gene_metadata = {}, {}
    for row in rows:
        tid, gid = row["transcript_id"], row["gene_id"]
        require(tid and gid, "Empty transcript/gene ID")
        require(tid not in transcripts, f"Duplicated transcript mapping: {tid}")
        require(int(row["sequence_length"]) > 0, f"Invalid reference length: {tid}")
        transcripts[tid] = row
        meta = tuple(row[k] for k in ("gene_version", "gene_name", "gene_biotype"))
        require(gid not in gene_metadata or gene_metadata[gid] == meta,
                f"Inconsistent gene metadata: {gid}")
        gene_metadata[gid] = meta
    # This order is the reference table order. Quantification rows are explicitly
    # joined by the complete, versioned transcript ID, never by row position.
    ids = list(transcripts)
    genes = sorted(gene_metadata)
    gene_index = {gid: i for i, gid in enumerate(genes)}
    tx_gene_index = [gene_index[transcripts[tid]["gene_id"]] for tid in ids]
    return transcripts, ids, genes, gene_metadata, tx_gene_index


def read_abundance(path, tx_ids, mapping):
    index = {tid: i for i, tid in enumerate(tx_ids)}
    values = {metric: array("d", [0.0]) * len(tx_ids) for metric in METRICS}
    seen, input_order = set(), []
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require({"target_id", "length", *METRICS}.issubset(reader.fieldnames or []),
                f"Unexpected abundance columns: {path}")
        for row in reader:
            tid = row["target_id"]
            require(tid in index, f"Unexpected transcript in {path}: {tid}")
            require(tid not in seen, f"Duplicated abundance transcript in {path}: {tid}")
            seen.add(tid)
            input_order.append(tid)
            require(float(row["length"]) == int(mapping[tid]["sequence_length"]),
                    f"Sequence length disagreement for {tid} in {path}")
            for metric in METRICS:
                values[metric][index[tid]] = finite_nonnegative(row[metric], f"{metric}/{tid}")
            require(values["eff_length"][index[tid]] > 0 or values["tpm"][index[tid]] == 0,
                    f"Positive TPM with zero effective length: {tid}")
    missing = set(tx_ids) - seen
    require(not missing, f"Missing {len(missing)} transcripts in {path}: {sorted(missing)[:3]}")
    total_tpm = math.fsum(values["tpm"])
    total_counts = math.fsum(values["est_counts"])
    require(abs(total_tpm - 1_000_000) <= TPM_SUM_ABS_TOLERANCE,
            f"TPM sum outside rounding tolerance in {path}: {total_tpm}")
    require(total_counts > 0, f"Zero estimated fragment count: {path}")
    return values, {"transcripts": len(tx_ids), "transcript_order_reindexed": input_order != tx_ids,
                    "transcript_tpm_sum": total_tpm, "transcript_est_counts_sum": total_counts,
                    "fractional_est_count_entries": sum(x != int(x) for x in values["est_counts"])}


def aggregate_gene(values, tx_gene_index, gene_count):
    result = {}
    for metric in ("est_counts", "tpm"):
        # Collect per-gene terms and use math.fsum to reduce accumulation error.
        terms = [[] for _ in range(gene_count)]
        for j, value in zip(tx_gene_index, values[metric], strict=True):
            terms[j].append(value)
        sums = array("d", (math.fsum(group) for group in terms))
        require(math.isclose(math.fsum(sums), math.fsum(values[metric]),
                             rel_tol=1e-12, abs_tol=1e-8), f"Gene mass conservation failed: {metric}")
        result[metric] = sums
    return result


def display_value(value):
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        require(math.isfinite(value), "Nonfinite output value")
        return format(value, ".17g")
    return value


def write_table(path, columns, rows):
    # Reproducible gzip headers avoid timestamps changing artifact hashes.
    path = Path(path)
    with path.open("wb") as binary:
        if path.suffix == ".gz":
            zipped = gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0)
            stream = io.TextIOWrapper(zipped, encoding="utf-8", newline="")
        else:
            stream = io.TextIOWrapper(binary, encoding="utf-8", newline="")
        with stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            for row in rows:
                writer.writerow([display_value(value) for value in row])


def direction(value):
    if value is None:
        return "NA"
    return "up" if value > 0 else "down" if value < 0 else "unchanged"


def log_ratio(a, b, pseudocount):
    if a + pseudocount == 0 or b + pseudocount == 0:
        return None
    return math.log2((a + pseudocount) / (b + pseudocount))


def direction_agreement(mouse, human):
    if mouse is None or human is None:
        return "NA"
    if mouse == 0 or human == 0:
        return "unchanged_in_at_least_one"
    return "same" if (mouse > 0) == (human > 0) else "opposite"


def read_human(path, candidates):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        records = list(csv.DictReader(stream))
    require(len(records) == len(candidates) == 25, "Human candidate denominator changed")
    human = {row["gene_id"]: row for row in records}
    require(len(human) == 25, "Duplicated human candidate")
    require(set(human) == {c["original_human_id"] for c in candidates}, "Human candidate set mismatch")
    for row in human.values():
        for column in ("log2_OA_vs_control", "log2_sesamin_vs_OA"):
            row[column] = float(row[column])
            require(math.isfinite(row[column]), f"Nonfinite original human contrast: {column}")
    return human


def candidate_reports(candidates, runs, genes, gene_values, human):
    gene_index = {gid: i for i, gid in enumerate(genes)}
    coverage, descriptive = [], []
    for candidate in candidates:
        hid = candidate["original_human_id"]
        homologues = candidate["mouse_orthologues"]
        primary = candidate["primary_comparison_eligible"]
        represented = [h["mouse_id"] in gene_index for h in homologues]
        # Primary reference coverage was frozen before any expression analysis.
        require(not primary or (len(homologues) == 1 and all(represented)),
                f"Missing/ambiguous primary candidate coverage: {hid}")
        coverage.append({"original_human_id": hid, "human_symbol": candidate["input_symbol"],
                         "orthology_status": candidate["status"], "primary_comparison_eligible": primary,
                         "mouse_orthologue_count": len(homologues),
                         "mouse_gene_ids": ";".join(h["mouse_id"] for h in homologues) or None,
                         "mouse_symbols": ";".join(h["mouse_symbol"] for h in homologues) or None,
                         "reference_gene_count": sum(represented),
                         "coverage_status": "no_mouse_orthologue" if not homologues else
                            "all_targets_in_reference" if all(represented) else "target_missing_from_reference"})
        # A one-to-many candidate gets one row per homologue, both retained. A
        # candidate without an orthologue gets one NA row, not a fabricated zero.
        for homologue in homologues or [None]:
            gid = homologue["mouse_id"] if homologue else None
            available = gid is not None and gid in gene_index
            row = {"original_human_id": hid, "human_symbol": candidate["input_symbol"],
                   "mouse_gene_id": gid, "mouse_symbol": homologue["mouse_symbol"] if homologue else None,
                   "orthology_status": candidate["status"], "primary_comparison_eligible": primary,
                   "comparison_role": "primary_one_to_one" if primary else
                         "secondary_one_to_many" if homologue else "not_comparable_no_orthologue",
                   "coverage_status": "quantified" if available else "no_mouse_orthologue" if gid is None
                         else "target_missing_from_reference"}
            values = [gene_values[r["run"]]["tpm"][gene_index[gid]] for r in runs] if available else [None] * len(runs)
            for run, value in zip(runs, values, strict=True):
                row[run["run"] + "_TPM"] = value
            means = {g: math.fsum(v for v, r in zip(values, runs, strict=True)
                                  if r["provisional_group"] == g) / 3 if available else None for g in GROUPS}
            row.update({f"mean_{g}_TPM_provisional": means[g] for g in GROUPS})
            row["all_libraries_TPM_lt1"] = all(v < 1 for v in values) if available else None
            for label, pc in PSEUDOCOUNTS:
                row[f"log2_HFD_vs_ctrl_pc{label}"] = log_ratio(means["HFD"], means["ctrl"], pc) if available else None
                row[f"log2_sess_vs_HFD_pc{label}"] = log_ratio(means["sess"], means["HFD"], pc) if available else None
            mouse_d, mouse_a = row["log2_HFD_vs_ctrl_pc0p1"], row["log2_sess_vs_HFD_pc0p1"]
            human_d, human_a = human[hid]["log2_OA_vs_control"], human[hid]["log2_sesamin_vs_OA"]
            row["mouse_HFD_vs_ctrl_direction"] = direction(mouse_d)
            row["mouse_sess_vs_HFD_direction"] = direction(mouse_a)
            row["human_OA_vs_ctrl_direction"] = direction(human_d)
            row["human_sesamin_vs_OA_direction"] = direction(human_a)
            row["cross_species_disease_direction"] = direction_agreement(mouse_d, human_d)
            row["cross_species_sesamin_direction"] = direction_agreement(mouse_a, human_a)
            row["both_cross_species_directions_same"] = (
                row["cross_species_disease_direction"] == row["cross_species_sesamin_direction"] == "same"
            ) if available else None
            row["mouse_two_contrasts_opposite_direction"] = (mouse_d * mouse_a < 0) if available else None
            descriptive.append(row)
    return coverage, descriptive


def execute(data_root, output_dir, manifest_path, orthology_path, plan_path, human_path, synthetic=False):
    begin = datetime.now(timezone.utc).isoformat()
    require(not output_dir.exists() or not any(output_dir.iterdir()),
            f"Refusing to overwrite non-empty output directory: {output_dir}")
    runs = read_run_manifest(manifest_path)
    plan, orthology = load_json(plan_path), load_json(orthology_path)
    candidates = orthology["candidates"]
    require(len(candidates) == 25, "Expected complete 25-candidate set")
    require(sum(c["primary_comparison_eligible"] for c in candidates) == 15, "Expected 15 primary candidates")
    require({c["original_human_id"] for c in candidates if c["primary_comparison_eligible"]}
            == set(plan["primary_human_gene_ids"]), "Frozen primary candidate set changed")
    mapping_path = data_root / "reference/tx2gene.tsv"
    inputs = [manifest_path, orthology_path, plan_path, human_path, mapping_path]
    expected_paths = [data_root / "quant" / r["run"] / "abundance.tsv" for r in runs]
    require(all(path.is_file() for path in expected_paths), "All nine abundance.tsv inputs must exist")
    inputs.extend(expected_paths)
    input_hashes = {str(p.resolve()): digest(p) for p in inputs}
    if not synthetic:
        for relative in ("fastq_manifest.json", "orthology/orthology.json"):
            selected = manifest_path if relative.startswith("fastq") else orthology_path
            require(digest(selected) == plan["input_snapshot_sha256"][relative], f"Frozen input changed: {relative}")
        require(digest(human_path) == orthology["input_sha256"], "Original human candidate source changed")
    mapping, tx_ids, genes, gene_metadata, tx_gene_index = read_mapping(mapping_path)
    require(set(plan["primary_mouse_gene_ids"]).issubset(genes), "Primary mouse gene coverage incomplete")
    human = read_human(human_path, candidates)
    tx_values, gene_values, qc = {}, {}, []
    for run, path in zip(runs, expected_paths, strict=True):
        values, metrics = read_abundance(path, tx_ids, mapping)
        aggregated = aggregate_gene(values, tx_gene_index, len(genes))
        tx_values[run["run"]], gene_values[run["run"]] = values, aggregated
        metrics.update(run=run["run"], original_label=run["original_label"],
                       gene_tpm_sum=math.fsum(aggregated["tpm"]),
                       gene_est_counts_sum=math.fsum(aggregated["est_counts"]),
                       mass_conservation_passed=True)
        qc.append(metrics)
    coverage, descriptive = candidate_reports(candidates, runs, genes, gene_values, human)
    # Detect input changes before emitting any products. No raw input is modified.
    require(all(digest(path) == fingerprint for path, fingerprint in input_hashes.items()),
            "An aggregation input changed during processing")
    output_dir.mkdir(parents=True, exist_ok=True)
    run_ids = [r["run"] for r in runs]
    outputs = []
    for level, ids, values, metrics in (("transcript", tx_ids, tx_values, METRICS),
                                      ("gene", genes, gene_values, ("est_counts", "tpm"))):
        for metric in metrics:
            filename = f"{level}_{metric}.tsv.gz"
            write_table(output_dir / filename, [level + "_id", *run_ids],
                        ([identifier, *(values[run][metric][i] for run in run_ids)]
                         for i, identifier in enumerate(ids)))
            outputs.append(filename)
    write_table(output_dir / "gene_metadata.tsv.gz",
                ["gene_id", "gene_version", "gene_name", "gene_biotype", "transcript_count"],
                ([gid, *gene_metadata[gid], n] for gid, n in sorted(Counter(row["gene_id"] for row in mapping.values()).items())))
    outputs.append("gene_metadata.tsv.gz")
    for filename, rows in (("library_manifest.tsv", runs), ("abundance_qc.tsv", qc),
                           ("candidate_coverage25.tsv", coverage), ("candidate_mouse_descriptive.tsv", descriptive)):
        write_table(output_dir / filename, list(rows[0]), ([r[key] for key in rows[0]] for r in rows))
        outputs.append(filename)
    primary_rows = [r for r in descriptive if r["primary_comparison_eligible"]]
    summary = {
        "status": "passed", "synthetic": synthetic, "started_utc": begin,
        "ended_utc": datetime.now(timezone.utc).isoformat(),
        "python": {"executable": sys.executable, "version": platform.python_version()},
        "command_argv": sys.argv, "script_sha256": digest(__file__), "input_sha256": input_hashes,
        "output_sha256": {name: digest(output_dir / name) for name in outputs},
        "library_count": len(runs), "transcript_count": len(tx_ids), "gene_count": len(genes),
        "human_candidate_count": len(coverage), "candidate_homologue_rows_including_NA": len(descriptive),
        "primary_comparison_gene_count": len(primary_rows),
        "coverage_counts": dict(Counter(row["coverage_status"] for row in coverage)),
        "primary_both_cross_species_directions_same_count": sum(r["both_cross_species_directions_same"] is True for r in primary_rows),
        "primary_all_libraries_TPM_lt1_count": sum(r["all_libraries_TPM_lt1"] is True for r in primary_rows),
        "primary_low_expression_gene_ids": [r["mouse_gene_id"] for r in primary_rows if r["all_libraries_TPM_lt1"] is True],
        "per_run": qc, "tpm_sum_absolute_tolerance": TPM_SUM_ABS_TOLERANCE,
        "validation": {"all_nine_libraries_retained": True, "exact_versioned_transcript_set": True,
                       "nonnegative_finite_abundance": True, "gene_mass_conservation": True,
                       "complete_human_candidate_denominator": True, "no_zero_imputation_for_absent_orthologues": True,
                       "input_hashes_unchanged_during_run": True},
        "interpretation": [
            "Descriptive exploratory estimates only; group labels ctrl/HFD/sess remain inferred.",
            "Nine submitted sequencing libraries do not establish nine independent animals or culture preparations.",
            "Gene estimated counts are sums of fractional transcript estimates; not raw integer counts.",
            "Gene TPM sums transcript TPM, retaining every biotype and candidate row, including all-below-1 TPM flags.",
            "Cross-species comparison uses direction only; human normalized counts and mouse TPM magnitudes are not pooled.",
            "ACAT2 homologues remain separate secondary rows; absent orthologues are NA and are not replaced with analogous genes.",
            "Positive directional agreement is not independent biological validation, causal evidence, a new mechanism, or clinical efficacy.",
            "No biological p values, FDR, confidence intervals, or inferential differential expression were calculated.",
            "FASTQ integrity and per-run quantification QC are performed upstream; this program validates abundance aggregation only."
        ]}
    write_json(output_dir / "aggregation_summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--manifest", type=Path, default=PREV / "fastq_manifest.json")
    parser.add_argument("--orthology", type=Path, default=PREV / "orthology/orthology.json")
    parser.add_argument("--plan", type=Path, default=PREV / "processing_plan.json")
    parser.add_argument("--human", type=Path, default=STEPS / "sesamin_exploration_04/priority_all.csv")
    args = parser.parse_args()
    summary = execute(args.data_root, args.output_dir or args.data_root / "matrices",
                      args.manifest, args.orthology, args.plan, args.human)
    print(json.dumps({k: summary[k] for k in ("status", "synthetic", "library_count", "transcript_count", "gene_count", "human_candidate_count", "primary_comparison_gene_count")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
