"""Meaningful synthetic checks; never uses actual mouse expression estimates."""
import csv
from datetime import datetime, timezone
import gzip
import json
import math
from pathlib import Path
import sys

import aggregate_quantification as agg


HERE = Path(__file__).resolve().parent
FIXTURE = HERE / "syntheticfixtures" / datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%S_%fZ")
FIXTURE.mkdir(parents=True, exist_ok=False)
checks = []


def check(name, condition):
    checks.append({"name": name, "passed": bool(condition)})
    if not condition:
        raise AssertionError(name)


def raises_valueerror(name, fn):
    try:
        fn()
    except ValueError:
        check(name, True)
    else:
        check(name, False)


def table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    agg.write_table(path, list(rows[0]), ([r[k] for k in rows[0]] for r in rows))


def read_table(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    # Real frozen identifiers/labels only. All abundances below are manufactured.
    plan = agg.load_json(agg.PREV / "processing_plan.json")
    orthology = agg.load_json(agg.PREV / "orthology/orthology.json")
    manifest = agg.load_json(agg.PREV / "fastq_manifest.json")
    for name, value in (("plan.json", plan), ("orthology.json", orthology), ("manifest.json", manifest)):
        agg.write_json(FIXTURE / name, value)
    runs = agg.read_run_manifest(FIXTURE / "manifest.json")
    genes = [h["mouse_id"] for c in orthology["candidates"] for h in c["mouse_orthologues"]]
    vnn1 = next(c["mouse_orthologues"][0]["mouse_id"] for c in orthology["candidates"] if c["input_symbol"] == "VNN1")
    genes.append("SYNTHETIC_BACKGROUND_GENE")
    mapping = []
    for i, gid in enumerate(genes):
        for isoform in range(2 if i == 0 else 1):
            mapping.append({"transcript_id": f"SYNTHETIC_TX_{i}_{isoform}.1", "gene_id": gid,
                            "gene_version": "1", "gene_name": f"SYNTHETIC_{i}",
                            "gene_biotype": "synthetic_biotype", "transcript_biotype": "synthetic_biotype",
                            "sequence_length": "100", "source_fasta": "SYNTHETIC_ONLY"})
    table(FIXTURE / "reference/tx2gene.tsv", mapping)
    abundance_by_run = {}
    for i, run in enumerate(runs):
        group = run["provisional_group"]
        multiplier = {"ctrl": 1, "HFD": 2, "sess": 0.5}[group]
        gene_tpm = {g: (j + 1) * 0.1 * multiplier for j, g in enumerate(genes[:-1])}
        gene_tpm[vnn1] = {"ctrl": 0.0, "HFD": 0.2, "sess": 0.1}[group]
        gene_tpm[genes[-1]] = 1_000_000 - math.fsum(gene_tpm.values())
        rows = []
        for m in mapping:
            tpm = gene_tpm[m["gene_id"]] / (2 if m["gene_id"] == genes[0] else 1)
            rows.append({"target_id": m["transcript_id"], "length": 100, "eff_length": 50,
                         "est_counts": tpm * 0.125, "tpm": tpm})
        if i % 2:
            rows.reverse()
        abundance_by_run[run["run"]] = rows
        table(FIXTURE / "quant" / run["run"] / "abundance.tsv", rows)
    with (FIXTURE / "human.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["gene_id", "log2_OA_vs_control", "log2_sesamin_vs_OA"])
        for candidate in orthology["candidates"]:
            writer.writerow([candidate["original_human_id"], 1, -1])
    summary = agg.execute(FIXTURE, FIXTURE / "output", FIXTURE / "manifest.json",
                          FIXTURE / "orthology.json", FIXTURE / "plan.json", FIXTURE / "human.csv", synthetic=True)
    check("synthetic_flag_and_all_nine", summary["synthetic"] and summary["library_count"] == 9)
    check("all_25_human_candidates_retained", summary["human_candidate_count"] == 25)
    check("26_rows_with_both_ACAT2_homologues", summary["candidate_homologue_rows_including_NA"] == 26)
    check("15_primary_pairs", summary["primary_comparison_gene_count"] == 15)
    check("nonpositional_join_reindexes_shuffled_rows", sum(r["transcript_order_reindexed"] for r in summary["per_run"]) == 4)
    check("fractional_counts_preserved_in_inputs", all(r["fractional_est_count_entries"] > 0 for r in summary["per_run"]))
    gene_counts = {r["gene_id"]: r for r in read_table(FIXTURE / "output/gene_est_counts.tsv.gz")}
    check("isoform_sum_fractional_counts_exact", math.isclose(float(gene_counts[genes[0]][runs[0]["run"]]), 0.0125, abs_tol=1e-15))
    descriptive = read_table(FIXTURE / "output/candidate_mouse_descriptive.tsv")
    vrow = next(r for r in descriptive if r["human_symbol"] == "VNN1")
    check("zero_pseudocount_undefined_ratio_NA", vrow["log2_HFD_vs_ctrl_pc0"] == "NA")
    check("zero_pseudocount_positive_means_valid", float(vrow["log2_sess_vs_HFD_pc0"]) == -1)
    check("primary_pseudocount_formula", math.isclose(float(vrow["log2_HFD_vs_ctrl_pc0p1"]), math.log2(3), abs_tol=1e-14))
    check("low_expression_retained_and_flagged", vrow["all_libraries_TPM_lt1"] == "true")
    cxcl8 = next(r for r in descriptive if r["human_symbol"] == "CXCL8")
    check("no_orthologue_not_fabricated_zero", cxcl8["mouse_gene_id"] == "NA" and all(cxcl8[r["run"] + "_TPM"] == "NA" for r in runs))
    acat = [r for r in descriptive if r["human_symbol"] == "ACAT2"]
    check("one_to_many_both_secondary", len(acat) == 2 and all(r["comparison_role"] == "secondary_one_to_many" for r in acat))
    first_run = runs[0]["run"]
    data = abundance_by_run[first_run]
    txmap, txids, _, _, _ = agg.read_mapping(FIXTURE / "reference/tx2gene.tsv")
    mutations = {
        "duplicate_transcript_rejected": data + [data[0]],
        "missing_transcript_rejected": data[:-1],
        "unknown_transcript_rejected": [{**data[0], "target_id": "UNKNOWN.1"}] + data[1:],
        "wrong_version_rejected": [{**data[0], "target_id": data[0]["target_id"] + "9"}] + data[1:],
        "negative_abundance_rejected": [{**data[0], "est_counts": -0.1}] + data[1:],
        "nan_abundance_rejected": [{**data[0], "tpm": "NaN"}] + data[1:],
        "wrong_reference_length_rejected": [{**data[0], "length": 101}] + data[1:],
        "TPM_sum_error_rejected": [{**r, "tpm": r["tpm"] * 2} for r in data],
    }
    for name, rows in mutations.items():
        path = FIXTURE / "negative_cases" / (name + ".tsv")
        table(path, rows)
        raises_valueerror(name, lambda p=path: agg.read_abundance(p, txids, txmap))
    raises_valueerror("completed_outputs_not_overwritten", lambda: agg.execute(
        FIXTURE, FIXTURE / "output", FIXTURE / "manifest.json", FIXTURE / "orthology.json",
        FIXTURE / "plan.json", FIXTURE / "human.csv", synthetic=True))
    result = {"synthetic_only": True, "actual_mouse_expression_used": False,
              "status": "passed", "checks": checks, "check_count": len(checks),
              "fixtures": str(FIXTURE), "command_argv": sys.argv,
              "test_script_sha256": agg.digest(__file__), "aggregation_script_sha256": agg.digest(agg.__file__)}
    agg.write_json(FIXTURE / "test_results.json", result)
    agg.write_json(HERE / "synthetic_test_results.json", result)
    print(json.dumps({"synthetic_only": True, "passed": len(checks), "fixtures": str(FIXTURE)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
