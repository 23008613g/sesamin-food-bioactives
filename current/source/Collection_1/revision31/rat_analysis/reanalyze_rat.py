"""Recompute the frozen rat comparison from deposited matrix and annotations.

Run with the project's .venv Python. Only this script's directory is written.
This is an explicit stage31 audit of known results, not preregistration.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
from itertools import combinations
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from scipy import stats

HERE = Path(__file__).resolve().parent
NEXT = HERE.parents[1]
PROJECT = NEXT.parent
OLD = NEXT / "analysis_strengthening_15/rat_validation"
OUTPUT = HERE / "results"
OUTPUT.mkdir(exist_ok=True)
PROJECT = HERE
INPUTS = {role: HERE / path for role, path in {'matrix': 'inputs/GSE42693_series_matrix.txt.gz', 'annotation': 'inputs/GPL1355.annot.gz', 'candidates': 'inputs/human_candidates_all25.csv', 'orthology': 'inputs/human_rat_orthology_ensembl116.tsv', 'entrez_crossreference': 'inputs/rat_entrez_ensembl116.tsv', 'old_mapping': 'inputs/candidate_mapping_all25.csv', 'old_results': 'inputs/candidate_results_all25_two_normalizations.csv', 'old_exact': 'inputs/posthoc_exact_sensitivity_all25.csv', 'old_samples': 'inputs/candidate_sample_values.csv', 'protocol': 'PROTOCOL.md'}.items()}

URLS = {
    "matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42693/matrix/GSE42693_series_matrix.txt.gz",
    "annotation": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1355/annot/GPL1355.annot.gz",
    "orthology": "https://jun2026.archive.ensembl.org/biomart/martservice",
    "entrez_crossreference": "https://jun2026.archive.ensembl.org/biomart/martservice",
    "source_paper": "https://doi.org/10.3177/jnsv.55.31",
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def save_csv(name, frame):
    frame.to_csv(OUTPUT / name, index=False, na_rep="NA", float_format="%.15g")


def table_block(text, kind):
    return text.split(f"!{kind}_table_begin\n", 1)[1].split(f"!{kind}_table_end", 1)[0]


def split_ids(value):
    return tuple(sorted({x.strip() for x in value.split("///") if x.strip() and x.strip() != "---"}))


checks = []


def check(name, condition, detail=""):
    checks.append({"check": name, "passed": bool(condition), "detail": detail})
    if not condition:
        raise AssertionError(f"{name}: {detail}")


manifest = [
    {"role": role, "path_relative_to_project": str(path.relative_to(PROJECT)), "sha256": sha(path), "bytes": path.stat().st_size,
     "source_url": URLS.get(role, "local frozen prior analysis"), "read_only": True}
    for role, path in INPUTS.items()
]
save_json("input_manifest.json", manifest)

matrix_text = gzip.open(INPUTS["matrix"], "rt", encoding="utf-8").read()
matrix = pd.read_csv(io.StringIO(table_block(matrix_text, "series_matrix")), sep="\t", index_col=0)
expected_samples = [f"GSM{n}" for n in range(1048036, 1048056)]
check("matrix_dimensions", matrix.shape == (31099, 20), str(matrix.shape))
check("matrix_probe_and_sample_unique", matrix.index.is_unique and matrix.columns.is_unique)
check("all_twenty_samples_in_deposited_order", matrix.columns.tolist() == expected_samples)
raw = matrix.to_numpy(dtype=float)
check("finite_nonnegative_signal", np.isfinite(raw).all() and (raw >= 0).all())

metadata = {}
for line in matrix_text.splitlines():
    if line.startswith("!Sample_") or line.startswith("!Series_"):
        row = next(csv.reader([line], delimiter="\t"))
        metadata.setdefault(row[0], []).append(row[1:])
processing = sorted(set(metadata["!Sample_data_processing"][0]))
check("deposited_processing_GCOS_global_scaling", len(processing) == 1 and "global scaling" in processing[0] and "1.4" in processing[0])
quantiles = np.quantile(raw, [0, .01, .25, .5, .75, .99, 1])
# No universal numerical cutoff identifies a log matrix. Interpret the observed
# distribution jointly with deposited GCOS processing and the source Methods.
scale_audit = {
    "decision": "retain_log2_deposited_signal_plus_1",
    "matrix_scale": "unlogged positive globally scaled Affymetrix signal; not raw CEL and not RMA",
    "scale_evidence": "Embedded GEO processing explicitly states GCOS 1.4 and global scaling; positive deposited signals extend to 18018.4, with Vnn1 control values 210.585 to 392.899. This is consistent with unlogged scaled signal, not a pre-logged expression matrix.",
    "deposited_processing": processing,
    "distribution_quantiles": {str(p): float(q) for p, q in zip([0, .01, .25, .5, .75, .99, 1], quantiles)},
    "all_values_above_100": int((raw > 100).sum()),
    "control_probe_means_close_to_one": int(np.isclose(raw[:, :5].mean(axis=1), 1).sum()),
    "original_paper_processing": "Ide et al. 2009, printed p.32, reports subsequent GeneSpring GX 7.3 per-chip 50th-percentile and per-gene control-mean normalization, followed by detection/intensity filtering. That pipeline is not reproduced or attributed to this deposited matrix.",
    "original_paper_animal_selection": "Experiment 1: seven male Sprague-Dawley rats per group; highest and lowest terminal-body-weight animals excluded from microarrays, leaving five individual livers per exposure. Groups balanced on initial mean body weight; random assignment not established from the inspected Methods.",
    "new_sample_or_probe_filter": False,
    "CEL_reprocessed": False,
    "already_log_or_double_log_error_found": False,
}
save_json("scale_audit.json", scale_audit)

sample_rows = []
group_order = ["control", "sesamin", "episesamin", "sesamolin"]
for j, sample in enumerate(expected_samples):
    sample_rows.append({
        "sample": sample, "source_title": metadata["!Sample_title"][0][j], "group": group_order[j // 5],
        "primary_comparison": j < 10, "biological_unit": "individual_rat_liver",
        "sex": "male", "strain": "Sprague-Dawley", "dietary_lignan_percent": 0 if j < 5 else .2,
        "duration_days": 15, "selected_per_original_group": 5, "original_animals_per_group": 7,
        "original_selection": "highest_and_lowest_terminal_body_weight_excluded",
        "current_exclusion": "none", "source_processing": processing[0],
        "raw_min": raw[:, j].min(), "raw_median": np.median(raw[:, j]), "raw_max": raw[:, j].max(),
    })
save_csv("sample_metadata_and_scale.csv", pd.DataFrame(sample_rows))

annotation_text = gzip.open(INPUTS["annotation"], "rt", encoding="utf-8").read()
ann = pd.read_csv(io.StringIO(table_block(annotation_text, "platform")), sep="\t", dtype=str, keep_default_na=False)
ann["entrez_set"] = ann["Gene ID"].map(split_ids)
ann["symbol_set"] = ann["Gene symbol"].map(split_ids)
check("annotation_dimensions_and_unique_ID", len(ann) == 31099 and ann.ID.is_unique)
check("annotation_exact_probe_set", set(ann.ID) == set(matrix.index))
platform_audit = pd.DataFrame({
    "probe_id": ann.ID, "deposited_entrez_ids": ann["Gene ID"], "deposited_symbols": ann["Gene symbol"],
    "n_unique_entrez_ids": ann.entrez_set.map(len), "n_unique_symbols": ann.symbol_set.map(len),
    "uniquely_annotated_by_entrez": ann.entrez_set.map(len) == 1,
})
platform_audit.to_csv(OUTPUT / "all_platform_probe_annotation_audit.csv.gz", index=False, compression="gzip")
orth = pd.read_csv(INPUTS["orthology"], sep="\t", dtype=str, keep_default_na=False)
orth.columns = ["human_id", "current_human_symbol", "rat_id", "rat_symbol", "orthology_type", "orthology_confidence", "rat_percent_id", "human_percent_id"]
xref = pd.read_csv(INPUTS["entrez_crossreference"], sep="\t", dtype=str, keep_default_na=False)
xref.columns = ["rat_id", "rat_symbol", "rat_entrez"]
candidates = pd.read_csv(INPUTS["candidates"], dtype={"gene_id": str})
check("fixed_candidate_family_size_and_uniqueness", len(candidates) == 25 and candidates.gene_id.is_unique)
mapping_rows, probe_decisions = [], []
for c in candidates.itertuples(index=False):
    targets = orth.loc[(orth.human_id == c.gene_id) & (orth.rat_id != "")]
    row = dict(human_id=c.gene_id, human_symbol=c.symbol, rat_id="", rat_symbol="", rat_entrez="",
               orthology_type="", orthology_confidence="", probe_ids="", probe_n=0, status="")
    if len(targets) != 1 or targets.iloc[0].orthology_type != "ortholog_one2one":
        row["status"] = "no_rat_orthologue" if targets.empty else "non_one_to_one_orthology"
        row["orthology_type"] = ";".join(sorted(set(targets.orthology_type)))
    else:
        target = targets.iloc[0]
        row.update(target[["rat_id", "rat_symbol", "orthology_type", "orthology_confidence"]].to_dict())
        entrez = sorted(set(xref.loc[xref.rat_id == row["rat_id"], "rat_entrez"]) - {""})
        if len(entrez) != 1:
            row["status"] = "ambiguous_or_missing_rat_entrez"
        else:
            row["rat_entrez"] = entrez[0]
            matches = ann.loc[ann.entrez_set.map(lambda ids: entrez[0] in ids)]
            eligible = matches.loc[matches.entrez_set.map(len) == 1]
            for p in matches.itertuples(index=False):
                probe_decisions.append({"human_id": c.gene_id, "human_symbol": c.symbol, "probe_id": p.ID,
                                        "target_entrez": entrez[0], "all_deposited_entrez": ";".join(p.entrez_set),
                                        "all_deposited_symbols": ";".join(p.symbol_set), "unique_entrez_count": len(p.entrez_set),
                                        "included": len(p.entrez_set) == 1,
                                        "reason": "unique_target_entrez" if len(p.entrez_set) == 1 else "multiple_gene_IDs"})
            row["probe_ids"] = ";".join(sorted(eligible.ID))
            row["probe_n"] = len(eligible)
            row["status"] = "mapped" if len(eligible) else "no_unambiguous_probe"
    mapping_rows.append(row)
mapping = pd.DataFrame(mapping_rows)
decisions = pd.DataFrame(probe_decisions)
save_csv("candidate_mapping_all25.csv", mapping)
save_csv("candidate_probe_mapping_decisions.csv", decisions)
save_csv("excluded_ambiguous_probes.csv", decisions.loc[~decisions.included])
old_map = pd.read_csv(INPUTS["old_mapping"], dtype=str, keep_default_na=False)
map_matches = mapping.astype(str).equals(old_map[mapping.columns].astype(str))
check("reconstructed_mapping_matches_frozen_stage15", map_matches)

log_signal = np.log2(raw + 1)


def quantile_normalize(a):
    """All-array mean order statistics, with exact tie blocks averaged."""
    rank_target = np.sort(a, axis=0).mean(axis=1)
    z = np.empty_like(a)
    for col in range(a.shape[1]):
        _, inverse, counts = np.unique(a[:, col], return_inverse=True, return_counts=True)
        starts = np.r_[0, np.cumsum(counts)[:-1]]
        group_targets = np.add.reduceat(rank_target, starts) / counts
        z[:, col] = group_targets[inverse]
    return z


qn = quantile_normalize(log_signal)
allocations = list(combinations(range(10), 5))
check("all_exact_allocations", len(allocations) == 252 and len(set(allocations)) == 252)


def comparison(x):
    control, treatment = x[:5], x[5:10]
    difference = float(treatment.mean() - control.mean())
    part_c = control.var(ddof=1) / len(control)
    part_t = treatment.var(ddof=1) / len(treatment)
    se = float(np.sqrt(part_c + part_t))
    df = float((part_c + part_t) ** 2 / (part_c**2 / 4 + part_t**2 / 4))
    t_stat = difference / se
    p = float(2 * stats.t.sf(abs(t_stat), df))
    scipy_test = stats.ttest_ind(treatment, control, equal_var=False)
    check("manual_vs_scipy_welch", np.isclose(p, scipy_test.pvalue, atol=1e-14, rtol=1e-10) and np.isclose(t_stat, scipy_test.statistic, atol=1e-12))
    radius = float(stats.t.ppf(.975, df) * se)
    finite10 = x[:10]
    perm_diff = []
    for chosen in allocations:
        mask = np.zeros(10, dtype=bool)
        mask[list(chosen)] = True
        perm_diff.append(float(finite10[mask].mean() - finite10[~mask].mean()))
    extreme = int(np.sum(np.abs(perm_diff) >= abs(difference) - 1e-12))
    check("exact_tail_includes_complements", extreme >= 2 and extreme % 2 == 0)
    return {"effect_log2_signal": difference, "control_mean_log2": float(control.mean()), "sesamin_mean_log2": float(treatment.mean()),
            "control_sd_log2": float(control.std(ddof=1)), "sesamin_sd_log2": float(treatment.std(ddof=1)),
            "welch_se": se, "welch_t": t_stat, "welch_df": df, "ci95_low": difference - radius, "ci95_high": difference + radius,
            "welch_p": p, "exact_p": extreme / 252, "exact_exceedances": extreme, "exact_allocations": 252,
            "complete_group_separation": bool(control.max() < treatment.min() or treatment.max() < control.min()),
            "median_group_difference": float(np.median(treatment) - np.median(control))}


results, animal_values, probe_values = [], [], []
for method, transformed in [("deposited_log2", log_signal), ("quantile_sensitivity", qn)]:
    for mapped in mapping.to_dict("records"):
        result = {**mapped, "method": method}
        if mapped["status"] != "mapped":
            results.append(result)
            # Keep explicit 25 x 20 x 2 candidate/animal grid; absent is NA, never zero.
            for j, sample in enumerate(expected_samples):
                animal_values.append({"human_id": mapped["human_id"], "human_symbol": mapped["human_symbol"], "method": method,
                                      "sample": sample, "group": group_order[j // 5], "primary_comparison": j < 10,
                                      "mapping_status": mapped["status"], "probe_n": 0, "aggregate_log2_signal": np.nan})
            continue
        probe_ids = mapped["probe_ids"].split(";")
        probe_ix = matrix.index.get_indexer(probe_ids)
        check("eligible_probes_exist_in_matrix", (probe_ix >= 0).all())
        x = np.median(transformed[probe_ix, :], axis=0)
        effects = transformed[probe_ix, 5:10].mean(axis=1) - transformed[probe_ix, :5].mean(axis=1)
        result.update(comparison(x))
        result.update({"probe_effect_min": float(effects.min()), "probe_effect_max": float(effects.max()),
                       "probe_directions_mixed": bool(effects.min() < 0 < effects.max())})
        results.append(result)
        for j, sample in enumerate(expected_samples):
            animal_values.append({"human_id": mapped["human_id"], "human_symbol": mapped["human_symbol"], "method": method,
                                  "sample": sample, "group": group_order[j // 5], "primary_comparison": j < 10,
                                  "mapping_status": mapped["status"], "probe_n": mapped["probe_n"], "aggregate_log2_signal": float(x[j])})
            for probe_id, ix in zip(probe_ids, probe_ix):
                probe_values.append({"human_id": mapped["human_id"], "human_symbol": mapped["human_symbol"], "method": method,
                                     "probe_id": probe_id, "sample": sample, "group": group_order[j // 5], "primary_comparison": j < 10,
                                     "deposited_unlogged_signal": float(raw[ix, j]), "transformed_signal": float(transformed[ix, j])})
results = pd.DataFrame(results)


def bh(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p, kind="stable")
    q_order = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.minimum(q_order, 1)
    check("BH_vs_scipy", np.allclose(q, stats.false_discovery_control(p, method="bh"), atol=1e-14))
    return q


for method in results.method.unique():
    ix = results.index[results.method == method]
    check("all25_per_method", len(ix) == 25)
    for test in ["welch", "exact"]:
        p = results.loc[ix, f"{test}_p"]
        available = p.notna()
        results.loc[ix, f"{test}_q_bh25"] = bh(p.fillna(1))
        results.loc[ix[~available], f"{test}_q_bh25"] = np.nan
        results.loc[ix[available], f"{test}_q_bh_measurable_secondary"] = bh(p.loc[available])
        results.loc[ix, f"{test}_measurable_family_n"] = int(available.sum())
        check("missing_reported_P_and_q_stay_NA", results.loc[ix[~available], [f"{test}_p", f"{test}_q_bh25", f"{test}_q_bh_measurable_secondary"]].isna().all().all())
results["analysis_status"] = "existing_analysis_recomputed"
results["measurable_BH_status"] = "new_posthoc_secondary_sensitivity"
save_csv("candidate_results_all25_welch_exact_two_normalizations.csv", results)
animals = pd.DataFrame(animal_values)
save_csv("candidate_sample_values_all25.csv", animals)
save_csv("candidate_probe_values_long.csv", pd.DataFrame(probe_values))

# Compare reconstructed values with saved source tables, never use old values in calculation.
old_result = pd.read_csv(INPUTS["old_results"])
old_exact = pd.read_csv(INPUTS["old_exact"])
differences = []
for new in results.to_dict("records"):
    prior = old_result.loc[(old_result.human_id == new["human_id"]) & (old_result.method == new["method"])].iloc[0]
    prior_exact = old_exact.loc[(old_exact.human_id == new["human_id"]) & (old_exact.method == new["method"])].iloc[0]
    field_map = {"effect_log2_signal": "effect_log2_signal", "ci95_low": "ci95_low", "ci95_high": "ci95_high",
                 "welch_df": "welch_df", "welch_p": "p_value", "welch_q_bh25": "q_bh_fixed25"}
    for new_field, old_field in field_map.items():
        old_value, new_value = prior.get(old_field, np.nan), new.get(new_field, np.nan)
        differences.append({"human_id": new["human_id"], "human_symbol": new["human_symbol"], "method": new["method"],
                            "field": new_field, "old_value": old_value, "new_value": new_value, "new_minus_old": new_value - old_value,
                            "within_saved_precision_tolerance": bool(np.isclose(new_value, old_value, atol=1e-9, rtol=1e-10, equal_nan=True))})
    for new_field, old_field in {"exact_p": "permutation_p", "exact_q_bh25": "permutation_q_bh25"}.items():
        old_value, new_value = prior_exact.get(old_field, np.nan), new.get(new_field, np.nan)
        differences.append({"human_id": new["human_id"], "human_symbol": new["human_symbol"], "method": new["method"],
                            "field": new_field, "old_value": old_value, "new_value": new_value, "new_minus_old": new_value - old_value,
                            "within_saved_precision_tolerance": bool(np.isclose(new_value, old_value, atol=1e-9, rtol=1e-10, equal_nan=True))})
diff = pd.DataFrame(differences)
save_csv("old_new_statistic_differences.csv", diff)
old_animals = pd.read_csv(INPUTS["old_samples"])
animal_compare = animals.loc[animals.mapping_status == "mapped"].merge(old_animals[["human_id", "method", "sample", "aggregate_log2_signal"]],
                            on=["human_id", "method", "sample"], suffixes=("_new", "_old"), validate="one_to_one")
animal_compare["new_minus_old"] = animal_compare.aggregate_log2_signal_new - animal_compare.aggregate_log2_signal_old
save_csv("old_new_animal_value_differences.csv", animal_compare)
check("saved_statistic_precision_matches", diff.within_saved_precision_tolerance.all())
check("saved_animal_precision_matches", np.allclose(animal_compare.aggregate_log2_signal_new, animal_compare.aggregate_log2_signal_old, atol=1e-9, rtol=1e-10))
check("full_candidate_animal_grid_preserved", len(animals) == 25 * 20 * 2)
for item in manifest:
    check("input_hash_unchanged_" + item["role"], sha(PROJECT / item["path_relative_to_project"]) == item["sha256"])
save_csv("verification_checks.csv", pd.DataFrame(checks))

vnn = results.loc[results.human_symbol == "VNN1"]
summary = {
    "execution_status": "EXECUTED_AND_CHECKED",
    "executed_utc": datetime.now(timezone.utc).isoformat(),
    "script_sha256": sha(__file__), "protocol_sha256": sha(HERE / "PROTOCOL.md"),
    "environment": {"executable": sys.executable, "python": sys.version, "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "platform": platform.platform()},
    "dimensions": {"probes": 31099, "arrays": 20, "primary_control_animals": 5, "primary_sesamin_animals": 5, "fixed_candidate_family": 25,
                   "measurable_candidates": int((mapping.status == "mapped").sum()), "unavailable_candidates": int((mapping.status != "mapped").sum())},
    "mapping_status_counts": mapping.status.value_counts().to_dict(),
    "single_probe_mapped_candidates": int(((mapping.status == "mapped") & (mapping.probe_n == 1)).sum()),
    "mult_probe_mapped_candidates": int((mapping.probe_n > 1).sum()),
    "ambiguous_matching_probes_excluded": int((~decisions.included).sum()),
    "scale_decision": scale_audit["decision"], "double_log_error_found": False,
    "old_new_comparison": {"statistics_compared_including_NA": len(diff), "all_match_within_saved_precision": bool(diff.within_saved_precision_tolerance.all()),
                           "maximum_absolute_statistic_difference": float(diff.new_minus_old.abs().max()), "animal_values_compared": len(animal_compare),
                           "maximum_absolute_animal_difference": float(animal_compare.new_minus_old.abs().max())},
    "VNN1": vnn[["method", "probe_ids", "effect_log2_signal", "ci95_low", "ci95_high", "welch_p", "welch_q_bh25", "exact_p", "exact_exceedances", "exact_q_bh25", "welch_q_bh_measurable_secondary", "exact_q_bh_measurable_secondary"]].to_dict("records"),
    "fdr05_bh25": {method: {test: results.loc[(results.method == method) & (results[f"{test}_q_bh25"] < .05), "human_symbol"].tolist() for test in ["welch", "exact"]} for method in results.method.unique()},
    "checks_passed": len(checks), "checks_failed": 0,
    "bounded_scope": ["Existing deposited-matrix statistics and mapping recomputed from unchanged source inputs; prior conclusions were known.",
                      "New secondary measurable-family BH adjustment reported alongside, never replacing fixed BH25.",
                      "No CEL reprocessing, detection-call filtering, probe sequence remapping, raw read rerun or new experiment.",
                      "Intervals are nominal and not simultaneous; original selection and unverified random allocation remain limitations."]
}
save_json("analysis_summary.json", summary)
save_json("output_manifest.json", [{"path": p.relative_to(HERE).as_posix(), "bytes": p.stat().st_size, "sha256": sha(p)} for p in sorted(OUTPUT.glob("*"))])
print(json.dumps({"status": summary["execution_status"], "dimensions": summary["dimensions"], "scale": summary["scale_decision"], "checks_passed": len(checks), "VNN1": summary["VNN1"], "old_new_comparison": summary["old_new_comparison"]}, ensure_ascii=False, indent=2))
