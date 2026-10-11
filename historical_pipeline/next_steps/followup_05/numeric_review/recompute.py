"""Independent limited numerical review; does not import analysis scripts."""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
PREV = STAGE.parent / "sesamin_exploration_04"
ORG = STAGE.parent / "organoid_expression"
CHECKS = []
USED = set()


def digest(path):
    USED.add(path.resolve())
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    USED.add(path.resolve())
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path, **kwargs):
    USED.add(path.resolve())
    return pd.read_csv(path, **kwargs)


def check(name, condition, **details):
    CHECKS.append({"name": name, "passed": bool(condition), **details})


def close(name, actual, expected):
    a = np.asarray(actual, dtype=float)
    b = np.asarray(expected, dtype=float)
    check(name, a.shape == b.shape and np.allclose(a, b, rtol=1e-11, atol=1e-10, equal_nan=True),
          scalar_comparisons=int(a.size), max_absolute_error=float(np.max(np.abs(a-b)[np.isfinite(a-b)])) if np.isfinite(a-b).any() else None)


def ratio(numerator, denominator, pseudo):
    if numerator + pseudo <= 0 or denominator + pseudo <= 0:
        return float("nan")
    return math.log2(numerator + pseudo) - math.log2(denominator + pseudo)


def bulk():
    protocol = read_json(PREV / "protocol.json")
    freeze = read_json(STAGE / "robustness_freeze.json")
    check("robustness_protocol_freeze", digest(STAGE / "robustness_protocol.json") == freeze["protocol_sha256"])
    candidates = read_csv(PREV / "priority_all.csv", keep_default_na=False).set_index("gene_id")
    check("unchanged_candidate_set", digest(PREV / "priority_all.csv") == freeze["candidate_sha256"] and len(candidates) == 25)
    path = (PREV / protocol["input"]).resolve()
    check("raw_bulk_input_hash", digest(path) == protocol["input_sha256"])
    names = protocol["sample_order"]
    raw = read_csv(path, sep="\t", usecols=["name"] + names).set_index("name")[names]
    # Compute the full positive-reference normalization again, using log ratios.
    positive = raw.loc[(raw > 0).all(axis=1)].to_numpy(dtype=float)
    log_values = np.log(positive)
    sample_ratios = np.exp(log_values - log_values.mean(axis=1)[:, None])
    factors = np.median(sample_ratios, axis=0)
    factors /= math.exp(float(np.log(factors).mean()))
    qc = read_csv(PREV / "sample_qc.csv").set_index("sample").loc[names]
    close("all_nine_size_factors_from_raw", qc.median_ratio_factor, factors)
    close("all_nine_library_totals", qc.assigned_counts, raw.sum(axis=0))
    normalized = raw.loc[candidates.index].div(pd.Series(factors, index=names), axis="columns")
    sample_output = read_csv(STAGE / "bulk_candidate_samples.csv").set_index(["gene_id", "sample"])
    check("bulk_sample_table_keys", len(sample_output) == 225 and sample_output.index.is_unique)
    expected_raw, expected_normalized, actual_raw, actual_normalized = [], [], [], []
    for gene in candidates.index:
        for sample in names:
            row = sample_output.loc[(gene, sample)]
            actual_raw.append(row.raw_count)
            actual_normalized.append(row.normalized_count)
            expected_raw.append(raw.at[gene, sample])
            expected_normalized.append(normalized.at[gene, sample])
            check("sample_group_" + gene + "_" + sample, row.group == qc.at[sample, "group"])
    close("225_raw_counts", actual_raw, expected_raw)
    close("225_normalized_counts", actual_normalized, expected_normalized)
    groups = {group: qc.index[qc.group == group].tolist() for group in ["Control", "OA", "OA + sesamin"]}
    pseudo_output = read_csv(STAGE / "bulk_pseudocount_sensitivity.csv").set_index(["gene_id", "pseudocount"])
    ranges = read_csv(STAGE / "bulk_group_ranges.csv").set_index(["gene_id", "group"])
    separation = read_csv(STAGE / "bulk_group_separation.csv").set_index("gene_id")
    spotlight = {}
    expected_effects, actual_effects = [], []
    for gene in candidates.index:
        arrays = {group: normalized.loc[gene, samples].to_numpy() for group, samples in groups.items()}
        means = {group: float(values.mean()) for group, values in arrays.items()}
        for group, values in arrays.items():
            row = ranges.loc[(gene, group)]
            close("bulk_range_" + gene + "_" + group,
                  row[["mean_normalized", "median_normalized", "min_normalized", "max_normalized", "max_library_share_normalized"]],
                  [values.mean(), np.median(values), values.min(), values.max(), values.max() / values.sum() if values.sum() else np.nan])
            check("bulk_nonzero_" + gene + "_" + group, row.raw_nonzero_libraries == int((raw.loc[gene, groups[group]] > 0).sum()))
        gene_effects = []
        for pseudo in [0, 0.5, 1, 5]:
            expected = [ratio(means["OA"], means["Control"], pseudo), ratio(means["OA + sesamin"], means["OA"], pseudo), ratio(means["OA + sesamin"], means["Control"], pseudo)]
            actual = pseudo_output.loc[(gene, pseudo), ["OA_vs_control", "sesamin_vs_OA", "sesamin_vs_control"]].to_numpy(dtype=float)
            expected_effects.extend(expected)
            actual_effects.extend(actual)
            gene_effects.append({"pseudocount": pseudo, "OA_vs_control": expected[0], "sesamin_vs_OA": expected[1], "sesamin_vs_control": expected[2]})
        counts = []
        for numerator, denominator, direction in [("OA", "Control", candidates.at[gene, "log2_OA_vs_control"]), ("OA + sesamin", "OA", candidates.at[gene, "log2_sesamin_vs_OA"])]:
            n = sum((a > b if direction > 0 else a < b) for a in arrays[numerator] for b in arrays[denominator])
            counts.append(n)
        close("dependent_pair_counts_" + gene, separation.loc[gene, ["OA_control_pairs_expected_direction_of9", "sesamin_OA_pairs_expected_direction_of9"]], counts)
        if candidates.at[gene, "symbol"] in ["CXCL8", "VNN1"]:
            spotlight[candidates.at[gene, "symbol"]] = {"gene_id": gene, "raw_counts_in_protocol_order": [int(x) for x in raw.loc[gene, names]], "group_means_normalized": means, "pseudo_effects": gene_effects, "dependent_pair_counts_of9": [int(x) for x in counts]}
    close("all_100_pseudocount_rows_300_effects", actual_effects, expected_effects)
    return candidates, {"candidate_n": len(candidates), "raw_gene_n": len(raw), "positive_reference_gene_n": len(positive), "sample_order": names, "size_factors": factors.tolist(), "spotlight": spotlight}


def organoid(candidates):
    matrix = read_csv(ORG / "all_gene_library_subtype_counts.tsv.gz", sep="\t", keep_default_na=False).set_index("gene_id")
    meta = read_csv(ORG / "library_subtype_summary.csv").set_index("count_column")
    libraries = meta.gsm.drop_duplicates().tolist()
    subtypes = [x for x in meta.subtype.drop_duplicates() if x != "ALL_HEPATOCYTE_LIKE"]
    all_label = "ALL_HEPATOCYTE_LIKE"
    counts = matrix[meta.index]
    totals = counts.sum(axis=0)
    close("organoid_24_total_UMIs_from_all_gene_matrix", meta.total_counts, totals)
    check("organoid_complete_design", len(libraries) == 4 and len(subtypes) == 5)
    library_weights = {}
    for lib in libraries:
        columns = [lib + "|" + subtype for subtype in subtypes]
        close("organoid_all_gene_subtype_sum_" + lib, counts[columns].sum(axis=1), counts[lib + "|" + all_label])
        library_weights[lib] = np.array([totals[lib + "|" + subtype] / totals[lib + "|" + all_label] for subtype in subtypes])
    fixed = np.vstack([library_weights[lib] for lib in libraries]).mean(axis=0)
    reported_weights = read_csv(STAGE / "organoid_reference_weights.csv").set_index("subtype").loc[subtypes]
    close("five_reference_weights", reported_weights.reference_UMI_weight, fixed)
    fractions = read_csv(STAGE / "organoid_library_UMI_fractions.csv").set_index("gsm")
    close("twenty_library_UMI_fractions", fractions.loc[libraries, subtypes], np.vstack(list(library_weights.values())))
    min_cells = [int(meta.loc[meta.subtype == subtype, "cells_aggregated"].min()) for subtype in subtypes]
    close("five_minimum_subtype_cell_counts", reported_weights.minimum_cells_across_four_libraries, min_cells)
    reweighted = read_csv(STAGE / "organoid_reweighted_samples.csv").set_index(["gene_id", "gsm"])
    contrasts = read_csv(STAGE / "organoid_composition_sensitivity.csv").set_index("gene_id")
    context = read_csv(PREV / "organoid_candidate_context.csv").set_index(["gene_id", "gsm", "subtype"])
    within = read_csv(STAGE / "organoid_within_subtype.csv").set_index(["gene_id", "subtype"])
    check("organoid_unique_output_keys", reweighted.index.is_unique and len(reweighted) == 100 and within.index.is_unique and len(within) == 125)
    covered, missing, reconstruction_errors, spotlight = [], [], [], {}
    source_cpm_actual, source_cpm_expected = [], []
    source_umi_actual, source_umi_expected = [], []
    fixed_actual, fixed_expected, observed_actual, observed_expected = [], [], [], []
    for gene in candidates.index:
        if gene not in matrix.index:
            missing.append(gene)
            row = reweighted.loc[gene]
            check("uncovered_preserved_" + gene, (~row.feature_covered).all() and row.observed_CPM.isna().all() and row.fixed_weight_CPM.isna().all())
            source = context.loc[gene]
            check("uncovered_context_" + gene, source.status.eq("feature_not_covered").all() and source.CPM.isna().all() and source.candidate_UMI.isna().all())
            continue
        covered.append(gene)
        per_lib = {}
        for lib in libraries:
            subtype_cpm = []
            for subtype in subtypes + [all_label]:
                column = lib + "|" + subtype
                umi = int(counts.at[gene, column])
                cpm = umi * 1000000 / int(totals[column])
                source = context.loc[(gene, lib, subtype)]
                source_umi_actual.append(source.candidate_UMI)
                source_umi_expected.append(umi)
                source_cpm_actual.append(source.CPM)
                source_cpm_expected.append(cpm)
                check("organoid_context_status_" + gene + "_" + column, source.status == ("observed_positive" if umi else "observed_zero"))
                if subtype != all_label:
                    subtype_cpm.append(cpm)
                    out = within.loc[(gene, subtype)]
                    close("within_library_" + gene + "_" + column, [out[lib + "_CPM"], out[lib + "_UMI"]], [cpm, umi])
            observed = int(counts.at[gene, lib + "|" + all_label]) * 1000000 / int(totals[lib + "|" + all_label])
            reconstructed = math.fsum(w * x for w, x in zip(library_weights[lib], subtype_cpm))
            standardized = math.fsum(w * x for w, x in zip(fixed, subtype_cpm))
            reconstruction_errors.append(abs(observed-reconstructed))
            row = reweighted.loc[(gene, lib)]
            check("organoid_reweighted_coverage_" + gene + "_" + lib, row.feature_covered)
            fixed_actual.append(row.fixed_weight_CPM)
            fixed_expected.append(standardized)
            observed_actual.append(row.observed_CPM)
            observed_expected.append(observed)
            per_lib[lib] = {"condition": meta.at[lib + "|" + all_label, "group"], "observed": observed, "standardized": standardized, "reconstructed": reconstructed}
        values = {}
        for prefix, metric in [("observed_CPM", "observed"), ("fixed_weight_CPM", "standardized")]:
            control = np.mean([x[metric] for x in per_lib.values() if x["condition"] == "CTRL-PA"])
            pa = np.mean([x[metric] for x in per_lib.values() if x["condition"] == "PA"])
            values.update({prefix + "_CTRL_mean": float(control), prefix + "_PA_mean": float(pa), prefix + "_log2_ratio_plus1": ratio(pa, control, 1)})
        close("organoid_contrast_" + gene, contrasts.loc[gene, list(values)], list(values.values()))
        for subtype in subtypes:
            row = within.loc[(gene, subtype)]
            controls = [context.at[(gene, lib, subtype), "CPM"] for lib in libraries if meta.at[lib + "|" + subtype, "group"] == "CTRL-PA"]
            pa = [context.at[(gene, lib, subtype), "CPM"] for lib in libraries if meta.at[lib + "|" + subtype, "group"] == "PA"]
            close("within_summary_" + gene + "_" + subtype, row[["control_mean_CPM", "PA_mean_CPM", "control_min_CPM", "control_max_CPM", "PA_min_CPM", "PA_max_CPM", "log2_CPM_ratio_plus1"]],
                  [np.mean(controls), np.mean(pa), min(controls), max(controls), min(pa), max(pa), ratio(np.mean(pa), np.mean(controls), 1)])
            check("within_low_cell_" + gene + "_" + subtype, row.minimum_cells == min_cells[subtypes.index(subtype)] and row.any_library_below50_cells == (min_cells[subtypes.index(subtype)] < 50))
        if candidates.at[gene, "symbol"] in ["CXCL8", "VNN1"]:
            spotlight[candidates.at[gene, "symbol"]] = {"by_library": per_lib, "contrasts": values}
    close("all_576_covered_source_counts", source_umi_actual, source_umi_expected)
    close("all_576_covered_source_CPM", source_cpm_actual, source_cpm_expected)
    close("all_96_covered_observed_CPM", observed_actual, observed_expected)
    close("all_96_covered_fixed_weight_CPM", fixed_actual, fixed_expected)
    check("all_96_observed_CPM_reconstructions", max(reconstruction_errors) < 1e-10, maximum_absolute_error=max(reconstruction_errors))
    check("24_covered_and_one_absent", len(covered) == 24 and len(missing) == 1)
    return {"covered_candidate_n": len(covered), "uncovered_gene_ids": missing, "libraries": libraries, "subtypes": subtypes, "minimum_cells": min_cells, "fixed_UMI_weights": fixed.tolist(), "reconstruction_max_abs_error_CPM": max(reconstruction_errors), "spotlight": spotlight,
            "boundary": "Started from existing all-gene UMI aggregation; did not reread large original single-cell MTX files or reassign cell annotations."}


def phh(candidates):
    spec = read_json(STAGE / "phh_protocol.json")
    freeze = read_json(STAGE / "phh_freeze.json")
    record = read_json(STAGE / "phh_summary.json")
    check("phh_protocol_freeze", digest(STAGE / "phh_protocol.json") == freeze["protocol_sha256"] == record["protocol_sha256"])
    check("phh_script_record_hash", digest(STAGE / "analyze_phh.py") == record["script_sha256"])
    check("phh_candidate_hash", digest(PREV / "priority_all.csv") == freeze["candidate_sha256"] == record["candidate_sha256"])
    raw_path = (STAGE / spec["input"]).resolve()
    check("phh_raw_hash", digest(raw_path) == freeze["raw_count_sha256"] == record["input_sha256"])
    manifest_path = STAGE / spec["manifest"]
    check("phh_manifest_execution_hash", digest(manifest_path) == record["manifest_sha256"])
    meta = read_csv(manifest_path, keep_default_na=False).set_index("matrix_column")
    soft = STAGE / "phh_design" / "GSE247407_family.soft.gz"
    USED.add(soft.resolve())
    source_records = {}
    current = None
    with gzip.open(soft, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("^SAMPLE = "):
                current = line.split(" = ", 1)[1]
                source_records[current] = {}
            elif current is not None and " = " in line:
                key, val = line.split(" = ", 1)
                source_records[current].setdefault(key, []).append(val)
    for column, row in meta.iterrows():
        source = source_records[row.gsm]
        check("phh_manifest_source_" + row.gsm,
              source["!Sample_title"] == [row.title]
              and "individual: " + row.donor_id in source["!Sample_characteristics_ch1"]
              and "treatment: " + row.condition in source["!Sample_characteristics_ch1"]
              and "Column in matrix: " + column in source["!Sample_description"])
    labels = meta.donor_id.drop_duplicates().tolist()
    clusters = {(label, condition): meta.index[(meta.donor_id == label) & (meta.condition == condition)].tolist() for label in labels for condition in ["Control", "FFA"]}
    check("phh_thirty_libraries_ten_clusters_five_labels", len(meta) == 30 and meta.index.is_unique and meta.gsm.is_unique and len(labels) == 5 and len(clusters) == 10 and all(len(v) == 3 for v in clusters.values()))
    raw_table = read_csv(raw_path, sep="\t", keep_default_na=False)
    check("phh_unique_gene_ids", raw_table.gene_id.is_unique)
    raw_table = raw_table.set_index("gene_id")
    raw = raw_table[meta.index].astype(float)
    check("phh_integer_nonnegative_finite_source", raw.shape == (58735, 30) and np.isfinite(raw.to_numpy()).all() and (raw >= 0).all().all() and (raw % 1 == 0).all().all())

    def independent_factors(matrix):
        x = matrix.to_numpy(dtype=float)
        reference = x[(x > 0).all(axis=1)]
        logged = np.log(reference)
        ratios = np.exp(logged - logged.mean(axis=1, keepdims=True))
        factors = np.median(ratios, axis=0)
        factors /= math.exp(float(np.log(factors).mean()))
        return pd.Series(factors, index=matrix.columns), len(reference)

    factors, ref_n = independent_factors(raw)
    totals = raw.sum(axis=0)
    norm = raw.div(factors, axis="columns")
    cpm = raw.div(totals, axis="columns") * 1000000
    qc = read_csv(STAGE / "phh_library_qc.csv").set_index("matrix_column").loc[meta.index]
    close("phh_thirty_size_factors", qc.median_ratio_factor, factors)
    close("phh_thirty_library_totals", qc.assigned_counts, totals)
    close("phh_thirty_detected_gene_counts", qc.detected_genes, (raw > 0).sum(axis=0))
    pool = pd.DataFrame({key: raw[columns].sum(axis=1) for key, columns in clusters.items()})
    pool_factors, pool_ref_n = independent_factors(pool)
    pool_norm = pool.div(pool_factors, axis="columns")
    pool_output = read_csv(STAGE / "phh_pooled_factors.csv").set_index(["donor_label", "condition"])
    close("phh_ten_pooled_factors", pool_output.loc[list(clusters), "pooled_factor"], pool_factors)
    close("phh_ten_pooled_totals", pool_output.loc[list(clusters), "pooled_total"], pool.sum(axis=0))
    check("phh_reference_gene_counts", ref_n == record["reference_genes_30_libraries"] and pool_ref_n == record["reference_genes_10_pooled_profiles"])
    check("phh_25_exact_gene_ids_matched", all(gene in raw.index for gene in candidates.index))
    sample_out = read_csv(STAGE / "phh_candidate_samples.csv").set_index(["gene_id", "sample"])
    effect_out = read_csv(STAGE / "phh_donor_effects.csv").set_index(["gene_id", "donor_label"])
    summary_out = read_csv(STAGE / "phh_candidate_summary.csv").set_index("gene_id")
    omission_out = read_csv(STAGE / "phh_label_omission.csv").set_index(["gene_id", "omitted_label"])
    check("phh_output_unique_keys", len(sample_out) == 750 and sample_out.index.is_unique and len(effect_out) == 125 and effect_out.index.is_unique and len(omission_out) == 125 and omission_out.index.is_unique)
    sample_actual, sample_expected = [], []
    all_effect_actual, all_effect_expected, spotlight = [], [], {}
    for gene in candidates.index:
        for sample in meta.index:
            row = sample_out.loc[(gene, sample)]
            sample_actual.extend([row.raw_count, row.normalized_count, row.CPM])
            sample_expected.extend([raw.at[gene, sample], norm.at[gene, sample], cpm.at[gene, sample]])
            check("phh_sample_labels_" + gene + "_" + sample, row.feature_covered and row.gsm == meta.at[sample, "gsm"] and row.donor_label == meta.at[sample, "donor_id"] and row.condition == meta.at[sample, "condition"] and row.reference_symbol == raw_table.at[gene, "gene_name"])
        gene_effects = {}
        for label in labels:
            actual = effect_out.loc[(gene, label)]
            values = {}
            for condition in ["Control", "FFA"]:
                columns = clusters[(label, condition)]
                nvalues, cvalues = norm.loc[gene, columns], cpm.loc[gene, columns]
                values.update({condition + "_mean_norm": float(nvalues.mean()), condition + "_min_norm": float(nvalues.min()), condition + "_max_norm": float(nvalues.max()), condition + "_mean_CPM": float(cvalues.mean()), condition + "_CPM_ge1_libraries": int((cvalues >= 1).sum()), condition + "_nonzero_libraries": int((raw.loc[gene, columns] > 0).sum())})
            main_effect = ratio(values["FFA_mean_norm"], values["Control_mean_norm"], 1)
            cpm_effect = ratio(values["FFA_mean_CPM"], values["Control_mean_CPM"], 1)
            pooled_effect = ratio(pool_norm.at[gene, (label, "FFA")], pool_norm.at[gene, (label, "Control")], 1)
            values.update({"log2_FFA_Control_mean_norm_plus1": main_effect, "log2_FFA_Control_mean_CPM_plus1": cpm_effect, "log2_FFA_Control_pooled_norm_plus1": pooled_effect})
            all_effect_actual.extend(actual[list(values)].to_numpy(dtype=float))
            all_effect_expected.extend(values.values())
            check("phh_effect_flags_" + gene + "_" + label, actual.feature_covered and actual.sign_matches_CPM == (np.sign(main_effect) == np.sign(cpm_effect)) and actual.sign_matches_pooled == (np.sign(main_effect) == np.sign(pooled_effect)) and actual.same_disease_direction_as_HepG2 == (np.sign(main_effect) == np.sign(candidates.at[gene, "log2_OA_vs_control"])))
            gene_effects[label] = {"primary": main_effect, "CPM": cpm_effect, "pooled": pooled_effect}
        effects = np.array([x["primary"] for x in gene_effects.values()])
        summary_values = {"label_n": 5, "median_log2_effect": float(np.median(effects)), "minimum_log2_effect": float(effects.min()), "maximum_log2_effect": float(effects.max()), "positive_labels": int((effects > 0).sum()), "negative_labels": int((effects < 0).sum()), "CPM_direction_agree_labels": sum(np.sign(x["primary"]) == np.sign(x["CPM"]) for x in gene_effects.values()), "pool_direction_agree_labels": sum(np.sign(x["primary"]) == np.sign(x["pooled"]) for x in gene_effects.values())}
        close("phh_summary_" + gene, summary_out.loc[gene, list(summary_values)], list(summary_values.values()))
        omit = {}
        for label in labels:
            omit[label] = float(np.median([x["primary"] for other, x in gene_effects.items() if other != label]))
            close("phh_omit_" + gene + "_" + label, [omission_out.at[(gene, label), "median_effect_other4"]], [omit[label]])
        if candidates.at[gene, "symbol"] in ["CXCL8", "VNN1"]:
            spotlight[candidates.at[gene, "symbol"]] = {"by_label": gene_effects, "summary": {k: int(v) if isinstance(v, np.integer) else v for k, v in summary_values.items()}, "median_omitting_label": omit}
    close("phh_750_sample_rows_2250_numbers", sample_actual, sample_expected)
    close("phh_125_effect_rows_1875_numbers", all_effect_actual, all_effect_expected)
    for filename, expected in record["artifact_hashes"].items():
        check("phh_artifact_hash_" + filename, digest(STAGE / filename) == expected)
    return {"status": "passed" if all(x["passed"] for x in CHECKS if x["name"].startswith("phh_")) else "needs_revision", "source_shape": list(raw.shape), "library_n": 30, "label_condition_profiles": 10, "author_label_n": 5, "labels": labels, "reference_gene_n_30": ref_n, "reference_gene_n_pooled10": pool_ref_n, "spotlight": spotlight,
            "boundary": "GEO SOFT sample-label mapping verified. Unique-person identity, F95/F95P relationship and replicate nature were not independently established; no claim of 30 independent people or independent sesamin validation."}


def main():
    candidates, bulk_result = bulk()
    organoid_result = organoid(candidates)
    phh_result = phh(candidates)
    final_summary = read_json(STAGE / "robustness_summary.json")
    check("final_robustness_script_hash", digest(STAGE / "check_robustness.py") == final_summary["script_sha256"])
    for filename, expected in final_summary["artifact_hashes"].items():
        check("final_robustness_artifact_" + filename, digest(STAGE / filename) == expected)
    review = {"executed_utc": datetime.now(timezone.utc).isoformat(), "status": "passed" if all(x["passed"] for x in CHECKS) else "needs_revision",
              "assessment": "Share with caveats: descriptive numerical identities checked, not biological validation.",
              "script_sha256": digest(Path(__file__)), "python": platform.python_version(), "executable": sys.executable,
              "versions": {"numpy": np.__version__, "pandas": pd.__version__}, "bulk": bulk_result, "organoid": organoid_result,
              "PHH": phh_result,
              "checks": CHECKS, "check_count": len(CHECKS), "failures": [x for x in CHECKS if not x["passed"]],
              "design_caveats": ["25 candidates and CXCL8/VNN1 were selected in prior analysis; this is post-selection follow-up, not held-out confirmation.",
                                 "Shared OA group couples the two bulk contrasts; no sample pairing or independence inferred from numeric suffixes.",
                                 "Each organoid condition has two technical libraries without established independent differentiation replication; there is no sesamin treatment.",
                                 "HB1 has only 4 or 9 cells in some libraries; fixed-weight values may be unstable.",
                                 "UMI weights depend on RNA content and capture as well as cell numbers; standardization is not pure cell-composition or causal adjustment.",
                                 "PHH source metadata supports five author labels with three libraries per label/condition, not 30 independent donors. F95/F95P identity relationship and replicate nature remain unresolved.",
                                 "PHH label-wise responses disagree in direction for both spotlight genes; positive median and normalization sensitivity do not establish universal disease response or sesamin effects.",
                                 "Pseudo-count ratios are regularized descriptive ratios, not fitted model coefficients, significance or effect probabilities."],
              "not_checked": ["New literature novelty claims", "Source clinical or dietary efficacy", "Plots", "Independent donor/culture provenance", "Original cell-level MTX aggregation in this review"],
              "execution_note": "Initial PHH-inclusive run computed checks but failed while formatting an input path outside next_steps; path inventory root corrected to project root and complete review rerun. No source or analysis outputs modified.",
              "source_sha256": {str(p.relative_to(STAGE.parent.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(USED)}}
    (HERE / "review.json").write_text(json.dumps(review, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": review["status"], "checks": len(CHECKS), "failures": review["failures"], "spotlight_PHH": phh_result["spotlight"]}))
    if review["failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
