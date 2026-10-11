"""Build publication-safe, deterministic display data from frozen analysis tables.

Run with the existing project scientific environment. This performs packaging and
reconciliation only; it does not fit models, recompute inferential tests, or select
candidates. Only DATA_DIR and AUDIT_DIR are written. Source artifacts are read only.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

APP = Path(__file__).resolve().parents[1]
WORKSPACE = APP.parents[1]
NEXT = WORKSPACE / "research_planning/hepatic_lipotoxicity_2026-10-07/next_steps"
STAGE = NEXT / "analysis_strengthening_15"
OLD = NEXT / "publication_package_14"
DATA_DIR = APP / "space/data"
AUDIT_DIR = APP / "audit"
DATA_DIR.mkdir(parents=True, exist_ok=True)
AUDIT_DIR.mkdir(parents=True, exist_ok=True)
INPUTS: dict[Path, str] = {}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> pd.DataFrame:
    INPUTS[path] = sha(path)
    return pd.read_csv(path)


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def dump(path: Path, obj, pretty=False):
    path.write_text(json.dumps(plain(obj), ensure_ascii=False, allow_nan=False,
                               indent=2 if pretty else None,
                               separators=None if pretty else (",", ":")),
                    encoding="utf-8")


def one(frame, **filters):
    mask = np.ones(len(frame), dtype=bool)
    for key, value in filters.items():
        mask &= frame[key].eq(value).to_numpy()
    selected = frame.loc[mask]
    assert len(selected) == 1, (filters, len(selected))
    return selected.iloc[0]


human = read(OLD / "data/human_candidates_all25.csv")
mouse = read(OLD / "data/mouse_candidates_all26.csv")
library = read(STAGE / "robustness/candidate_library_values.csv")
consensus = read(STAGE / "robustness/cross_species_configuration_consensus.csv")
ranges = read(STAGE / "robustness/configuration_effect_summary.csv")
rat = read(STAGE / "rat_validation/outputs/candidate_results_all25_two_normalizations.csv")
rat_exact = read(STAGE / "rat_validation/outputs/posthoc_exact_sensitivity_all25.csv")
rat_values = read(STAGE / "rat_validation/outputs/candidate_sample_values.csv")
sc = read(STAGE / "single_cell/library_category_candidates_all25.csv")
kitagawa = read(STAGE / "single_cell/kitagawa_summary_all25.csv")
cell_df = read(OLD / "figures/single_cell/cell_level_umap_qc_expression.csv.gz")

datasets = [
    {"id": "GSE248903", "name": "Human hepatocyte cell model", "species": "human",
     "role": "Candidate discovery", "design": "Nine HepG2 libraries: Control, OA model and sesamin treatment, three per group.",
     "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248903",
     "limitations": ["Candidates were selected using this dataset; stability is not independent validation.",
                     "Independent culture identities, matching of exposure time/media/vehicle and explicit source-paper linkage remain incompletely documented."]},
    {"id": "PRJCA021826", "name": "Mouse hepatocyte cell model", "species": "mouse",
     "role": "Cross-species directional comparison", "design": "Nine libraries: ctrl, HFD, sess; three per group. HFD is an original label for the cell fatty-acid model, not a whole-animal feeding experiment.",
     "url": "https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA021826",
     "aliases": ["CRA013797", "PRJDB40955", "DRP018055"],
     "limitations": ["Mirrored accessions represent one dataset, not multiple independent validations.",
                     "Animal/culture identities, timing and dose remain incompletely documented.",
                     "TPM comparisons are descriptive and the low-abundance flags must be retained."]},
    {"id": "GSE42693", "name": "Rat dietary sesamin study", "species": "rat",
     "role": "Independent intervention context", "design": "Healthy male rats received 0.2% dietary lignan for 15 days; primary comparison is five control versus five sesamin arrays.",
     "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42693",
     "publication": "https://doi.org/10.3177/jnsv.55.31",
     "limitations": ["Original study selected five of seven animals in each group after excluding the lowest and highest final body weights.",
                     "Whole healthy liver: there is no induced fatty-acid injury, so this is not a matching lipotoxicity replication.",
                     "Confidence intervals are nominal Welch 95% intervals conditional on the selected samples, not simultaneous intervals.",
                     "Welch and post hoc exact-label sensitivity have different multiple-testing conclusions: VNN1 q approximately 0.00257 versus 0.1984.",
                     "Exact-label results require exchangeability; random assignment was not verified. They are not randomized causal tests.",
                     "Two other lignans are separate groups and are not sesamin replicates."]},
    {"id": "GSE207889", "name": "Single-cell liver organoid reference", "species": "human",
     "role": "Cell localization and fatty-acid context", "design": "24,034 author-retained cells from four technical libraries, two control and two PA, from one differentiation.",
     "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207889",
     "limitations": ["There is no sesamin treatment in these libraries and no independent biological replication.",
                     "Author cell labels and saved UMAP coordinates from the full 82,467-cell atlas are reused; no new clustering or embedding was performed.",
                     "UMAP is not spatial transcriptomics, and distances or visual separation are not effect sizes or significance.",
                     "Detection and composition can vary with sampling, sequencing depth and technical library effects.",
                     "Measured zeros differ from missing features; AOC4P is unavailable and remains null."]}
]

genes = []
for _, h in human.iterrows():
    symbol = h.symbol
    mrows = mouse.loc[mouse.original_human_id.eq(h.gene_id)]
    primary_rows = mrows.loc[mrows.primary_comparison_eligible.eq(True)]
    assert len(primary_rows) <= 1
    m = primary_rows.iloc[0] if len(primary_rows) else None
    human_samples = library.loc[library.species.eq("human") & library.human_symbol.eq(symbol)]
    msamples = []
    if m is not None:
        for number in range(1010532, 1010541):
            sample = f"DRR{number}"
            msamples.append({"sample": sample, "group": "ctrl" if number < 1010535 else "HFD" if number < 1010538 else "sess", "tpm": m[f"{sample}_TPM"]})
    mouse_obj = {
        "status": "primary_one_to_one" if m is not None else "secondary_one_to_many" if len(mrows) == 2 else "no_comparable_orthologue",
        "primary": m is not None,
        "low_abundance": m.all_libraries_TPM_lt1 if m is not None else None,
        "model_log2": m.log2_HFD_vs_ctrl_pc0p1 if m is not None else None,
        "treatment_log2": m.log2_sess_vs_HFD_pc0p1 if m is not None else None,
        "means": {"control": m.mean_ctrl_TPM_provisional if m is not None else None,
                  "model": m.mean_HFD_TPM_provisional if m is not None else None,
                  "sesamin": m.mean_sess_TPM_provisional if m is not None else None},
        "samples": msamples,
        "orthologues": mrows.to_dict("records"),
    }
    stability = None
    if m is not None:
        c = one(consensus, human_symbol=symbol, aggregation="mean")
        detail = ranges.loc[ranges.species.eq("mouse") & ranges.metric.eq("TPM") & ranges.aggregation.eq("mean") & ranges.human_symbol.eq(symbol)]
        r = one(detail, contrast="treatment_vs_model")
        stability = {"both_n": int(c.both_model_and_treatment_concordant_n), "total": int(c.configuration_pairs),
                     "share": c.both_concordant_configuration_share,
                     "mouse_treatment_range_percent": [r.min_percent_difference, r.max_percent_difference],
                     "details": detail.to_dict("records")}
    rat_obj = {"status": one(rat, human_id=h.gene_id, method="deposited_log2").status}
    for alias, method in [("deposited", "deposited_log2"), ("quantile", "quantile_sensitivity")]:
        r = one(rat, human_id=h.gene_id, method=method)
        e = one(rat_exact, human_id=h.gene_id, method=method)
        points = rat_values.loc[rat_values.human_id.eq(h.gene_id) & rat_values.method.eq(method) & rat_values.primary_comparison.eq(True)]
        rat_obj[alias] = {"effect": r.effect_log2_signal, "ci_low": r.ci95_low, "ci_high": r.ci95_high,
                          "p_welch": r.p_value, "q_welch": r.q_bh_fixed25, "p_exact": e.permutation_p,
                          "q_exact": e.permutation_q_bh25, "samples": [
                              {"sample": p["sample"], "group": p["group"], "value": p["aggregate_log2_signal"]}
                              for p in points.to_dict("records")],
                          "rat_symbol": r.rat_symbol, "probe_ids": r.probe_ids,
                          "fixed_bh_family_size": 25, "exact_allocations": 252 if r.status == "mapped" else None}
    srows = sc.loc[sc.gene_id.eq(h.gene_id)]
    assert len(srows) == 16
    sc_obj = {"covered": bool(srows.covered.iloc[0]), "libraries": [
        {"library": s["gsm"], "condition": s["condition"], "category": s["celltype_category"],
         "n_cells": s["n_cells"], "n_detected": s["n_detected_cells"], "detection_pct": s["detected_pct"], "mean_log1p": s["mean_log1p_cp10k"]}
        for s in srows.to_dict("records")], "kitagawa": {}}
    for alias, scheme in [("equal_library", "equal_library_actual"), ("pooled", "pooled_cell_actual")]:
        k = one(kitagawa, gene_id=h.gene_id, scheme=scheme)
        sc_obj["kitagawa"][alias] = {"total": k.total_PA_minus_control, "composition": k.composition_total,
                                     "within": k.within_category_total, "control_mean": k.control_profile_mean, "pa_mean": k.PA_profile_mean}
    genes.append({"id": h.gene_id, "symbol": symbol,
                  "human": {"model_log2": h.log2_OA_vs_control, "treatment_log2": h.log2_sesamin_vs_OA,
                            "treatment_vs_control_log2": h.log2_sesamin_vs_control,
                            "means": {"control": h.mean_norm_control, "model": h.mean_norm_OA, "sesamin": h.mean_norm_sesamin_OA},
                            "samples": [{"sample": s["library"], "group": s["group"], "value": s["median_ratio_count"]} for s in human_samples.to_dict("records")]},
                  "mouse": mouse_obj, "stability": stability, "rat": rat_obj, "single_cell": sc_obj})

pca = {}
for species, n_features in [("human", 13014), ("mouse", 15340)]:
    scores = read(STAGE / f"robustness/{species}_pca_scores.csv")
    variance = read(STAGE / f"robustness/{species}_pca_variance.csv")
    pca[species] = {"samples": [{"sample": r["library"], "group": r["group"], "x": r["PC1"], "y": r["PC2"]} for r in scores.to_dict("records")],
                    "variance": variance.explained_variance_ratio.iloc[:2].tolist(), "feature_n": n_features}
scores = read(STAGE / "rat_validation/outputs/all20_array_pca.csv")
variance = read(STAGE / "rat_validation/outputs/pca_variance.csv")
pca["rat"] = {"samples": [{"sample": r["sample"], "group": r["group"], "x": r["PC1"], "y": r["PC2"]} for r in scores.loc[scores.method.eq("deposited_log2")].to_dict("records")],
              "variance": [variance.iloc[0].pc1, variance.iloc[0].pc2], "feature_n": 31099}

download_specs = [
    (OLD / "data/human_candidates_all25.csv", "Human candidates: all 25"),
    (OLD / "data/mouse_candidates_all26.csv", "Mouse mapping and expression: all candidates, including secondary/NA mappings"),
    (STAGE / "rat_validation/outputs/candidate_results_all25_two_normalizations.csv", "Rat results: all 25 candidates and both normalizations"),
    (STAGE / "rat_validation/outputs/posthoc_exact_sensitivity_all25.csv", "Rat post hoc exact-label sensitivity: all 25 candidates"),
    (STAGE / "rat_validation/outputs/candidate_sample_values.csv", "Rat candidate array values: all 20 arrays and both normalizations"),
    (STAGE / "robustness/cross_species_configuration_consensus.csv", "Cross-species stability: all 15 primary candidates, means and medians"),
    (STAGE / "robustness/configuration_effect_summary.csv", "All configuration sensitivity summaries"),
    (STAGE / "robustness/candidate_library_values.csv", "Human and mouse values for all 15 primary candidates"),
    (STAGE / "single_cell/library_category_candidates_all25.csv", "Single-cell broad categories: all 25 candidates and four libraries"),
    (STAGE / "single_cell/library_candidates_all25.csv", "Single-cell library summaries: all 25 candidates and four libraries"),
    (STAGE / "single_cell/kitagawa_summary_all25.csv", "Single-cell composition decomposition: all candidates and weighting schemes"),
    (STAGE / "single_cell/weighting_sensitivity_all25.csv", "Single-cell weighting sensitivity: all 25 candidates"),
]
downloads = []
for source, label in download_specs:
    frame = read(source)
    destination = DATA_DIR / source.name
    shutil.copyfile(source, destination)
    assert sha(source) == sha(destination)
    downloads.append({"label": label, "path": f"data/{destination.name}", "rows": len(frame), "sha256": sha(destination)})

round_errors = {}


def rounded(series: pd.Series, name: str):
    arr = series.to_numpy(dtype=float)
    rounded_arr = np.round(arr, 6)
    finite = np.isfinite(arr)
    round_errors[name] = float(np.max(np.abs(rounded_arr[finite] - arr[finite]))) if finite.any() else None
    return plain(rounded_arr)


cells = {"schema_version": "1.0.0", "dataset": "GSE207889", "n_cells": len(cell_df),
         "coordinates": {"x": rounded(cell_df.umap_1, "umap_x"), "y": rounded(cell_df.umap_2, "umap_y")},
         "categories": {}, "codes": {}, "genes": human.symbol.tolist(), "expression": {},
         "metadata": {"embedding": "Author saved full-atlas UMAP; four-library subset; not recomputed and not spatial coordinates",
                      "expression_unit": "Natural log1p(10000 * gene UMI / all 33538-gene UMI total per cell)",
                      "precision_decimal_places": 6, "row_identity": "Anonymous ordinal array index only; original barcodes omitted",
                      "missing_feature": "AOC4P is null for every cell; null is not measured zero"}}
for alias, column in [("library", "gsm"), ("condition", "condition"), ("subtype", "celltype_scTypedb"), ("broad", "celltype_category")]:
    labels = sorted(cell_df[column].unique().tolist())
    cells["categories"][alias] = labels
    cells["codes"][alias] = cell_df[column].map({label: i for i, label in enumerate(labels)}).tolist()
for symbol in human.symbol:
    cells["expression"][symbol] = rounded(cell_df[f"{symbol}__log1p_cp10k"], symbol)

evidence = {"schema_version": "1.0.0", "metadata": {
    "title": "Sesamin Evidence Explorer", "data_version": "2026-10-09-analysis-strengthening-15",
    "purpose": "Explore fixed research candidates across study contexts; not a treatment or product efficacy assessment.",
    "datasets": datasets,
    "units": {"human": "Log2((mean median-ratio normalized count + 1) ratio); sample dots use median-ratio normalized counts.",
              "mouse": "Log2((mean TPM + 0.1) ratio); individual dots and means are TPM.",
              "rat": "Difference in mean log2(deposited GCOS signal + 1); quantile sensitivity normalizes all 20 log arrays before comparison.",
              "sc": "Natural log1p(CP10k) per cell, including observed zeros; detection is percent with raw UMI > 0.",
              "stability": "Descriptive share among 64 x 64 = 4096 cross-species configuration pairs; not a probability or confidence level.",
              "kitagawa": "PA minus control mean per-cell log1p(CP10k), split exactly into broad-category composition and within-category terms; not causal shares."},
    "limitations": ["The 25 candidates and 15 primary one-to-one comparisons are fixed; absent mappings and low expression remain visible.",
                    "No new experiments, therapeutic benefit, new target validation or quantitatively benchmarked AI improvement is demonstrated.",
                    "VNN1 was previously reported using the human dataset; this app does not claim it is a newly discovered target.",
                    "Configuration ranges are not confidence intervals. Single-cell cells or technical libraries are not independent biological replicates.",
                    "Composition decomposition is descriptive: within-category terms can include finer subtype, sequencing-depth and technical differences.",
                    "Effect units, species, tissue contexts and exposure designs differ; magnitudes should not be pooled across studies."],
    "fixed_candidates": 25, "primary_mouse_candidates": 15, "single_cell_count": 24034,
    "inference": {"rat_primary": "Exploratory two-sided Welch test, nominal 95% CI and BH over fixed 25-candidate family.",
                  "rat_sensitivity": "Explicitly post hoc enumeration of all 252 five-versus-five label allocations, same fixed BH family; missing candidates padded with p=1 only in adjustment and displayed as null.",
                  "other_analyses": "Descriptive; no biological significance tests."}},
    "genes": genes, "pca": pca, "downloads": downloads,
    "provenance": [{"source_file": path.name, "sha256": digest} for path, digest in sorted(INPUTS.items(), key=lambda x: str(x[0]))]}

evidence = plain(evidence)
cells = plain(cells)
dump(DATA_DIR / "evidence.json", evidence)
dump(DATA_DIR / "cells.json", cells)

# Validate packaging against independent invariants and reconciled source values.
checks = {}
checks["all25_unique_candidates"] = len(genes) == 25 and len(set(g["id"] for g in genes)) == 25
checks["all15_primary_mouse"] = sum(g["mouse"]["primary"] for g in genes) == 15
checks["all15_rat_mapped"] = sum(g["rat"]["status"] == "mapped" for g in genes) == 15
checks["all24034_cells"] = cells["n_cells"] == 24034
checks["all_columnar_lengths"] = all(len(v) == 24034 for section in ["coordinates", "codes", "expression"] for v in cells[section].values())
checks["all24_measured_genes"] = sum(any(v is not None for v in arr) for arr in cells["expression"].values()) == 24
checks["AOC4P_missing_not_zero"] = all(v is None for v in cells["expression"]["AOC4P"])
checks["rounding_max_error_le_half_unit"] = all(v is None or v <= 5.00000001e-7 for v in round_errors.values())
checks["all400_library_category_rows"] = sum(len(g["single_cell"]["libraries"]) for g in genes) == 400
checks["all_sc_counts_reconcile"] = all(sum(row["n_cells"] for row in g["single_cell"]["libraries"]) == 24034 for g in genes)
checks["all_rat_primary_dots_5_plus5"] = all(len(g["rat"][method]["samples"]) == (10 if g["rat"]["status"] == "mapped" else 0) and all(sum(s["group"] == group for s in g["rat"][method]["samples"]) == (5 if g["rat"]["status"] == "mapped" else 0) for group in ["control", "sesamin"]) for g in genes for method in ["deposited", "quantile"])
lookup = {g["symbol"]: g for g in evidence["genes"]}
checks["VNN1_model_direction_contrast"] = lookup["VNN1"]["human"]["treatment_log2"] < 0 and lookup["VNN1"]["mouse"]["treatment_log2"] < 0 and lookup["VNN1"]["rat"]["deposited"]["effect"] > 0
checks["VNN1_4096_of4096"] = lookup["VNN1"]["stability"]["both_n"] == lookup["VNN1"]["stability"]["total"] == 4096
checks["AREG_3072_of4096"] = lookup["AREG"]["stability"]["both_n"] == 3072
checks["JAG2_2560_of4096_low"] = lookup["JAG2"]["stability"]["both_n"] == 2560 and lookup["JAG2"]["mouse"]["low_abundance"] is True
checks["VNN1_Welch_vs_exact_kept"] = abs(lookup["VNN1"]["rat"]["deposited"]["q_welch"] - .002569) < .000001 and abs(lookup["VNN1"]["rat"]["deposited"]["q_exact"] - 25/126) < 1e-10
checks["CXCL8_rat_NA_preserved"] = lookup["CXCL8"]["rat"]["deposited"]["effect"] is None
checks["ACAT2_secondary_two_mappings"] = len(lookup["ACAT2"]["mouse"]["orthologues"]) == 2 and not lookup["ACAT2"]["mouse"]["primary"]
checks["all_primary_samples_means_reconcile"] = all(abs(np.mean([s["tpm"] for s in g["mouse"]["samples"] if s["group"] == group]) - g["mouse"]["means"][alias]) < 1e-8 for g in evidence["genes"] if g["mouse"]["primary"] for group, alias in [("ctrl", "control"), ("HFD", "model"), ("sess", "sesamin")])
checks["kitagawa_identities_preserved"] = all(abs(k["composition"] + k["within"] - k["total"]) < 1e-12 for g in evidence["genes"] if g["single_cell"]["covered"] for k in g["single_cell"]["kitagawa"].values())
sc_reconciliation = []
for g in evidence["genes"]:
    if not g["single_cell"]["covered"]:
        continue
    values = np.asarray(cells["expression"][g["symbol"]])
    for row in g["single_cell"]["libraries"]:
        lib_code = cells["categories"]["library"].index(row["library"])
        cat_code = cells["categories"]["broad"].index(row["category"])
        mask = (np.asarray(cells["codes"]["library"]) == lib_code) & (np.asarray(cells["codes"]["broad"]) == cat_code)
        sc_reconciliation.append(abs(values[mask].mean() - row["mean_log1p"]) < 5.01e-7 and int(np.count_nonzero(values[mask] > 0)) == int(row["n_detected"]))
checks["384_rounded_sc_aggregate_reconciliations"] = len(sc_reconciliation) == 384 and all(sc_reconciliation)
checks["PCA_all_samples_retained"] = [len(pca[s]["samples"]) for s in ["human", "mouse", "rat"]] == [9, 9, 20]
checks["source_hashes_unchanged"] = all(sha(path) == digest for path, digest in INPUTS.items())
payload_files = [DATA_DIR / "evidence.json", DATA_DIR / "cells.json"] + [DATA_DIR / Path(d["path"]).name for d in downloads]
leak_regex = re.compile(r"(?:(?<![A-Za-z])[A-Za-z]:[\\/]|OneDrive|Users[\\/]|hf_[A-Za-z0-9]{20,}|sk-[A-Za-z0-9]{20,})")
leaks = {p.name: bool(leak_regex.search(p.read_text(encoding="utf-8-sig"))) for p in payload_files}
checks["no_local_path_or_token_patterns"] = not any(leaks.values())
audit = {"status": "PASS" if all(checks.values()) else "FAIL", "operation": "Packaging and reconciliation only; no new research analysis", "checks": checks,
         "row_counts": {"genes": len(genes), "cells": len(cell_df), "sc_library_category": len(sc), "rat_normalization_results": len(rat), "stability_mean_median": len(consensus)},
         "rounding_max_abs_error": round_errors, "source_hashes": evidence["provenance"],
         "outputs": [{"file": p.name, "bytes": p.stat().st_size, "sha256": sha(p)} for p in payload_files]}
dump(AUDIT_DIR / "data_build_checks.json", audit, pretty=True)
assert all(checks.values()), {key: value for key, value in checks.items() if not value}
print(json.dumps({"status": audit["status"], "checks": len(checks), "files": len(payload_files), "evidence_bytes": (DATA_DIR / "evidence.json").stat().st_size, "cells_bytes": (DATA_DIR / "cells.json").stat().st_size}))
