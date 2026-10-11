"""Independent, read-only review of the public stage-15 payload.

Uses no builder imports. Recomputes display conversions from retained counts,
per-array values, and configuration effects. Writes only the adjacent review
JSON/Markdown. This is a separate Codex calculation review, not expert review.
"""
from pathlib import Path
from itertools import combinations
from datetime import datetime, timezone
import hashlib
import json
import math
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import t

APP = Path(__file__).resolve().parents[1]
BASE = APP.parents[1]
NEXT = BASE / "research_planning/hepatic_lipotoxicity_2026-10-07/next_steps"
OLD = NEXT / "publication_package_14"
STAGE = NEXT / "analysis_strengthening_15"
DATA = APP / "space/data"
checks = []
sources = {}
notes = []


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path):
    sources[str(path.relative_to(BASE))] = digest(path)
    return pd.read_csv(path)


def check(name, condition, details=None):
    checks.append({"name": name, "pass": bool(condition), "details": details})


def close(a, b, tol=1e-9):
    if a is None or (isinstance(a, float) and math.isnan(a)):
        return b is None or (isinstance(b, float) and math.isnan(b))
    if b is None:
        return False
    return math.isclose(float(a), float(b), rel_tol=1e-8, abs_tol=tol)


def bh(values):
    # Independent rank-loop implementation; missing candidates supplied as 1.
    indexed = sorted(enumerate(values), key=lambda pair: pair[1])
    result = [None] * len(values)
    previous = 1.0
    for rank in range(len(values), 0, -1):
        i, p = indexed[rank - 1]
        previous = min(previous, p * len(values) / rank)
        result[i] = previous
    return result


payload_hashes = {p.name: digest(p) for p in sorted(DATA.iterdir()) if p.is_file()}
evidence = json.loads((DATA / "evidence.json").read_text(encoding="utf-8"))
cells = json.loads((DATA / "cells.json").read_text(encoding="utf-8"))
human = table(OLD / "data/human_candidates_all25.csv")
mouse = table(OLD / "data/mouse_candidates_all26.csv")
original_cells = table(OLD / "figures/single_cell/cell_level_umap_qc_expression.csv.gz")
libraries = table(STAGE / "robustness/candidate_library_values.csv")
rat_samples = table(STAGE / "rat_validation/outputs/candidate_sample_values.csv")
rat_results = table(STAGE / "rat_validation/outputs/candidate_results_all25_two_normalizations.csv")
config = table(STAGE / "robustness/all_configuration_effects.csv.gz")
gene_map = {g["symbol"]: g for g in evidence["genes"]}
check("fixed_25_identity_order_and_uniqueness", len(evidence["genes"]) == 25
      and [g["id"] for g in evidence["genes"]] == human.gene_id.tolist()
      and list(gene_map) == human.symbol.tolist() and len(gene_map) == 25)
check("primary_mouse_15_with_secondary_ACAT2_separate",
      sum(g["mouse"]["primary"] for g in gene_map.values()) == 15
      and not gene_map["ACAT2"]["mouse"]["primary"]
      and gene_map["ACAT2"]["mouse"]["status"] == "secondary_one_to_many"
      and len(gene_map["ACAT2"]["mouse"]["orthologues"]) == 2)

human_missing_dots = []
for row in human.to_dict("records"):
    name = row["symbol"]
    g = gene_map[name]
    h = g["human"]
    means = [row[k] for k in ["mean_norm_control", "mean_norm_OA", "mean_norm_sesamin_OA"]]
    recomputed = [math.log2((means[1] + 1) / (means[0] + 1)),
                  math.log2((means[2] + 1) / (means[1] + 1)),
                  math.log2((means[2] + 1) / (means[0] + 1))]
    check(f"human_three_effects_and_means_{name}",
          all(close(h[k], v) for k, v in zip(["model_log2", "treatment_log2", "treatment_vs_control_log2"], recomputed))
          and all(close(h["means"][k], v) for k, v in zip(["control", "model", "sesamin"], means)))
    if not h["samples"]:
        human_missing_dots.append(name)
    else:
        ref = libraries.query("species == 'human' and human_symbol == @name").set_index("library")
        check(f"human_sample_labels_values_{name}", len(h["samples"]) == 9 and all(
            s["sample"] in ref.index and s["group"] == ref.loc[s["sample"], "group"]
            and close(s["value"], ref.loc[s["sample"], "median_ratio_count"]) for s in h["samples"]))
    mr = mouse[mouse.original_human_id.eq(row["gene_id"])]
    pr = mr[mr.primary_comparison_eligible.eq(True)]
    mo = g["mouse"]
    check(f"mouse_mapping_{name}", len(mo["orthologues"]) == len(mr)
          and sorted(str(x["mouse_gene_id"]) for x in mo["orthologues"]) == sorted(
              "None" if pd.isna(x) else str(x) for x in mr.mouse_gene_id)
          and bool(len(pr)) == mo["primary"])
    if not len(pr):
        check(f"mouse_nonprimary_no_false_effects_{name}", all(mo[x] is None for x in ["model_log2", "treatment_log2", "low_abundance"])
              and not mo["samples"] and g["stability"] is None)
        continue
    r = pr.iloc[0]
    ordered = [r[f"DRR{i}_TPM"] for i in range(1010532, 1010541)]
    c, m, s = [float(np.mean(ordered[start:start + 3])) for start in [0, 3, 6]]
    check(f"mouse_sample_groups_effects_units_{name}",
          len(mo["samples"]) == 9 and all(close(x["tpm"], y) for x, y in zip(mo["samples"], ordered))
          and [x["group"] for x in mo["samples"]] == ["ctrl"] * 3 + ["HFD"] * 3 + ["sess"] * 3
          and [x["sample"] for x in mo["samples"]] == [f"DRR{i}" for i in range(1010532, 1010541)]
          and all(close(mo["means"][k], v) for k, v in zip(["control", "model", "sesamin"], [c, m, s]))
          and close(mo["model_log2"], math.log2((m + .1) / (c + .1)))
          and close(mo["treatment_log2"], math.log2((s + .1) / (m + .1)))
          and mo["low_abundance"] == all(v < 1 for v in ordered))
    part = config[config.human_symbol.eq(name) & config.aggregation.eq("mean")]
    by_species = {}
    for species, metric in [("human", "median_ratio_reestimated"), ("mouse", "TPM")]:
        cs = part[part.species.eq(species) & part.metric.eq(metric)]
        by_species[species] = cs.pivot(index="configuration", columns="contrast", values="log2_ratio").sort_index()
    hg, mg = by_species["human"], by_species["mouse"]
    n = 0
    for _, hv in hg.iterrows():
        for _, mv in mg.iterrows():
            n += all(np.sign(hv[k]) == np.sign(mv[k]) and np.sign(hv[k]) != 0
                     for k in ["model_vs_control", "treatment_vs_model"])
    rng = part[part.species.eq("mouse") & part.metric.eq("TPM") & part.contrast.eq("treatment_vs_model")].percent_difference_no_pseudocount
    st = g["stability"]
    check(f"independent_cartesian_configuration_count_{name}", len(hg) == len(mg) == 64
          and st["total"] == 4096 and st["both_n"] == n and close(st["share"], n / 4096)
          and all(close(a, b) for a, b in zip(st["mouse_treatment_range_percent"], [rng.min(), rng.max()])))

# Rat inference is recalculated from all frozen individual observations, not
# simply compared to the inferential-result table used by the builder.
rat_recomputed = {}
for alias, method in [("deposited", "deposited_log2"), ("quantile", "quantile_sensitivity")]:
    p_values, exact_values = [], []
    results = []
    for name, g in gene_map.items():
        r = g["rat"][alias]
        original = rat_samples[rat_samples.human_symbol.eq(name) & rat_samples.method.eq(method)
                               & rat_samples.primary_comparison.eq(True)]
        if g["rat"]["status"] != "mapped":
            check(f"rat_missing_remains_null_{alias}_{name}", len(original) == 0 and len(r["samples"]) == 0
                  and all(r[k] is None for k in ["effect", "ci_low", "ci_high", "p_welch", "q_welch", "p_exact", "q_exact"]))
            p_values.append(1.0); exact_values.append(1.0); results.append(None)
            continue
        ref = original.set_index("sample")
        check(f"rat_sample_join_{alias}_{name}", len(ref) == len(r["samples"]) == 10
              and ref.group.value_counts().to_dict() == {"control": 5, "sesamin": 5}
              and all(s["sample"] in ref.index and s["group"] == ref.loc[s["sample"], "group"]
                      and close(s["value"], ref.loc[s["sample"], "aggregate_log2_signal"]) for s in r["samples"]))
        c = original[original.group.eq("control")].aggregate_log2_signal.to_numpy()
        s = original[original.group.eq("sesamin")].aggregate_log2_signal.to_numpy()
        effect = s.mean() - c.mean()
        terms = [np.var(v, ddof=1) / len(v) for v in [c, s]]
        se = math.sqrt(sum(terms)); df = sum(terms) ** 2 / sum(x * x / 4 for x in terms)
        pw = 2 * t.sf(abs(effect / se), df)
        ci = [effect - t.ppf(.975, df) * se, effect + t.ppf(.975, df) * se]
        pooled = list(c) + list(s); total = sum(pooled)
        exceedances = sum(abs(sum(pooled[j] for j in allocation) / 5 -
                              (total - sum(pooled[j] for j in allocation)) / 5) >= abs(effect) - 1e-12
                          for allocation in combinations(range(10), 5))
        pe = exceedances / 252
        p_values.append(float(pw)); exact_values.append(float(pe))
        results.append({"symbol": name, "effect": float(effect), "p_welch": float(pw), "p_exact": pe})
        check(f"rat_independent_effect_CI_tests_{alias}_{name}", close(r["effect"], effect)
              and close(r["ci_low"], ci[0]) and close(r["ci_high"], ci[1])
              and close(r["p_welch"], pw) and close(r["p_exact"], pe)
              and r["fixed_bh_family_size"] == 25 and r["exact_allocations"] == 252)
    qw, qe = bh(p_values), bh(exact_values)
    for i, (name, g) in enumerate(gene_map.items()):
        if results[i] is not None:
            check(f"rat_independent_BH25_{alias}_{name}", close(g["rat"][alias]["q_welch"], qw[i])
                  and close(g["rat"][alias]["q_exact"], qe[i]))
            results[i].update(q_welch=qw[i], q_exact=qe[i])
    rat_recomputed[alias] = [r for r in results if r]
check("VNN1_opposite_direction_and_method_dependent_FDR_explicit",
      gene_map["VNN1"]["human"]["treatment_log2"] < 0 and gene_map["VNN1"]["mouse"]["treatment_log2"] < 0
      and all(gene_map["VNN1"]["rat"][k]["effect"] > 0 and gene_map["VNN1"]["rat"][k]["q_welch"] < .05
              and close(gene_map["VNN1"]["rat"][k]["q_exact"], 25 / 126) for k in ["deposited", "quantile"]))

check("cell_24034_retained_rows_complete", len(original_cells) == cells["n_cells"] == 24034
      and original_cells.cell_key.nunique() == 24034)
max_errors = {}
for source, key in [("umap_1", "x"), ("umap_2", "y")]:
    target = np.asarray(cells["coordinates"][key], dtype=float)
    delta = float(np.max(np.abs(target - original_cells[source].to_numpy())))
    max_errors[source] = delta
    check(f"all24034_author_UMAP_{key}", len(target) == 24034 and delta <= 5.01e-7)
decoded = {}
for key, column in [("library", "gsm"), ("condition", "condition"), ("subtype", "celltype_scTypedb"), ("broad", "celltype_category")]:
    labels = cells["categories"][key]
    codes = cells["codes"][key]
    decoded[key] = [labels[i] for i in codes]
    check(f"all24034_labels_{key}", len(codes) == 24034 and decoded[key] == original_cells[column].tolist()
          and set(codes) == set(range(len(labels))))
check("cell_fixed_gene_order", cells["genes"] == human.symbol.tolist()
      and list(cells["expression"]) == human.symbol.tolist())
measured = 0
for name, g in gene_map.items():
    exp = cells["expression"][name]
    if name == "AOC4P":
        check("AOC4P_not_measured_is_null_never_zero", len(exp) == 24034 and all(v is None for v in exp)
              and not g["single_cell"]["covered"] and original_cells[f"{name}__counts"].isna().all()
              and all(r["mean_log1p"] is None and r["detection_pct"] is None for r in g["single_cell"]["libraries"]))
        continue
    measured += 1
    original_counts = original_cells[f"{name}__counts"].to_numpy(dtype=float)
    recomputed = np.log1p(10000 * original_counts / original_cells.raw_total_counts.to_numpy())
    array = np.asarray(exp, dtype=float)
    error = float(np.max(np.abs(array - recomputed)))
    max_errors[name] = error
    check(f"all24034_count_to_expression_{name}", len(exp) == 24034 and error <= 5.01e-7
          and np.isfinite(array).all() and np.array_equal(array > 0, original_counts > 0))
    observed = original_cells[["gsm", "condition", "celltype_category"]].copy()
    observed["value"] = recomputed
    observed["detected"] = original_counts > 0
    groups = observed.groupby(["gsm", "condition", "celltype_category"], observed=True).agg(
        n=("value", "size"), mean=("value", "mean"), detected=("detected", "sum"))
    sc = g["single_cell"]
    check(f"independent_16_library_category_aggregates_{name}", len(sc["libraries"]) == 16 and all(
        (r["library"], r["condition"], r["category"]) in groups.index
        and r["n_cells"] == groups.loc[(r["library"], r["condition"], r["category"]), "n"]
        and r["n_detected"] == groups.loc[(r["library"], r["condition"], r["category"]), "detected"]
        and close(r["mean_log1p"], groups.loc[(r["library"], r["condition"], r["category"]), "mean"])
        and close(r["detection_pct"], 100 * r["n_detected"] / r["n_cells"])
        for r in sc["libraries"]))
    for alias in ["equal_library", "pooled"]:
        profiles = {}
        direct = {}
        for condition in ["CTRL", "PA"]:
            obs = observed[observed.condition.eq(condition)]
            if alias == "pooled":
                profile = obs.groupby("celltype_category", observed=True).agg(n=("value", "size"), mu=("value", "mean"))
                profile["p"] = profile.n / len(obs)
                direct[condition] = obs.value.mean()
            else:
                items = []
                for lib, values in obs.groupby("gsm"):
                    part = values.groupby("celltype_category", observed=True).agg(n=("value", "size"), mu=("value", "mean"))
                    part["p"] = part.n / len(values)
                    part["pmu"] = part.p * part.mu
                    items.append(part)
                joined = pd.concat(items)
                profile = joined.groupby(level=0).agg(p=("p", "mean"), pmu=("pmu", "mean"))
                profile["mu"] = profile.pmu / profile.p
                direct[condition] = obs.groupby("gsm").value.mean().mean()
            profiles[condition] = profile
        control, pa = profiles["CTRL"], profiles["PA"]
        comp = ((pa.p - control.p) * (pa.mu + control.mu) / 2).sum()
        within = ((pa.mu - control.mu) * (pa.p + control.p) / 2).sum()
        k = sc["kitagawa"][alias]
        check(f"independent_cell_to_Kitagawa_{alias}_{name}", close(k["total"], direct["PA"] - direct["CTRL"])
              and close(k["composition"], comp) and close(k["within"], within)
              and close(k["control_mean"], direct["CTRL"]) and close(k["pa_mean"], direct["PA"])
              and close(comp + within, direct["PA"] - direct["CTRL"]))
check("exactly24_measured_features", measured == 24)

for species in ["human", "mouse", "rat"]:
    p = evidence["pca"][species]
    if species == "rat":
        reference = table(STAGE / "rat_validation/outputs/all20_array_pca.csv").query("method == 'deposited_log2'")
        var = table(STAGE / "rat_validation/outputs/pca_variance.csv").iloc[0]
        variance = [var.pc1, var.pc2]; id_col = "sample"; expected_n = 20; features = 31099
    else:
        reference = table(STAGE / f"robustness/{species}_pca_scores.csv")
        var = table(STAGE / f"robustness/{species}_pca_variance.csv")
        variance = var.explained_variance_ratio.iloc[:2].tolist()
        id_col = "library"; expected_n = 9; features = 13014 if species == "human" else 15340
    indexed = reference.set_index(id_col)
    check(f"PCA_retained_scores_labels_variance_{species}", len(p["samples"]) == expected_n
          and p["feature_n"] == features and all(close(a, b) for a, b in zip(p["variance"], variance))
          and all(s["sample"] in indexed.index and s["group"] == indexed.loc[s["sample"], "group"]
                  and close(s["x"], indexed.loc[s["sample"], "PC1"])
                  and close(s["y"], indexed.loc[s["sample"], "PC2"]) for s in p["samples"]))

locations = [OLD / "data", STAGE / "robustness", STAGE / "single_cell", STAGE / "rat_validation/outputs"]
for dl in evidence["downloads"]:
    target = APP / "space" / dl["path"]
    candidates = [directory / target.name for directory in locations if (directory / target.name).is_file()]
    check(f"download_original_bytes_sha_{target.name}", len(candidates) == 1 and target.read_bytes() == candidates[0].read_bytes()
          and digest(target) == dl["sha256"] and len(pd.read_csv(target)) == dl["rows"])
    if len(candidates) == 1:
        sources[str(candidates[0].relative_to(BASE))] = digest(candidates[0])

all_public_text = "\n".join(p.read_text(encoding="utf-8-sig") for p in sorted(DATA.iterdir()) if p.is_file())
leak_patterns = {
    "absolute_windows_path": r"\b[A-Za-z]:[\\/]",
    "home_path": r"(?:OneDrive|/home/|/Users/|Users\\)",
    "credential_like": r"(?:hf_[A-Za-z0-9]{20,}|sk-[A-Za-z0-9]{20,}|Bearer\s+[A-Za-z0-9._-]{20,})",
    "10x_barcode": r"\b[ACGT]{16}-[0-9]+\b",
    "barcode_identity_field": r'"(?:cell_key|barcode|original_obs_index|author_obs_row|raw_barcode_column_0based)"',
}
for label, pattern in leak_patterns.items():
    check(f"no_public_{label}", re.search(pattern, all_public_text) is None)
check("source_hashes_unchanged_by_review", all(digest(BASE / path) == value for path, value in sources.items()))
check("payload_hashes_stable_during_review", all(digest(DATA / name) == value for name, value in payload_hashes.items()))

notes.append({"severity": "nonblocking", "issue": "Ten human candidates have aggregate expression/effects but no per-library dots in the public payload.",
              "genes": human_missing_dots,
              "action": "If individual samples are displayed, mark these as not packaged, not absent biological expression. Do not call these candidates unmeasured in the human study."})
notes.append({"severity": "required_caveat", "issue": "Human/mouse/organoid analyses are descriptive; rat Welch q and post hoc exact q lead to different fixed-family decisions. Four organoid libraries come from one differentiation."})
notes.append({"severity": "scope", "issue": "This review verifies data packaging against frozen study artifacts and independently recomputes the stated display transformations. It does not redo raw sequence quantification, orthology discovery, whole-transcriptome PCA, original atlas embedding, or all primary study methods; no expert/clinical review or browser/UI validation is claimed."})
failed = [c for c in checks if not c["pass"]]
report = {
    "status": "PASS_WITH_CAVEATS" if not failed else "NEEDS_REVISION",
    "reviewer": "Separate Codex agent; independent code path, not an independent human/domain expert",
    "reviewed_at_utc": datetime.now(timezone.utc).isoformat(),
    "python": sys.version,
    "checks_n": len(checks), "passed_n": len(checks) - len(failed),
    "blockers": failed, "checks": checks, "notes": notes,
    "reviewed_payload_sha256": payload_hashes,
    "reviewed_source_sha256": sources,
    "max_rounding_absolute_error": max_errors,
    "independently_recomputed_rat_statistics": rat_recomputed,
}
destination = Path(__file__).resolve().with_suffix(".json")
destination.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
summary = ["# Independent public-data review", "", f"Assessment: **{report['status']}**.", "",
           f"{report['passed_n']} / {report['checks_n']} checks passed. Review date: {report['reviewed_at_utc']}.", "",
           "This is a separate Codex agent's calculation and provenance review, not human scientific peer review. The builder was read but neither imported nor executed by the review script.", "",
           "## Verified", "",
           "- All 25 fixed candidate identities and missing mappings are retained; 15 primary mouse mappings remain separate from the two ACAT2 secondary mappings.",
           "- Human and mouse display effects were independently recomputed from frozen means/individual TPM values; sample groups and units were checked.",
           "- All 15 rat candidate effects, nominal Welch intervals, Welch P values, all 252 label allocations and both BH-25 adjustments were independently recalculated for both normalizations. VNN1 exact q = 25/126 = 0.1984126984 remains visible alongside Welch q < 0.05.",
           "- All 15 Cartesian configuration shares were recomputed from the 64 human and 64 mouse effect configurations; range endpoints were checked. These are descriptive configurations, not biological replicates or probabilities.",
           "- All 24,034 author-saved UMAP coordinates and every library/condition/subtype/broad-category code match retained source rows within six-decimal rounding precision.",
           "- Every value for 24 measured genes (576,816 cell-gene combinations) was independently recomputed as natural log1p(10,000 × raw candidate UMI / all-gene raw UMI total). AOC4P remains null, not zero.",
           "- All 384 library/category aggregates and 48 composition decompositions were independently derived from the retained cells. Author coordinates were not reconstructed or refitted.",
           "- All PCA scores, sample labels, variance fractions, feature counts and all 12 downloadable CSV byte hashes match frozen sources.",
           "- No Windows/home paths, credential-like strings, original 10x barcodes or original cell-identity field names were found in public data files.", "",
           "## Blockers", ""]
summary += [f"- {c['name']}: {c['details']}" for c in failed] if failed else ["None in the reviewed data payload."]
summary += ["", "## Nonblocking limitations", "",
            "- Ten human candidates have aggregate effects but no individual sample dots packaged: " + ", ".join(human_missing_dots) + ". The interface must distinguish unavailable packaged dots from missing human measurements.",
            "- Different units, species, tissue models and exposure designs must remain apparent. Rat Welch/exact disagreement and organoid technical replication limits must be visible to readers.",
            "- Source-data provenance and old scientific reviews are carried forward; this review did not rerun raw sequencing, discover orthology, refit PCA, rerun UMAP or certify clinical benefit.",
            "- Browser rendering and front-end labels are outside this review. A later payload change invalidates the corresponding hash match and needs focused rechecking.", "",
            "The adjacent JSON records every check, source/data SHA-256 and independently recalculated rat statistics. Re-run this script using the project's existing scientific Python environment."]
Path(__file__).resolve().with_suffix(".md").write_text("\n".join(summary) + "\n", encoding="utf-8")
print(json.dumps({"status": report["status"], "checks": len(checks), "failed": [f["name"] for f in failed],
                  "report": str(destination.relative_to(APP)), "human_missing_dots_n": len(human_missing_dots)}))
if failed:
    raise SystemExit(1)
