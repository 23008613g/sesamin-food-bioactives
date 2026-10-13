"""Descriptive full-universe rat/HepG2 comparison under Stage 39 protocol.

No human inferential test, gene-independence p value, correlation CI, candidate
selection, pathway activation claim, or equivalence of physical effect scales.
Only reads inputs; writes results/cross_species and its own execution record.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata


ROOT = Path(__file__).resolve().parents[1]
PROJECT = Path("C:/Users/AI/OneDrive/Desktop/公司文章/research_planning/hepatic_lipotoxicity_2026-10-07")
OUT = ROOT / "results" / "cross_species"
REVIEW = ROOT / "review" / "cross_species_execution.json"
HUMAN_ORDER = [f"{group}_{i}" for group in ("Control", "Model", "Treat") for i in range(1, 4)]
RAT_ORDER = [f"GSM{i}" for i in range(1048036, 1048046)]
VARIANTS = ("deposited", "qn10", "qn20")
HUMAN_NORMALIZATIONS = ("mr", "cpm")
CONTRASTS = ("treat_vs_model", "model_vs_control")
HUMAN_SOURCE_SHA = "645e62b32be158a113ad31d1e35e0cb7431dc6efbcc4e7d5b83a46f5a63f18c0"


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def save(frame, name):
    path = OUT / name
    if name.endswith(".gz"):
        frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})
    else:
        frame.to_csv(path, index=False)


def close_report(actual, expected, name, atol=1e-10, rtol=1e-10):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    if a.shape != b.shape or not np.allclose(a, b, atol=atol, rtol=rtol, equal_nan=False):
        raise ValueError(f"Reproduction failed: {name}, shapes {a.shape} vs {b.shape}")
    return {"name": name, "passed": True, "max_absolute_error": float(np.max(np.abs(a - b))),
            "atol": atol, "rtol": rtol, "n": int(a.size)}


def directions(values):
    # Exact numerical zeros are a separate category, not assigned to positive.
    return np.where(values > 0, "positive", np.where(values < 0, "negative", "zero"))


def correlation(x, y):
    # Average tied ranks; deliberately calculate no inferential p value.
    rx, ry = rankdata(x, method="average"), rankdata(y, method="average")
    if np.std(rx) == 0 or np.std(ry) == 0:
        return None
    return float(np.corrcoef(rx, ry)[0, 1])


def describe_pair(x, y):
    xd, yd = directions(x), directions(y)
    same = (xd == yd) & (xd != "zero")
    opposite = (xd != yd) & (xd != "zero") & (yd != "zero")
    both_nonzero = (xd != "zero") & (yd != "zero")
    n = len(x)
    return {"n": n, "spearman_descriptive": correlation(x, y),
            "same_nonzero_direction_n": int(same.sum()),
            "opposite_nonzero_direction_n": int(opposite.sum()),
            "either_zero_n": int((~both_nonzero).sum()),
            "both_zero_n": int(((xd == "zero") & (yd == "zero")).sum()),
            "same_nonzero_direction_fraction_all": float(same.mean()) if n else None,
            "opposite_nonzero_direction_fraction_all": float(opposite.mean()) if n else None,
            "same_direction_fraction_among_both_nonzero": float(same.sum() / both_nonzero.sum()) if both_nonzero.any() else None}


def gene_z(matrix):
    if not np.isfinite(matrix).all():
        raise ValueError("Nonfinite expression value")
    sd = matrix.std(axis=1, ddof=1)
    z = np.zeros_like(matrix, dtype=float)
    positive = sd > 0
    z[positive] = (matrix[positive] - matrix[positive].mean(axis=1, keepdims=True)) / sd[positive, None]
    if positive.any():
        assert np.allclose(z[positive].mean(axis=1), 0, atol=1e-10)
        assert np.allclose(z[positive].std(axis=1, ddof=1), 1, atol=1e-10)
    return z, int((~positive).sum())


def main(project):
    global PROJECT
    PROJECT = project.resolve()
    started = utc()
    source = PROJECT / "next_steps/candidate_screen_03/sesamin_inputs/GSE248903_total_count.anno.txt.gz"
    archive = PROJECT / "next_steps/sesamin_exploration_04/all_gene_descriptive.csv.gz"
    old_qc = PROJECT / "next_steps/sesamin_exploration_04/sample_qc.csv"
    mapping_path = ROOT / "sources/orthology/human_rat_one2one_highconfidence_unique_entrez.tsv"
    map_audit_path = ROOT / "sources/orthology/human_rat_mapping_audit.tsv"
    rat_ann_path = ROOT / "results/rat_gene_annotation.csv"
    go_members_path = ROOT / "results/go_membership.csv.gz"
    go_results_path = ROOT / "results/go_results_deposited.csv"
    protocol = ROOT / "PROTOCOL.txt"
    inputs = [source, archive, old_qc, mapping_path, map_audit_path, rat_ann_path, go_members_path,
              go_results_path, protocol, ROOT / "PROTOCOL_ADDENDUM_1.txt", ROOT / "protocol_freeze.json",
              ROOT / "PROTOCOL_ADDENDUM_1.freeze.json", ROOT / "sources/orthology/orthology_freeze_report.json"]
    for variant in VARIANTS:
        inputs += [ROOT / f"results/rat_gene_expression_{variant}.csv.gz",
                   ROOT / f"results/rat_gene_results_{variant}.csv.gz"]
    missing = [str(p) for p in inputs if not p.is_file()]
    if missing:
        raise FileNotFoundError("Upstream inputs not ready; analysis not started: " + "; ".join(missing))
    if sha(source) != HUMAN_SOURCE_SHA:
        raise ValueError("Human count input differs from frozen historical input")
    protocol_freeze = json.loads((ROOT / "protocol_freeze.json").read_text(encoding="utf-8"))
    addendum_freeze = json.loads((ROOT / "PROTOCOL_ADDENDUM_1.freeze.json").read_text(encoding="utf-8"))
    if sha(protocol) != protocol_freeze["protocol_sha256"]:
        raise ValueError("Protocol hash differs from frozen value")
    if sha(ROOT / "PROTOCOL_ADDENDUM_1.txt") != addendum_freeze["sha256"]:
        raise ValueError("Addendum hash differs from frozen value")
    OUT.mkdir(parents=True, exist_ok=True)
    record = {"status": "running", "started_utc": started, "command": sys.argv,
              "script": str(Path(__file__).resolve()), "script_sha256": sha(Path(__file__)),
              "protocol_sha256": sha(protocol), "protocol_addendum_sha256": sha(ROOT / "PROTOCOL_ADDENDUM_1.txt"),
              "python": sys.version, "executable": sys.executable,
              "platform": platform.platform(), "versions": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
              "inputs": [{"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size} for p in inputs],
              "interpretation": "Post-hoc exploratory cross-study description; rat animals and HepG2 libraries differ in biological model, species, challenge, exposure, assay and experimental-unit certainty. Neither gene correlations nor matched-member scores validate normal-human responses or causal rescue.",
              "inference_policy": "No human p-values, no gene-independence correlation tests, no binomial intervals. Human libraries are not asserted to be independent cultures. GO members are not independent biological replicates.",
              "effect_scales": "Rat: difference of group means of log2 microarray signals. Human: log2 ratio of arithmetic group means of normalized counts plus 1. These are not identical physical effect scales.",
              "matched_score_scales": "Rat gene matrices and human log2(normalized_count+1) matrices separately standardized over all 10 animals / all 9 libraries using ddof=1; identical matched gene members averaged equally. Standardized score differences are study-relative descriptors, not quantitative equivalence of intervention effects.",
              "selection_policy": "All eligible rat genes and CPM-retained human genes with frozen strict orthology; no candidate, effect, significance, stability or direction filtering."}
    write_json(REVIEW, record)
    table = pd.read_csv(source, sep="\t", usecols=["name", "GeneName"] + HUMAN_ORDER,
                        dtype={"name": str, "GeneName": str}, keep_default_na=False)
    counts = table[HUMAN_ORDER].to_numpy(dtype=float)
    if counts.shape != (60605, 9) or not table["name"].is_unique:
        raise ValueError("Unexpected human full count matrix or duplicate gene key")
    if not np.isfinite(counts).all() or (counts < 0).any() or (counts % 1 != 0).any():
        raise ValueError("Invalid human counts")
    totals = counts.sum(axis=0)
    cpm = counts / totals[None, :] * 1e6
    retained = (cpm >= 1).sum(axis=1) >= 3
    if retained.sum() != 13014:
        raise ValueError("Historical CPM filter not reproduced")
    positive = (counts > 0).all(axis=1)
    if positive.sum() != 15412:
        raise ValueError("Historical median-ratio reference set not reproduced")
    geometric = np.exp(np.log(counts[positive]).mean(axis=1))
    size_factors = np.median(counts[positive] / geometric[:, None], axis=0)
    size_factors /= np.exp(np.log(size_factors).mean())
    normalized = {"mr": counts / size_factors, "cpm": cpm}
    effects = {}
    human_result = pd.DataFrame({"human_ensembl_gene_id": table["name"], "human_original_symbol": table["GeneName"],
                                 "retained_cpm_filter": retained, "cpm_ge1_libraries": (cpm >= 1).sum(axis=1)})
    for norm, matrix in normalized.items():
        groupmeans = np.column_stack([matrix[:, start:start+3].mean(axis=1) for start in (0, 3, 6)])
        c, m, t = groupmeans.T
        effects[norm] = {"treat_vs_model": np.log2((t + 1) / (m + 1)),
                         "model_vs_control": np.log2((m + 1) / (c + 1)),
                         "treat_vs_control": np.log2((t + 1) / (c + 1))}
        for contrast, values in effects[norm].items():
            human_result[f"{norm}_{contrast}"] = values
        for group, values in zip(("control", "model", "treat"), groupmeans.T):
            human_result[f"{norm}_mean_{group}"] = values
    archived = pd.read_csv(archive, dtype={"gene_id": str}, keep_default_na=False).set_index("gene_id")
    if not archived.index.is_unique or set(archived.index) != set(table["name"]):
        raise ValueError("Archived full-gene table does not match source keys")
    archived = archived.loc[table["name"]]
    old_keep = archived["retained_cpm_filter"].astype(str).str.lower().eq("true").to_numpy()
    if not np.array_equal(old_keep, retained):
        raise ValueError("Archived CPM filter differs")
    checks = [close_report((cpm >= 1).sum(axis=1), archived["cpm_ge1_libraries"], "human_cpm_filter_counts")]
    old_names = {"treat_vs_model": "log2_sesamin_vs_OA", "model_vs_control": "log2_OA_vs_control", "treat_vs_control": "log2_sesamin_vs_control"}
    for norm in HUMAN_NORMALIZATIONS:
        for contrast, old_name in old_names.items():
            old_name = ("cpm_" if norm == "cpm" else "") + old_name
            checks.append(close_report(effects[norm][contrast], archived[old_name], f"human_{norm}_{contrast}"))
    old_sample = pd.read_csv(old_qc).set_index("sample").loc[HUMAN_ORDER]
    checks.append(close_report(size_factors, old_sample["median_ratio_factor"], "human_median_ratio_size_factors"))
    checks.append(close_report(totals, old_sample["assigned_counts"], "human_all_record_library_totals"))
    for group, oldgroup in (("control", "control"), ("model", "OA"), ("treat", "sesamin_OA")):
        checks.append(close_report(human_result[f"mr_mean_{group}"], archived[f"mean_norm_{oldgroup}"], f"human_mr_mean_{group}"))
    save(human_result, "human_all_gene_recomputed.csv.gz")
    save(pd.DataFrame({"sample": HUMAN_ORDER, "group": np.repeat(["Control", "Model", "Treat"], 3),
                       "library_total_all_60605": totals, "median_ratio_size_factor": size_factors}), "human_normalization_factors.csv")
    rat_ann = pd.read_csv(rat_ann_path, dtype={"rat_entrez": str}, keep_default_na=False)
    if not rat_ann["rat_entrez"].is_unique:
        raise ValueError("Rat annotation Entrez keys are not unique")
    orth = pd.read_csv(mapping_path, sep="\t", dtype=str, keep_default_na=False)
    for key in ("human_ensembl_gene_id", "rat_ensembl_gene_id", "rat_entrez_id"):
        if not orth[key].is_unique or (orth[key] == "").any():
            raise ValueError("Strict orthology map key invalid: " + key)
    if not orth["orthology_type"].eq("ortholog_one2one").all() or not orth["orthology_confidence"].eq("1").all():
        raise ValueError("Non-strict orthology row")
    humans = set(table["name"])
    kept = set(table.loc[retained, "name"])
    rats = set(rat_ann["rat_entrez"])
    mapped = orth.copy()
    mapped["human_in_source"] = mapped["human_ensembl_gene_id"].isin(humans)
    mapped["human_retained_cpm"] = mapped["human_ensembl_gene_id"].isin(kept)
    mapped["rat_in_eligible_universe"] = mapped["rat_entrez_id"].isin(rats)
    mapped["included_cross_species"] = mapped["human_retained_cpm"] & mapped["rat_in_eligible_universe"]
    mapped["exclusion_reason"] = ["included" if included else ";".join(filter(None, [
        "human_absent_from_counts" if not hsource else "human_fails_cpm_filter" if not hkeep else "",
        "rat_not_in_eligible_universe" if not rkeep else ""])) for hsource, hkeep, rkeep, included in
        zip(mapped.human_in_source, mapped.human_retained_cpm, mapped.rat_in_eligible_universe, mapped.included_cross_species)]
    save(mapped, "strict_mapping_inclusion_audit.csv")
    full_audit = pd.read_csv(map_audit_path, sep="\t", dtype=str, keep_default_na=False)
    strict_pairs = set(zip(orth.human_ensembl_gene_id, orth.rat_ensembl_gene_id))
    included_pairs = set(zip(mapped.loc[mapped.included_cross_species, "human_ensembl_gene_id"], mapped.loc[mapped.included_cross_species, "rat_ensembl_gene_id"]))
    reason_lookup = dict(zip(zip(mapped.human_ensembl_gene_id, mapped.rat_ensembl_gene_id), mapped.exclusion_reason))
    full_audit["strict_mapping_eligible"] = [(h, r) in strict_pairs for h, r in zip(full_audit.human_ensembl_gene_id, full_audit.rat_ensembl_gene_id)]
    full_audit["included_cross_species"] = [(h, r) in included_pairs for h, r in zip(full_audit.human_ensembl_gene_id, full_audit.rat_ensembl_gene_id)]
    full_audit["cross_species_reason"] = [reason_lookup.get((h, r), ";".join(filter(None, [
        status if status != "included_one2one" else "", "not_highconfidence" if conf != "1" else "",
        entrez if entrez != "unique_bidirectional" else ""]))) for h, r, status, conf, entrez in
        zip(full_audit.human_ensembl_gene_id, full_audit.rat_ensembl_gene_id, full_audit.orthology_mapping_status,
            full_audit.orthology_confidence, full_audit.entrez_mapping_status)]
    save(full_audit, "all_orthology_rows_inclusion_audit.csv.gz")
    shared = mapped.loc[mapped.included_cross_species, orth.columns].sort_values("human_ensembl_gene_id").reset_index(drop=True)
    if len(shared) < 1000:
        raise ValueError("Unexpectedly small shared universe; investigate rather than proceed")
    shared = shared.merge(rat_ann.rename(columns={"rat_entrez": "rat_entrez_id", "gene_symbol": "rat_platform_symbol"}),
                          on="rat_entrez_id", how="left", validate="one_to_one")
    save(shared, "shared_gene_mapping.csv")
    human_audit = human_result[["human_ensembl_gene_id", "human_original_symbol", "retained_cpm_filter"]].copy()
    human_audit["strict_orthology_available"] = human_audit.human_ensembl_gene_id.isin(orth.human_ensembl_gene_id)
    human_audit["included_cross_species"] = human_audit.human_ensembl_gene_id.isin(shared.human_ensembl_gene_id)
    save(human_audit, "human_gene_universe_inclusion_audit.csv.gz")
    rat_audit = rat_ann.copy()
    rat_audit["strict_orthology_available"] = rat_audit.rat_entrez.isin(orth.rat_entrez_id)
    rat_audit["included_cross_species"] = rat_audit.rat_entrez.isin(shared.rat_entrez_id)
    save(rat_audit, "rat_gene_universe_inclusion_audit.csv")
    human_position = pd.Series(np.arange(len(table)), index=table["name"]).loc[shared.human_ensembl_gene_id].to_numpy()
    shared_rat_ids = shared.rat_entrez_id.tolist()
    rat_matrices, rat_effects = {}, {}
    for variant in VARIANTS:
        frame = pd.read_csv(ROOT / f"results/rat_gene_expression_{variant}.csv.gz", dtype={"rat_entrez": str}, keep_default_na=False).set_index("rat_entrez")
        res = pd.read_csv(ROOT / f"results/rat_gene_results_{variant}.csv.gz", dtype={"rat_entrez": str}).set_index("rat_entrez")
        if not frame.index.is_unique or not res.index.is_unique or set(frame.index) != rats or set(res.index) != rats:
            raise ValueError("Rat gene universe differs by input or normalization")
        if frame.columns.tolist() != RAT_ORDER:
            raise ValueError("Rat sample columns do not match explicit ordered controls/sesamin")
        vals = frame.loc[shared_rat_ids, RAT_ORDER].to_numpy(dtype=float)
        if not np.isfinite(vals).all():
            raise ValueError("Nonfinite rat values")
        delta = vals[:, 5:].mean(axis=1) - vals[:, :5].mean(axis=1)
        checks.append(close_report(delta, res.loc[shared_rat_ids, "difference"], f"rat_{variant}_shared_effect_recomputation"))
        rat_matrices[variant], rat_effects[variant] = vals, delta
        save(pd.concat([shared[["human_ensembl_gene_id", "rat_ensembl_gene_id", "rat_entrez_id"]], pd.DataFrame(vals, columns=RAT_ORDER)], axis=1), f"shared_rat_expression_{variant}.csv.gz")
    gene_summaries, gene_contingencies = [], []
    for norm in HUMAN_NORMALIZATIONS:
        loghuman = np.log2(normalized[norm][human_position] + 1)
        save(pd.concat([shared[["human_ensembl_gene_id", "rat_ensembl_gene_id", "rat_entrez_id"]], pd.DataFrame(loghuman, columns=HUMAN_ORDER)], axis=1), f"shared_human_log2_expression_{norm}.csv.gz")
        for variant in VARIANTS:
            pair = shared.copy()
            pair["rat_normalization"], pair["human_normalization"] = variant, norm
            pair["rat_sesamin_minus_control"] = rat_effects[variant]
            pair["rat_direction"] = directions(rat_effects[variant])
            for contrast in CONTRASTS:
                human_effect = effects[norm][contrast][human_position]
                pair[f"human_{contrast}"] = human_effect
                pair[f"human_{contrast}_direction"] = directions(human_effect)
                context = {"rat_normalization": variant, "human_normalization": norm, "human_contrast": contrast}
                gene_summaries.append(context | describe_pair(rat_effects[variant], human_effect))
                for rd in ("negative", "zero", "positive"):
                    for hd in ("negative", "zero", "positive"):
                        n = int(((pair.rat_direction == rd) & (pair[f"human_{contrast}_direction"] == hd)).sum())
                        gene_contingencies.append(context | {"rat_direction": rd, "human_direction": hd, "n": n, "fraction_all_shared": n / len(shared)})
            save(pair, f"shared_gene_effects_rat_{variant}_human_{norm}.csv.gz")
    save(pd.DataFrame(gene_summaries), "gene_comparison_summary.csv")
    save(pd.DataFrame(gene_contingencies), "gene_direction_contingency.csv")
    membership = pd.read_csv(go_members_path, dtype={"rat_entrez": str, "go_id": str}, keep_default_na=False)
    metadata = pd.read_csv(go_results_path, dtype={"go_id": str}, keep_default_na=False)
    if not metadata.go_id.is_unique or membership.duplicated(["go_id", "rat_entrez"]).any():
        raise ValueError("GO ID / membership keys are not unique")
    if not set(membership.go_id).issubset(set(metadata.go_id)) or not set(membership.rat_entrez).issubset(rats):
        raise ValueError("GO membership outside frozen gene/term universe")
    orig_n = membership.groupby("go_id").rat_entrez.nunique()
    if set(orig_n.index) != set(metadata.go_id) or not orig_n.between(15, 500).all():
        raise ValueError("Primary GO universe not restricted to 15-500 genes")
    checks.append(close_report(orig_n.loc[metadata.go_id], metadata.n_genes, "primary_GO_membership_counts"))
    go_map = shared[["human_ensembl_gene_id", "human_symbol", "rat_ensembl_gene_id", "rat_entrez_id", "rat_symbol"]].rename(columns={"rat_entrez_id": "rat_entrez"})
    matched = membership.merge(go_map, on="rat_entrez", how="inner", validate="many_to_one")
    match_n = matched.groupby("go_id").rat_entrez.nunique()
    go_audit = metadata[["go_id", "go_name", "n_genes"]].rename(columns={"n_genes": "n_original_rat_genes"}).copy()
    go_audit["n_matched_genes"] = go_audit.go_id.map(match_n).fillna(0).astype(int)
    go_audit["included_matched_score"] = go_audit.n_matched_genes.between(15, 500)
    go_audit["matched_fraction"] = go_audit.n_matched_genes / go_audit.n_original_rat_genes
    go_audit["reason"] = np.where(go_audit.included_matched_score, "included", "fewer_than_15_shared_genes")
    save(go_audit, "matched_go_inclusion_audit.csv")
    go_ids = sorted(go_audit.loc[go_audit.included_matched_score, "go_id"])
    matched = matched.loc[matched.go_id.isin(go_ids)].sort_values(["go_id", "rat_entrez"]).reset_index(drop=True)
    if len(go_ids) == 0:
        raise ValueError("No GO sets meet matched-member range")
    save(matched, "matched_go_membership.csv.gz")
    pos = {key: i for i, key in enumerate(shared_rat_ids)}
    members = {go: np.array([pos[key] for key in frame.rat_entrez], dtype=int) for go, frame in matched.groupby("go_id", sort=True)}
    go_name = metadata.set_index("go_id").go_name.to_dict()
    gene_set_sha = {go: hashlib.sha256("\n".join(sorted(shared_rat_ids[i] for i in members[go])).encode()).hexdigest() for go in go_ids}
    scores, score_rows, zero_variance = {}, [], {}
    for study, norms, matrices, order, groups in [
        ("rat", VARIANTS, rat_matrices, RAT_ORDER, ["Control"] * 5 + ["Sesamin"] * 5),
        ("human", HUMAN_NORMALIZATIONS, {n: np.log2(normalized[n][human_position] + 1) for n in HUMAN_NORMALIZATIONS}, HUMAN_ORDER, ["Control"] * 3 + ["Model"] * 3 + ["Treat"] * 3)]:
        for norm in norms:
            z, nzero = gene_z(matrices[norm])
            zero_variance[f"{study}_{norm}"] = nzero
            scored = np.vstack([z[members[go]].mean(axis=0) for go in go_ids])
            scores[(study, norm)] = scored
            for i, go in enumerate(go_ids):
                for j, sample in enumerate(order):
                    score_rows.append({"go_id": go, "go_name": go_name[go], "n_matched_genes": len(members[go]),
                                       "member_entrez_sha256": gene_set_sha[go], "study": study, "normalization": norm,
                                       "sample": sample, "group": groups[j], "mean_member_gene_z": float(scored[i, j])})
            wide = pd.DataFrame(scored, columns=order)
            wide.insert(0, "go_id", go_ids)
            save(wide, f"matched_go_scores_{study}_{norm}.csv")
    save(pd.DataFrame(score_rows), "matched_go_sample_scores.csv.gz")
    go_pairs, go_summary, go_contingency = [], [], []
    for variant in VARIANTS:
        rs = scores[("rat", variant)]
        rd = rs[:, 5:].mean(axis=1) - rs[:, :5].mean(axis=1)
        for norm in HUMAN_NORMALIZATIONS:
            hs = scores[("human", norm)]
            hd = {"treat_vs_model": hs[:, 6:].mean(axis=1) - hs[:, 3:6].mean(axis=1),
                  "model_vs_control": hs[:, 3:6].mean(axis=1) - hs[:, :3].mean(axis=1)}
            for contrast in CONTRASTS:
                context = {"rat_normalization": variant, "human_normalization": norm, "human_contrast": contrast}
                rdir, hdir = directions(rd), directions(hd[contrast])
                go_summary.append(context | describe_pair(rd, hd[contrast]))
                for r in ("negative", "zero", "positive"):
                    for h in ("negative", "zero", "positive"):
                        n = int(((rdir == r) & (hdir == h)).sum())
                        go_contingency.append(context | {"rat_direction": r, "human_direction": h, "n": n, "fraction_all_matched_sets": n / len(go_ids)})
                for i, go in enumerate(go_ids):
                    go_pairs.append(context | {"go_id": go, "go_name": go_name[go], "n_matched_genes": len(members[go]),
                        "n_original_rat_genes": int(orig_n[go]), "member_entrez_sha256": gene_set_sha[go],
                        "rat_mean_control_score": float(rs[i, :5].mean()), "rat_mean_sesamin_score": float(rs[i, 5:].mean()),
                        "rat_score_difference": float(rd[i]), "human_mean_control_score": float(hs[i, :3].mean()),
                        "human_mean_model_score": float(hs[i, 3:6].mean()), "human_mean_treat_score": float(hs[i, 6:].mean()),
                        "human_score_difference": float(hd[contrast][i]), "rat_direction": rdir[i], "human_direction": hdir[i]})
    save(pd.DataFrame(go_pairs), "matched_go_cross_species_effects.csv.gz")
    save(pd.DataFrame(go_summary), "matched_go_comparison_summary.csv")
    save(pd.DataFrame(go_contingency), "matched_go_direction_contingency.csv")
    readme = """# Descriptive cross-study expression comparison

All outputs use the full retained human and eligible rat gene universes; no 25-candidate selection was applied.
`shared_gene_mapping.csv` defines the exact shared universe; full inclusion/exclusion audits are retained.
`gene_comparison_summary.csv` reports all 3 rat normalizations x 2 human normalizations x 2 human contrasts.
`gene_direction_contingency.csv` includes all nine negative/zero/positive combinations for every comparison.
`shared_gene_effects_rat_*_human_*.csv.gz` retain every mapped gene without effect or significance thresholds.

Rat effects are differences of group means of log microarray signals; human effects are log ratios of mean normalized counts plus 1. They are not the same physical effect scale. Spearman coefficients are descriptive and no correlation p values or confidence intervals are computed.

`matched_go_membership.csv.gz` defines identical human/rat members for each included GO set. Original sets must retain 15-500 shared genes. Human log2(normalized counts+1) and rat log-signal gene matrices are standardized separately across their nine libraries and ten animals (ddof=1), then member gene z scores are averaged equally. Constant genes contribute zero, if present. Complete individual sample scores, contrasts, set-membership hashes, exclusion counts and normalization sensitivities are retained.

Human libraries have unresolved independence of culture units and exposure timing. No human inferential tests were performed. Gene members and overlapping GO sets are not independent biological replicates. No result establishes normal human hepatocyte validation, treatment rescue, pathway activation, metabolic flux or causality.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    record.update({"status": "complete", "completed_utc": utc(), "checks": checks,
        "counts": {"human_raw_genes": len(table), "human_retained_cpm": int(retained.sum()),
            "human_median_ratio_reference_genes": int(positive.sum()), "rat_eligible_genes": len(rat_ann),
            "strict_orthology_pairs": len(orth), "shared_genes": len(shared), "primary_GO_sets": len(metadata),
            "matched_GO_sets": len(go_ids), "matched_GO_membership_rows": len(matched)},
        "zero_variance_shared_genes": zero_variance,
        "strict_mapping_exclusion_counts": {str(k): int(v) for k, v in mapped.exclusion_reason.value_counts().items()},
        "gene_comparison_summary": gene_summaries, "matched_GO_comparison_summary": go_summary,
        "outputs": [{"path": str(p.relative_to(ROOT)), "sha256": sha(p), "bytes": p.stat().st_size} for p in sorted(OUT.iterdir()) if p.is_file()]})
    write_json(REVIEW, record)
    print(json.dumps({"status": "complete", "counts": record["counts"], "max_reproduction_error": max(c["max_absolute_error"] for c in checks)}, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=PROJECT)
    args = parser.parse_args()
    try:
        main(args.project)
    except Exception as exc:
        failure = {"status": "failed", "utc": utc(), "error": repr(exc), "traceback": traceback.format_exc(),
                   "command": sys.argv, "script_sha256": sha(Path(__file__))}
        if REVIEW.exists():
            failure["execution_record_before_failure"] = json.loads(REVIEW.read_text(encoding="utf-8"))
        write_json(REVIEW, failure)
        raise
