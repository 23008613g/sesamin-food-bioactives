"""Descriptive post-hoc diagnostics of frozen limma outputs; no model rerun.

All original 14,604 genes are retained. BH q<0.05 is the existing reporting
threshold. Empirical signal quartiles are descriptive strata, not detection
thresholds or new filters. The same ten animals underlie every normalization.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

VARIANTS = ["deposited", "qn10", "qn20"]
LABELS = {"deposited": "Deposited", "qn10": "QN10", "qn20": "QN20"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    def scalar(v):
        if isinstance(v, np.generic):
            return v.item()
        raise TypeError(type(v).__name__)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False, default=scalar) + "\n", encoding="utf-8")


def save_csv(frame, path):
    opts = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    frame.to_csv(path, index=False, compression=opts)


def quantiles(x):
    q = np.quantile(np.asarray(x, float), [0, .25, .5, .75, 1], method="linear")
    return dict(zip(["minimum", "q25", "median", "q75", "maximum"], q.tolist()))


def rho(x, y):
    a, b = rankdata(x, method="average"), rankdata(y, method="average")
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage39", type=Path, default=Path("D:/ResearchData/pfhn_genome39_2026-10-13"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    stage, out = args.stage39.resolve(), args.output_dir.resolve()
    assert stage != out and stage not in out.parents, "Diagnostic output must be outside the read-only source stage."
    out.mkdir(parents=True, exist_ok=True)
    records, summaries, quartiles, correlations = {}, [], [], []
    checks, inputs, hit_sets = [], {}, {}
    reference_ids = None
    for v in VARIANTS:
        lp = stage / "results/limma" / f"rat_limma_{v}.csv.gz"
        gp = stage / "results" / f"rat_gene_results_{v}.csv.gz"
        inputs[str(lp)] = sha(lp); inputs[str(gp)] = sha(gp)
        l = pd.read_csv(lp, dtype={"gene_id": str})
        g = pd.read_csv(gp, dtype={"rat_entrez": str})
        assert len(l) == len(g) == 14604
        assert l.gene_id.is_unique and g.rat_entrez.is_unique
        assert l.gene_id.tolist() == g.rat_entrez.tolist()
        if reference_ids is None:
            reference_ids = l.gene_id.tolist()
        assert l.gene_id.tolist() == reference_ids
        assert (l.n_control == 5).all() and (l.n_sesamin == 5).all()
        assert np.isfinite(l[["AveExpr", "logFC", "adj.P.Val", "mean_control", "mean_sesamin"]]).all().all()
        signal_error = float(np.max(np.abs(l.AveExpr - g.mean_sample_signal)))
        signal_mean_error = float(np.max(np.abs(l.AveExpr - (l.mean_control + l.mean_sesamin) / 2)))
        effect_error = float(np.max(np.abs(l.logFC - g.difference)))
        assert max(signal_error, signal_mean_error, effect_error) < 1e-10
        d = l.rename(columns={"gene_id": "rat_entrez"}).copy()
        d.insert(1, "gene_symbol", g.gene_symbol.fillna("").to_numpy())
        d.insert(2, "n_probes", g.n_probes.to_numpy())
        d.insert(0, "normalization", v)
        d["signal_rank_right"] = np.searchsorted(np.sort(d.AveExpr), d.AveExpr, side="right")
        d["signal_empirical_percentile"] = 100 * d.signal_rank_right / len(d)
        d["signal_descriptive_quartile"] = np.minimum(4, np.ceil(d.signal_empirical_percentile / 25)).astype(int)
        d["limma_BH_hit"] = d["adj.P.Val"] < .05
        d["effect_sign"] = np.sign(d.logFC).astype(int)
        d["absolute_logFC"] = d.logFC.abs()
        assert ((d.signal_empirical_percentile > 0) & (d.signal_empirical_percentile <= 100)).all()
        hit_sets[v] = set(d.loc[d.limma_BH_hit, "rat_entrez"])
        records[v] = d
        save_csv(d, out / f"limma_all_genes_signal_context_{v}.csv.gz")
        save_csv(d.loc[d.limma_BH_hit].sort_values(["signal_empirical_percentile", "rat_entrez"]), out / f"limma_hits_{v}.csv")
        for pop, mask in [("all_eligible_genes", np.ones(len(d), bool)), ("BH_hits", d.limma_BH_hit), ("not_BH_hits", ~d.limma_BH_hit)]:
            sub = d.loc[mask]
            row = {"normalization": v, "population": pop, "n": len(sub)}
            for name, vals in [("mean_log2_signal", sub.AveExpr), ("signal_percentile", sub.signal_empirical_percentile), ("logFC", sub.logFC), ("absolute_logFC", sub.absolute_logFC)]:
                row.update({name + "_" + k: val for k, val in quantiles(vals).items()})
            row.update({"positive_effect_n": int((sub.effect_sign > 0).sum()), "negative_effect_n": int((sub.effect_sign < 0).sum()), "zero_effect_n": int((sub.effect_sign == 0).sum())})
            summaries.append(row)
            correlations.append({"normalization": v, "population": pop, "n": len(sub),
                                 "spearman_mean_signal_vs_signed_logFC": rho(sub.AveExpr, sub.logFC),
                                 "spearman_mean_signal_vs_absolute_logFC": rho(sub.AveExpr, sub.absolute_logFC),
                                 "inference": "descriptive rank correlation only; no p value or interval"})
        for q in range(1, 5):
            mask = d.signal_descriptive_quartile == q
            n_hits = int((mask & d.limma_BH_hit).sum())
            quartiles.append({"normalization": v, "empirical_signal_quartile": q,
                              "percentile_lower_exclusive": 25 * (q - 1), "percentile_upper_inclusive": 25 * q,
                              "all_gene_n": int(mask.sum()), "BH_hit_n": n_hits,
                              "percent_of_all_hits": 100 * n_hits / len(hit_sets[v]),
                              "interpretation": "exploratory descriptive rank stratum; not an assay detection threshold or filtering rule"})
        checks.append({"normalization": v, "genes": len(d), "BH_hits": len(hit_sets[v]),
                       "maximum_AveExpr_vs_frozen_all_gene_mean_error": signal_error,
                       "maximum_AveExpr_vs_two_group_mean_average_error": signal_mean_error,
                       "maximum_logFC_vs_frozen_difference_error": effect_error,
                       "all_gene_ids_and_order_retained": True})

    union = set.union(*hit_sets.values()); intersection = set.intersection(*hit_sets.values())
    wide = pd.DataFrame({"rat_entrez": reference_ids, "gene_symbol": records["deposited"].gene_symbol})
    for v in VARIANTS:
        for c in ["limma_BH_hit", "logFC", "AveExpr", "signal_empirical_percentile", "effect_sign", "adj.P.Val", "CI95_moderated_low", "CI95_moderated_high"]:
            wide[f"{v}_{c}"] = records[v][c].to_numpy()
    wide["hit_normalization_count"] = wide[[f"{v}_limma_BH_hit" for v in VARIANTS]].sum(axis=1)
    wide["same_nonzero_effect_direction_all_three"] = ((wide.deposited_effect_sign != 0) &
        (wide.deposited_effect_sign == wide.qn10_effect_sign) & (wide.deposited_effect_sign == wide.qn20_effect_sign))
    union_frame = wide.loc[wide.rat_entrez.isin(union)].copy()
    assert len(union_frame) == len(union)
    save_csv(union_frame, out / "limma_union_hits_stability.csv")
    save_csv(wide, out / "limma_all_gene_normalization_stability.csv.gz")
    patterns = []
    for bits in itertools.product([0, 1], repeat=3):
        mask = np.ones(len(wide), dtype=bool)
        for v, flag in zip(VARIANTS, bits):
            mask &= wide[f"{v}_limma_BH_hit"].to_numpy() == bool(flag)
        patterns.append({"deposited_hit": bits[0], "qn10_hit": bits[1], "qn20_hit": bits[2],
                         "gene_count": int(mask.sum()),
                         "same_nonzero_direction_all_three_n": int((mask & wide.same_nonzero_effect_direction_all_three).sum())})
    pairs = []
    for a, b in itertools.combinations(VARIANTS, 2):
        common, combined = hit_sets[a] & hit_sets[b], hit_sets[a] | hit_sets[b]
        both = wide.loc[wide.rat_entrez.isin(common)]
        either = wide.loc[wide.rat_entrez.isin(combined)]
        delta = either[f"{b}_logFC"] - either[f"{a}_logFC"]
        pairs.append({"normalization_a": a, "normalization_b": b, "hits_a_n": len(hit_sets[a]), "hits_b_n": len(hit_sets[b]),
                      "intersection_n": len(common), "union_n": len(combined), "jaccard": len(common) / len(combined),
                      "fraction_a_retained_in_b": len(common) / len(hit_sets[a]),
                      "fraction_b_retained_in_a": len(common) / len(hit_sets[b]),
                      "intersection_same_nonzero_sign_n": int(((both[f"{a}_effect_sign"] == both[f"{b}_effect_sign"]) & (both[f"{a}_effect_sign"] != 0)).sum()),
                      "union_same_nonzero_sign_n": int(((either[f"{a}_effect_sign"] == either[f"{b}_effect_sign"]) & (either[f"{a}_effect_sign"] != 0)).sum()),
                      "union_median_absolute_effect_change": float(delta.abs().median()), "union_maximum_absolute_effect_change": float(delta.abs().max()),
                      "union_effect_spearman": rho(either[f"{a}_logFC"], either[f"{b}_logFC"]),
                      "all_gene_effect_spearman": rho(wide[f"{a}_logFC"], wide[f"{b}_logFC"])})

    summary_frame, quartile_frame = pd.DataFrame(summaries), pd.DataFrame(quartiles)
    save_csv(summary_frame, out / "signal_distribution_summary.csv")
    save_csv(quartile_frame, out / "signal_quartile_counts.csv")
    save_csv(pd.DataFrame(correlations), out / "effect_signal_correlations.csv")
    save_csv(pd.DataFrame(patterns), out / "hit_membership_patterns.csv")
    save_csv(pd.DataFrame(pairs), out / "normalization_pairwise_stability.csv")
    s20 = []
    for v in VARIANTS:
        full = summary_frame.loc[(summary_frame.normalization == v) & (summary_frame.population == "all_eligible_genes")].iloc[0]
        hits = summary_frame.loc[(summary_frame.normalization == v) & (summary_frame.population == "BH_hits")].iloc[0]
        qs = quartile_frame.loc[quartile_frame.normalization == v].sort_values("empirical_signal_quartile").BH_hit_n.astype(int).tolist()
        v_union = wide.loc[wide[f"{v}_limma_BH_hit"]]
        s20.append({"Normalization": LABELS[v], "All genes (n)": int(full.n), "BH hits (n)": int(hits.n),
                    "Median mean log2 signal, all genes": full.mean_log2_signal_median,
                    "Median mean log2 signal, hits": hits.mean_log2_signal_median,
                    "Median hit signal percentile": hits.signal_percentile_median,
                    "Hit counts by signal quartile Q1 / Q2 / Q3 / Q4": " / ".join(map(str, qs)),
                    "Hits shared by all three normalizations (n)": len(intersection),
                    "Hits with same nonzero effect sign in all three (n)": int(v_union.same_nonzero_effect_direction_all_three.sum())})
    save_csv(pd.DataFrame(s20), out / "Table_S20_limma_signal_summary.csv")

    caveats = [
        "This is a post-hoc descriptive audit of existing limma results; no model, normalization, p value, multiple-testing family or gene filter was changed.",
        "Mean signal is AveExpr: the arithmetic mean of gene-level log2 signals across the same ten arrays, not log2 of the arithmetic mean untransformed intensity.",
        "Empirical percentile is 100 times the number of all 14,604 same-normalization gene means <= the gene's mean, divided by 14,604; tied means share the upper empirical rank.",
        "Quartiles (0,25], (25,50], (50,75], (75,100] are exploratory descriptive strata. They are not assay detection thresholds and no gene was excluded with them.",
        "Hits in the lower signal ranks require context, but relative ranks do not establish absolute detectability; background noise and gene-specific probe performance were not calibrated.",
        "The same five control and five sesamin animals underlie all normalizations. Retention and sign concordance are analytical stability, not independent biological replication.",
        "Signal/effect correlations are descriptive, with no independence-based p values or intervals. Conditioning on BH hits is outcome selection and is not an independent test of a signal-effect relationship.",
        "Mostly high-ranked signal among hits does not establish that individual hits are noise-free, causal, normal-human responses, or clinically beneficial.",
    ]
    (out / "Table_S20_legend.txt").write_text(
        "Supplementary Table S20. Signal-intensity distribution and normalization stability of existing limma BH hits. "
        "BH q <0.05 is applied separately to the same 14,604-gene family in each previously executed normalization. "
        "Mean log2 signal is the arithmetic mean across the ten primary arrays. Percentiles use the right-continuous empirical distribution of all gene means within each normalization; "
        "Q1-Q4 denote (0,25], (25,50], (50,75] and (75,100] percentiles. These exploratory descriptive strata are not detection thresholds or filtering rules. "
        "The three analyses use the same animals. The shared hit intersection contains 24 genes; the union contains 41. Complete hit lists, all-gene ranks and effect/interval values are supplied separately.\n",
        encoding="utf-8")
    d_hits = summary_frame.loc[(summary_frame.normalization == "deposited") & (summary_frame.population == "BH_hits")].iloc[0]
    d_all = summary_frame.loc[(summary_frame.normalization == "deposited") & (summary_frame.population == "all_eligible_genes")].iloc[0]
    all_union_same_sign = bool(union_frame.same_nonzero_effect_direction_all_three.all())
    assert {v: len(hit_sets[v]) for v in VARIANTS} == {"deposited": 24, "qn10": 41, "qn20": 41}
    assert len(union) == 41 and len(intersection) == 24
    assert hit_sets["qn10"] == hit_sets["qn20"] and all_union_same_sign
    assert quartile_frame.BH_hit_n.tolist() == [1, 2, 5, 16, 2, 2, 6, 31, 2, 1, 7, 31]
    methods = (
        "As a post-hoc diagnostic, we characterized the mean-signal distribution of the existing limma BH q <0.05 results without refitting models or changing the retained gene universe. "
        "For each normalization, AveExpr was checked against the mean gene-level log2 signal across the ten primary arrays and converted to an empirical percentile among all 14,604 genes "
        "using the proportion of gene means less than or equal to that value. Descriptive quartiles summarize relative signal ranks; they are neither detection thresholds nor additional filters. "
        "We enumerated complete hit-set intersections and compared the signs and magnitudes of the existing sesamin-minus-control coefficients across normalizations. "
        "Signal-effect Spearman coefficients were descriptive and were not assigned inferential P values or confidence intervals. All normalizations use the same animals."
    )
    results = (
        f"The 24 deposited-signal limma hits were not concentrated at the low end of the retained signal distribution. Their median mean log2 signal was {d_hits.mean_log2_signal_median:.3f} "
        f"(IQR {d_hits.mean_log2_signal_q25:.3f}-{d_hits.mean_log2_signal_q75:.3f}), compared with {d_all.mean_log2_signal_median:.3f} "
        f"(IQR {d_all.mean_log2_signal_q25:.3f}-{d_all.mean_log2_signal_q75:.3f}) across all 14,604 genes. "
        "The median empirical signal percentile of the hits was 89.0%; 16 of 24 lay in the upper quartile and one in the lower quartile. "
        "QN10 and QN20 each yielded 41 hits, with median signal percentiles of 92.1%; each list contained 31 upper-quartile and two lower-quartile genes. "
        "All 24 deposited-signal hits were retained in both quantile-normalized analyses, whose 41-gene hit lists were identical. "
        f"The union contained {len(union)} genes, and {int(union_frame.same_nonzero_effect_direction_all_three.sum())} retained the same nonzero effect direction in all three analyses. "
        "These findings indicate stability of the reported hit membership and direction under the examined normalizations, while neither relative signal rank nor analytical stability establishes detectability or biological validation."
    )
    all_gene_correlations = [next(c for c in correlations if c["normalization"] == v and c["population"] == "all_eligible_genes") for v in VARIANTS]
    signed_rhos = ", ".join(f"{c['spearman_mean_signal_vs_signed_logFC']:.3f}" for c in all_gene_correlations)
    absolute_rhos = ", ".join(f"{c['spearman_mean_signal_vs_absolute_logFC']:.3f}" for c in all_gene_correlations)
    optional_effect_text = (
        "The absence of a general low-signal concentration among BH hits did not eliminate intensity-related effect patterns. "
        f"Across all 14,604 genes, descriptive Spearman correlations between mean signal and signed log2-signal differences were {signed_rhos} "
        f"for deposited, QN10 and QN20 results, respectively; correlations with absolute differences were {absolute_rhos}. "
        "These post-hoc associations were reported without inferential tests and do not identify the source of intensity dependence. "
        f"Across the 41-gene hit union, the largest absolute coefficient changes from deposited signals to QN10 and QN20 were {pairs[0]['union_maximum_absolute_effect_change']:.3f} "
        f"and {pairs[1]['union_maximum_absolute_effect_change']:.3f} log2-signal units, respectively."
    )
    (out / "proposed_methods_results_en.txt").write_text("Methods\n" + methods + "\n\nResults\n" + results + "\n\nInterpretation boundaries\n" + "\n".join("- " + c for c in caveats) + "\n", encoding="utf-8")
    (out / "optional_effect_signal_results_en.txt").write_text(optional_effect_text + "\n", encoding="utf-8")
    low_context = pd.concat([records[v].loc[records[v].limma_BH_hit & (records[v].signal_descriptive_quartile == 1),
                              ["normalization", "rat_entrez", "gene_symbol", "AveExpr", "signal_empirical_percentile", "logFC", "adj.P.Val"]]
                             for v in VARIANTS], ignore_index=True)
    save_csv(low_context, out / "lower_quartile_hits_context_not_detection_calls.csv")
    param = stage / "results/limma/parameters.json"
    inputs[str(param)] = sha(param)
    for path, digest in inputs.items():
        assert sha(path) == digest
    output_paths = sorted(p for p in out.iterdir() if p.is_file() and (p.suffix in [".csv", ".gz", ".txt"]))
    report = {
        "status": "executed and checked", "utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv, "python_executable": sys.executable, "python_version": platform.python_version(),
        "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
        "script_sha256": sha(Path(__file__)), "input_sha256": inputs,
        "source_files_unchanged_after_run": True,
        "no_model_rerun_or_new_filtering": True, "original_threshold": "BH q <0.05 separately in each existing 14,604-gene limma family",
        "checks": checks, "hit_union_n": len(union), "all_three_hit_intersection_n": len(intersection),
        "all_union_hits_same_nonzero_direction_all_three": all_union_same_sign,
        "QN10_QN20_hit_sets_identical": hit_sets["qn10"] == hit_sets["qn20"],
        "deposited_hits_subset_of_both_QN": hit_sets["deposited"] <= hit_sets["qn10"] and hit_sets["deposited"] <= hit_sets["qn20"],
        "Table_S20": s20, "signal_distribution_summaries": summaries,
        "signal_quartile_counts": quartiles, "effect_signal_correlations": correlations,
        "membership_patterns": patterns, "pairwise_stability": pairs,
        "lower_quartile_hits_descriptive_only": low_context.to_dict(orient="records"),
        "interpretation_boundaries": caveats,
        "output_sha256": {p.name: sha(p) for p in output_paths},
        "proposed_methods": methods, "proposed_results": results,
        "optional_effect_signal_results": optional_effect_text,
    }
    dump(out / "limma_signal_diagnostics.json", report)
    print(json.dumps({"status": report["status"], "hits": {v: len(hit_sets[v]) for v in VARIANTS},
                      "union": len(union), "intersection": len(intersection), "all_union_signs_match": all_union_same_sign}, ensure_ascii=True))


if __name__ == "__main__":
    main()
