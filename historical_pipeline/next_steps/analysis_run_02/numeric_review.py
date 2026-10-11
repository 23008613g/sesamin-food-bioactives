"""Read-only source validation; writes only this review's numeric_review.json.

Run using the project .venv. No downloads, installs, raw-source edits, or models.
Quantization scenarios are hypothetical numerical checks, not observed data.
"""
from datetime import datetime, timezone
from fractions import Fraction
from itertools import product
from pathlib import Path
import csv
import gzip
import hashlib
import importlib.metadata
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bh(p):
    """Independent BH implementation; all submitted tests retain their slots."""
    order = np.argsort(p)
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q = np.empty(len(p))
    q[order] = np.minimum(adjusted, 1.0)
    return q


def rational_welch(a, b):
    """Exact rational means/variances/df, with SciPy t CDF/quantile at the end."""
    n, m = len(a), len(b)
    ma, mb = sum(a) / n, sum(b) / m
    va = sum((z - ma) ** 2 for z in a) / (n - 1)
    vb = sum((z - mb) ** 2 for z in b) / (m - 1)
    se2 = va / n + vb / m
    effect = ma - mb
    df = se2 ** 2 / ((va / n) ** 2 / (n - 1) + (vb / m) ** 2 / (m - 1))
    se = float(se2) ** 0.5
    tv = float(effect) / se
    p = float(2 * stats.t.sf(abs(tv), float(df)))
    width = float(stats.t.ppf(.975, float(df))) * se
    return {"effect": float(effect), "variance_treated": float(va),
            "variance_reference": float(vb), "SE": se, "df": float(df),
            "t": tv, "p": p,
            "CI95": [float(effect) - width, float(effect) + width]}


def main():
    matrix = BASE / "downloads/GSE200765_series_matrix.txt.gz"
    lines = gzip.open(matrix, "rt", encoding="utf-8").read().splitlines()
    begin = lines.index("!series_matrix_table_begin")
    end = lines.index("!series_matrix_table_end")
    records = list(csv.reader(lines[begin + 1:end], delimiter="\t"))
    samples = records[0][1:]
    raw = [r for r in records[1:] if r[0].startswith("TC") and ".hg." in r[0]]
    ids = [r[0] for r in raw]
    x = np.array([[float(v) for v in r[1:]] for r in raw])
    manifest_path = BASE / "sample_manifest.csv"
    manifest = pd.read_csv(manifest_path).set_index("gsm")
    groups = np.array([manifest.loc[s, "group"] for s in samples])
    titles = next(list(csv.reader([line], delimiter="\t"))[0][1:]
                  for line in lines[:begin] if line.startswith("!Sample_title"))
    contract = json.loads((HERE / "analysis_contract.json").read_text())
    original_summary = json.loads((HERE / "analysis_summary.json").read_text())
    assert x.shape == (19911, 14)
    assert np.isfinite(x).all() and len(set(ids)) == len(ids)
    evidence = []
    for a, b in contract["contrasts"]:
        name = a + "_vs_" + b
        ix = np.flatnonzero(groups == a)
        iy = np.flatnonzero(groups == b)
        xx, yy = x[:, ix], x[:, iy]
        n, m = xx.shape[1], yy.shape[1]
        table_path = HERE / (name + ".csv")
        out = pd.read_csv(table_path, keep_default_na=False)
        assert out.feature_id.tolist() == ids
        assert (out.n_treated == n).all() and (out.n_reference == m).all()
        effect = xx.mean(1) - yy.mean(1)
        # Centering first at a real observation makes constant rows exactly zero.
        dx, dy = xx - xx[:, [0]], yy - yy[:, [0]]
        vx = np.sum((dx - dx.mean(1, keepdims=True)) ** 2, axis=1) / (n - 1)
        vy = np.sum((dy - dy.mean(1, keepdims=True)) ** 2, axis=1) / (m - 1)
        ax, ay = vx / n, vy / m
        se = np.sqrt(ax + ay)
        with np.errstate(divide="ignore", invalid="ignore"):
            df = (ax + ay) ** 2 / (ax ** 2 / (n - 1) + ay ** 2 / (m - 1))
            tv = effect / se
            p = 2 * stats.t.sf(np.abs(tv), df)
        valid = np.isfinite(p) & (se > 0)
        p_for_bh = np.where(valid, p, 1.0)
        q = bh(p_for_bh)
        width = stats.t.ppf(.975, df) * se
        constant_a, constant_b = np.ptp(xx, axis=1) == 0, np.ptp(yy, axis=1) == 0
        constant_either = constant_a | constant_b
        reported_hits = out.q_BH_this_contrast.to_numpy(float) < .05
        def maxdiff(col, calculated):
            actual = out[col].to_numpy(float)
            finite = np.isfinite(actual) & np.isfinite(calculated)
            return float(np.max(np.abs(actual[finite] - calculated[finite])))
        differences = {col: maxdiff(col, calc) for col, calc in [
            ("mean_treated", xx.mean(1)), ("mean_reference", yy.mean(1)),
            ("mean_signal_difference", effect), ("SE_assuming_independent_samples", se),
            ("welch_df", df), ("welch_t", tv), ("p_two_sided", p),
            ("q_BH_this_contrast", q), ("CI95_low_pointwise", effect - width),
            ("CI95_high_pointwise", effect + width)]}
        assert max(differences.values()) < 1e-10
        sensitivity = []
        all_hit_max_p = p_for_bh.copy()
        for j in np.flatnonzero(reported_hits):
            aa = [Fraction(raw[j][int(k) + 1]) for k in ix]
            bb = [Fraction(raw[j][int(k) + 1]) for k in iy]
            exact = rational_welch(aa, bb)
            scenarios = []
            # 64 corners for each 3-versus-3 hit; all fit inside +/-0.005.
            delta = Fraction(49, 10000)
            for signs in product((-1, 1), repeat=n + m):
                ta = [v + sign * delta for v, sign in zip(aa, signs[:n])]
                tb = [v + sign * delta for v, sign in zip(bb, signs[n:])]
                value = rational_welch(ta, tb)
                altered_p = p_for_bh.copy()
                altered_p[j] = value["p"]
                scenarios.append({"p": value["p"], "q": float(bh(altered_p)[j]),
                                  "signs": signs, "treated": [float(z) for z in ta],
                                  "reference": [float(z) for z in tb]})
            high_p = max(scenarios, key=lambda z: z["p"])
            all_hit_max_p[j] = high_p["p"]
            sensitivity.append({"feature_id": ids[j], "symbol": str(out.iloc[j].gene_symbols),
                "submitted_treated": [float(z) for z in aa],
                "submitted_reference": [float(z) for z in bb],
                "exact_rational_moments": exact, "original_q": float(q[j]),
                "constant_either_group": bool(constant_either[j]),
                "hypothetical_scenarios": len(scenarios),
                "corner_p_range": [min(z["p"] for z in scenarios), max(z["p"] for z in scenarios)],
                "corner_q_range_one_feature_changed": [min(z["q"] for z in scenarios), max(z["q"] for z in scenarios)],
                "max_p_scenario": high_p})
        q_joint = bh(all_hit_max_p)
        for hit, j in zip(sensitivity, np.flatnonzero(reported_hits)):
            hit["q_if_each_of_four_reported_hits_uses_its_max_p_corner"] = float(q_joint[j])
        evidence.append({"contrast": name, "csv_sha256": sha(table_path),
            "n_treated": n, "n_reference": m, "rows": len(out), "BH_family_size": len(p_for_bh),
            "FDR_scope": "This contrast only; no experiment-wide four-contrast control claimed",
            "max_absolute_differences": differences,
            "valid_rows_stable_recalculation": int(valid.sum()),
            "reported_valid_but_stably_undefined": int(np.count_nonzero(out.test_defined.to_numpy(bool) & ~valid)),
            "constant_treated_rows": int(constant_a.sum()), "constant_reference_rows": int(constant_b.sum()),
            "constant_either_rows": int(constant_either.sum()),
            "constant_both_rows": int((constant_a & constant_b).sum()),
            "constant_either_feature_ids": [ids[j] for j in np.flatnonzero(constant_either)],
            "smallest_q_constant_either": float(q[constant_either].min()),
            "q_lt_005_reported": int(reported_hits.sum()),
            "q_lt_005_stable_recalculation": int(np.count_nonzero(valid & (q < .05))),
            "significant_hits_and_quantization_sensitivity": sensitivity})
    review = {"checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "skill_used": "validate-data", "interpreter": sys.executable,
        "execution_note": "First launch stopped at a syntax error in this review script before any data processing. The extra bracket was corrected, then the complete review was executed.",
        "package_versions": {m: importlib.metadata.version(m) for m in ["numpy", "scipy", "pandas"]},
        "review_script_sha256": sha(Path(__file__)),
        "reviewed_main_script_sha256": sha(HERE / "run_bulk_analysis.py"),
        "main_summary_sha256": sha(HERE / "analysis_summary.json"),
        "source_sha256": sha(matrix), "manifest_sha256": sha(manifest_path),
        "source_hash_matches_original_run": sha(matrix) == original_summary["source_sha256"],
        "contract_hash_matches_original_run": sha(HERE / "analysis_contract.json") == original_summary["contract_sha256"],
        "main_script_hash_matches_original_run": sha(HERE / "run_bulk_analysis.py") == original_summary["script_sha256"],
        "main_script_hash_note": "Parent reports an in-progress change to figure-export/QA calls after the original run; final main rerun must refresh provenance.",
        "shape": list(x.shape),
        "all_submitted_values_representable_at_2_decimal_places": bool(np.allclose(x * 100, np.rint(x * 100), atol=1e-10, rtol=0)),
        "sample_mapping": [{"gsm": s, "group": str(g), "source_title": t} for s, g, t in zip(samples, groups, titles)],
        "sample_mapping_review": "All 14 titles and manifest assignments agree; active-vitamin-D arm has n=2; other arms n=3.",
        "method": "Read original gzip with independent CSV parser. Independently recompute stable group variances, Welch-Satterthwaite df, t, and hand-coded BH. SciPy supplies t CDF/quantiles. Exact rational moments independently checked for the four reported q<.05 rows.",
        "quantization_method": "Hypothetical +/-0.0049 perturbations at every corner for the six values of each reported 3-vs-3 hit, compatible with values rounded to two decimals. Other features are unchanged for single-feature q ranges. A separate joint scenario changes all four reported hits to their own max-p corner. These are deterministic sensitivity scenarios, not added/observed samples, probabilities, corrected values, or exhaustive global extrema over continuous rounding intervals.",
        "contrasts": evidence,
        "limited_findings": [
            "No material mismatch in mapping, Welch df, pointwise CI, p, or within-contrast BH. This does not validate independent cultures, normality, the source limma analysis, or biological effects.",
            "Precision-loss warnings need localized constant/near-constant-group flags for traceability. Constant rows have one constant group only; all four q<.05 rows have nonzero variance in both groups. Recalculation does not change the submitted-value result.",
            "HTRA1, ALG1 and DEFB115 each cross q=.05 under at least one checked, hypothetical within-rounding scenario. HP remains below .05 in the checked scenarios only; this is not a proof of continuous global robustness or biological replication. Preserve the submitted-value results and label this sensitivity; do not switch methods or exclude features to seek significance."
        ],
        "scope_limits": ["No search or download", "No raw-source modification", "No model/API calls", "No package install", "No biological/independent validation", "No evidence from protein assays treated as RNA replication", "No rendered-figure review", "No experiment-wide FDR guarantee across four contrasts"]}
    (HERE / "numeric_review.json").write_text(json.dumps(review, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"file": str(HERE / "numeric_review.json"),
          "max_abs_difference": max(max(c["max_absolute_differences"].values()) for c in evidence),
          "constant_either": {c["contrast"]: c["constant_either_rows"] for c in evidence},
          "hits": [h for c in evidence for h in c["significant_hits_and_quantization_sensitivity"]]}, ensure_ascii=True))


if __name__ == "__main__":
    main()
