"""Reproduce the frozen descriptive matched-scale analysis from copied summaries.

Usage: python run.py [--output-dir PATH]
Requires Python >=3.10, numpy and pandas. No raw-cell renormalization is run.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from decimal import Decimal, getcontext
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TOL = 1e-12
LIBRARIES = {"GSM6322577": "CTRL", "GSM6322578": "CTRL", "GSM6322579": "PA", "GSM6322580": "PA"}
CATEGORIES = {"Hepatocytes", "Cholangiocytes", "Hepatic stellate cells", "Fibroblasts"}
POPULATIONS = ("Whole organoid", "Hepatocytes")
NUMERIC = ("n_detected_cells", "detected_pct", "mean_log1p_cp10k", "sum_raw_counts")
CHECKS = []


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(name, passed, **details):
    CHECKS.append({"check": name, "passed": bool(passed), **details})
    if not passed:
        raise AssertionError(f"{name}: {details}")


def close(name, observed, expected):
    observed = np.asarray(observed, dtype=float)
    expected = np.asarray(expected, dtype=float)
    valid = np.isfinite(observed) & np.isfinite(expected)
    missing_match = np.array_equal(np.isnan(observed), np.isnan(expected))
    error = float(np.max(np.abs(observed[valid] - expected[valid]))) if valid.any() else 0.0
    check(name, observed.shape == expected.shape and missing_match and error <= TOL,
          n_compared=int(valid.sum()), max_absolute_error=error, tolerance=TOL)


def read(name):
    return pd.read_csv(HERE / "inputs" / name, float_precision="round_trip")


def sign(value):
    return "unavailable" if pd.isna(value) else "positive" if value > 0 else "negative" if value < 0 else "zero"


def json_write(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def execute(output):
    started = datetime.now(timezone.utc).isoformat()
    freeze = json.loads((HERE / "protocol_freeze.json").read_text(encoding="utf-8"))
    check("frozen_protocol_hash", sha(HERE / "protocol.md") == freeze["protocol_sha256"])
    for item in freeze["inputs"]:
        check(f"frozen_input_hash:{item['copy']}", sha(HERE / item["copy"]) == item["sha256"])
    data = read("library_category_candidates_all25.csv")
    direct = read("library_candidates_all25.csv")
    counts = read("category_sampling_counts.csv")
    kit = read("kitagawa_summary_all25.csv")
    sensitivity = read("weighting_sensitivity_all25.csv")
    symbols = data.symbol.drop_duplicates().tolist()
    check("complete_25_candidates", len(symbols) == 25 and "AREG" in symbols and "AOC4P" in symbols)
    check("unique_category_gene_keys", not data.duplicated(["gsm", "celltype_category", "symbol"]).any())
    check("complete_category_grid", len(data) == 400 and set(data.celltype_category) == CATEGORIES
          and all(len(g) == 25 for _, g in data.groupby(["gsm", "celltype_category"])))
    check("fixed_library_condition_map", data[["gsm", "condition"]].drop_duplicates().set_index("gsm").condition.to_dict() == LIBRARIES)
    check("consistent_feature_and_coverage", all(g.gene_id.nunique() == 1 and g.covered.nunique() == 1 for _, g in data.groupby("symbol")))
    check("only_AOC4P_unavailable", set(data.loc[~data.covered, "symbol"]) == {"AOC4P"})
    check("missing_feature_not_zero", data.loc[~data.covered, list(NUMERIC)].isna().all().all())
    measured = data[data.covered]
    check("valid_measured_values", np.isfinite(measured[list(NUMERIC)]).all().all()
          and (measured[list(NUMERIC)] >= 0).all().all()
          and (measured.n_detected_cells <= measured.n_cells).all())
    close("source_detection_percent", measured.detected_pct, 100 * measured.n_detected_cells / measured.n_cells)
    check("consistent_category_denominators", all(g.n_cells.nunique() == 1 and g.n_cells.iloc[0] > 0
          for _, g in data.groupby(["gsm", "celltype_category"])))
    denominator = data[["gsm", "condition", "celltype_category", "n_cells"]].drop_duplicates()
    counted = denominator.merge(counts, on=["gsm", "condition", "celltype_category"], validate="one_to_one", suffixes=("", "_reference"))
    check("category_count_reference", len(counted) == 16 and (counted.n_cells == counted.n_cells_reference).all())
    check("retained_cell_totals", denominator.n_cells.sum() == 24034
          and denominator.loc[denominator.celltype_category == "Hepatocytes", "n_cells"].sum() == 18097,
          all_cells=int(denominator.n_cells.sum()), hepatocyte_cells=int(denominator.loc[denominator.celltype_category == "Hepatocytes", "n_cells"].sum()))
    close("category_fraction_counts", counted.fraction, counted.n_cells / counted.denominator_n_cells)

    lib_records = []
    for symbol in symbols:
        for gsm, condition in LIBRARIES.items():
            g = data[(data.symbol == symbol) & (data.gsm == gsm)]
            for population in POPULATIONS:
                selected = g if population == "Whole organoid" else g[g.celltype_category == "Hepatocytes"]
                covered = bool(selected.covered.iloc[0])
                n = int(selected.n_cells.sum())
                record = dict(gsm=gsm, condition=condition, population=population, gene_id=selected.gene_id.iloc[0],
                              symbol=symbol, covered=covered, n_cells=n)
                if covered:
                    detected = int(selected.n_detected_cells.sum())
                    record.update(n_detected_cells=detected, detected_pct=100 * detected / n,
                                  mean_log1p_cp10k=float(np.dot(selected.n_cells, selected.mean_log1p_cp10k) / n),
                                  sum_raw_counts=int(selected.sum_raw_counts.sum()))
                else:
                    record.update({name: np.nan for name in NUMERIC})
                lib_records.append(record)
    library = pd.DataFrame(lib_records)
    check("library_output_grid", len(library) == 200 and not library.duplicated(["gsm", "symbol", "population"]).any())
    whole = library[library.population == "Whole organoid"].merge(direct, on=["gsm", "condition", "gene_id", "symbol", "covered"], validate="one_to_one", suffixes=("", "_reference"))
    check("direct_whole_library_join", len(whole) == 100)
    for name in ["n_cells", *NUMERIC]:
        close(f"independent_source_whole_library:{name}", whole[name], whole[name + "_reference"])

    summary_records = []
    for symbol in symbols:
        for population in POPULATIONS:
            p = library[(library.symbol == symbol) & (library.population == population)]
            for weighting in ("equal_library", "pooled_cell"):
                r = dict(gene_id=p.gene_id.iloc[0], symbol=symbol, covered=bool(p.covered.iloc[0]), population=population, weighting=weighting)
                for condition, prefix in [("CTRL", "control"), ("PA", "PA")]:
                    g = p[p.condition == condition]
                    check(f"two_libraries:{symbol}:{population}:{weighting}:{condition}", len(g) == 2)
                    weights = np.full(len(g), 1 / len(g)) if weighting == "equal_library" else g.n_cells.to_numpy() / g.n_cells.sum()
                    close(f"weights_sum:{symbol}:{population}:{weighting}:{condition}", weights.sum(), 1.0)
                    r[prefix + "_mean"] = float(np.dot(weights, g.mean_log1p_cp10k)) if r["covered"] else np.nan
                    r[prefix + "_n_libraries"] = len(g)
                    r[prefix + "_n_cells"] = int(g.n_cells.sum())
                r["PA_minus_control"] = r["PA_mean"] - r["control_mean"]
                r["direction"] = sign(r["PA_minus_control"])
                summary_records.append(r)
    summary = pd.DataFrame(summary_records)
    check("summary_output_grid", len(summary) == 100 and not summary.duplicated(["symbol", "population", "weighting"]).any())
    for weighting, scheme, refcol in [("equal_library", "equal_library_actual", "equal_library_total"), ("pooled_cell", "pooled_cell_actual", "pooled_cell_total")]:
        g = summary[(summary.population == "Whole organoid") & (summary.weighting == weighting)]
        k = g.merge(kit[kit.scheme == scheme], on=["symbol", "gene_id", "covered"], validate="one_to_one")
        check(f"kitagawa_all25_join:{weighting}", len(k) == 25)
        for output_col, input_col in [("control_mean", "control_profile_mean"), ("PA_mean", "PA_profile_mean"), ("PA_minus_control", "total_PA_minus_control")]:
            close(f"prior_kitagawa:{weighting}:{output_col}", k[output_col], k[input_col])
        s = g.merge(sensitivity, on=["symbol", "covered"], validate="one_to_one")
        check(f"weighting_reference_all25_join:{weighting}", len(s) == 25)
        close(f"prior_weighting_sensitivity:{weighting}", s.PA_minus_control, s[refcol])

    # Independent CSV/Decimal implementation bypasses pandas grouping and numpy dot products.
    getcontext().prec = 50
    with (HERE / "inputs" / "library_category_candidates_all25.csv").open(encoding="utf-8", newline="") as f:
        original_rows = list(csv.DictReader(f))
    decimal_library = {}
    decimal_errors = []
    for row in library.itertuples(index=False):
        selected = [r for r in original_rows if r["gsm"] == row.gsm and r["symbol"] == row.symbol
                    and (row.population == "Whole organoid" or r["celltype_category"] == "Hepatocytes")]
        dn = sum(Decimal(r["n_cells"]) for r in selected)
        check(f"decimal_library_count:{row.symbol}:{row.gsm}:{row.population}", int(dn) == row.n_cells)
        dm = sum(Decimal(r["n_cells"]) * Decimal(r["mean_log1p_cp10k"]) for r in selected) / dn if row.covered else None
        decimal_library[(row.symbol, row.gsm, row.population)] = (dn, dm)
        if dm is not None:
            decimal_errors.append(abs(float(dm) - row.mean_log1p_cp10k))
    close("independent_decimal_library", decimal_errors, np.zeros(len(decimal_errors)))
    decimal_summary_errors = []
    for row in summary.itertuples(index=False):
        if not row.covered:
            continue
        dm = {}
        for condition in ("CTRL", "PA"):
            vals = [decimal_library[(row.symbol, gsm, row.population)] for gsm, t in LIBRARIES.items() if t == condition]
            dm[condition] = sum(v[1] for v in vals) / Decimal(len(vals)) if row.weighting == "equal_library" else sum(n * m for n, m in vals) / sum(n for n, _ in vals)
        decimal_summary_errors.extend([abs(float(dm["CTRL"]) - row.control_mean), abs(float(dm["PA"]) - row.PA_mean), abs(float(dm["PA"] - dm["CTRL"]) - row.PA_minus_control)])
    close("independent_decimal_summary", decimal_summary_errors, np.zeros(len(decimal_summary_errors)))

    paired = library[library.population == "Whole organoid"].merge(library[library.population == "Hepatocytes"],
                on=["gsm", "condition", "gene_id", "symbol", "covered"], suffixes=("_whole", "_hep"), validate="one_to_one")
    paired["whole_minus_hepatocyte"] = paired.mean_log1p_cp10k_whole - paired.mean_log1p_cp10k_hep
    check("paired_nested_population_rows", len(paired) == 100 and (paired.n_cells_whole >= paired.n_cells_hep).all())
    pairing_errors = []
    comparisons = []
    for symbol in symbols:
        for weighting in ("equal_library", "pooled_cell"):
            g = summary[(summary.symbol == symbol) & (summary.weighting == weighting)].set_index("population")
            w, h = g.loc["Whole organoid"], g.loc["Hepatocytes"]
            comparisons.append(dict(gene_id=w.gene_id, symbol=symbol, covered=bool(w.covered), weighting=weighting,
                     whole_PA_minus_control=w.PA_minus_control, hepatocyte_PA_minus_control=h.PA_minus_control,
                     whole_minus_hepatocyte_delta=w.PA_minus_control-h.PA_minus_control,
                     whole_direction=w.direction, hepatocyte_direction=h.direction,
                     same_direction=(w.direction == h.direction) if w.covered else pd.NA))
            if w.covered and weighting == "equal_library":
                pair_means = {}
                for condition, col in [("CTRL", "control_mean"), ("PA", "PA_mean")]:
                    p = paired[(paired.symbol == symbol) & (paired.condition == condition)]
                    pair_means[condition] = float(p.whole_minus_hepatocyte.mean())
                    pairing_errors.append(abs(pair_means[condition] - (w[col] - h[col])))
                pairing_errors.append(abs((pair_means["PA"] - pair_means["CTRL"]) - (w.PA_minus_control-h.PA_minus_control)))
    close("paired_equal_library_identities", pairing_errors, np.zeros(len(pairing_errors)))
    comparison = pd.DataFrame(comparisons)
    areg = summary[summary.symbol == "AREG"]
    areg_primary = areg[(areg.population == "Whole organoid") & (areg.weighting == "equal_library")].iloc[0]
    check("reported_AREG_whole_reconstruction", f"{areg_primary.PA_minus_control:+.6f}" == "+0.054898", actual=float(areg_primary.PA_minus_control))
    for name, frame in [("library", library), ("summary", summary)]:
        check(f"AOC4P_preserved:{name}", len(frame[frame.symbol == "AOC4P"]) == (8 if name == "library" else 4)
              and frame.loc[frame.symbol == "AOC4P", "mean_log1p_cp10k" if name == "library" else "PA_minus_control"].isna().all())
    check("measured_source_zeros_retained", (measured.mean_log1p_cp10k == 0).sum() > 0
          and (library[library.covered].mean_log1p_cp10k == 0).sum() > 0,
          measured_source_category_zeros=int((measured.mean_log1p_cp10k == 0).sum()),
          measured_output_library_zeros=int((library[library.covered].mean_log1p_cp10k == 0).sum()))

    output.mkdir(parents=True, exist_ok=True)
    tables = {"matched_scale_library_all25.csv": library, "matched_scale_summary_all25.csv": summary,
              "paired_population_library_all25.csv": paired, "population_contrast_all25.csv": comparison,
              "AREG_library.csv": library[library.symbol == "AREG"], "AREG_summary.csv": areg,
              "AREG_population_contrast.csv": comparison[comparison.symbol == "AREG"]}
    for name, table in tables.items():
        table.to_csv(output / name, index=False, na_rep="NA", float_format="%.17g")
    report = {"assessment": "Share with caveats", "all_checks_passed": all(c["passed"] for c in CHECKS),
              "number_of_checks": len(CHECKS), "checks": CHECKS,
              "raw_count_rerun_performed": False, "scope": "Aggregation of previously computed category summaries; independent Decimal arithmetic and reconciliation to prior direct-library outputs.",
              "required_caveats": ["One differentiation; four technical libraries and cells are not biological replicates.",
                "No sesamin intervention, P values or biological confidence intervals.",
                "Normalization provenance inherited from source; raw matrices and retained-cell identifiers were not reprocessed or revalidated.",
                "Nested populations are summarized within the same library; no cross-condition library pairing or isolated causal composition claim.",
                "AOC4P remains unavailable rather than zero; all 25 candidates and contrary results retained."]}
    json_write(output / "validation.json", report)
    finished = datetime.now(timezone.utc).isoformat()
    execution = {"started_utc": started, "finished_utc": finished, "protocol_frozen_utc": freeze["frozen_utc"],
                 "command": [sys.executable, *sys.argv], "working_directory": str(Path.cwd()),
                 "python_executable": sys.executable, "python_version": platform.python_version(),
                 "platform": platform.platform(), "dependencies_used": {"numpy": np.__version__, "pandas": pd.__version__},
                 "randomness": "None; deterministic aggregation, no random seed needed", "raw_count_rerun_performed": False,
                 "outputs_directory": str(output), "validation_passed": True}
    json_write(output / "execution.json", execution)
    assets = [HERE / "protocol.md", HERE / "protocol_freeze.json", HERE / "freeze_protocol.py", HERE / "run.py"]
    assets += sorted((HERE / "inputs").glob("*.csv")) + sorted(output.glob("*.csv")) + [output / "validation.json", output / "execution.json"]
    for optional in ["data_dictionary.md", "README.md"]:
        if (HERE / optional).exists():
            assets.append(HERE / optional)
    json_write(output / "hash_manifest.json", {"algorithm": "SHA-256", "created_utc": finished,
            "files": [{"path": str(p.relative_to(HERE)) if p.is_relative_to(HERE) else str(p), "bytes": p.stat().st_size, "sha256": sha(p)} for p in assets]})
    print(areg.to_string(index=False))
    print(json.dumps({"all_checks_passed": True, "checks": len(CHECKS), "outputs": str(output)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE / "outputs")
    args = parser.parse_args()
    execute(args.output_dir.resolve())
