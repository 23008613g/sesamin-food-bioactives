# PFHN P6 matched-scale analysis

The completed matched-scale comparison gives positive AREG PA-minus-control differences in both the whole organoid (+0.054898) and author-labelled Hepatocytes (+0.039797), using equal library weights and mean per-cell log1p(CP10k). Pooled-cell sensitivity remains positive for both (+0.063003 and +0.037042). This result must be retained when interpreting the older negative hepatocyte pooled-UMI CPM ratios.

## Reproduction

From this package directory, run `python run.py` using an environment with numpy and pandas. Paths are anchored to the script and only local `inputs/` copies are required. `--output-dir PATH` optionally writes a separate result directory. The executable used in this project was `research_planning/hepatic_lipotoxicity_2026-10-07/.venv/Scripts/python.exe`, Python 3.13.9, numpy 2.4.6 and pandas 3.0.3. No packages were installed. The run is deterministic; no random seed is needed.

The executed workspace command was:

```powershell
& 'research_planning/hepatic_lipotoxicity_2026-10-07/.venv/Scripts/python.exe' 'research_planning/hepatic_lipotoxicity_2026-10-07/next_steps/pfhn_revision_26/analysis/matched_scale_areg/run.py'
```

`freeze_protocol.py` documents the initial copy/hash step and is not needed to reproduce results. It deliberately refuses to overwrite the existing freeze. `protocol.md` and `protocol_freeze.json` are immutable historical records; `run.py` verifies their hashes and all five input-copy hashes before calculating. Freeze occurred at 2026-10-10 05:03:04.968716 UTC, before the new matched-scale expression calculation.

## Delivered analysis and validation

The two main CSVs are `outputs/matched_scale_library_all25.csv` (200 rows) and `outputs/matched_scale_summary_all25.csv` (100 rows). They retain all 25 candidates, including NA AOC4P, all four technical libraries, measured zeros and contrary directions. Compact AREG tables and paired nested-population accounting are also provided. `data_dictionary.md` defines every table and unit; `methods_results_en_zh.md` supplies manuscript-ready wording for author review.

All implemented assertions passed. Cell counts reconcile exactly to 24,034 overall and 18,097 Hepatocytes. Reconstructed whole-library candidate means agree with the pre-existing directly aggregated whole-library table to at most 1.11e-16; both whole-organoid weightings agree with prior Kitagawa and weighting-sensitivity totals to at most 3.34e-16. Independently coded Decimal arithmetic (50-digit precision, CSV input strings) agrees with the floating-point library and condition calculations within 2.36e-16. Paired equal-library identities agree within 1.67e-16. The check tolerance is 1e-12. Original files were read-only and their copied hashes match.

The `validate-data` skill guided the question, denominator, missingness and aggregation checks. These are internal arithmetic and source-output reconciliations, not independent experimental validation. No raw matrices, per-cell all-gene denominator calculation, retained-cell barcode matching, normalization or annotation were rerun. Four technical libraries derive from one differentiation and cannot supply biological-replicate inference. The data contain no sesamin intervention. No P values, confidence intervals, plots or manuscript edits were produced in this analysis package.

## Interpretation boundary

Matching per-cell expression scale and equal library weighting yields the same positive AREG direction in whole organoids and Hepatocytes. Therefore, the old positive-whole versus negative-hepatocyte-CPM contrast cannot be treated as a population-only sign reversal. The present result does not apportion the old discrepancy uniquely among normalization, transformation and UMI weighting, nor isolate a causal composition contribution. Same-library whole/hepatocyte values are nested; CTRL/PA libraries are not paired biological samples.
