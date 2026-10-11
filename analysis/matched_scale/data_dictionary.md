# Matched-scale output dictionary

All CSVs are UTF-8, comma-separated with a header. `NA` denotes unavailable data, never zero. Measured zero expression is numeric `0`. Numeric expression outputs are written with 17 significant digits; reporting to six decimals is generally sufficient. Rows retain the source order of all 25 candidates. There are no P values, confidence intervals or biological-replicate estimates.

## Library table

`outputs/matched_scale_library_all25.csv` contains 200 rows: 25 candidates x 4 technical libraries x 2 nested populations. `outputs/AREG_library.csv` contains its eight AREG rows.

| Column | Meaning |
| --- | --- |
| gsm | Original technical-library accession: GSM6322577, GSM6322578, GSM6322579 or GSM6322580. |
| condition | Original condition grouping: CTRL for 577/578, PA for 579/580. No cross-condition paired-library design is implied. |
| population | Whole organoid = all four broad categories; Hepatocytes = the author category containing five hepatocyte-like subtypes. |
| gene_id | Fixed source Ensembl identifier. |
| symbol | Fixed candidate symbol, retaining unresolved source names. |
| covered | True when the feature is measured; False for unavailable AOC4P. |
| n_cells | Retained cells in this library and population, regardless of gene coverage. This is an aggregation denominator, not biological n. |
| n_detected_cells | Count with raw UMI >0, summed from source strata; NA if unavailable. |
| detected_pct | 100 x n_detected_cells / n_cells. |
| mean_log1p_cp10k | Mean of per-cell ln(1 + 10,000 x / all-gene UMI total), including measured zeros. Whole means are reconstructed by category cell-count weighting; Hepatocytes means are copied from the source category. |
| sum_raw_counts | Candidate UMI total summed across selected source strata; audit context only, not the normalization denominator. |

The four exact source category strings are `Hepatocytes`, `Cholangiocytes`, `Fibroblasts`, and `Hepatic stellate cells` (HSC in the protocol/manuscript shorthand). The source denominator has 33,538 features; this revision reuses the prior per-cell scale without accessing raw matrices.

## Condition summary table

`outputs/matched_scale_summary_all25.csv` has 100 rows: 25 candidates x 2 populations x 2 weighting schemes. `outputs/AREG_summary.csv` has its four AREG rows.

| Column | Meaning |
| --- | --- |
| gene_id, symbol, covered, population | Same definitions as above. |
| weighting | `equal_library` is primary: each of two condition libraries has weight 1/2. `pooled_cell` is a sensitivity: a library has weight n_cells / sum of condition cells for that population. |
| control_mean, PA_mean | Weighted arithmetic means of per-library mean_log1p_cp10k. They are not logarithms of pooled counts or of average CP10k. |
| control_n_libraries, PA_n_libraries | Number of technical libraries retained in each condition (always 2); not independent biological repeats. |
| control_n_cells, PA_n_cells | Retained-cell totals for the relevant population and condition; whole 8,603/15,431 and Hepatocytes 6,283/11,814. |
| PA_minus_control | PA_mean - control_mean, in mean log1p(CP10k) units. |
| direction | positive, negative, zero, or unavailable, based on the unrounded numeric difference. |

The whole-cell library denominators are 3,093; 5,510; 7,399; 8,032. The corresponding Hepatocytes denominators are 2,405; 3,878; 5,958; 5,856. These denominators make the pooled-cell weighting explicit and reconstructible. Each population uses the same four library identifiers; its pooled weights differ because the sampled populations differ.

## Paired population accounting

`outputs/paired_population_library_all25.csv` contains 100 same-library pairs. Columns shared by the merge are gsm, condition, gene_id, symbol and covered; the suffixes `_whole` and `_hep` mark the two populations for the library-table fields. `whole_minus_hepatocyte` is the same-library whole mean minus Hepatocytes mean. The populations are nested, so this difference has no independent-samples interpretation. This table does not pair CTRL and PA libraries.

`outputs/population_contrast_all25.csv` contains 50 rows (25 x 2 weighting schemes). `outputs/AREG_population_contrast.csv` has two AREG rows. `whole_PA_minus_control` and `hepatocyte_PA_minus_control` are the two population-specific contrasts. `whole_minus_hepatocyte_delta` is their arithmetic difference. `whole_direction` and `hepatocyte_direction` retain their individual signs; `same_direction` is True or False for covered candidates and NA for unavailable candidates. This difference does not isolate a causal composition effect.

## Validation and execution files

- `protocol_freeze.json`: immutable first freeze time, protocol hash, source and input-copy hashes; created before new expression arithmetic.
- `outputs/validation.json`: all programmatic checks, tolerances, counts, maximum residuals and required caveats. Multiple assertions reflect table entries and combinations, not independent validation experiments.
- `outputs/execution.json`: actual UTC start/end, command, interpreter, platform, dependency versions, output directory and no-raw-rerun statement.
- `outputs/hash_manifest.json`: hashes of executed code, frozen material, copied data, numerical outputs and execution/validation records. The manifest does not hash itself.
- `package_manifest.json`: final package hashes including explanatory documents; excludes itself. Re-running the analysis will change time-stamped execution files, so refresh this final-package snapshot if redistributing a rerun.

All original inputs remain unchanged. The numerical scope is aggregation and reconciliation of existing summaries; raw per-cell normalization and barcode linkage are inherited provenance, not newly validated operations.
