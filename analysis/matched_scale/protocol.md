# Frozen protocol: matched-scale organoid population comparison

Protocol version: 1.0. Date: 2026-10-10. The freeze time and SHA-256 of this protocol and its input copies are recorded in `protocol_freeze.json` before running `run.py`. This is a local prospective analysis specification, not a registry preregistration. P6 of the PFHN targeted review and the prior manuscript supplied the question and existing reference values; new matched-scale hepatocyte results have not been computed at protocol freeze.

## Question and scope

Using the same author-retained GSE207889 cells, compare the PA-minus-control AREG expression difference for the whole organoid and author category `Hepatocytes` on the same per-cell scale and with the same condition-library weights. AREG is highlighted; all 25 previously fixed candidates, including missing AOC4P, remain in every output. No new selection, annotation, normalization, cell filtering, gene filtering, raw-count processing, inferential tests or confidence intervals are permitted.

The four technical libraries are GSM6322577/GSM6322578 (CTRL) and GSM6322579/GSM6322580 (PA), from one differentiation. They contain 24,034 author-retained cells overall and 18,097 cells in the author-labelled Hepatocytes category (five hepatocyte-like subtypes). Cells and technical libraries are not independent biological replicates. No sesamin intervention is present.

## Inputs and provenance

All original inputs are read-only. Byte-identical minimal copies in `inputs/` are used for the portable calculation:

- `library_category_candidates_all25.csv`: existing per-library, four-category candidate means and counts; main calculation input.
- `library_candidates_all25.csv`: pre-existing directly aggregated whole-library summaries; independent aggregation-path reconciliation target.
- `category_sampling_counts.csv`: pre-existing category and library cell-count reconciliation target.
- `kitagawa_summary_all25.csv`: pre-existing whole-organoid equal-library and pooled-cell totals and condition means; reconciliation target, not a source for the new result.
- `weighting_sensitivity_all25.csv`: pre-existing whole-organoid differences and direction comparisons; reconciliation target.

Original paths and SHA-256 values, byte-identical input-copy hashes, source manuscript hash and review-plan hash are recorded at freeze. No expression result calculated by this revision is used to decide the specification. Existing known whole-organoid AREG reference: +0.054898 mean log1p(CP10k) units when printed to six decimals.

## Definitions and arithmetic

For gene g and retained cell c, the already computed source scale is y_gc = ln(1 + 10,000 x_gc / T_c), where x_gc is raw gene UMI count and T_c is that cell's raw total across all 33,538 source features. Means include measured zero-count cells. This revision consumes existing stratum means; it does not re-execute or independently verify raw-cell normalization. No candidate-only denominator is introduced.

Let n_lk be cells in library l and one of the four disjoint author categories k: Hepatocytes, Cholangiocytes, HSC, Fibroblasts. Let m_glk be the supplied mean y. For each library:

- Whole organoid: N_l = sum_k n_lk; M_g,l,whole = sum_k(n_lk m_glk) / N_l.
- Hepatocytes: N_l,hep = n_l,hep; M_g,l,hep = m_g,l,hep.

For population p and condition t in {CTRL, PA}, with L_t = 2 libraries:

- Primary equal-library mean: E_gpt = (1/L_t) sum_(l in t) M_glp.
- Explicit pooled-cell sensitivity: W_gpt = sum_(l in t) N_lp M_glp / sum_(l in t) N_lp.
- Difference, for either weighting: Delta_gp = mean_gp,PA - mean_gp,CTRL.

Both populations use identical library sets and weight 1/2 for the primary comparison. Pooled-cell weights are population-specific recovered-cell fractions and therefore deliberately assess a different weighting choice. Within a library, whole and hepatocyte means are paired summaries of nested populations, not independent samples. No pairing across CTRL and PA libraries is assumed. Optional accounting output reports within-library whole-minus-hepatocyte differences and the difference between the two population-specific condition contrasts; this is descriptive arithmetic, not an isolated causal composition effect.

## Missingness and retention

Assert unique library/category/gene keys, the complete 4 x 4 x 25 category grid, matching candidate/gene identifiers and coverage, exactly four named libraries and two condition labels, positive stratum cell counts, and per-gene invariant cell denominators. Retain zero means and all negative or positive differences. Covered candidates must have finite, nonnegative means and valid detection counts. AOC4P must be uncovered in all strata and retain NA expression and detection values; never impute unavailable genes as zero. No missing-library/category replacement or renormalization of available library weights is allowed: input incompleteness fails the run.

## Verification before interpretation

1. Reconcile cell counts to category_sampling_counts and independent library_candidates counts, and to 24,034 overall / 18,097 Hepatocytes.
2. Reconstruct every covered whole-library mean from the disjoint category summaries and compare with the pre-existing direct whole-library mean (absolute tolerance 1e-12).
3. Reconstruct every covered whole-organoid condition mean and difference for each weighting; compare with both existing Kitagawa totals and weighting-sensitivity results (absolute tolerance 1e-12), including six-decimal AREG +0.054898.
4. Independently recompute all library and condition arithmetic from CSV strings using standard-library Decimal (precision 50); compare with the numpy/pandas calculation (absolute tolerance 1e-12). This validates the aggregation arithmetic, not raw-cell source generation.
5. Verify equal-library mean of same-library whole-minus-hepatocyte differences equals the difference of population means within each condition; verify the corresponding PA-minus-CTRL contrast identity. Check the population-specific pooled weights sum to one and reconstruct both pooled means.
6. Assert all 25 candidates, missing AOC4P, measured zeros and any contrary results are retained, without inferential claims. Save results only when checks pass.

## Outputs and reporting

`run.py` must work using only its local input copies and Python/numpy/pandas; relative paths anchored to the script make it portable. Save library and condition tables for all25, separate compact AREG tables, paired-population accounting, validation records, dependency/interpreter and execution records, source/code/output hashes, a data dictionary and concise bilingual Methods/Results. No plot or manuscript is edited by this analysis task. No stochastic computation is used and no random seed is needed.

Interpretation is descriptive and limited to this differentiation. New concordance or discordance cannot establish a mechanism, separate all subtype/depth/technical influences, restore missing exposure metadata, provide independent replication, or establish a sesamin effect. Distinguish the new same-scale comparison from the older hepatocyte pooled-UMI CPM ratios, which use different transformations and weighting. Report unfavorable or same-direction results faithfully.
