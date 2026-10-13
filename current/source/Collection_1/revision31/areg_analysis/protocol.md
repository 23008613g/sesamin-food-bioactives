# AREG raw-source verification protocol, PFHN revision 31

Recorded 2026-10-11, before executing this revision's raw-expression recalculation.
This is a post-hoc audit. The analyst has read the old manuscript, prior AREG
results and the revision strategy. Neither the analysis nor its thresholds are
claimed to be prospectively preregistered or independent biological validation.

## Fixed scope

1. Read the original local deposited `GSE207889_os.h5ad` only for `obs` labels;
   do not use processed X for raw expression and do not reannotate cells.
   Read the four original MTX matrices, feature lists and barcode lists for
   GSM6322577--GSM6322580. Join by library accession and uniquely parsed barcode.
   Retain every author-annotated cell in these four libraries without a new QC
   threshold, downsampling or a forced target cell count. Report unmatched and
   duplicate keys. Preserve raw broad/subtype labels separately from labels
   defined for this analysis.
2. Verify whether broad `Cholangiocytes` is deposited in the author's obs column
   or introduced by the earlier local scripts. Report its complete source-label
   crosswalk. Keep the broad category and all constituent subtypes. Add the
   explicitly post-hoc sensitivity group Cholangiocytes + Ductal cells excluding
   Smooth muscle cells, without changing the deposited labels.
3. For AREG, x is its exact Ensembl raw UMI count and T is the sum across every
   measured feature in each original matrix (expected 33,538, to be verified).
   Compute ln(1 + 10,000*x/T) in float64 and include measured zeros in means.
   Absent features or absent groups produce NA, not an invented zero. Report
   each library and subtype's cell count, raw nonzero count/fraction, raw UMI
   total, and mean log1p(CP10k). Keep all four technical libraries.
4. Recompute whole-organoid and hepatocyte-like condition means with equal
   library weighting (two libraries per condition), plus pooled-cell weighting
   as a sensitivity, keeping transformation and zeros fixed. Also retain
   pooled-UMI CPM as a different estimand, not a population-only comparison.
   No cell-level P values, biological confidence intervals, or causal claims.
5. For GSE115469, reread the deposited author-normalized expression matrix and
   cell annotation. Preserve author labels; collapse only Hepatocyte_1..6 for
   the declared broad comparison. Keep all five donors and report paired
   Cholangiocytes-minus-Hepatocytes means and source-value nonzero fractions.
   Do not call those fractions raw-UMI detection. Retain P1's contrary result.
   Both groups >=20 cells is a post-hoc coverage sensitivity, not an exclusion
   from the primary analysis. Report every broad cell type including immune
   populations, without cell-type specificity claims.
6. Hash all source files and this protocol before execution. Preserve a run
   log, environment versions, script hash and output hashes. Compare recomputed
   values to archived summaries only after raw recomputation. Reconcile any
   differences rather than altering computation to reproduce old numbers.

## Verification

Check unique accession/barcode joins, conditions, raw matrix shape, nonnegative
integer counts and all-feature denominators. Confirm the label counts sum to
the actual retained population. Confirm fractions are bounded, all intended
library/subtype combinations exist (NA where empty), subtype raw sums reconcile
with totals, and direct cell weighting agrees with weighted library summaries.
Use scalar sums as an independent arithmetic check of key means. This is source
and arithmetic validation only; deposited annotations are not independently
validated cell identities and these references contain no sesamin intervention.
