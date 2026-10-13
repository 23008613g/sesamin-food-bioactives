Post-hoc signal and stability diagnostic of existing limma results
================================================================

Scope
-----
This directory contains descriptive diagnostics of already-executed Stage 39
results. No model, normalization, hypothesis test or BH family was rerun or
changed, and no gene was filtered. Five control and five sesamin animals are
the same under all three normalization variants. Stability across variants
is not independent biological replication.

Important definitions
---------------------
AveExpr is the arithmetic mean of the ten gene-level log2-signal values,
not log2 of their mean untransformed intensity. It was numerically checked
against both existing group means and the frozen all-gene result table.

Signal percentile = 100 * count(all 14,604 means <= this mean) / 14,604,
calculated separately within each normalization. Ties share the upper
empirical rank. The Q1-Q4 intervals are (0,25], (25,50], (50,75], (75,100].
These newly described quartile strata are exploratory summaries, not
detection thresholds, inferential tests or filtering rules.

BH hit membership uses the original adj.P.Val < 0.05, unchanged. The
complete 24/41/41 hit lists, all-gene context and all membership patterns,
including empty patterns, are retained. Sign and size comparisons use
existing sesamin-minus-control coefficients. Rank correlations are
descriptive only; no correlation P values or confidence intervals are used.

Files to integrate
------------------
Table_S20_limma_signal_summary.csv : compact three-row supplementary table
Table_S20_legend.txt               : complete table definition and caveats
FigS19.pdf / .svg / .png           : supplementary two-panel diagnostic
FigS19_legend.txt                  : complete figure caption
Methods and Results, including the signed/absolute effect-versus-signal
results, are integrated in ESM_1.pdf rather than distributed as drafting files.

Detailed reproducibility tables
-------------------------------
limma_hits_{deposited,qn10,qn20}.csv             : all hits, ranked by signal
limma_all_genes_signal_context_{variant}.csv.gz : all 14,604 genes per variant
limma_union_hits_stability.csv                  : all 41 union genes
limma_all_gene_normalization_stability.csv.gz  : full-universe comparison
signal_distribution_summary.csv                : min/Q1/median/Q3/max
signal_quartile_counts.csv                      : denominators and counts
effect_signal_correlations.csv                 : all/hit/non-hit rank effects
hit_membership_patterns.csv                    : all eight Boolean patterns
normalization_pairwise_stability.csv           : overlap, retention and signs
lower_quartile_hits_context_not_detection_calls.csv : descriptive subset only

Executed environment
--------------------
Existing project Python 3.13.9; NumPy 2.4.6, pandas 3.0.3, SciPy 1.17.1,
Matplotlib 3.10.9. Figure PDF audit uses the separately existing PyMuPDF
environment. Nothing is installed or upgraded; PATH is not modified.
The exact executable paths and hashes are recorded in the JSON reports.

Reproduction
------------
Use a prepared environment with the recorded dependencies. Extract ESM_3.zip
to a source directory, then copy revision43/diagnostics/ to a NEW working
directory OUTSIDE that source directory. This keeps the extracted source and
its recorded outputs unchanged. Run the copied scripts from that working
directory, using the absolute root containing results/ as --stage39:

  python diagnose_limma_hits.py --stage39 /path/to/extracted_ESM_3_root

The output defaults to the copied script's directory. The plotting script
reads diagnostic tables beside itself. Rebuild the figure using the supplied
auditors and an existing interpreter containing PyMuPDF:

  python make_limma_diagnostic_figure.py --qa-tools /path/to/extracted_ESM_3_root/figures/qa_tools --audit-python /path/to/pdf-audit/python

Default local paths merely record this execution environment; pass the above
arguments on another machine. The diagnostic script asserts the known frozen
24/41/41 lists and observed membership/quartile counts so changed inputs cannot
silently generate the same narrative. A changed source version requires review.

QA records
----------
limma_signal_diagnostics.json : executed values, boundaries, input/output hashes
figure_qa/figure_execution.json : actual render commands and format checks
figure_qa/visual_review.json : actual final PNG inspection, not author approval
final_qa.json : table roundtrip, source preservation and final integration checks

The 600-dpi PNG is a direct render. PDF/SVG remain vector graphics with embedded
PDF fonts. The plot includes every hit and all nonempty intersection patterns;
the diagonal is a uniform percentile reference, not a fitted null model.
