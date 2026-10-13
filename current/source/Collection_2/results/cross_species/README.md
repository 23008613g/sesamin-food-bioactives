# Descriptive cross-study expression comparison

All outputs use the full retained human and eligible rat gene universes; no 25-candidate selection was applied.
`shared_gene_mapping.csv` defines the exact shared universe; full inclusion/exclusion audits are retained.
`gene_comparison_summary.csv` reports all 3 rat normalizations x 2 human normalizations x 2 human contrasts.
`gene_direction_contingency.csv` includes all nine negative/zero/positive combinations for every comparison.
`shared_gene_effects_rat_*_human_*.csv.gz` retain every mapped gene without effect or significance thresholds.

Rat effects are differences of group means of log microarray signals; human effects are log ratios of mean normalized counts plus 1. They are not the same physical effect scale. Spearman coefficients are descriptive and no correlation p values or confidence intervals are computed.

`matched_go_membership.csv.gz` defines identical human/rat members for each included GO set. Original sets must retain 15-500 shared genes. Human log2(normalized counts+1) and rat log-signal gene matrices are standardized separately across their nine libraries and ten animals (ddof=1), then member gene z scores are averaged equally. Constant genes contribute zero, if present. Complete individual sample scores, contrasts, set-membership hashes, exclusion counts and normalization sensitivities are retained.

Human libraries have unresolved independence of culture units and exposure timing. No human inferential tests were performed. Gene members and overlapping GO sets are not independent biological replicates. No result establishes normal human hepatocyte validation, treatment rescue, pathway activation, metabolic flux or causality.
