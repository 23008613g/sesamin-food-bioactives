# Stage 15: library-aware single-cell strengthening protocol

Frozen scope before executing this extension; 2026-10-09 (client date).

## Question and inputs

Evaluate whether the observed localization of candidate expression is consistent across the four source technical libraries, and quantify how control/PA differences in mean log1p(CP10k) separate algebraically into cell-composition and within-category components. Only the already verified cell-level and summary tables supplied in the stage 14 package are read. No raw matrix, h5ad, embedding, cell label, candidate selection or original result is changed.

All 25 fixed human candidates are retained: 24 measured candidates plus AOC4P, absent from the feature universe and preserved as NA. All 24,034 author-retained cells from the two control and two PA libraries remain included. Ten author subtype labels and four author broad categories are retained. The four technical libraries derive from one differentiation and do not provide independent biological replication. No sesamin is present in this reference subset.

## Per-library localization

For each library and every author subtype or broad category, recompute candidate detection (raw count > 0) and mean log1p(CP10k), including measured zero-count cells. The denominator for CP10k is the existing all-gene raw UMI total per cell. Report the full candidate grid, sample sizes and NA entries. For every candidate, compare the author's Cholangiocytes and Hepatocytes broad categories within each of the four libraries. Report all four signed detection-percentage-point and mean-expression differences, together with pooled-cell and equally weighted library summaries. AREG is a prespecified highlighted candidate; it does not determine inclusion or the direction of any comparison. Direction agreement is descriptive, not a significance test.

## Symmetric composition decomposition

Let t=0 denote control and t=1 denote PA, with two technical libraries l per condition and four mutually exclusive broad categories k. Let N_tl be the retained-cell count in a library, N_tlk the category count, p_tlk=N_tlk/N_tl, and μ_tlk the arithmetic mean of individual-cell y=ln(1+10000×raw_gene_count/raw_total_UMI) in that library/category. Let m_tl=Σ_k p_tlk μ_tlk. All categories must have positive counts in all four libraries; if a category is empty, its mean remains unavailable and the affected decomposition is not imputed.

For a specified condition-level profile (p_tk,μ_tk), define M_t=Σ_k p_tk μ_tk and Δ=M_1−M_0. The symmetric Kitagawa terms are:

- Composition contribution C_k=(p_1k−p_0k)(μ_1k+μ_0k)/2.
- Within-category contribution W_k=(μ_1k−μ_0k)(p_1k+p_0k)/2.
- Exact identity: Δ=Σ_k(C_k+W_k).

These terms are in mean log1p(CP10k) units, not log fold change, percentage protection or causal effects. Negative and positive contributions are retained. Contributions are not divided by small total differences to create potentially unstable attribution percentages.

### Primary: equal-library actual mixture

Give each of the two libraries weight 1/2. Define p_tk=(p_t1k+p_t2k)/2 and μ_tk=(p_t1k μ_t1k+p_t2k μ_t2k)/(p_t1k+p_t2k). Then M_t=(m_t1+m_t2)/2 exactly. The within-category mean is weighted by each category's contribution to the equal-library mixture; it is not the unweighted average of the two category means. This construction reduces dominance by the larger library while preserving the actual equal-library overall mean.

### Sensitivity: pooled-cell actual mixture

Define p_tk=Σ_lN_tlk/Σ_lN_tl and μ_tk=Σ_lN_tlk μ_tlk/Σ_lN_tlk. Then M_t equals the direct mean over all retained cells of that condition. Larger libraries contribute more cells by construction.

### Diagnostic only: separate unweighted profile averages

Define pbar_tk=mean_l(p_tlk), mubar_tk=mean_l(μ_tlk), and R_t=Σ_k pbar_tk mubar_tk. This is a reference mixture and generally differs from the actual equal-library M_t. Preserve the exact population-covariance identity M_t−R_t=Σ_k mean_l[(p_tlk−pbar_tk)(μ_tlk−mubar_tk)]. Thus the actual equal-library contrast equals the diagnostic reference-mixture decomposition plus the PA-minus-control covariance residual. This diagnostic is reported separately and is not a primary observed total.

## Depth and technical-library context

Report each library's retained-cell count, mean and median raw UMIs, detected-gene counts, mitochondrial percentages and broad-category cell fractions alongside its candidate results. Do not add four-point correlation statistics, cell-level or technical-library P values, confidence intervals, automatic batch correction or new clustering. Technical-library and cell-composition differences may contribute to expression/detection differences, but this design cannot separate them causally.

## Verification and provenance

Check unique library/barcode keys, expected cells and all candidate/group combinations; all four category fractions must sum to one. Reconstruct normalization from raw candidate counts and raw total UMIs, and crosscheck the regenerated 25×4×10 table against the stage 14 summaries. Verify decomposition identities both algebraically by independent scalar evaluation and numerically against direct condition means, including the diagnostic covariance residual. Verify that swapping control and PA reverses the total and both components. Store code, exact formulas, inputs and SHA256 checksums, environment versions, commands, execution log, output tables and numerical validation. Preserve failed attempts truthfully. No new figure is drawn until the numerical outputs are returned to the root agent.
