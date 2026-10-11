# Source eligibility and normalization branch decision

Decision made before candidate-expression inspection on 2026-10-09. The original frozen protocol already prespecified retaining author-normalized data when raw counts are unavailable. This records the applicable branch and an explicit qualification to detection terminology; it does not alter the candidate set, comparison direction, donors, or sparse-coverage rule.

GEO family SOFT for all five GSE115469 samples states that the deposited CSV contains log2CPM expression, normalized using scran. The original article and its BaderLab/singleLiverCells/bareRawpipeline.R code describe filtering then scran normalization. The downloaded author code calls normalize(eb1F) (line 104) and passes do.logNormalize=F to Seurat (line 128). The precise scran version/pseudocount is not stated in these files. Therefore the analysis will use the deposited expression directly, label means as author-normalized log2CPM (GEO label), and report the percentage of entries strictly above zero as **nonzero deposited expression**, not verified raw-UMI detection. We will numerically check all entries for finite/nonnegative values and report zero and fractional-value counts. The expression means are interpretable on the deposited scale only, not as log fold changes or comparable magnitudes to organoid log1p(CP10k). No inverse transformation or re-normalization will be performed.

If negative values are present, the source-zero metric will not be presented as detection; means remain the primary expression outcome. Even if all values are nonnegative with many exact zeros, this numerical behavior alone does not prove raw-count equivalence. No SRA/FASTQ alignment or raw-count reconstruction is performed in this bounded analysis.

Public donor identity in the cell annotations is P1TLH through P5TLH, corresponding to GEO Patient 1 through Patient 5 (GSM3178782–GSM3178786). Both target categories occur in all five donors. Cholangiocyte counts are 8, 23, 41, 13, 34; hepatocyte counts are 96, 629, 2621, 43, 112. Retain all five donor pairs, and report the prespecified >=20-cells-per-type sensitivity using P2/P3/P5. These choices are made from metadata alone. All 20 author labels are retained in coverage and subtype tables; six Hepatocyte_# labels are combined only for the prespecified broad comparison.

GSE207889 family SOFT identifies W01 (WiCell) hPSC-derived HLOs, including control-for-PA and PA-500 uM day-25 samples. GSE115469 comprises independent primary adult liver transplant-donor specimens in a separate study. Independence is supported by documented study/sample provenance, not genotype matching. Healthy-atlas donors had normal histology but the article notes mild inflammation associated with neurological death. This dataset cannot confirm sesamin treatment response, protectiveness, clinical benefit, or identical adult-tissue and organoid localization.

Sources accessed 2026-10-09:

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3178784
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207889
- https://doi.org/10.1038/s41467-018-06318-7
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6197289/fullTextXML
- https://github.com/BaderLab/singleLiverCells/blob/master/bareRawpipeline.R (download resolved to commit in source_manifest.json)
