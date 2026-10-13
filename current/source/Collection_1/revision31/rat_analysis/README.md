# GSE42693 rat comparison: completed revision audit

Status: **EXECUTED AND CHECKED**. The stage31 revision reran the deposited-matrix analysis; it did not merely copy the earlier result tables. Numerical conclusions are unchanged. The deposited matrix is unlogged, globally scaled signal, so the existing log2(signal + 1) transformation is retained. No double-log error was found.

## Scope and sources

- Source expression: GSE42693, 31,099 probes × 20 arrays. Five control and five sesamin animals form the primary contrast. The remaining five episesamin and five sesamolin animals are distinct groups used only in the all-array quantile-normalization sensitivity and source diagnostics.
- Original study: Ide et al., Journal of Nutritional Science and Vitaminology 55:31–43 (2009), DOI [10.3177/jnsv.55.31](https://doi.org/10.3177/jnsv.55.31), printed page 32. The paper links to this accession via GEO's PMID 19352061. It reports seven male rats per diet, with the highest and lowest terminal-body-weight animals excluded from microarray analysis, leaving five individual livers per group. We exclude no additional animal.
- [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42693), [matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42693/matrix/GSE42693_series_matrix.txt.gz), [GPL1355 annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1355/annot/GPL1355.annot.gz).
- GEO's embedded processing states GCOS 1.4 with global scaling. The matrix's minimum/median/maximum are 0.0221933 / 27.04825 / 18,018.4. The source paper describes a subsequent GeneSpring workflow, including percentile/per-gene normalization and filtering. We do not conflate that published workflow with the deposited signal used here, and did not rerun CEL processing or detection-call filtering.
- Probe mapping uses the unchanged GPL1355 annotation dated 2016-08-09 and archived Ensembl 116 orthology/Entrez cross-references. It is identifier-based, not a newly performed sequence-remapping analysis.

## Results and interpretation

All 25 fixed human candidates remain in the results: 15 measurable, nine without a rat orthologue and one without one-to-one orthology. Each of the 15 measurable candidates has exactly one eligible probe under the fixed Entrez rule. Thus the specified within-animal median-across-probes rule reduces to the single probe here; there was no selection among multiple probes based on effect. VNN1 maps to rat Vnn1, Entrez 29142, probe 1389253_at, with a single annotated Entrez identifier. All 31,099 probe annotation records are supplied for audit.

| VNN1 statistic | Deposited log2(signal + 1) | All-20-array quantile sensitivity |
|---|---:|---:|
| Sesamin minus control mean | 1.67014258237446 | 1.77477683493553 |
| Nominal 95% Welch interval | 1.13384867687983–2.20643648786909 | 1.23933814076375–2.31021552910732 |
| Two-sided Welch P | 0.000102764047213 | 0.0000942400040482 |
| Welch BH25 adjusted P | 0.00256910118032 | 0.00235600010120 |
| Exact-label P | 2/252 = 0.00793650793651 | 2/252 = 0.00793650793651 |
| Exact-label BH25 adjusted P | 0.198412698413 | 0.198412698413 |
| New secondary Welch BH15 | 0.00154146070819 | 0.00141360006072 |
| New secondary exact-label BH15 | 0.119047619048 | 0.119047619048 |

The observed difference is positive and all five sesamin values exceed all five controls. Welch and exact-label inferential results differ. The exact reference is conditional on exchangeability; random allocation was not verified, and neither analysis removes the original terminal-weight selection. The confidence interval is nominal, not a multiplicity-adjusted simultaneous interval. The finite 252-allocation distribution makes the two-sided P discrete; this does not repair design limitations.

The fixed multiplicity family remains 25. Missing hypotheses contribute P=1 only inside adjustment and remain NA in the report. The additional BH15 adjustment is explicitly a **new post-hoc secondary sensitivity** over the measurable family. It does not replace BH25; no candidate passes exact-label adjusted P <0.05 with either family or either normalization.

The newly calculated values agree with 400 saved statistics (including NA states) and 600 saved animal values to their stored precision. The largest absolute statistic difference is 4.64×10^-12; the largest animal-value difference is 4.98×10^-11. These are numerical/storage-precision differences, not a change in scientific results. There are 128 passed implementation/data-integrity checks, including manual versus SciPy Welch calculation, independently implemented quantile normalization, BH agreement with SciPy, complete exact allocation enumeration, and unchanged source hashes. This is computational self-checking, not an independent expert review.

## Ready-to-use manuscript text

### Methods

For GSE42693, we analyzed the deposited GCOS 1.4 globally scaled Affymetrix signal rather than raw CEL files or the original paper's downstream GeneSpring-normalized values. The signal was unlogged; we applied log2(signal + 1). The primary contrast comprised five control and five sesamin-fed rats (0.2% dietary sesamin, 15 days). The original investigators had selected these animals from seven per group by excluding the highest and lowest terminal-body-weight animals. We retained all deposited animals. Fixed Ensembl 116 one-to-one orthologues were linked to GPL1355 probes with a single annotated Entrez identifier. Fifteen of the 25 candidates were measurable, each with one eligible probe; the general within-animal median-probe aggregation rule therefore introduced no multi-probe choice. Two-sided Welch tests and nominal 95% confidence intervals were supplemented by the previously specified post-hoc exact-label reference, enumerating all 252 five-versus-five allocations using the absolute mean difference. Fixed-family Benjamini–Hochberg adjustment included 25 hypotheses; unavailable hypotheses contributed P = 1 only during adjustment and remained NA in reported results. Sensitivities used quantile normalization across all 20 arrays and, separately, a newly added measurable-family adjustment over 15 candidates. Exact-label inference assumes exchangeability and does not establish randomized causal validity or correct the original animal selection.

### Results

Reanalysis confirmed the positive VNN1 difference in the dietary rat dataset (+1.670 log2 signal units; nominal 95% Welch interval 1.134–2.206; Welch BH25-adjusted P = 0.00257). The quantile-normalization sensitivity retained the direction (+1.775). The exact-label reference yielded P = 2/252 and BH25-adjusted P = 0.1984; no candidate met exact-label adjusted P <0.05. Restricting correction to the 15 measurable candidates as a secondary sensitivity gave VNN1 exact-label adjusted P = 0.1190, without changing that conclusion. The observed direction was positive, whereas the strength of inferential support depended on the test assumptions. These findings describe selected healthy-rat whole-liver samples and do not establish a fatty-liver treatment mechanism.

## Files

- `PROTOCOL.md`: bounded plan recorded before this rerun, explicitly acknowledging earlier results were known.
- `reanalyze_rat.py`: computation from unchanged deposited matrix and annotation inputs.
- `input_manifest.json`: source roles, project-relative paths, hashes, sizes and source URLs.
- `scale_audit.json`: scale decision, deposited metadata, distribution and source-paper processing distinction.
- `analysis_summary.json`: actual environment, dimensions, numerical results, comparisons, scope and limitations.
- `source_evidence_log.md`: source verification and access outcomes.
- `results/sample_metadata_and_scale.csv`: all 20 arrays, source labels, exposure groups, experimental units and raw-scale diagnostics.
- `results/all_platform_probe_annotation_audit.csv.gz`: every GPL1355 probe's identifiers and mapping uniqueness.
- `results/candidate_mapping_all25.csv` and `candidate_probe_mapping_decisions.csv`: complete candidate and eligible-probe mapping.
- `results/excluded_ambiguous_probes.csv`: exclusion ledger (zero excluded matching probes in this fixed set; header retained).
- `results/candidate_results_all25_welch_exact_two_normalizations.csv`: all 50 candidate/normalization records, both tests and both correction families.
- `results/candidate_sample_values_all25.csv`: complete 25 × 20 × 2 grid including explicit unavailable values.
- `results/candidate_probe_values_long.csv`: each contributing probe/animal value, before and after transformation.
- `results/old_new_statistic_differences.csv` and `old_new_animal_value_differences.csv`: full old/new comparison.
- `results/verification_checks.csv`: 128 checks with outcomes.
- `output_manifest.json`: result-file hashes.

## Rerun

From the company-article workspace, use the existing project environment; do not install or upgrade packages:

```powershell
& 'research_planning/hepatic_lipotoxicity_2026-10-07/.venv/Scripts/python.exe' -X utf8 'research_planning/hepatic_lipotoxicity_2026-10-07/next_steps/pfhn_strategy_revision_31/rat_analysis/reanalyze_rat.py'
```

The script uses the documented project-relative frozen inputs and writes only this module's own results. A recipient needs the listed inputs; the manuscript or a script alone is not a complete executable archive. Python 3.13.9, NumPy 2.4.6, pandas 3.0.3 and SciPy 1.17.1 were actually loaded for this execution. The original paper was read with the pre-existing document Python/PyMuPDF environment solely for Methods verification; it is not a statistical dependency.

Not performed: raw CEL reprocessing, sequence-level probe remapping, new assays, post-hoc animal deletion, a new randomization procedure, external upload, contact with source authors, or a claim that the healthy-rat dataset validates fatty-liver efficacy.
