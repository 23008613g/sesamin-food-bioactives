# Scientific module and figure map

The repository includes the current processed-input execution paths and earlier source code. Historical source files remain separated because the required directory layout and acquisition history differ from the current portable package.

| Scientific stage | Source |
|---|---|
| Input metadata, sample manifests, processed-table inspection | `historical_pipeline/audit_metadata.py`, `build_sample_manifest.py`, `inspect_processed_tables.py` |
| Initial human candidate selection and independent numerical review | `historical_pipeline/next_steps/sesamin_exploration_04/run_analysis.py`, `numeric_review/recompute.py` |
| Primary human hepatocyte analysis/checks | `historical_pipeline/next_steps/followup_05/analyze_phh.py`, `check_robustness.py`, `numeric_review/recompute.py` |
| Mouse metadata and frozen orthology | `historical_pipeline/next_steps/mouse_readiness_07/` |
| Mouse FASTQ download/recovery, read QC, reference/index, kallisto orchestration and gene aggregation | `historical_pipeline/next_steps/mouse_quantification_08/` |
| Original organoid preparation, barcode/feature linkage and count extraction | `historical_pipeline/next_steps/organoid_prepare.py`, `organoid_expression/` |
| Original single-cell derived features/checks | `historical_pipeline/next_steps/single_cell_extension_14/` |
| Full candidate robustness, rat and organoid summary reanalysis | `reproduction/code/`; original source/check variants under `historical_pipeline/next_steps/analysis_strengthening_15/` |
| Independent adult-liver context | `reproduction/code/external/` |
| Matched-scale whole/Hepatocyte comparison | `analysis/matched_scale/run.py` |
| Current main plots | `figures/build_figures.py --only all` |
| Display data and frontend | `app/`; historical builders/checks under `historical_pipeline/app_tools/` |

Historical raw-processing/setup scripts may download data or require external software if explicitly invoked. They are not called by any current reproduction entry point. Raw reads, matrices, software executables and historical workspace contents are not bundled. Redacted `<WORKSPACE>`/`<LOCAL_PATH>` strings identify an archived dependency rather than a usable default path.

## Final main figure routing

| Main | Question/source function |
|---|---|
| 1 | Six dataset roles; `build_stage26.figure1` |
| 2 | Full-gene cached PCA; `build_figures.pca` |
| 3 | All25 human candidates and15 primary categorical comparisons; `build_figures.candidates` |
| 4 | Individual VNN1 observations in human/mouse/rat; `build_stage26.figure2` |
| 5 | Corrected robustness, including five constant-zero cells; `build_stage26.figureS16` |
| 6 | One saved author UMAP shown by subtype/AREG; `build_stage26.figure3` |
| 7 | AREG library detection and matched-scale means; `build_stage26.figure4` |
| 8 | Five adult-liver donors and immune-cell context; `build_stage26.figure5` |

## Fourteen supplementary figures

Historical figure generators are retained under `historical_pipeline/publication_package_12/figures`, `publication_package_13/figures`, `publication_package_14/figures/single_cell`, `analysis_strengthening_15/figures` and `submission_finalization_25/figure_revision`.

| Final supplement | Source content |
|---|---|
| S1-S2 | Library quality/PCA and normalization/omission sensitivity; stage12 support, corrected display labels in stage25 |
| S3-S4 | Full PHH candidate and mouse abundance matrices; stage13 new-support |
| S5-S7 | Per-library saved UMAP, raw-count QC and cell composition; stage14 single-cell figure source |
| S8 | PHH/background and fixed UMI weighting; stage12 main `fig4_background` |
| S9 | Subtype contrasts/UMI shares; stage13 new-main `fig6_organoid_subtypes` |
| S10 | Four candidate expression UMAPs; stage14 single-cell `fig8_sc_features` |
| S11 | Original software-interface screenshot; app source and historical capture script, not a new quantitative plot |
| S12 | Three focal candidates in all human/mouse libraries; stage12 main `fig2_expression` |
| S13 | Complete rat candidate results and sensitivities; stage15 main generator |
| S14 | All25 candidates by10 author subtypes; stage14 single-cell `fig9_sc_dotplot` |

The final manuscript retains the rendered archival supplementary figures. This source repository includes their generators but does not claim that the earlier multi-stage layouts can all be rebuilt with only the current main-figure input subset. The full old-workspace dependency is explicit; the portable v24 numerical package and current main figures are separate executable paths. Duplicate copied generators, Word/manuscript packaging, private reviews, literature full texts, credentials, deployment tooling and external helper implementations are excluded, with per-file reasons in the inventory.
