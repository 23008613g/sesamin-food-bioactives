# Bounded quantitative reproduction package

This local release candidate rebuilds the manuscript's human/mouse candidate robustness, rat intervention analysis, library-specific organoid summaries, and external adult-liver atlas results from packaged processed numerical inputs. It contains real processed data, not simulated examples. It does not reproduce all work from raw sequencing or array files. The package itself does not create a public repository, licence, DOI, or publication.

The entry point uses paths relative to this directory. The original project, its Windows drive letters, network access, and data downloads are not needed at runtime. `provenance/source_manifest.json` records dataset accessions, public source URLs, relative artifact names, byte counts, and SHA-256. Workstation paths, raw local logs, author full text, and third-party helper source files are excluded. Execute `run.py`, which uses the scripts in `code/`.

## Run

Use Python 3.13.9 with the versions in `requirements.txt`; see `ENVIRONMENT.md`. The following command writes a new `results/` folder relative to the package, even when run from a different working directory:

```text
python /path/to/portable_package/run.py
```

For another run, choose a new output directory:

```text
python /path/to/portable_package/run.py --output results_second_run
```

Existing output directories are refused. No files are downloaded, installed, uploaded, or deleted by this command. It executes, in order:

1. `code/robustness.py`: human full-gene median-ratio normalization and mouse processed matrices; all 64 configurations per species, all 4,096 cross-species pairs, all fixed 15 primary orthologues, five measurement/normalization specifications, means and medians, and full-transcriptome QC.
2. `code/robustness_check.py`: a separate original calculation path checks 5,760 mouse and 270 human contrasts plus 30 complete cross-species candidate summaries, and recomputes abundance flags.
3. `code/rat.py`: map the fixed 25 human candidates through the frozen rat annotation/orthology inputs, calculate `log2(signal+1)`, median probe aggregation, Welch effects/intervals/tests, and BH over a family of 25. Quantile normalization includes all 20 deposited arrays; only the 5 control and 5 sesamin animals enter the primary comparison.
4. `code/rat_exact.py`: enumerate all 252 allocations of the 10 selected animals, then apply BH over 25 candidates. This is the original post hoc exchangeability sensitivity, not a verified randomized experiment.
5. `code/organoid.py`: reconstruct candidate `log1p(CP10k)` from cell counts and stored whole-cell totals; summarize 24,034 cells in four libraries, all 25 candidates (24 measured, AOC4P missing), four broad categories and 10 subtypes; rebuild library contrasts and symmetric Kitagawa decompositions under all three weighting schemes.
6. `code/external/analyze_external.py`: read the full deposited GSE115469 author-normalized matrix, preserve five donors and all fixed candidates, and recompute donor-level descriptive cellular-context results. See `EXTERNAL_MODULE.md`.
7. `verify.py`: compare all four modules to independently retained reference tables, cell by cell. The three earlier modules comprise 48 tables, 198,755 rows, and 624,747 finite numeric values; external tables are added to the same verification, with totals recorded in `verification/observed_verification.json`. Labels and missingness are also checked. Integer/Boolean values are exact; other numbers use `rtol=1e-9`, `atol=1e-10`. A global sign reversal of rat PCA axes may be aligned because PCA axis signs are arbitrary.

Run logs are in `results/logs/`. `results/execution.json` records commands, runtime versions, timestamps, and input integrity. `results/verification.json` records all comparisons and main reproduced conclusions. To compare an already computed result folder again, use `python verify.py --results /path/to/results`.

## Expected quantitative findings

| Result | Expected result |
|---|---|
| Human/mouse primary candidates | 15, fixed before this rerun; no candidate reselection |
| Configurations | 64 per species; 4,096 Cartesian product pairs, not paired biological samples |
| Both mean-based contrast directions concordant | VNN1: 4,096/4,096; AREG: 3,072/4,096; JAG2: 2,560/4,096 |
| Mouse VNN1 treatment-versus-model mean TPM difference range | -48.438929% to -27.183572% |
| Mouse AREG range | -23.166687% to +9.761990% |
| Mouse JAG2 range | -76.626442% to +78.109024% |
| Rat mapped candidates | 15/25; retain the other 10 as missing, not as null effects |
| Rat Welch BH <0.05 | VNN1 only, under each normalization |
| Rat VNN1 exact sensitivity | 2/252, p = 0.00793650793651; BH q = 0.198412698413 under each normalization |
| Rat exact BH <0.05 | None |
| Organoid AREG cholangiocyte-minus-hepatocyte detection differences | +48.874102, +40.547119, +44.091607, +36.384743 percentage points in GSM6322577–GSM6322580 order |
| Organoid weighting sign agreement | 21/24 measured candidates; disagreements: PSTPIP2, BCAS1, C2orf49-DT |

Full precision and every candidate, including negative and unstable findings, remain in the output tables. Enumeration fractions are not replication probabilities or confidence intervals. Organoid observations are four technical libraries from one differentiation, without sesamin treatment; neither cell numbers nor library counts establish biological replication. The decomposition is an algebraic description, not causal adjustment.

## Scope and input boundary

The 18 core numerical/mapping/reference inputs total **29,344,241 bytes**; five external-atlas inputs add **108,159,417 bytes** before adding protocols, expected tables, scripts, and documentation. All required inputs for the stated rerun are present. `DATA_DICTIONARY.md` defines their roles; `provenance/source_manifest.json` gives source accessions/URLs, relative original artifact names, byte counts, and SHA-256. Original protocols and freeze records are preserved unchanged; this rerun is not a new prospective protocol registration. `provenance/script_adaptations.json` records original and portable script hashes and the scoped changes: runtime paths/imports/invocation metadata and omission of one unused private upstream metadata log. Original source files and detailed patches remain in the local audit. Numerical algorithms and statistical parameters are unchanged.

The following earlier stages are deliberately outside this package: raw FASTQ quality control/alignment/quantification; mouse transcript-to-gene extraction and orthology selection; rat CEL processing and sequence-level probe remapping; original human candidate selection; raw organoid MTX/h5ad barcode reconciliation, all-gene UMI extraction, original cell labels/UMAP generation, and cell exclusion decisions. Stored cell-wide totals and labels are input data, not regenerated observations. The `__log1p_cp10k` columns are reference values; candidate normalization is genuinely recalculated from counts and totals. All-gene matrices needed to regenerate cell totals are not included. Figure drawing remains in the parent manuscript delivery. The external adult-liver atlas is a separately documented fourth module, explicitly included in the rerun and numerical comparison.

The supplied matrices and derived data retain source-dataset provenance. No blanket licence has been assigned to third-party datasets; see `LICENSE_PROVENANCE.md`. Construction and verification were performed locally; hosting status is controlled by the parent delivery, not by this program. The independent-directory check uses the same existing scientific Python environment: portability across locations was executed; a freshly installed environment, Linux/macOS run, and download-from-public-repository workflow were not tested.
