OFFLINE REPRODUCTION OF THE STAGE 39 ANALYSIS
Client date: 2026-10-13 (Asia/Taipei)

Scope
-----
scripts/reproduce_analysis.py executes the original rat genome-wide analysis,
the existing limma sensitivity analysis, and the full cross-species comparison
without changing their statistical code. It reproduces all result tables made
by those three main functions and compares them numerically with this bundle.
The wrapper does not download data, update annotations, install packages,
rewrite manuscript files, build publication figures, or repeat the separate
independent-audit scripts. The limma main function still generates its own
plotSA diagnostic PDFs as part of its original execution.

Required existing software
--------------------------
The completed validation used Python 3.13.9, NumPy 2.4.6, pandas 3.0.3,
SciPy 1.17.1, R 4.6.1, and limma 3.68.4. The R script explicitly requires
the recorded R/limma versions and stops if they differ. Python versions are
recorded in the execution report; numeric comparisons detect result drift.
No package manager, PATH alteration, authentication, or external service is
used by this wrapper. Use an existing appropriate Python environment.

Example command (PowerShell; run from the extracted bundle directory)
-------------------------------------------------------------------
& "C:/path/to/existing/python.exe" "./scripts/reproduce_analysis.py" --output-dir "D:/ResearchData/pfhn_stage39_fresh_run" --rscript "C:/Program Files/R/R-4.6.1/bin/Rscript.exe"

Pass --rscript explicitly on another computer. Its default is the known
Windows installation above; the wrapper does not search for or install R.
The output directory MUST NOT EXIST, even as an empty directory, and must be
outside the source bundle. If execution fails, its logs are retained; choose
another new output directory for a later run. Existing results are never
silently overwritten or deleted.

Optional: --reference-dir "path/to/another/complete/reference_bundle" changes
the reference result location. It defaults to the directory containing this
wrapper's scripts/ folder. Keep the original results/ folder with the bundle
when using the default numeric-comparison step.

Bundled inputs and frozen annotations
------------------------------------
sources/project_inputs/ contains verified read-only copies of five inputs,
preserving their original next_steps/ relative hierarchy:
  - GSE42693_series_matrix.txt.gz (deposited rat expression signals)
  - GPL1355.annot.gz (full historical platform annotation)
  - GSE248903_total_count.anno.txt.gz (original human count matrix)
  - all_gene_descriptive.csv.gz (archived full human descriptive results)
  - sample_qc.csv (archived human normalization factors and QC)
Their relative paths, original locations and SHA-256 hashes are retained in
sources/project_inputs/source_manifest.json. The original project is no
longer required for a reproduction run.

sources/go/ retains the fixed ontology and download provenance. The ontology
SHA-256 is checked against its frozen source manifest.
sources/orthology/ retains the full Ensembl 116 source tables, audit tables,
strict one-to-one analysis map and freeze report. Every output identified in
that freeze report is hash-checked. No BioMart or GO request is made.
Both protocol files and both freeze records are copied and verified.

How execution is isolated
-------------------------
The wrapper copies the complete frozen inputs and three original analysis
scripts into the new output directory and confirms copy hashes. It imports
the unchanged Python analysis modules from that new directory, redirects only
their ROOT/OUT/PROJECT/NEXT/INPUTS/REVIEW locations, calls their original main
functions, and invokes Rscript as an explicit subprocess. It does not execute
the scripts/fetch_* files. Source input hashes are checked again after the run.

Execution records
-----------------
The new output directory contains:
  review/portable_reproduction.json -- commands, runtimes, input/code hashes,
      phase status and per-table comparisons
  analysis_summary.json -- reproduced rat/GO analysis summary
  review/limma_execution.json -- R model and synthetic software checks
  review/cross_species_execution.json -- cross-species execution and checks
  logs/ -- separate stdout/stderr for each of the three phases
  results/ -- regenerated rat, limma and cross-species result tables
  results/limma/sessionInfo.txt -- actual R and loaded package versions

The table comparison checks shapes, column order, row order, identifiers and
non-numeric content. Numerical measures must match with atol=1e-10 and
rtol=1e-10. It compares decompressed CSV values, not gzip timestamps. PDF/RDS
binary bytes and timestamp/path-containing JSON records need not be identical.
Tables added later for visualization or independent audits are not themselves
products of the three original main functions and are outside this comparison.

Actual validation completed
---------------------------
On 2026-10-13 the wrapper was executed end-to-end in the new directory
D:/ResearchData/pfhn_genome39_reproductioncheck_2026-10-13 using the existing
project Python and R installations. After a comparator-only missing-value
parsing improvement, the final wrapper version was executed in another new
directory, D:/ResearchData/pfhn_genome39_reproductioncheck_final_2026-10-13.
Rat, limma and cross-species phases all completed in both runs. All 62
generated CSV/CSV.GZ tables matched the reference in each run, including
identifiers and row order; maximum numerical difference was 0.0 at the written
CSV precision. The existing-output refusal was also executed and passed,
without changing that successful execution record.

The corresponding summary is review/portable_reproduction.json in this bundle.
This is an actual fresh-directory, same-computer validation. A different
operating system, CPU/BLAS implementation or dependency set has not been
executed here; strict comparisons can flag platform-related numeric or PCA
sign differences. Reproduction confirms computation, not biological validity,
independent experimental replication, normal-human response or clinical effect.
