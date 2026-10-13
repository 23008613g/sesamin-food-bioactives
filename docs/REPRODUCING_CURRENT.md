# Reproducing the current numerical analyses

## Acquire and inspect

The release attachment `Supplementary_Data.zip` is the same 102,653,038-byte archive as the current manuscript supplement. Its SHA-256 is `f369474b2856389de42671112c9fdc2eea1e85c1dbb0b05feb79ff34368e57bd`. `tools/acquire_current.py` checks these values before safely extracting to a new directory. No extraction runs analysis code. The top-level archive README supersedes distribution descriptions in historical nested READMEs.

After extraction, run `python .local/current_data/restore_exact_duplicates.py --restore`. This restores 182 original duplicate paths from byte-identical canonical copies, refuses incompatible existing content and makes no network calls. The archive has 626 scientific/provenance files physically stored plus its packaging files. It omits 57 visual exports and three regenerable R fit caches; original numerical results and analysis code are retained. Omitted artwork is not necessarily regenerable from plotting code. Current figure numbering and illustrations belong to the manuscript, not all historical scripts.

## Whole-transcriptome numerical comparison

The archived validation used Python 3.13.9, NumPy 2.4.6, pandas 3.0.3, SciPy 1.17.1, R 4.6.1 and limma 3.68.4. The R script enforces the recorded R/limma versions. Use a suitable existing environment; this repository does not install or upgrade software automatically.

```sh
python .local/current_data/Collection_2/scripts/reproduce_analysis.py --output-dir /path/to/new_output_directory --rscript /path/to/Rscript
```

The output directory must not exist and must be outside the input bundle. This wrapper uses the frozen deposited rat matrix, historical probe annotation, human input matrix/descriptive results and normalization factors, fixed GO ontology and Ensembl 116 orthology. It does not download newer annotations. It runs the rat, limma and cross-species numerical modules and compares 62 CSV/CSV.GZ tables. Its documented tolerances are atol=1e-10 and rtol=1e-10. Publication figure rendering and the later signal-intensity diagnostic are separate workflows.

Read `Collection_2/revision43/diagnostics/README_reproduction.txt` for the latter diagnostic. Read `Collection_1/revision31/PORTABILITY.md` for candidate, rat-candidate and AREG modules: their input requirements and portability differ, and some raw-data stages are outside the supplied-input scope. Revision31 scripts write beside their own module directories. Before invoking them, create a separate fresh working copy of the complete extracted archive, and preserve the original extraction and reference results. Do not run scientific scripts in the browseable Git mirror. This is not a single command that regenerates every manuscript element.

## Adult-liver external dependency

The older processed-input modules require `GSE115469_Data.csv.gz` (108,065,293 bytes; SHA-256 `9d20630b98dba970ef3ec616d48a18857d94aa38e29886ab358a6dbdb10cbdde`). It was already external to the processed-input distribution. See `Collection_1/legacy_processed_inputs_v24/EXTERNAL_INPUT.json` and `fetch_external_matrix.py`, which support a local copy of the recorded original v24 archive or an explicit download. The pinned archive is:

https://huggingface.co/spaces/fc28/sesamin-atlas/resolve/7eaf4e31c5c59be5af8da4b7f5604c90f7249662/research/v24/research_reproduction_v24.zip

Recorded archive size: 141,075,041 bytes. Recorded SHA-256: `7496138851da9bb509bb7bc13a1077b7f8081abbd08d3e3a68b9ff598652e946`. The URL was inherited from provenance; this GitHub update did not freshly download or validate that remote dependency. The entire set is not self-contained offline until it is restored.

## What has and has not been checked

Archived reports distinguish executed numerical comparisons, independent formula checks and visual checks. The GitHub packaging checks only distribution integrity and helper behavior. No new model fit, raw FASTQ quantification, raw CEL processing, UMAP fit or scientific result was generated for this release. Original execution records do not establish a fresh-environment, Linux/macOS, independent wet-lab or clinical validation. Local workstation paths in preserved records describe earlier execution, not required installation locations.
