# Food bioactives in fatty liver research: assessing molecular evidence through a sesamin case study

**食源性活性成分在脂肪肝研究中的分子证据评估：以芝麻素为例**

Research code for an exploratory secondary analysis of public transcriptomic studies. This repository brings together the processed-input reproduction modules, the later matched-scale AREG analysis, all eight main-figure generators, historical processing and supplementary-figure source, and the static research explorer. It contains no new biological experiment, clinical efficacy assessment or validated AI prediction.

The manuscript uses sesamin as a case study for evaluating molecular evidence about food bioactives. The studies answer different questions: healthy rat dietary liver and an adult-liver atlas are contextual references, and the organoid/primary-hepatocyte datasets have no sesamin intervention. Opposite, missing, zero and unstable findings are retained.

## Start here

Python 3.13.9 with the exact existing versions in `requirements.txt` was used in the project. No new environment installation or full scientific rerun was performed while preparing this repository. The static app needs only an HTTP server and a browser.

```sh
# Complete eight-main-figure redraw from included small processed inputs
python figures/build_figures.py --only all

# Bounded matched-scale update from the five included derived summary tables
python analysis/matched_scale/run.py --output-dir analysis/matched_scale/results_local

# Serve the existing explorer locally; no remote publishing occurs
python -m http.server 8000 --directory app
```

For the four-module quantitative reproduction, obtain source inputs separately. The repository intentionally excludes the original count/expression matrices, raw FASTQ/CEL/MTX/h5ad files and third-party tool binaries.

```sh
# Inspect which of the 23 separately acquired numerical inputs are missing
python data/acquire_inputs.py

# Import from a local copy of the existing checksummed v24 archive
python data/acquire_inputs.py --archive /path/to/research_reproduction_v24.zip

# Alternatively, explicitly download that existing public artifact
python data/acquire_inputs.py --download

# Run the processed-input analyses and compare against retained reference tables
python reproduction/run.py --output results_local
```

The acquisition commands do not grant data rights. Consult the source terms and `LICENSE_PROVENANCE.md`. Imported inputs are ignored by Git. The full reproduction refuses to overwrite an existing output directory. It rebuilds human/mouse robustness, rat Welch/BH and exact-label sensitivity, organoid library summaries, and the independent adult-liver analysis; it does not claim a complete rerun from raw sequencing.

## Code map

| Location | Included work | Execution status and scope |
|---|---|---|
| `reproduction/code/` | Human/mouse configurations, independent formula checks, rat primary/exact analysis, organoid aggregation/decomposition, adult-liver analysis/checks | Existing v24 numerical implementation; relative paths; source inputs acquired separately |
| `analysis/matched_scale/` | Whole-organoid versus Hepatocyte means on the same per-cell log1p(CP10k) scale, equal-library and pooled-cell weights | Current portable update; all 25 candidates including missing AOC4P |
| `figures/` | Eight main figures, including complete-candidate and cached-PCA redraws | Portable code and small processed inputs; no PCA refit, annotation or inference |
| `historical_pipeline/` | Earlier candidate selection, PHH processing, kallisto orchestration/reference construction, FASTQ QC, aggregation, organoid barcode/UMI linkage, independent checks and supplementary plotting | Source archive requiring the original directory layout and original inputs; not claimed newly runnable end to end |
| `app/` | Static HTML/CSS/JavaScript explorer and previously published derived display tables | Display snapshot v1.1.3; predates the matched-scale revision and does not replace current manuscript results |
| `provenance/` | Included-source hashes, complete inclusion/exclusion inventory and v24 archive inventory | Distinguishes original source hashes, distributed modifications and omitted files |

The earlier figure-tool helper library is not vendored. Current plots use the small project-local rectangle/label check in `figures/alignment.py`. This is narrower than the historical full glyph/collision audits and is identified as such. Supplementary source scripts are preserved with their historical dependencies rather than silently presented as a freshly tested pipeline. See `MODULE_MAP.md` for the figure routing and scientific processing chain.

## Data sources and scientific boundaries

| Dataset | Role | Boundary |
|---|---|---|
| [GSE248903](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248903) | Human HepG2 count-based candidate context | Treatment timing and independent culture replication remain unresolved |
| [CRA013797](https://ngdc.cncb.ac.cn/gsa/browse/CRA013797), PRJCA021826; mirror PRJDB40955/DRP018055 | Mouse hepatocyte comparison | The mirror is the same study; treatment grouping is provisional |
| [GSE42693](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42693), GPL1355 | Independent dietary rat liver | Healthy-liver context, five selected animals per group; distinct lignans are not extra sesamin replicates |
| [GSE247407](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247407) | Primary-hepatocyte fatty-acid context | No sesamin treatment; donor-label relationships unresolved |
| [GSE207889](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207889) | Organoid localization and library-aware expression | Four technical libraries from one differentiation; no sesamin treatment; original labels and saved coordinates retained |
| [GSE115469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469) | Independent adult-liver donor/immune context | Deposited normalized expression, five donors; contrary donor retained; no intervention |

Frozen orthology is Ensembl release 116; exact accession, URL, archive-member and checksum records are in `data/acquisition_manifest.json` and the historical v24 provenance. Public accessibility does not imply a repository-wide data licence.

The matched-scale AREG update is positive in both whole organoids (+0.054898) and Hepatocytes (+0.039797) using equal-library mean log1p(CP10k). The older pooled-UMI ratio must not be interpreted as a population-only sign reversal. Human–mouse configuration stability is descriptive and not a probability of independent replication. Rat VNN1 direction differs from the cell datasets, and the post-hoc exact-label BH sensitivity finds no candidate at q<0.05.

## Versions and checks

This local source preparation is dated **2026-10-11**. It combines numerical reproduction v24, the matched-scale update from revision26, and the eight-main-figure presentation from revision27. No GitHub publication, DOI, accepted-paper status or new licence is implied by this date.

Historical v24 records report a copied-directory comparison of 62 tables, 240,855 rows and 845,526 finite numerical values. Those are earlier execution records, not a new result produced by repository construction. Current bounded checks, their exact scope and unexecuted paths are listed in `EXECUTION_HISTORY.md` and `provenance/local_checks.json`.

Source protocols and scientific input hashes are preserved. Any distributed metadata with redacted local paths is a separately identified public derivative; it is not described as the original immutable freeze. Originals remain in the private project audit. No manuscript Word/PDF files, email addresses, keys, virtual environments, third-party tool source or company image assets are included.

## Licence status

No repository-wide open-source licence has been selected. The pre-existing code did not provide one, so this preparation does not silently relicense it. Existing dependency and original-data rights remain separate. See `LICENSE_PROVENANCE.md`; obtain the relevant permission before reusing material beyond the rights already applicable to it.
