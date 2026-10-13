# Food-derived sesamin and hepatic lipid metabolism: a computational assessment of cross-study transcriptomic evidence

Analysis code and numerical evidence for a secondary analysis of public transcriptomic studies of food-derived sesamin. The current manuscript finds no stable cross-model correspondence of transcriptional effects and distinguishes candidate readouts from their cellular sources. HepG2 is a tumour-derived model; these results do not establish a protective mechanism or a sesamin response in normal human hepatocytes.

**Current manuscript snapshot: 2026-10-13 — [v2026.10.13](https://github.com/23008613g/sesamin-food-bioactives/releases/tag/v2026.10.13).** This is a research-code version, not an accepted-publication claim. Authors: Tingting ZHENG and Honglie ZHANG, Infinia Intelligent Health (Hong Kong) Limited.

## Start here

- [Current code and result index](current/README.md): all 53 Python/R/shell scripts from the current data archive, run documentation, and selected small numerical tables.
- [Reproduction guide](docs/REPRODUCING_CURRENT.md): exact input archive, environments, commands and verification boundaries.
- [Dataset access and roles](docs/DATA_ACCESS.md): public accessions, model limitations and external inputs.
- [Version history](CHANGELOG.md): which analyses supersede the earlier repository snapshot.

To verify the browseable files using Python's standard library:

```sh
python tools/check_current_sources.py
```

Download `Supplementary_Data.zip` from the [versioned release](https://github.com/23008613g/sesamin-food-bioactives/releases/tag/v2026.10.13), or use the identical attachment supplied with the manuscript. Then import it into a new directory:

```sh
python tools/acquire_current.py --archive /path/to/Supplementary_Data.zip --dest .local/current_data
python .local/current_data/restore_exact_duplicates.py --restore
```

An explicit `--download` option is available once the release is accessible to you. A private repository requires authorized access; the helper does not read credentials or change authentication. Neither acquisition command installs dependencies or executes scientific analyses. The archive includes frozen inputs for the whole-transcriptome comparison, but the adult-liver GSE115469 matrix is an existing separately acquired dependency. Read the reproduction guide before running an analysis.

## Current analytical modules

| Module | Location inside the complete extracted archive | Scope |
|---|---|---|
| Candidate comparison and revised AREG context | `Collection_1/revision31/` | Fixed candidates; raw-count and annotation sensitivity; technical-library boundaries |
| Rat whole-transcriptome analysis | `Collection_2/scripts/analyze_rat_genome.py` | 14,604 eligible genes; gene and GO tests; normalization sensitivity |
| limma sensitivity | `Collection_2/scripts/limma_sensitivity.R` | Separate inferential specification with recorded R/limma versions |
| Full human–rat comparison | `Collection_2/scripts/compare_cross_species.py` | 7,912 strict one-to-one shared genes; normalization-dependent, weak GO correspondence |
| Signal-intensity diagnostic | `Collection_2/revision43/diagnostics/` | Signal distribution of limma hits and stability across normalizations |
| Main numerical reproduction | `Collection_2/scripts/reproduce_analysis.py` | Runs three original functions and compares 62 result tables with the frozen reference |

Single-cell and organoid analyses provide cellular context. Their datasets do not include sesamin exposure; four organoid libraries come from one differentiation. Rat Areg has low signal and a wide interval, so an estimate near zero is not evidence of no biological response. No new wet-lab experiment, prospective validation, clinical efficacy study or validated AI predictor is supplied.

## Evidence and verification

Existing execution records document a same-computer, fresh-directory reproduction of 62 whole-transcriptome result tables, with maximum numerical difference 0.0 at written CSV precision. This refers to the archived original execution, not a new scientific rerun during GitHub preparation. Repository checks cover source/archive hashes, extraction safeguards and file mapping. They do not establish biological validity or cross-platform reproducibility. See [execution history](EXECUTION_HISTORY.md).

The earlier `reproduction/`, `analysis/matched_scale/`, `figures/`, `historical_pipeline/` and `app/` directories remain for traceability. Their code and numbering describe the 2026-10-11 snapshot and are not the current complete manuscript pipeline. The retained app is an older display snapshot. The [old overview](docs/README_2026-10-11.md) is explicitly historical.

## Citation and rights

Use [CITATION.cff](CITATION.cff) and cite the version/commit used. No DOI or published-paper identifier has been assigned here. A later archival DOI can be added after it actually exists.

No repository-wide open-source licence has been selected. Public visibility alone would not grant a blanket reuse licence. Source datasets and dependencies retain their own rights; see [LICENSE_PROVENANCE.md](LICENSE_PROVENANCE.md). This repository does not distribute the manuscript, cover letter, business-card assets or private editorial correspondence.
