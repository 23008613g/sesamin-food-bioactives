# Executed checks and reproducibility limits

## Earlier project history

The retained v24 records describe an earlier copied-directory quantitative rerun: 62 compared tables, 240,855 rows and 845,526 finite numerical values, with external-atlas bounded/focused checks. Those records belong to the earlier versioned artifact; they are not a new full reanalysis performed during this repository preparation.

Revision26 added the matched-scale whole-organoid/Hepatocyte comparison. Revision27 expanded the presentation to eight main figures, reusing cached PCA, fixed candidate data, saved UMAP coordinates and existing analysis results. Source protocols/input hashes and the distinction between observed and missing/zero results remain explicit.

## Actually executed on 2026-10-11

Using the existing project Python 3.13.9 environment (NumPy 2.4.6, pandas 3.0.3, SciPy 1.17.1, Matplotlib 3.10.9, Pillow 12.3.0), the local repository preparation performed:

- Syntax parsing of every included Python source file without executing historical processing modules.
- v24 `run.py --help` to check its command-line entry point only.
- The input presence/checksum command with no download; it reported all 23 intentionally omitted source inputs and returned its documented missing-input status.
- The current matched-scale analysis in a new private audit output directory: all 647 implemented checks passed. Every generated CSV was byte-identical to its retained revision27 counterpart.
- A complete eight-main-figure redraw from included small processed plotting inputs: all 30 source/transform checks passed. All eight PNG files were byte-identical to the final revision27 manuscript artworks.
- Fifty-four explicit bounded packaging comparisons, including source script hashes, copied plotting-input hashes, output identity and acquisition status. The machine-readable record is `provenance/local_checks.json`.

The generated smoke-test plots and private runtime metadata are not distributed in Git. The original immutable project sources were read only. Runtime inputs and paths were relative for the current analysis/figure entry points; the private audit retains the exact local command context.

## Not executed in this preparation

No dependency was installed or upgraded. No full v24 numerical rerun, raw FASTQ/CEL processing, all-gene organoid extraction, historical directory-layout pipeline, fresh-environment installation, Linux/macOS test, browser automation, data download, GitHub creation or push was executed by the repository-preparation agent. The input-acquisition download/import path is supplied but has not been executed in this preparation. The static frontend was inspected as source/data and was not newly browser-certified.

The local alignment helper checks rectangle/label geometry. It is narrower than the earlier separate full glyph/collision tool, which is not redistributed. Exact agreement of newly generated PNGs with previously reviewed manuscript artwork supports visual preservation; it does not establish biological validity.

Metadata containing workstation locations was sanitized only in public derivatives. Original/distributed hashes and notes identify those changes. A sanitized freeze record is not claimed to be the original file; scientific protocol/input hashes and historical timestamps are retained.
