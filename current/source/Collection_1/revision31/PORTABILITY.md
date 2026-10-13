# Reproduction scope

Use Python 3.13.9 with the recorded existing versions; these instructions do not install anything.

1. `python candidate_analysis/run.py` recomputes selection, adjacent thresholds, all 25 human candidates and 15 mouse directions. Human deposited counts are included; mouse gene TPM is previously quantified, not a raw-read rerun.
2. `python rat_analysis/reanalyze_rat.py` recomputes the deposited-matrix statistics, probe audit, both normalizations, all 252 label allocations and both correction families. All required inputs are included. No CEL reprocessing occurs.
3. `python areg_analysis/reproduce_from_cell_extracts.py` verifies the four-population by four-library organoid summaries and the five-donor cholangiocyte/hepatocyte summaries (not every subtype or immune-category output) from the included per-cell extracts. This is a processed-input check, not the raw-matrix extraction step.
4. `areg_analysis/recompute_areg.py` is the unchanged script actually used for the raw organoid and adult matrix audit. It expects the project-relative layout documented in `external_raw_inputs.json`, plus the auxiliary inputs in `raw_input_layout/`. Copy the contents of `raw_input_layout/` into the `next_steps/` root, preserving subdirectories. The original author H5AD and large expression matrices are external and must be retrieved from their public sources and hash-checked. To rerun it, create `next_steps/pfhn_strategy_revision_31/areg_analysis/`, copy the module there, and place the listed sources relative to `next_steps/`. Preserve the raw files read-only. The 3.24 GB H5AD is uncompressed; GEO may distribute a compressed version that must be decompressed before comparing this hash.

Original source-level scripts and their execution records are preserved. Portable versions of candidate and rat scripts change only input routing; their executed result checks are recorded separately. The included cell extracts were generated in this revision from raw organoid counts/original author annotation and the adult author-normalized matrix. They are not independent validation data. No external download is needed for steps 1–3.

`legacy_processed_inputs_v24.zip` at the supplement root preserves the prior processed-input workflow and retrieval instructions for the original adult matrix. It predates this revision. Its results must not be described as new raw-read or raw-UMI analyses.

The source rat article is cited by DOI and checked locally, but its complete PDF is not redistributed. Portable rat execution checks 127 computational/input conditions; the original source-level execution additionally checked the source-PDF hash (128).
