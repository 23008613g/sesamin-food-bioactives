# Code, data and licence provenance

No repository-wide open-source licence has been selected. This document records provenance; it does not grant new rights, assert ownership of third-party material or place public data under a new licence. The absence of a licence is deliberate pending an authorized choice by the rights holder.

| Component | Origin and treatment |
|---|---|
| `reproduction/` analysis code | Existing project v24 portable numerical implementation, retaining original script hashes and adaptation records. No separate source-code licence was supplied in the project archive. |
| `analysis/matched_scale/` | Project-authored revision26 analysis of five already derived numerical summaries. Scientific inputs and protocol content retained. Distributed freeze/execution metadata with local paths is explicitly a sanitized derivative; immutable originals remain private. |
| `figures/` | Project-authored revision26/27 plots and processed plotting inputs. Only import/output-directory plumbing changes for the local alignment implementation; no scientific values or plotting algorithms changed by repository preparation. |
| `figures/alignment.py`, `data/acquire_inputs.py` and repository documentation | Written for this local preparation. No new open-source licence is asserted for these additions. |
| `historical_pipeline/` | Earlier project processing, verification and plotting source. Workstation paths/contact strings are redacted in distributed copies. These files retain the original algorithmic code but are not claimed to be a newly portable full workflow. |
| `app/` | Existing static research explorer and its previously published derived display tables. Original public developer attribution is retained; contact addresses and company image assets are omitted from this source release. No source licence was present in the inspected app files. |
| Small processed tables and expected results | Numerical derivatives of the credited source datasets. Their inclusion does not relicense original source records or constitute permission for unrelated reuse. Dataset provenance remains attached. |
| Separately acquired source inputs | All v24 `inputs/` are excluded from Git. The accession/URL/hash manifest and explicit acquisition command point to the pre-existing versioned research artifact; consult each original source's terms. Public availability is not treated as a blanket licence. |
| Dependencies and earlier helper tools | Python, NumPy, pandas, SciPy, Matplotlib, Pillow, kallisto and historical optional libraries retain their own terms. Their code/binaries and the earlier figure-skill helper source are not vendored. |
| Original atlas-author code/articles | Not included. Article or author-pipeline licences are not substituted for expression-matrix terms. |

`provenance/construction_sources.json` records original/distributed hashes and scoped modifications. `provenance/code_inventory.json` records relevant inclusions, exclusions and duplicates. Existing v24 `reproduction/LICENSE_PROVENANCE.md` is retained as a historical component record; this top-level file states the boundary of the combined source release.

Before making this repository public or assigning an open-source licence, the responsible rights holder should resolve the project-code licence and confirm the intended data redistribution boundary. A private repository does not itself resolve those rights questions. No code or data is presented as licensed merely because it was previously hosted online.
