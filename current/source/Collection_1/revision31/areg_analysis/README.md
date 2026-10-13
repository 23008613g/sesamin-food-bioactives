# AREG source-level audit and recomputation

Status: executed and checked on 2026-10-11. This is a **post-hoc verification**
with prior results known, not a prospective validation experiment. The inputs
were read-only, and the 95 checks in `checks.json` passed. All outputs come from
the project Python 3.13.9 environment; see `execution.json`. No installations,
new UMAP, inferential P values, biological CIs or external uploads occurred.

## What changed versus the old PFHN manuscript

The old matched-scale analysis reweighted archived category means. This version
returns to all four original 33,538-feature MTX count matrices and the original
deposited author's H5AD annotations. It recomputes cell totals and AREG counts,
joins library accession plus barcode, and then derives each cell's
`ln(1 + 10000 * AREG_UMI / all_feature_UMI)`. It retains measured zeros and all
author-retained cells: 3,093 / 5,510 / 7,399 / 8,032, total 24,034. There were no
unmatched or duplicate retained keys and no new QC threshold. The raw barcode
matrices contained 4,051 / 6,561 / 8,634 / 9,073 cells; cells outside the author's
retained set were not relabelled as new QC failures.

**The broad Cholangiocytes category really is deposited in the author's
`obs/celltype_category` field.** The finer labels are the corresponding
`obs/celltype_scTypedb` values: Cholangiocytes 213, Ductal cells 2,965, Smooth
muscle cells 241. The earlier local extraction read both columns without
creating the broad merge. This establishes label provenance, not the validity
of mature cell identities. `cell_label_mapping.csv` retains original and
analysis labels separately, library/barcode and denominator data.

## Verified organoid results

The raw recalculation reproduces the archived matched-scale AREG differences:

| Population | Equal-library PA minus control | Pooled-cell sensitivity |
|---|---:|---:|
| Whole organoid | +0.0548982433584 | +0.0630033324451 |
| Hepatocyte-like (deposited Hepatocytes) | +0.0397972911031 | +0.0370423750846 |
| Deposited broad Cholangiocytes | -0.129474259490 | -0.043813228010 |
| Cholangiocytes + Ductal cells, excluding SMC | -0.124006656389 | -0.0484243385734 |

Values are arithmetic mean individual-cell log1p(CP10k), including zero-count
cells. The original whole-versus-hepatocyte finding is unchanged: changing only
the cell population, with transformation and weighting held fixed, does not
reverse the positive PA-minus-control direction. The original broad category
and the new SMC-exclusion sensitivity both show negative within-population
differences. Thus this result must not be written as all populations sharing a
single response direction. All four libraries are technical observations from
one differentiation, and these differences are descriptive.

| Deposited subtype | Cells | Raw AREG nonzero | Raw detection | Mean log1p(CP10k) |
|---|---:|---:|---:|---:|
| Cholangiocytes | 213 | 188 | 88.26% | 2.308258 |
| Ductal cells | 2,965 | 2,453 | 82.73% | 2.362759 |
| Smooth muscle cells | 241 | 180 | 74.69% | 2.050846 |

These pooled values give equal weight to observed cells. All ten subtypes and
all 40 library/subtype combinations are retained in the machine-readable
tables. The new SMC-exclusion group has 3,178 cells and retains high AREG raw
detection in all four libraries (80.70%, 83.70%, 76.75%, 87.96%). Excluding SMC
does not remove the descriptive difference relative to Hepatocytes. It also
does not establish cholangiocyte specificity: expression occurs in the SMC
subtype and in immune populations in the separate adult atlas.

The prior pooled-UMI CPM estimate is a different estimand and was independently
recomputed: hepatocyte-like log2[(PA mean library CPM + 1)/(control mean library
CPM + 1)] = -0.111612246449. It should not be compared to mean log1p(CP10k) as
though only cell identity changed. `areg_pooled_umi_different_estimand.csv`
retains that distinction, without interpreting it as an erroneous statistic or
as a causal composition effect.

## Verified adult-liver results

The deposited GSE115469 normalized CSV was reread directly, extracting the
unique AREG row and joining all 8,444 original cell IDs to annotation. The matrix
is author-normalized log2CPM as described by GEO. It was not normalized again.
`pct_nonzero_source` is a property of deposited values and is **not verified
raw-UMI detection**.

| Donor | Cholangiocytes | Hepatocytes | Chol minus hep mean source log2CPM |
|---|---:|---:|---:|
| P1TLH | 8 | 96 | -0.001592189 |
| P2TLH | 23 | 629 | +0.129370849 |
| P3TLH | 41 | 2,621 | +0.124695575 |
| P4TLH | 13 | 43 | +0.084744310 |
| P5TLH | 34 | 112 | +0.019541414 |

All five donors remain in the primary summary: four positive, one negative;
equal-donor mean difference +0.071351991646 source log2CPM and +8.279730470702
percentage points in source nonzero fraction. The >=20-cells-in-both-types
coverage sensitivity includes P2/P3/P5 only, all positive, with means +0.091203
and +9.018568 percentage points. P1 and P4 remain visible rather than excluded
from the primary result. This audit repeats a known coverage sensitivity and
does not claim a new prospective test.

The all-cell-type table preserves stronger expression in several immune labels.
For example, pooled means are 1.1800 in NK-like cells and 0.7625 in inflammatory
macrophages, compared with 0.10094 in Cholangiocytes and 0.01278 in Hepatocytes.
These are descriptive source-normalized values, not between-cell significance
tests. Donor coverage is supplied (NK-like cells occur in 3 donors, for example).
Neither atlas includes sesamin, so neither validates a sesamin response.

## Key deliverables

- `protocol.md`, `input_hashes_before_run.json`, `execution.json`, `run.log`:
  post-hoc plan, input provenance, environment and execution.
- `label_origin_audit.json`, `label_crosswalk.csv`, `cell_label_mapping.csv`:
  direct author-label provenance and accession/barcode mapping.
- `organoid_areg_cells.csv.gz`: recomputed individual-cell raw counts, totals,
  values and labels for bounded arithmetic reproduction.
- `areg_subtype_by_library.csv`, `areg_subtype_pooled.csv`: all ten subtypes,
  zero-inclusive means and raw nonzero counts/fractions.
- `areg_excluding_smc.csv`, `areg_broad_exclusion_pooled.csv`,
  `areg_broad_vs_hepatocyte_by_library.csv`: complete original-versus-exclusion
  comparisons, including all libraries and their hepatocyte comparisons.
- `areg_matched_scale.csv`, `areg_population_by_library.csv`: primary equal
  library and pooled-cell sensitivity estimates, with denominators.
- `areg_adult_donor_context.csv`, `areg_adult_summary.csv`,
  `areg_adult_all_celltypes_by_donor.csv`, `areg_adult_immune_context.csv`:
  all five donors, coverage sensitivity and complete cellular context.
- `organoid_old_new_comparison.csv`, `adult_old_new_comparison.csv`:
  agreement with previous estimates (organoid <=1e-12; adult <=1e-9, allowing
  archived rounding). Results were not forced to agree.
- `checks.json`, `additional_source_validation.json`, `output_hashes.json`:
  arithmetic/identity checks and archived original-object hash reconciliation.

## Run

Use the project `.venv/Scripts/python.exe` to run `recompute_areg.py`, followed
by `summarize_verified_results.py`. They resolve input files relative to the
original project directory and write only to this folder. There is no random
step. Do not confuse this bounded AREG reanalysis with full raw-read processing,
independent reannotation, or a rerun of all 25 candidates or all project modules.

## Suggested manuscript wording

Methods: "For this post-hoc audit, we reread the four deposited raw count
matrices and the author's H5AD annotations, joined cells by library accession
and barcode, and recalculated AREG log1p(CP10k) using each cell's total across
all 33,538 measured features. The deposited broad Cholangiocytes category was
retained, while its three source subtypes were reported separately; exclusion
of Smooth muscle cells was examined as an additional sensitivity analysis. All
technical libraries and zero-count cells were retained, with no cell-level
significance tests or biological confidence intervals."

Results: "AREG raw counts were nonzero in 188/213 source-labelled
Cholangiocytes, 2,453/2,965 Ductal cells and 180/241 Smooth muscle cells. The
first two labels remained higher than Hepatocytes in each library when pooled
after exclusion of Smooth muscle cells. Raw-source recalculation reproduced
positive matched-scale PA-minus-control differences for the whole organoid
(+0.05490) and hepatocyte-like subset (+0.03980). These descriptive findings
support careful specification of cellular source and expression summary,
rather than a cholangiocyte-specific or sesamin-mediated mechanism."
