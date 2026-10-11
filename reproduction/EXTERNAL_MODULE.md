# Independent adult-liver cellular context: GSE115469

This fourth module executes `code/external/analyze_external.py` unchanged from the new manuscript analysis, using five package-local input files. The frozen protocol and scale-decision records are beside the script. `run.py` calls this module after the three earlier analyses; a separate invocation is also available:

```text
python code/external/analyze_external.py --input-dir inputs/external --output-dir external_results
```

Use a new output folder. `run.py` additionally compares every external CSV output against the original reference analysis; the standalone module executes its own 35 bounded data/calculation checks.

| Input | Meaning |
|---|---|
| `GSE115469_Data.csv.gz` | 20,007 deposited feature rows × 8,444 cells. The source is author-normalized **log2CPM (GEO label)**, not raw UMI. The full matrix is read, checked, and used to extract the frozen candidates. |
| `GSE115469_CellClusterType.txt.gz` | Unique cell IDs, public donor labels P1TLH–P5TLH, barcodes, 20 author clusters and cell labels. Six hepatocyte subtype labels are combined for the prespecified broad comparison. |
| `GSE115469_family.soft.gz` | Public GEO sample/protocol metadata used to verify donor titles and the deposited scale. |
| `human_candidates_all25.csv` | Fixed original human candidate family; no reselection on this atlas. |
| `original_candidate_gene_coverage.csv` | Frozen prior gene-ID/symbol/raw-feature-alias mapping, used only for exact allowed aliases. |

There are 19 matched candidates and six absent features, retained as missing. The primary comparison calculates cholangiocyte-minus-hepatocyte differences within each of five donors, giving donors equal weight in descriptive summaries. The prespecified coverage sensitivity retains donors with at least 20 cells in both types: P2TLH, P3TLH, and P5TLH. Full donor, subtype, broad-type, all-cell context, missingness, and zero results remain in outputs. No hypothesis test, confidence interval, or de novo cell annotation is computed.

The nonzero percentage counts deposited values >0. It is not verified raw-UMI detection. No new normalization, raw-count inverse transformation, FASTQ processing, or original scran workflow is executed. Means are on the source log2CPM scale and are not directly comparable magnitudes to organoid log1p(CP10k). This separate healthy adult-liver dataset tests the generality of cellular localization; it cannot validate a sesamin response or hepatoprotection. Contrary and sparse results must be retained.

The original article and author R pipeline were inspected by the external-analysis module's source audit, but their full text and source code are **not** redistributed here. Source references remain in `code/external/source_scale_decision.md`; the numerical input inventory records hashes and provenance. The author pipeline commit reviewed locally was `9bff516c9427626f8ece79da4370d5be35012e5c`; it is not a dependency of this rerun.
