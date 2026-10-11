# Code, data, and licence provenance

This manifest documents origins. It does not assign a blanket licence to third-party material, claim copyright ownership, or create a new permission grant.

| Component | Origin and treatment |
|---|---|
| `code/robustness.py`, `code/robustness_check.py`, `code/rat.py`, `code/rat_exact.py`, `code/organoid.py` | Existing project stage 15 scientific scripts, adapted for package-local paths. Original and distributed hashes and change scopes are in `provenance/script_adaptations.json`. The source project did not supply a separate licence declaration for these helper scripts during construction; no third-party licence is invented. |
| `code/external/analyze_external.py` | Project stage 24 external-analysis script, copied unchanged. It is newly written analysis code, not the original atlas authors' R pipeline. The same project licence status applies. |
| `run.py`, `verify.py`, documentation and manifests | Written for this package's construction and numerical verification. No repository-wide open-source licence has been selected or asserted here. |
| NumPy, pandas, SciPy, Python | External runtime dependencies, identified by version in `requirements.txt` and `ENVIRONMENT.md`. Their source/binaries are not vendored in this artifact. Their own distributions retain their own licence notices. |
| Human counts and rat signals/platform annotation | Publicly deposited GSE248903, GSE42693 and GPL1355 numerical records. Source accessions and URLs are included. Public availability is recorded without replacing the source's terms with a package licence. |
| Mouse gene TPM and estimated counts | Previously generated project numerical derivatives of PRJCA021826/CRA013797, mirrored as PRJDB40955/DRP018055. The mirror denotes the same study, not an independent dataset. Raw reads and third-party quantification executables are not included. |
| Organoid candidate counts, QC totals, annotations and references | Previously generated project numerical derivatives of GSE207889. Cell labels and saved coordinates originate from the authors' deposit. No author manuscript or full raw h5ad/MTX is included. |
| External adult-liver matrix and cell annotations | Deposited GSE115469 author-normalized numerical inputs. The original article has DOI 10.1038/s41467-018-06318-7. The article XML/PDF and BaderLab R code are not included, so no author code licence is represented as a licence for this project. |
| Ensembl mappings | Frozen Ensembl 116 orthology and rat Entrez cross-reference tables from the earlier analysis. Their provenance is retained as mapping inputs; no live download or invented version update occurs. |

Local private paths, private execution logs, full-text articles, manuscript Word files, reviewer materials, credentials, and service configuration files are excluded from this subset. A separate local construction audit retains unredacted original locations and command evidence. Any later hosting or repository-level licence choice must preserve these component distinctions; this document itself is not a representation that all source materials share one licence.
