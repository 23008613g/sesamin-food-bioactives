CURRENT SUPPLEMENTARY DATA DISTRIBUTION
Food-derived sesamin and hepatic lipid metabolism: a computational assessment of cross-study transcriptomic evidence

This ordinary ZIP contains the numerical results, analysis code, frozen annotations,
metadata, source inputs and reproduction records used in this study. Collection_1
contains candidate/cellular-context analyses; Collection_2 contains whole-transcriptome
rat and cross-study analyses plus the signal-intensity diagnostic. These map to the
previous Supplementary Data 1 and Supplementary Data 2, respectively.

CONTENTS AND SIZE REDUCTION
The original five archives contain 747 outer members. Their historical nested
legacy_processed_inputs_v24.zip is expanded into
Collection_1/legacy_processed_inputs_v24/ (122 extracted members). All unique
analysis code, numerical tables, annotations and archived source input contents
are retained. No expression matrix, result CSV/TSV or analysis script was removed.
The source maps distinguish stored files, exact aliases, omitted visual outputs,
three omitted fit-cache outputs, and the replaced historical ZIP container.
There are 626 physically stored scientific/provenance files,
182 exact duplicate aliases, and 60 omitted visual/fit files.
The archive is a curated distribution, not a byte-complete copy of the old ZIPs.

START HERE: RESTORE EXACT DUPLICATE PATHS LOCALLY
Extract this ZIP into a NEW directory. The two Collection directories must remain
separate. The cited supplementary result paths are physically present under their
respective Collection prefix. Before running the original analysis scripts, run:

  python restore_exact_duplicates.py --restore

The standard-library helper verifies canonical SHA-256 hashes, copies identical
content to missing original paths, and refuses to overwrite different existing
content. It makes no network requests, installs nothing and runs no analyses.
Running it without --restore only checks the aliases. Do not run analyses until
the aliases required by the original directory layout have been restored.
SOURCE_MEMBER_MAP.csv lists every original outer member and extracted nested
member, original SHA-256, current path, status, and alias/omission reason.
PACKAGE_MANIFEST.csv describes the physically supplied current archive files.

EXISTING EXTERNAL ADULT-LIVER INPUT
The historical processed-input package already excludes one source matrix:
research_reproduction_v24/inputs/external/GSE115469_Data.csv.gz
Original bytes: 108065293
SHA-256: 9d20630b98dba970ef3ec616d48a18857d94aa38e29886ab358a6dbdb10cbdde
This packaging does not introduce a new omission of that matrix. The retained
Collection_1/legacy_processed_inputs_v24/EXTERNAL_INPUT.json gives its original
commit-pinned archive URL, size and checksum. The retained fetch_external_matrix.py
supports retrieval from a supplied local original archive or a checked download:

  python Collection_1/legacy_processed_inputs_v24/fetch_external_matrix.py --archive /path/to/original_v24_archive.zip

For the network option, consult that script and EXTERNAL_INPUT.json first. The
recorded archive URL is:
https://huggingface.co/spaces/fc28/sesamin-atlas/resolve/7eaf4e31c5c59be5af8da4b7f5604c90f7249662/research/v24/research_reproduction_v24.zip
Archive bytes: 141075041
Archive SHA-256: 7496138851da9bb509bb7bc13a1077b7f8081abbd08d3e3a68b9ff598652e946
No remote download was performed or newly verified during this packaging. The
entire distribution is not self-contained offline until that dependency is restored.
Other raw-data dependencies outside the documented processed-input reproduction
scope remain subject to their original protocols and accessions. UMAP fitting,
raw FASTQ quantification and raw CEL preprocessing were not reproduced here.

REPRODUCTION ENTRY POINTS
After restoring aliases, Collection_2/README_reproduction.txt describes the
whole-transcriptome 62-table reproduction with its existing Python/R dependencies.
Collection_1/revision31/PORTABILITY.md describes the revision31 modules.
For the legacy processed-input modules, first restore the adult-liver matrix,
then read Collection_1/legacy_processed_inputs_v24/research_reproduction_v24/README.md
and ENVIRONMENT.md. Existing source README files, historical journal names,
figures and container filenames are retained as provenance. This CURRENT guide
supersedes their distribution descriptions, including claims that all inputs are
embedded and instructions referring to an unexpanded legacy ZIP.

OMITTED VISUAL AND FIT OUTPUTS
Archived PNG, TIFF, JPEG, SVG, EPS and PDF renditions are not supplied in this
data ZIP. Current figures are supplied with the manuscript and Supplementary
Information. Plotting scripts and numerical figure-source tables are retained.
Many quantitative plots, including limma plotSA diagnostics, can be regenerated;
not all historical artwork or the two source-publication reference JPEGs are
reconstructible from analysis scripts. Their provenance records remain supplied.
Three rat_limma_*_fit.rds files are omitted as regenerable model-fit caches:
scripts/limma_sensitivity.R writes these fits, while full CSV results remain supplied.
OMITTED_VISUALS_AND_FIT_OBJECTS.csv lists every omitted rendition/cache precisely.

VERIFICATION BOUNDARIES
Packaging verifies member hashes, archive integrity, cited relative paths and
local exact-alias restoration. It does not rerun the statistical analyses or
establish independent biological validation. Existing execution records retain
their original dated scope; new package checks do not imply new experiments.
