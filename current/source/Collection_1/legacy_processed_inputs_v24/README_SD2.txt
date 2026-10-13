Scientific Reports - Supplementary Data 2
Processed-input reproducibility package v24; source bytes preserved.

The research_reproduction_v24/ directory contains the original code, expected
tables, protocols, dependency record and observed verification records, extracted
from the original v24 ZIP. Exactly one original input is not embedded:
research_reproduction_v24/inputs/external/GSE115469_Data.csv.gz (108,065,293 bytes).
It remains publicly available inside the commit-pinned archive in EXTERNAL_INPUT.json.
All 118 other original v24 members are preserved byte-for-byte. This is not a
self-contained offline distribution until the external matrix has been restored.

After extracting both Supplementary_Data_1.zip and Supplementary_Data_2.zip into
one chosen directory, restore the single dependency without changing its bytes:
  python fetch_external_matrix.py --archive PATH_TO_ORIGINAL_V24_ZIP
or, for a checked network download:
  python fetch_external_matrix.py
The script verifies the original ZIP and the extracted matrix by size and SHA-256.
It does not overwrite a different existing destination. No dependencies are needed
for this retrieval step beyond Python's standard library. The local-archive path
was tested during packaging; no new remote download was performed in this step.

Then read research_reproduction_v24/README.md and ENVIRONMENT.md before running
its existing run.py and verify.py in a suitable scientific environment. Those
programs and historical records are unchanged. This packaging step did not rerun
the research, install packages or reprocess original FASTQ/CEL data. The matched-
scale module is supplied separately in Supplementary Data 1.

FILE_MANIFEST_SD2.csv lists this ZIP's members excluding itself. Source lineage is
in Supplementary Data 1's SOURCE_MEMBER_MAP.csv. Original provenance statements
are historical; this wrapper describes the current two-ZIP distribution.
