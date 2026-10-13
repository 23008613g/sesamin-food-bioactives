# Sesamin hepatic-expression research reproduction, version 24

This versioned numerical/code archive accompanies the cross-model expression and cellular-context research. It is separate from the interactive explorer's software display.

- [Download the reproduction archive](https://huggingface.co/spaces/fc28/sesamin-atlas/resolve/main/research/v24/research_reproduction_v24.zip)
- [SHA-256 checksums](https://huggingface.co/spaces/fc28/sesamin-atlas/resolve/main/research/v24/SHA256SUMS)

The archive contains relative-path Python programs, 23 numerical/mapping inputs, data dictionaries, frozen protocols, environment requirements, provenance/hashes, and observed verification records. It rebuilds human/mouse candidate robustness, rat Welch/BH and exact sensitivity, library-aware organoid summaries, and the independent GSE115469 adult-liver cellular-context analysis.

A copied-directory execution compared 62 tables, 240,855 rows, and 845,526 finite numeric values against retained references. The external module also passed 35 bounded checks and 24 focused checks. Inputs total 137,503,658 bytes. The tested existing environment was Python 3.13.9, NumPy 2.4.6, pandas 3.0.3, and SciPy 1.17.1; no fresh installation or cross-platform run was performed.

After extracting, read README.md and run `python run.py` from an environment with those dependencies. This reproduces calculations from processed matrices and candidate cell counts; it does not rerun FASTQ/CEL processing, original all-gene UMI extraction, embedding, or upstream candidate selection. It does not establish clinical benefit or validated hepatoprotection. All negative and missing findings remain in the numerical tables.

Dataset accessions, exact input hashes and the archived Ensembl 116 query provenance are inside the package. The archive contains no article full text, manuscript Word files, private review materials, credentials, or workstation paths. No blanket licence is assigned to third-party data; consult LICENSE_PROVENANCE.md. This download location is a versioned research artifact, not a DOI or a peer-review certification.
