# Ensembl 116 full-genome orthology freeze

All source data queried from the June 2026 archive without a gene-list filter.

Strict analysis map: `human_rat_one2one_highconfidence_unique_entrez.tsv`.
Use stable Ensembl and Entrez IDs as keys; symbols are labels only.
All rejected/ambiguous/unmapped rows remain in raw data and mapping audit tables.

## Counts

- human_rat_raw_rows: 91806
- rat_entrez_raw_rows: 57760
- human_rat_rows: 91806
- human_ids: 86411
- human_ids_with_rat_ortholog: 18606
- rat_ids_with_human_ortholog: 20037
- human_ids_with_multiple_rat_ids: 1282
- rat_ids_with_multiple_human_ids: 1112
- rat_entrez_rows: 57760
- rat_ids_total: 43360
- rat_ids_with_entrez: 23497
- entrez_ids_total: 23793
- rat_ids_with_multiple_entrez: 522
- entrez_ids_with_multiple_rat_ids: 312
- rat_entrez_unique_bidirectional: 22905
- one2one: 16449
- one2one_highconfidence: 13733
- strict_analysis_map: 13531

## Provenance

Queries, raw TSV (including completion stamp), source attribute definitions, BioMart registry, dataset definitions, HTTP metadata, retrieval UTC times, and SHA256 hashes accompany this report.
No expression statistics or candidate-gene selection were performed by this script.
