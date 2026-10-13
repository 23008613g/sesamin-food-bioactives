# Sample metadata evidence log

Audit date: 11 October 2026. This is a bounded source review, not an author reply, independent experimental validation or new statistical analysis. Source retrieval, including failures, is recorded in `source_manifest_validated.json`; saved original source fields are in `sample_metadata_original_fields.json`. `sample_metadata_audit.csv` contains 77 rows: human 9, mouse 9, rat 20, PHH 30, selected organoid 4, adult liver 5.

## Human HepG2: retained with explicit limits

Fresh [GSE248903](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248903) and all nine GSM records were retrieved through the official GEO text endpoint. All nine sample snapshots are byte-identical to the earlier stage 25 snapshots. `Sample_source_name_ch1` and `Sample_description` explicitly distinguish the Control, Model and Treat libraries. The common `Sample_treatment_protocol_ch1` describes initial adhesion for 24 h, OA treatment for 24 h and another 36 h of sesamin for Treat. It deposits OA as `400  M` and sesamin as `16   M`: missing unit prefixes are not silently repaired to micromolar. The exact protocol is preserved in JSON.

The record does not establish matched total collection times in Control/Model/Treat, whether OA continues during the additional period, washing or medium replacement, vehicle matching, or independent culture identities. The study-level summary reports experiments by the depositors, but these are not new experiments or independently validated findings in our reanalysis. No linked publication is supplied by the fresh GEO record; similarity of author names or topic is insufficient to import dose details from another paper.

**Decision:** retain as exploratory library-level candidate-selection data. Describe a deposited treatment-associated contrast, not a controlled causal estimate of sesamin benefit. Do not add a biological-replicate DESeq2 claim while units remain unresolved.

## Mouse: project intervention known; sample-to-condition mapping unresolved

Fresh official [NGDC project PRJCA021826](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA021826) and [ENA study mirror DRP018055](https://www.ebi.ac.uk/ena/browser/view/DRP018055) both describe primary mouse hepatocytes and normal, PA/OA and Ses-treated PA/OA conditions. This is real project-level support for the study's intervention subject. However, the nine NGDC BioSamples and all nine ENA sample and experiment records give labels `ctrl1–3`, `HFD1–3` and `sess1–3` without per-sample treatment fields or an explicit correspondence table. Names are suggestive but cannot independently establish the intervention assignment.

All nine official ENA runs link to original SAMC identifiers in sample titles/attributes and CRX identifiers in experiment titles. The BioSamples also explicitly link `ngdc_project_id=PRJCA021826`. The mirror therefore supplies access to the same nine libraries, not a second cohort or additional biological replicates. Source material is hepatocyte, male, age 8 weeks. Animal IDs, isolation IDs, cell-culture units, pooling, concentrations, duration and vehicle remain absent. `HFD` must not be expanded into a demonstrated in-vivo high-fat-diet group.

[Bai et al., 2024](https://doi.org/10.1016/j.bcp.2024.116240) is biologically compatible with this project, but the accessible abstract does not establish the accession-to-library linkage; the publisher full-text endpoint returned 403 in this audit. No concentration or group mapping is borrowed from it.

**Decision:** retain source-labelled numerical exploration and complete orthologue tables, preferentially in supplementary material. Exclude this mouse comparison from claims of confirmed sesamin-response replication and from the abstract's biological result statement. The manuscript may state that a mouse accession describing a sesamin study was explored using deposited labels, with the mapping unresolved.

## Rat: dietary comparison retained

The fresh [GSE42693](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42693) series record and the original [Ide et al. article](https://doi.org/10.3177/jnsv.55.31) establish male Sprague-Dawley rats, lignan-free or 0.2% sesamin/episesamin/sesamolin diets, 15 days and 7 rats per group. Five individual liver arrays per group were selected after excluding the highest and lowest terminal body weights. The 20 sample rows come from the existing original GEO series-matrix headers, with their source-file hash recorded; fresh article and series metadata corroborate the design. No array-level treatment ambiguity was found.

**Decision:** retain sesamin versus control as the dietary comparison; keep animal selection explicit. The other two lignan groups may remain in the declared normalization sensitivity analysis, but are not sesamin replicates. Matrix-scale and probe-statistic checks are separate analysis work and are not claimed as executed by this metadata audit. Physiological context is documented in `nutrition_context_evidence.md`.

## PHH: all 30 headers retained as context

The fresh [GSE247407](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247407) sample list exactly matches the 30 complete official sample headers archived on 10 October. The source describes five patient labels, three repeats per condition and seven-day Control/FFA culture. Cached original fields have been copied, not represented as 30 new requests. F95 and F95P remain separate labels whose patient relationship is unresolved. No sesamin arm exists in the selected matrix.

**Decision:** context-only descriptive comparison, not sesamin-response validation. Do not treat the labels as five independently established donors or all 30 libraries as independent biological units.

## Organoid: correction and technical units explicitly handled

Four fresh GSM records confirm W01-derived orbital-shaker HLOs, two CTRL-PA and two PA 500 micromolar libraries, harvested day 25. The GEO protocol states five treatment days, while the [corrected source article methods](https://doi.org/10.1038/s44318-024-00308-w) describe four days from day 21 to day 25. This discrepancy is retained rather than silently reconciled. The original article's Statistical analysis section defines replicated scRNA-seq samples as suspensions arising from the same differentiation experiment. A general library-preparation sentence uses “biological replicates”, but it does not override the explicit experimental-unit definition or establish separate differentiation batches for these four libraries.

**Decision:** cellular reference with technical-library descriptive summaries only. No sesamin arm, no independent sesamin validation, no cell-level replication claim. `correction_note.md` documents the corrected preparation protocol and access limits.

## Adult liver: five donor-labelled references

Fresh [GSE115469](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115469) and all five GSMs establish Patient 1–5 liver specimens. The submitted processing description names scran normalization and author-supplied log2CPM values. All donors remain in the reference analysis; their sparse cell subgroups must remain visible.

**Decision:** independent donor context for AREG expression, not a drug-response test. Nonzero author-normalized values must not be labelled raw-UMI detection.

## Access and execution limits

The normal web viewer returned challenges for some GEO pages; their public GEO text endpoint was accessible. GSA CRA013797 returned a browser security challenge; no challenge was bypassed. The DDBJ page returned a client shell; ENA's official mirror provided the inspectable XML. An ENA run-table timeout succeeded on one retry. The Bai publisher text was not obtained. Hess correction XML and publisher HTML were read, but image/PDF requests returned HTML instead of the requested media, and have been classified as unusable media responses. Their filenames do not imply successful image/PDF verification.

The project's existing Python 3.13.9 was used for standard-library retrieval, parsing and coverage/hash checks. An existing document interpreter with PyMuPDF extracted and rendered the original rat PDF for visual verification of Table 2. No dependencies, environment changes, raw-data edits, author contact or uploads were performed.
