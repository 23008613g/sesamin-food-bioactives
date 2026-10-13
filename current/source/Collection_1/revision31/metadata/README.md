# Stage 31 bounded metadata and nutritional-context audit

Read `metadata_evidence_log.md` and `sample_metadata_audit.csv` first. This folder contains source evidence, a 77-row sample table, source-field extracts, a corrected-method note and original rat nutritional context. It does not perform expression quantification or statistical reanalysis.

Rebuild the audit from the workspace root with the existing project environment:

```powershell
& 'research_planning/hepatic_lipotoxicity_2026-10-07/.venv/Scripts/python.exe' -X utf8 'research_planning/hepatic_lipotoxicity_2026-10-07/next_steps/pfhn_strategy_revision_31/metadata/build_sample_audit31.py'
```

The builder expects the preserved prior-stage rat matrix metadata and PHH header extract at the recorded project-relative locations. Cached header copies are included in `sources/`, and original file SHA256 values are recorded in `validation.json`. This is the project-local rebuilding command, not a claim of a fresh end-to-end portable raw-data analysis.

`fetch_metadata31.py` and `fetch_followup31.py` document bounded retrieval of official endpoints. Re-running them intentionally fetches then-current records. `source_manifest_validated.json` preserves the original HTTP outcome and adds content classification; HTTP 200 alone is not evidence that the requested content was obtained. Five publisher image/PDF responses were HTML and are flagged as `html_error_not_requested_media`; they should not be treated as figures or PDFs for a submission package. The valid rat PDF begins with a PDF signature and its Table 2 was actually rendered and checked.

Sample-design details not established by the source records remain unresolved, as documented in the sample audit and evidence log. No email correspondence was used as evidence. Original data were not modified.
