# Current manuscript snapshot

This folder supersedes the 2026-10-11 repository overview for the current manuscript.
`source/` is a byte-identical browseable selection from the current supplementary
archive: all Python/R/shell scripts, selected run/protocol documentation and small
result tables. It is not a complete runnable directory. Use the complete extracted
archive for execution; omitted inputs are not silently replaced with demo data.

`source_manifest.json` identifies every mirrored member and its SHA-256.
`archive_inventory.csv` lists all physically stored ZIP members and whether each
has a Git mirror. `archive.json` pins the full release attachment's size and hash.
The archive-root maps inside `source/` refer to the complete extracted archive,
not to this intentionally smaller browseable selection.

The original member text is preserved for provenance, including historical
journal names, figure numbers and local execution paths. The top-level
README and docs/REPRODUCING_CURRENT.md govern current usage. A historical
script or log is not a claim that the complete pipeline was rerun for this release.
