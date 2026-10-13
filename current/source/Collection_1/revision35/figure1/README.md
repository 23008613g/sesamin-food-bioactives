# Figure 1: dataset evidence map

The figure is an original Python vector schematic, not an AI-generated image or quantitative plot. All six dataset roles, accession identifiers, library/cell counts and interpretation limits are retained from the preceding figure and manuscript. Original pictograms are navigation aids and do not depict observed cells or encode sample sizes.

## Rebuild

From this directory, use an environment with matplotlib:

```text
python draw_figure1.py --output-dir ./exports
```

The output directory receives editable PDF/SVG, EPS and 600-dpi PNG/TIFF. QA files are written beside the script. The alignment helper is included in `qa_tools`. The source was executed with Python 3.13.9 and matplotlib 3.10.9.

Use an environment with PyMuPDF and Pillow to audit:

```text
python check_figure1.py --output-dir ./exports
```

The final audit was executed with PyMuPDF 1.27.2.2. It checks source preflight, actual PDF lettering, rendered PDF collisions, all 31 required text assertions, dimensions, editable text and raster resolution. The helpers in `qa_tools` are copied unchanged from the installed figure QA toolkit; no credentials or external services are required.

`figure_contract.md` documents the scientific and visual contract. `figure1_render_record.json` records text geometry and render versions. `figure1_export_checks.json` records final checks and output hashes. The alignment SVG and collision PDF are QA-only diagnostics, not submission artwork.

No scientific analyses or data values were changed. The existing manuscript Figure 1 caption remains applicable. Final placement and page flow must still be checked after the figure is inserted into Word.
