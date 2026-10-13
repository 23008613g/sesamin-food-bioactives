Reproducing the Stage 39 figures
================================

Keep the package directory structure intact:
  scripts/make_genome_figures.py
  results/                              executed full-gene/GO/cross-study results
  figures/qa_tools/                     bundled unmodified audit helpers

The plotting script resolves all inputs and outputs relative to its own parent
directory. It does not require an installed nature-figure skill or any fixed
user directory. The four actual audit helpers are bundled under qa_tools, with
byte-identical source hashes, the source repository/commit, original headers,
and the original repository LICENSE retained in qa_tools/source_manifest.json
and LICENSE.nature-skills.txt.

Dependencies used in this execution
----------------------------------
Plotting interpreter: Python 3.13.9, matplotlib 3.10.9, numpy 2.4.6, pandas 3.0.3.
PDF collision auditor: PyMuPDF 1.27.2.2 in an existing Python 3.13.9 interpreter.
Pillow 12.3.0 was used separately to verify PNG pixels and resolution.
The exact auditor interpreter and runtime versions are recorded by each build
in qa/figure_build.json. Arial was available for the original rendering; the
script permits DejaVu Sans as a fallback. Different fonts or package versions
can alter geometry, so rerun the QA checks and inspect changed output images.

The script does not install dependencies, modify PATH, access credentials,
download data, or call an external service. Use a separately prepared Python
environment containing the stated dependencies. The source/text/alignment
auditor cores otherwise require only the Python standard library.

Run from the package root
-------------------------
If the active Python contains both plotting dependencies and PyMuPDF:

  python scripts/make_genome_figures.py

The default PDF-audit interpreter is sys.executable, namely the same Python
that is running the plotting script. A second existing interpreter containing
PyMuPDF can be selected explicitly. For PowerShell:

  $env:PFHN_FIGURE_AUDIT_PYTHON = 'C:/path/to/existing/audit-environment/python.exe'
  & 'C:/path/to/existing/plot-environment/python.exe' scripts/make_genome_figures.py

For a POSIX shell with existing environments:

  PFHN_FIGURE_AUDIT_PYTHON=/path/to/audit/python /path/to/plot/python scripts/make_genome_figures.py

The setting is an executable path, not a shell command or a package installer.
If the selected interpreter cannot import PyMuPDF, the script stops before
rendering and reports the missing PDF-audit dependency; it does not silently
omit QA or attempt installation.

Outputs and checks
------------------
The run regenerates four 183 mm-wide figures as editable-text PDF/SVG and
300 dpi PNG, figure_legends.txt, figure_contract.json, per-panel source_data,
and machine-readable QA files. It leaves upstream results untouched.

Every run executes source validation, final rendered panel-alignment checking
(1.5 pt tolerance), PDF text-size checking (5 pt floor), and PDF collision
checking. Source warnings are reviewed in qa/visual_review.md. The plotting
script does not claim to perform human or agent visual inspection itself.
Inspect new/changed final PNGs and panels after rerendering.

For the portability-only change, qa/portability_before.json records the already
visually inspected PNG bytes/pixels; qa/portability_verification.json records
the actual equality check after rerunning with the bundled audit scripts.
If all PNGs are byte- and pixel-identical, the preceding documented visual
inspection remains applicable; this is not represented as a new viewing event.
PDF/SVG file hashes may change due to export metadata even when PNGs remain
identical. The final exported files have updated hashes in the delivery audit.
