# UMAP title edit

The sole visible change is Author cell subtypes -> Cell subtypes in main Figure4. The original study annotation remains acknowledged in the manuscript caption. No data, embedding, category, color or geometry changes were made. figures/legacy_support.py now carries the shorter label; its only source change is that string. Its --only3 entry still exports the historic Figure_03 name, which maps to current main Figure4. All source tables and scientific analyses remain unchanged.

For the self-contained current plotting entry, run python revision33/figure/plot_umap.py in the recorded project environment. It uses the included frozen source table and outputs main/Fig4 in four formats. Required directories are included as empty placeholder paths. No dependency installation is performed. Current label/Word provenance and pixel checks are bundled here; revision32 asset hashes describe the preceding label version.
