Food-derived sesamin and hepatic lipid metabolism: a computational assessment of cross-study transcriptomic evidence
Plant Foods for Human Nutrition
Tingting ZHENG; Honglie ZHANG
Infinia Intelligent Health (Hong Kong) Limited, Hong Kong
Correspondence: Honglie ZHANG; hz@infinia-hk.com

Online Resource 3: Whole-transcriptome analyses and current figure revisions

This resource accompanies the final main manuscript and ESM_1.pdf. ESM_2.zip
contains the earlier fixed-candidate and cellular-context analyses. Earlier
figure files in ESM_2 are archival: use the final manuscript and its separately
numbered figures for the current presentation. No new raw-read quantification,
CEL preprocessing, UMAP fitting or biological experiment was performed here.

Start with README_reproduction.txt. scripts/reproduce_analysis.py accepts an
explicit existing Python/R environment and a NEW output directory. It uses
the included frozen inputs without network downloads. Its scope is the three
whole-genome rat, limma and cross-study analysis functions, with 62 result
tables compared against the references. These are distinct from the earlier
candidate/cellular-context reproduction tables described in Online Resource 2.
The archived successful run was performed on the same Windows computer;
another operating system or dependency version has not been executed here.

Figure numbering:
  figures/Fig_G1_sample_structure = final Fig. S17
  figures/Fig_G2_method_sensitivity = final Fig. 4
  figures/Fig_G3_GO_score_evidence = final Fig. S18
  figures/Fig_G4_cross_study_alignment = final Fig. 5
  revision41/Fig1 = current author-supplied Fig. 1 (raster)
  revision40/figures_new/Fig3 = final Fig. 3
  revision40/figures_new/Fig1 = previous programmatic workflow, archival only
  final Fig. 6 is the retained organoid UMAP; Figs. 7 and 8 retain its
  organoid quantification and adult-liver context from the earlier analyses.

revision40/scripts/update_main_figures.py reuses its verified source_data
copies. Pass --qa-tools ../../figures/qa_tools as an absolute path and
--audit-python with an existing interpreter containing PyMuPDF when running
from a different location. No new gene selection or statistical test is added
by the earlier programmatic figure update. Its Fig1 output reproduces the
previous layout, not the current author-supplied artwork in revision41.
The Fig4/Fig5 export wrapper additionally documents
its --stage39 argument, which should point to the extracted ESM_3 root,
and exports existing figure objects without rerunning the statistical analyses.

Measurement boundaries: HepG2 is tumour-derived; the dietary rat experiment
is not normal-human intervention validation. Organoid and adult-liver analyses
give cellular context without sesamin exposure. Full-family adjusted results
and non-significant results are retained. Source data and reused tools retain
their original provenance and applicable terms; inclusion does not transfer
ownership of the original datasets or publications.

Figure 1 replacement on 13 October 2026
--------------------------------------
revision41 contains the supplied PNG, a lossless TIFF conversion and a raster
PDF. All three use the same original 1448 x 1086 pixels. At the manuscript
width of 6.85 inches the effective density is 211.39 pixels per inch; these
are not upscaled 600-dpi or editable-vector exports. The creation tool and
individual illustration sources were not provided. The manuscript caption
retains tumour-derived HepG2, unconfirmed mouse exposure mapping and the
distinction between an analytical illustration and an intervention study.
All original scientific data and analysis code remain unchanged.


Descriptive signal-intensity and hit-stability diagnostic
--------------------------------------------------------
The English manuscript and ESM_1.pdf include Table S20 and Fig. S19.
revision43/diagnostics/ contains a post-hoc diagnostic of the existing limma
results. It does not refit models, replace the original multiple-testing
families, or add expression filters. It uses the same ten primary rat arrays.
All existing scientific result and code members are retained unchanged.

Run diagnose_limma_hits.py with an existing Python environment containing
NumPy, pandas and SciPy. Pass --stage39 <extracted_ESM_3_root> and
--output-dir <new_directory_outside_the_ESM_3_root>. This rebuilds the
descriptive diagnostic outputs from the frozen results/ tables. Consult
revision43/diagnostics/README_reproduction.txt for the recorded environment and figure
command. Existing full-pipeline reproduction and its 62-table comparison
remain separate from this diagnostic. No fresh raw-data pipeline run was
performed for this document revision.

Current supplemental source paths:
  S16: figures/source_data/G2a_full_genome_BH_counts.csv and results/
  S17: figures/source_data/G2b_preexisting_examples_complete_results.csv
       plus results/limma/preexisting_VNN1_AREG_examples.csv
  S18: figures/source_data/G3_all31_lipid_name_GO_terms_not_selected_for_plot.csv
       plus results/go_results_*.csv
  S20 and S19 figure: revision43/diagnostics/
The candidate archive manifest is ESM_2/PACKAGE_MANIFEST.csv.
