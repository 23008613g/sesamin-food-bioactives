"""Verify source tables, captions and final figure QA; record completed visual review."""
from pathlib import Path
import hashlib
import json
import pandas as pd
import numpy as np

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[3]
INPUT=PROJECT/'next_steps/single_cell_extension_14/analysis'
STEMS=['fig7_sc_umap','fig8_sc_features','fig9_sc_dotplot','figS5_sc_library_umap','figS6_sc_quality','figS7_sc_composition']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,obj): (HERE/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

numeric=json.loads((HERE/'numerical_qa.json').read_text(encoding='utf-8'))
composition=pd.read_csv(HERE/'figS7_cell_composition_source.csv')
refs={
    'a':('cell_composition_library_category4.csv','celltype_category'),
    'b':('cell_composition_library_hepatocyte_subtype5.csv','celltype_scTypedb')
}
compchecks=[]
for panel,(filename,category) in refs.items():
    original=pd.read_csv(INPUT/filename)
    displayed=composition[composition.panel==panel]
    assert len(displayed)==len(original)
    for row in displayed.itertuples():
        ref=original[(original.gsm==row.library)&(original[category]==row.category)]
        assert len(ref)==1
        value=ref.iloc[0]
        assert row.cells==value.n_cells and row.denominator_cells==value.denominator_n_cells
        assert np.isclose(row.percent,value.cell_fraction_pct,rtol=1e-9,atol=1e-10)
    compchecks.append({'panel':panel,'rows_checked':len(displayed),'matches_independent_analysis_table':True,
                       'source_file':filename,'source_sha256':sha(INPUT/filename)})
numeric['composition_independent_table_checks']=compchecks
write('numerical_qa.json',numeric)

legends=json.loads((HERE/'figure_legends.json').read_text(encoding='utf-8'))
wordcounts={k:len(v.split()) for k,v in legends.items()}
assert all(n<=200 for n in wordcounts.values())
qa=json.loads((HERE/'qa_results.json').read_text(encoding='utf-8'))
for stem in STEMS:
    c=json.loads((HERE/(stem+'.collision.json')).read_text(encoding='utf-8'))
    assert c['verdict']=='PASS' and c['summary']['fail']==c['summary']['warn']==0
    assert qa['figures'][stem]['text_returncode']==0
    assert qa['figures'][stem]['alignment'] in ['PASS','NOT APPLICABLE']

panels={
 'fig7_sc_umap':{
  'a':'All ten subtype colors, two coordinate axes, all cells shown, legend readable.',
  'b':'Condition overlay uses same extent and equal aspect; overlapping colors do not imply statistical separation.'},
 'fig8_sc_features':{
  'a':'VNN1 sparse 105/24034 (0.44%) detected; sparse marks remain sparse on common range.',
  'b':'AREG 10479/24034 (43.60%) detected; full shared color range preserved.',
  'c':'JAG2 518/24034 (2.16%) detected; no gene-specific color amplification.',
  'd':'CXCL8 2667/24034 (11.10%) detected; same limits, shared colorbar and no clipped observations.'},
 'fig9_sc_dotplot':{
  'single':'All 25 candidate labels and ten subtypes visible; AOC4P crosses differ from measured zero open circles; complete size/color legends.'},
 'figS5_sc_library_umap':{
  'a':'Control 1 includes 3093 cells; label, count, legend and coordinates legible.',
  'b':'Control 2 includes 5510 cells; identical coordinate limits.',
  'c':'PA 1 includes 7399 cells; revised row gutter separates title from upper x-axis label.',
  'd':'PA 2 includes 8032 cells; revised row gutter separates title from upper x-axis label.'},
 'figS6_sc_quality':{
  'a':'All-gene UMI log10 density, full observed range and median/IQR visible in four libraries.',
  'b':'Detected-gene log10 density, full range and median/IQR visible in four libraries.',
  'c':'Raw mitochondrial percent (not fraction), four complete ranges shown; no significance annotations.'},
 'figS7_sc_composition':{
  'a':'Four broad-class cell fractions and all-retained-cell denominators clear; bars sum to 100%.',
  'b':'Five hepatocyte subtype fractions use separate hep-like denominator; sparse HB1 retained; legends clear.'}
}
record={'status':'PASS','date':'2026-10-09','review_type':'Codex internal visual and technical check, not expert biological validation',
    'method':'All six final PNG exports actually opened with view_image; every panel and complete figure inspected. Figure S5 reopened after spacing repair.',
    'panels_total':sum(len(v) for v in panels.values()),'figures':[],
    'caption_word_counts':wordcounts,
    'reviewed_source_warnings':{
      'FINAL-WIDTH':'Static parser misses WIDTH variable; actual PDF width measured as 180 mm for all six.',
      'DEMO-DATA':'np.random.default_rng is only used to permute complete real-cell draw order; no synthetic data generation or row sampling.',
      'UNCERTAINTY-ENCODING':'Seed applies only to draw order; cell summaries are descriptive, not multi-seed/fold estimates. No independent biological replicate interval can be inferred.'},
    'resolved_attempts':[
      'Initial S5 export had two text-text collisions between upper x-axis labels and lower panel titles (21.4% overlap); hspace changed from 0.32 to 0.53 and only S5 rerendered. Final alignment/collision/text gates pass and final S5 was viewed.',
      'An auxiliary stdin command containing a non-ASCII absolute path was corrupted by Windows pipe encoding before changing files. Repeated using a relative working-directory path; captions were successfully saved.'
    ]}
for stem in STEMS:
    record['figures'].append({'stem':stem,'panels':[{'panel':k,'inspection':v,'status':'PASS'} for k,v in panels[stem].items()],
       'png_sha256':sha(HERE/(stem+'.png')),'pdf_sha256':sha(HERE/(stem+'.pdf')),
       'alignment':qa['figures'][stem]['alignment'],'collision':'PASS; 0 FAIL, 0 WARN','minimum_pdf_glyph_pt':7,
       'dimensions':qa['figures'][stem]['geometry'],'all_points_within_axis_limits':True})
write('visual_review.json',record)
manifest=[{'file':str(p.relative_to(HERE)).replace('\\','/'),'bytes':p.stat().st_size,'sha256':sha(p)}
          for p in sorted(HERE.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json' and '__pycache__' not in p.parts]
write('artifact_manifest.json',manifest)
print(json.dumps({'status':'PASS','figures':6,'panels':record['panels_total'],'composition_rows_checked':len(composition),'caption_words':wordcounts}))
