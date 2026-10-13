from fetch_metadata31 import fetch, ROOT
from concurrent.futures import ThreadPoolExecutor
import json

targets=[
('rat_original_pdf','https://www.jstage.jst.go.jp/article/jnsv/55/1/55_1_31/_pdf','Ide_2009_jnsv_55_31.pdf'),
('Hess_corrected_figure','https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs44318-024-00308-w/MediaObjects/44318_2024_308_Figb_HTML.png','Hess_corrected_Fig4.png'),
('Hess_original_figure','https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs44318-024-00308-w/MediaObjects/44318_2024_308_Figc_HTML.png','Hess_original_Fig4.png'),
('PRJDB40955_runs_retry','https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB40955&result=read_run&fields=run_accession,experiment_accession,sample_accession,secondary_sample_accession,sample_title,experiment_title,library_name,library_source,library_strategy&format=json','PRJDB40955_runs.json'),
('Bai_2024_publisher','https://www.sciencedirect.com/science/article/pii/S0006295224002235','Bai_2024_publisher.html'),
]
for i in range(987140,987149):
    a=f'DRX{i}'
    targets.append((a,f'https://www.ebi.ac.uk/ena/browser/api/xml/{a}',f'{a}.xml'))
for i in range(3178782,3178787):
    a=f'GSM{i}'
    targets.append((a,f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={a}&targ=self&form=text&view=full',f'{a}.soft.txt'))
with ThreadPoolExecutor(max_workers=3) as pool:
    manifest=list(pool.map(fetch,targets))
(ROOT/'source_manifest_followup.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
