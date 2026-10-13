"""Build a 77-row audit from source fields; unknowns are not inferred as facts."""
from pathlib import Path
from collections import defaultdict, Counter
import csv, gzip, hashlib, json, re, shutil, xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parent
STAGE=ROOT.parent
NEXT=STAGE.parent
OLD=NEXT/'submission_finalization_25'/'metadata_audit'
SRC=ROOT/'sources'
raw=[]
rows=[]
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def soft(path):
    d=defaultdict(list)
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('!') and ' = ' in line:
            k,v=line[1:].split(' = ',1); d[k].append(v)
    return dict(d)
def value(d,key): return ' | '.join(d.get(key,[])) or 'unknown'
def rel(d,regex):
    s=value(d,'Sample_relation'); m=re.findall(regex,s)
    return ';'.join(m) or 'unknown'
def item(dataset,acc,d,path):
    raw.append(dict(dataset=dataset,accession=acc,path=str(path),source_sha256=sha(path),fields=d))
    return dict(dataset=dataset,accession=acc,original_sample_name=value(d,'Sample_title'),run_accession='unknown',experiment_accession=rel(d,r'SRX\d+'),biosample=rel(d,r'SAM[NDC]\d+'),mirror_accessions='not_applicable',biological_material=value(d,'Sample_source_name_ch1'),source_label='unknown',analysis_label='unknown',exposure='unknown',concentration='unknown',duration='unknown',sequence_and_washout='unknown',control='unknown',harvest='unknown',experimental_unit='unknown',independence_status='unknown',pooling='unknown',sesamin_intervention='unknown',evidence_status='partly_confirmed',analysis_decision='unknown',unresolved_and_impact='unknown',source_url=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}',source_snapshot=path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path),source_date='2026-10-11',evidence_fields='unknown')

# Human discovery libraries: biological labels are in explicit source and description fields.
for n in range(7922593,7922602):
    acc=f'GSM{n}'; p=SRC/f'{acc}.soft.txt'; d=soft(p); row=item('GSE248903',acc,d,p)
    row.update(biological_material='HepG2 human hepatocellular carcinoma cell line',source_label=value(d,'Sample_source_name_ch1'),analysis_label=value(d,'Sample_description'),exposure='Control / OA-model / post-OA sesamin treatment as deposited',concentration='OA deposited as 400 M; sesamin as 16 M; missing concentration prefixes, not reconstructed',duration='OA 24 h; Treat another 36 h; initial adhesion 24 h',sequence_and_washout='OA precedes sesamin for Treat; OA continuation and washout unknown',control='Control source label; matched vehicle and elapsed harvest time unknown',harvest='matched collection across Control/Model/Treat unknown',experimental_unit='one RNA-seq library; independent culture identifier unknown',independence_status='three libraries per label; biological independence unknown',sesamin_intervention='series and sample protocol explicit; missing implementation details',analysis_decision='retain exploratory library-level candidate selection; no biological-replicate inferential claim',unresolved_and_impact='dose prefixes; elapsed collection time; OA continuation; vehicle; independent culture IDs unresolved. Cannot isolate sesamin from possible timing/protocol differences.',evidence_fields='Sample_source_name_ch1; Sample_description; Sample_treatment_protocol_ch1; Sample_relation')
    row['exposure']={'control':'Control-labelled culture; absence of OA/sesamin per group design, vehicle unknown','model':'OA-model-labelled culture','sesamin treated':'OA-model followed by sesamin as deposited'}[d['Sample_source_name_ch1'][0].lower()]
    rows.append(row)

# Mouse: project-level design exists, but no sample-to-treatment protocol is deposited.
runs=json.loads((SRC/'PRJDB40955_runs.json').read_text())
for run in sorted(runs,key=lambda r:r['run_accession']):
    samd=run['sample_accession']; p=SRC/f'{samd}.xml'; x=ET.parse(p)
    attrs={z.findtext('TAG'):z.findtext('VALUE') for z in x.findall('.//SAMPLE_ATTRIBUTE')}
    samc=attrs['ngdc_sample_id']; title=x.findtext('.//SAMPLE/TITLE')
    label=title.split(' ')[0]; group=re.sub(r'\d+$','',label)
    d={'Sample_title':[title],'Sample_source_name_ch1':[attrs.get('tissue','unknown')],'Sample_relation':[],'ENA_attributes':attrs}
    row=item('CRA013797 / PRJCA021826',samc,d,p)
    row.update(run_accession=run['run_accession'],experiment_accession=run['experiment_accession'],biosample=samc,mirror_accessions=f"{samd};{run['secondary_sample_accession']};DRP018055;PRJDB40955",biological_material='mouse primary hepatocytes at project level; male source animals aged 8 weeks',source_label=label,analysis_label=group,exposure='project: normal hepatocytes, PA/OA hepatocytes, Ses+PA/OA hepatocytes; sample label mapping unconfirmed',concentration='unknown',duration='unknown',sequence_and_washout='unknown',control='ctrl label only; vehicle and protocol unknown',harvest='unknown',experimental_unit='one RNA-seq library; animal/isolation/culture unit unknown',independence_status='three libraries per source label; independent animals or isolations not established',sesamin_intervention='project-level yes; sample-specific assignment not independently established',evidence_status='source labels and mirrors confirmed; biological mapping unresolved',analysis_decision='label-level supplementary exploration only; do not support a confirmed mouse sesamin-response claim',unresolved_and_impact='HFD/sess are names, not sample protocols. Do not interpret HFD as dietary treatment or treat mouse dataset as confirmed intervention validation. Dose, timing, pooling and biological units unknown.',source_url=f'https://ngdc.cncb.ac.cn/biosample/browse/{samc}; https://www.ebi.ac.uk/ena/browser/view/{samd}',evidence_fields='NGDC sample name; ENA TITLE, ngdc_sample_id, ngdc_project_id, tissue, sex, age; project DESCRIPTION')
    row['mirror_accessions']+=';'+re.search(r'CRX\d+',run['experiment_title']).group(0)
    rows.append(row)

# Rat: extract all 20 sample headers without reprocessing expression data.
rat_path=NEXT/'evidence_closure_06'/'intervention_access'/'sources'/'GSE42693_series_matrix.txt.gz'
rat_fields=defaultdict(list)
with gzip.open(rat_path,'rt',encoding='utf-8') as h:
    for cells in csv.reader(h,delimiter='\t'):
        if cells and cells[0]=='!series_matrix_table_begin': break
        if cells and cells[0].startswith('!Sample_'): rat_fields[cells[0][1:]].append(cells[1:])
rat_extract=[]
for i,acc in enumerate(rat_fields['Sample_geo_accession'][0]):
    d={k:[v[i] for v in vs] for k,vs in rat_fields.items()}; rat_extract.append(d)
rat_out=SRC/'GSE42693_cached_sample_headers.json'
rat_out.write_text(json.dumps(rat_extract,ensure_ascii=False,indent=2),encoding='utf-8')
for d in rat_extract:
    acc=d['Sample_geo_accession'][0]; name=d['Sample_title'][0]; group=name.removeprefix('liver_').split('_rep')[0]
    row=item('GSE42693',acc,d,rat_out)
    main=group in ('control diet','sesamin')
    row.update(run_accession='not_applicable_microarray',experiment_accession='not_applicable_microarray',biosample='unknown_not_in_source_header',source_label=group,analysis_label={'control diet':'Control','sesamin':'Sesamin','episesamin':'Episesamin','sesamolin':'Sesamolin'}[group],biological_material='whole liver; male Sprague-Dawley rat',exposure='lignan-free diet' if group=='control diet' else f'dietary {group}',concentration='0%' if group=='control diet' else '0.2% dietary lignan',duration='15 days',sequence_and_washout='dietary exposure throughout experiment; no washout described',control='lignan-free diet',harvest='end of 15-day feeding period',experimental_unit='individual animal; one array per selected rat',independence_status='5 animals per group selected from 7; highest and lowest terminal body weights excluded',pooling='individual-liver arrays; no pooling reported',sesamin_intervention='yes' if group=='sesamin' else 'control comparator' if group=='control diet' else 'no; different lignan',evidence_status='group/dose/duration/animal unit confirmed',analysis_decision='retain dietary contrast with selection caveat' if main else 'normalization sensitivity context only; never count as sesamin replicates',unresolved_and_impact='microarray animals lack a one-to-one link to original group-level phenotype data; 5/7 selection limits generalization',source_date='sample headers cached from prior audit; series and original article freshly checked 2026-10-11',evidence_fields='Sample_title; Sample_source_name_ch1; Sample_treatment_protocol_ch1; Series_overall_design; original article p32')
    rows.append(row)

# PHH: reuse all 30 preserved official headers; freshly verify the series sample list.
phh_cached=OLD/'phh_all_thirty_metadata.json'
phh=json.loads(phh_cached.read_text(encoding='utf-8'))
phh_out=SRC/'GSE247407_cached_sample_headers.json'; shutil.copyfile(phh_cached,phh_out)
fresh_ids=set(soft(SRC/'GSE247407.soft.txt')['Series_sample_id'])
assert fresh_ids=={d['accession'] for d in phh}
for d in phh:
    acc=d['accession']; row=item('GSE247407',acc,d,phh_out)
    ch={v.split(':',1)[0].strip():v.split(':',1)[1].strip() for v in d['Sample_characteristics_ch1']}
    individual=ch.get('individual','unknown'); treatment=ch.get('treatment','unknown')
    row.update(biological_material='primary human hepatocytes in 3D collagen sandwich culture',source_label=f'{individual}; {treatment}',analysis_label=value(d,'Sample_description').split(' | ')[0].replace('Column in matrix: ',''),exposure=treatment,concentration='unknown_from_accessed_sample_fields',duration='7 days',sequence_and_washout='unknown_from_accessed_sample_fields',control='without FFA',harvest='7-day culture',experimental_unit=f'library from donor label {individual}; within-label culture replication unknown',independence_status='five depositor patient labels; F95/F95P relationship and within-label replicates unresolved',sesamin_intervention='no',analysis_decision='context-only descriptive comparison; not sesamin-response validation',unresolved_and_impact='F95 and F95P cannot automatically be treated as independent patients; no sesamin treatment',source_date='sample fields from official snapshot 2026-10-10; complete series sample list checked 2026-10-11',evidence_fields='Sample_characteristics_ch1; Sample_description; Sample_relation; fresh Series_overall_design')
    rows.append(row)

# Four selected organoid libraries: paper defines scRNA-seq replicates from one differentiation.
for n in range(6322577,6322581):
    acc=f'GSM{n}'; p=SRC/f'{acc}.soft.txt'; d=soft(p); row=item('GSE207889',acc,d,p)
    pa=n>=6322579; label='PA' if pa else 'CTRL-PA'
    row.update(biological_material='W01/WiCell hPSC-derived orbital-shaker liver organoids',source_label=value(d,'Sample_title'),analysis_label=f'{label}_{1+(n-6322577)%2}',exposure='palmitic acid' if pa else 'PA-matched carrier control',concentration='PA 500 micromolar' if pa else 'no PA',duration='GEO 5 days; corrected paper 4 days (day 21 to day 25); conflict retained',sequence_and_washout='corrected paper: medium change day 21; harvest day 25',control='carrier solution without active fatty acid; corrected method',harvest='day 25',experimental_unit='technical single-cell-suspension/culture-well library from same differentiation experiment; independent differentiation not established',independence_status='two technical libraries per condition; four libraries do not establish four independent biological experiments',pooling='unknown at individual-organoid level',sesamin_intervention='no',analysis_decision='cellular context and descriptive library summaries only; no cell-level inferential replication',unresolved_and_impact='duration differs between GEO and corrected article; do not upgrade technical libraries or cells to independent differentiations',evidence_fields='Sample_title; characteristics; treatment/growth protocol; Hess Statistical analysis and 2024 correction')
    rows.append(row)

for n in range(3178782,3178787):
    acc=f'GSM{n}'; p=SRC/f'{acc}.soft.txt'; d=soft(p); row=item('GSE115469',acc,d,p)
    donor='P'+re.search(r'Patient (\d+)',d['Sample_title'][0]).group(1)
    row.update(biological_material='adult human liver, dissociated caudate lobe',source_label=value(d,'Sample_title'),analysis_label=donor,exposure='no sesamin intervention',concentration='not_applicable',duration='not_applicable',sequence_and_washout='not_applicable',control='cellular reference; no intervention comparator',harvest='after tissue dissociation',experimental_unit='individual donor',independence_status='five source patients/donors per series',pooling='separate donor-labelled libraries',sesamin_intervention='no',evidence_status='donor identity labels and reference role confirmed',analysis_decision='all-donor descriptive cellular context; no sesamin-response validation',unresolved_and_impact='author matrix described as log2CPM after scran normalization; nonzero fraction is not equated with raw-UMI detection',evidence_fields='Sample_title; extract_protocol; data_processing; Series_overall_design')
    rows.append(row)

with (ROOT/'sample_metadata_audit.csv').open('w',newline='',encoding='utf-8-sig') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
(ROOT/'sample_metadata_original_fields.json').write_text(json.dumps(raw,ensure_ascii=False,indent=2),encoding='utf-8')
assert len(rows)==77
assert len({(r['dataset'],r['accession']) for r in rows})==77
# Preserve HTTP success separately from actual document usability.
responses=[]
for p in sorted(ROOT.glob('source_manifest*.json')):
    if p.name=='source_manifest_validated.json': continue
    for m in json.loads(p.read_text(encoding='utf-8')):
        m['manifest_file']=p.name
        if 'path' in m:
            q=ROOT/m['path']; b=q.read_bytes(); m['hash_verified']=sha(q)==m['sha256']
            t=b[:5000].decode('utf-8',errors='replace')
            if '正在进行安全检查' in t: m['content_classification']='security_challenge_not_metadata'
            elif q.suffix in ('.png','.jpg','.pdf') and (b.lstrip().startswith(b'<!DOCTYPE html') or b.lstrip().startswith(b'<html')): m['content_classification']='html_error_not_requested_media'
            elif m['id']=='PRJDB40955_DDBJ': m['content_classification']='client_application_shell_not_bioproject_fields'
            else: m['content_classification']='source_document'
        responses.append(m)
(ROOT/'source_manifest_validated.json').write_text(json.dumps(responses,ensure_ascii=False,indent=2),encoding='utf-8')
validation=dict(rows=len(rows),counts_by_dataset=dict(Counter(r['dataset'] for r in rows)),duplicate_sample_rows=0,fresh_human_headers_identical_to_stage25=all(sha(SRC/f'GSM{i}.soft.txt')==sha(OLD/'sources'/f'GSM{i}.soft.txt') for i in range(7922593,7922602)),fresh_phh_series_matches_all_30_cached_samples=True,all_downloaded_hashes_match=all(m.get('hash_verified',True) for m in responses),response_classification_counts=dict(Counter(m.get('content_classification',m['status']) for m in responses)),cached_input_hashes={str(p):sha(p) for p in [rat_path,phh_cached]},interpretation='metadata verification, not new biological validation')
(ROOT/'validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(validation,ensure_ascii=False,indent=2))
