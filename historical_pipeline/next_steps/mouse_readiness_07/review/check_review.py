"""Read-only stage07 evidence checks; writes only review/review.json. No network."""
from pathlib import Path
import csv
import hashlib
import json
import re
import tarfile
import xml.etree.ElementTree as ET
from collections import Counter

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
STEPS = BASE.parent
checks = []
read_paths = set()

def read(path):
    read_paths.add(path)
    return path.read_bytes()

def js(path):
    return json.loads(read(path).decode('utf-8-sig'))

def rows(path, delimiter=','):
    return list(csv.DictReader(read(path).decode('utf-8-sig').splitlines(), delimiter=delimiter))

def check(name, passed, details=None):
    checks.append({'check': name, 'passed': bool(passed), 'details': details})

o = js(BASE / 'orthology/orthology.json')
cs = o['candidates']
old = rows(STEPS / 'sesamin_exploration_04/priority_all.csv')
check('candidate_identity_and_order_unchanged', [(c['original_human_id'], c['input_symbol']) for c in cs] == [(c['gene_id'], c['symbol']) for c in old])
check('candidate_input_hash', hashlib.sha256(read(STEPS / 'sesamin_exploration_04/priority_all.csv')).hexdigest() == o['input_sha256'])
counts = Counter(c['status'] for c in cs)
check('15_plus_1_plus_9_equals_original_25', dict(counts) == {'one2one_high_confidence': 15, 'non_one2one_mapping': 1, 'current_human_id_no_mouse_orthologue_in_ensembl116': 9} and len(cs) == 25)
bm = rows(BASE / 'orthology/raw/biomart_25.txt', '\t')
for c in cs:
    gid = c['original_human_id']
    lookup = js(BASE / f'orthology/raw/lookup_{gid}.json')
    h = js(BASE / f'orthology/raw/homology_{gid}.json')['data']
    hits = h[0]['homologies']
    bm_hits = [r for r in bm if r['Gene stable ID'] == gid and r['Mouse gene stable ID']]
    output_hits = c['mouse_orthologues']
    check(f'{gid}_lookup_exact_human', lookup['id'] == gid and lookup['species'] == 'homo_sapiens')
    check(f'{gid}_REST_target_set', {r['target']['id'] for r in hits} == {r['mouse_id'] for r in output_hits})
    check(f'{gid}_BioMart_target_set', {r['Mouse gene stable ID'] for r in bm_hits} == {r['mouse_id'] for r in output_hits})
    for r in output_hits:
        match = [x for x in bm_hits if x['Mouse gene stable ID'] == r['mouse_id']][0]
        check(f'{gid}_{r["mouse_id"]}_confidence_type', int(match['Mouse orthology confidence [0 low, 1 high]']) == r['official_orthology_confidence'] and match['Mouse homology type'] == r['orthology_type'])
    eligible = len(hits) == 1 and hits[0]['type'] == 'ortholog_one2one' and len(bm_hits) == 1 and bm_hits[0]['Mouse orthology confidence [0 low, 1 high]'] == '1'
    check(f'{gid}_primary_rule', c['primary_comparison_eligible'] == eligible)
for rec in js(BASE / 'orthology/source_records.json'):
    check(f'orthology_source_hash_{rec["key"]}', hashlib.sha256(read(BASE / 'orthology' / rec['file'])).hexdigest() == rec['sha256'])

ena = rows(BASE / 'sources/ena_files.tsv', '\t')
manifest = js(BASE / 'fastq_manifest.json')
headers = js(BASE / 'fastq_head_checks.json')
prior = js(STEPS / 'evidence_closure_06/intervention_access/raw_read_manifest.json')
prior_by_run = {r['mirror_run']: r for r in prior}
entries = manifest['entries']
check('ENA_runs_9_unique', len(ena) == len({r['run_accession'] for r in ena}) == 9)
check('FASTQ_files_18_unique', len(entries) == len({r['filename'] for r in entries}) == 18)
check('ENA_total_bytes', sum(int(n) for r in ena for n in r['fastq_bytes'].split(';')) == manifest['total_bytes'] == 31230644073)
check('HEAD_count18', len(headers) == 18)
prior_xml = STEPS / 'evidence_closure_06/intervention_access/sources/sra_sesamin_records.xml'
xml = ET.fromstring(read(prior_xml))
for r in ena:
    run = r['run_accession']
    old = prior_by_run[run]
    check(f'{run}_mirror_identity', r['study_accession'] == old['mirror_project'] and r['sample_accession'] == old['biosample_mirror'] and r['library_layout'] == 'PAIRED')
    packages = [p for p in xml.findall('.//EXPERIMENT_PACKAGE') if p.find(f'.//RUN[@accession="{run}"]') is not None]
    check(f'{run}_original_XML_identity', len(packages) == 1 and packages[0].find('SAMPLE').get('accession') == old['mirror_sample'] and packages[0].find('EXPERIMENT').get('accession') == old['mirror_experiment'])
    expected = list(zip(r['fastq_ftp'].split(';'), r['fastq_md5'].split(';'), map(int, r['fastq_bytes'].split(';'))))
    ee = [e for e in entries if e['run'] == run]
    check(f'{run}_paired_mates', len(ee) == 2 and sorted(e['mate'] for e in ee) == [1, 2])
    for ftp, md5, nbytes in expected:
        e = [e for e in ee if e['url'] == 'https://' + ftp][0]
        h = [h for h in headers if h['url'] == e['url']][0]
        check(f'{e["filename"]}_metadata_and_HEAD', e['bytes'] == nbytes == int(h['content_length']) and e['expected_md5'] == md5 and h['status'] == 200 and h['method'] == 'HEAD' and h['payload_bytes_downloaded'] == 0)

design = js(BASE / 'design/sample_design.json')
check('design_9_records', len(design['samples']) == 9)
design_sources = js(BASE / 'design/source_records.json')['sources']
for rec in design_sources:
    path = STEPS / rec['local_path']
    check(f'design_source_hash_{rec["id"]}', hashlib.sha256(read(path)).hexdigest() == rec['sha256'])
for sample in design['samples']:
    f = sample['fields']
    check(f'{sample["record_key"]}_treatment_inferred_and_biology_unknown', f['treatment_group']['status'] == 'inferred' and all(f[k]['status'] == 'unknown' for k in ['independent_biological_replicate', 'pooled_animals', 'donor_count', 'paired_donors_across_conditions']))
    rec = next(r for r in design_sources if r['id'] == sample['record_key'])
    html = read(STEPS / rec['local_path']).decode('utf-8-sig')
    treatments = re.findall(r'<th>\s*处理方法\s*</th>\s*<td>(.*?)</td>', html, flags=re.S)
    check(f'{sample["record_key"]}_source_treatment_field_empty', len(treatments) == 1 and not treatments[0].strip())

plan = js(BASE / 'processing_plan.json')
commands = js(BASE / 'planned_commands.json')
tool = js(BASE / 'tool_package_plan.json')
env = js(BASE / 'environment_preflight.json')
follow = js(BASE / 'environment_followup.json')
latest = js(BASE / 'sources/kallisto_latest.json')
release = js(BASE / 'sources/kallisto_v0.51.1.json')
check('pilot_lexicographic_selection', plan['pilot_run'] == min(r['run_accession'] for r in ena))
check('pilot_bytes', sum(e['bytes'] for e in entries if e['run'] == plan['pilot_run']) == plan['pilot_fastq_bytes'] == 3610761909)
check('reference_total_bytes', sum(int(r['bytes']) for r in plan['reference']['files']) == plan['reference_download_bytes'] == 264151884)
check('reference_same_release_assembly_unfetched', all('release-116' in r['url'] and 'GRCm39' in r['url'] and r['payload_bytes_downloaded'] == 0 and r['status'] == 200 and r['content_validated'] is False for r in plan['reference']['files']))
check('pin_not_latest_windows_asset', latest['tag_name'] == 'v0.52.0' and not any('windows' in a['name'].lower() for a in latest['assets']) and any(a['name'] == 'kallisto_windows-v0.51.1.tar.gz' for a in release['assets']))
archive = BASE / 'sources/kallisto_windows-v0.51.1.tar.gz'
check('archive_hash_bytes', len(read(archive)) == tool['archive_bytes'] == 3956678 and hashlib.sha256(read(archive)).hexdigest() == tool['sha256'] == plan['tool']['archive_sha256'])
with tarfile.open(archive, 'r:gz') as tar:
    members = tar.getmembers()
    check('archive_regular_or_dirs_only', all(m.isfile() or m.isdir() for m in members))
    check('archive_exe_path_exists', any(m.name == 'kallisto/kallisto.exe' for m in members))
    check('archive_saved_member_sizes', [(m.name,m.size) for m in members] == [(m['name'],m['bytes']) for m in tool['archive_members']])
check('tool_unexecuted_in_records', tool['extracted'] is False and tool['executable_run'] is False and plan['tool']['runtime_tested'] is False)
check('data_root_still_absent', not Path(plan['data_root_proposed']).exists())
check('commands_all_unexecuted', all(c['executed'] is False for c in commands))
quants = [c for c in commands if 'quant' in c['argv']]
check('quant_all9_paired_plaintext8', len(quants) == 9 and all(c['argv'][c['argv'].index('-t')+1] == '8' and '--plaintext' in c['argv'] and '--single' not in c['argv'] and sum(a.endswith('.fastq.gz') for a in c['argv']) == 2 for c in quants))
disks = {d['DeviceID']: d['FreeSpace'] for d in env['fixed_disks']['data']}
check('C_insufficient_and_D_reserve_available_at_check', disks['C:'] < plan['fastq_bytes'] and disks['D:'] > plan['proposed_free_space_gate_bytes'] == 100000000000)
check('environment_followup_beyond_PATH', env['path_tools']['kallisto'] is None and follow['registered_environments'] == ['<LOCAL_PATH>'] and follow['Rscript_exists'] and follow['sequence_packages_in_base'] == [])
check('15_primary_gene_ID_lists_match', plan['primary_human_gene_ids'] == [c['original_human_id'] for c in cs if c['primary_comparison_eligible']] and plan['primary_mouse_gene_ids'] == [c['mouse_orthologues'][0]['mouse_id'] for c in cs if c['primary_comparison_eligible']])
for name, expected in plan['input_snapshot_sha256'].items():
    check(f'plan_input_hash_{name}', hashlib.sha256(read(BASE/name)).hexdigest() == expected)

for name in ['report.md','processing_plan.md','build_processing_plan.py','assistance_record.json']:
    read(BASE/name)
result = {
    'reviewer': 'Codex internal subagent readiness_review07; not Claude or expert biological validation',
    'user_facing_date': '2026-10-08',
    'skill_loaded': 'validate-data',
    'scope': 'Saved-source consistency and processing-readiness review only; no network, installation, archive extraction, sequence download or quantification',
    'overall_assessment': 'Share with caveats' if all(c['passed'] for c in checks) else 'Needs revision',
    'passed': all(c['passed'] for c in checks),
    'check_count': len(checks),
    'checks': checks,
    'issues_requiring_current_document_revision': [],
    'execution_blockers': ['Project-local portable tool configuration requires explicit authorization under supplied AGENTS policy as interpreted in current report.', 'Windows binary compatibility and runtime dependencies not tested.', 'Reference and FASTQ payloads absent; their content/integrity and actual compute needs remain unknown.'],
    'scientific_caveats': ['Treatment mapping is inferred; donor count, independence, pooling and pairing remain unknown.', '15 orthologue pairs are mapping eligibility, not novel targets or validated intervention effects.', 'CXCL8 has no returned mouse orthologue; functional analog substitution is prohibited; ACAT2 retains both mouse targets as secondary.', 'TPM and estimated counts do not supply biological replication; no inferential test or direct cross-study effect magnitude pooling is supported.'],
    'not_verified': ['Current network availability (saved responses checked only)', 'Payload MD5, gzip and paired-read integrity', 'Actual binary version/help behavior and official toy test', 'No writes anywhere on D drive historically (only proposed project path absence and saved logs checked)', 'Live free disk space and peak runtime resource use', 'Experimental biology beyond archived public records'],
    'reviewed_input_sha256': {str(p.relative_to(STEPS)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(read_paths)},
    'reproduction': 'Run existing project .venv/Scripts/python.exe -X utf8 next_steps/mouse_readiness_07/review/check_review.py; output only review/review.json',
}
(HERE/'review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'passed': result['passed'], 'check_count':len(checks),'failures':[c for c in checks if not c['passed']]}))
