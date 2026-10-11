"""Bounded public metadata audit; no reads download, expression analysis, or installs.

Default: reparse archived official metadata. --fetch: fetch only six missing
NGDC BioSample pages with a 250 KB per-page cap and no challenge handling.
"""
from pathlib import Path
from datetime import datetime, timezone
from html.parser import HTMLParser
import argparse
import hashlib
import json
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
PRIOR = ROOT.parent.parent / 'evidence_closure_06' / 'intervention_access'
EXTRA_IDS = ['SAMC3214403', 'SAMC3214404', 'SAMC3214406', 'SAMC3214407', 'SAMC3214409', 'SAMC3214410']
OLD_IDS = ['SAMC3214402', 'SAMC3214405', 'SAMC3214408']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(name, obj):
    (ROOT / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


class Rows(HTMLParser):
    def __init__(self):
        super().__init__()
        self.stack = []
        self.rows = []
    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.stack.append({'cells': [], 'cell': None})
        elif tag in ('td', 'th') and self.stack:
            self.stack[-1]['cell'] = []
    def handle_data(self, data):
        if self.stack and self.stack[-1]['cell'] is not None:
            self.stack[-1]['cell'].append(data)
    def handle_endtag(self, tag):
        if tag in ('td', 'th') and self.stack:
            row = self.stack[-1]
            if row['cell'] is not None:
                row['cells'].append(' '.join(' '.join(row['cell']).split()))
            row['cell'] = None
        elif tag == 'tr' and self.stack:
            self.rows.append(self.stack.pop()['cells'])


def parse_html(path):
    parser = Rows()
    parser.feed(path.read_text(encoding='utf-8'))
    return {r[0]: r[1] for r in parser.rows if len(r) == 2}


def evidence(source, locator):
    return {'source_id': source, 'locator': locator}


def field(value, status, scope, sources, note=None, raw_value=None):
    result = {'value': value, 'status': status, 'scope': scope, 'evidence': sources}
    if note:
        result['note'] = note
    if raw_value is not None:
        result['raw_value'] = raw_value
    return result


def main():
    arg = argparse.ArgumentParser()
    arg.add_argument('--fetch', action='store_true')
    args = arg.parse_args()
    (ROOT / 'sources').mkdir(parents=True, exist_ok=True)
    attempts_path = ROOT / 'fetch_log.json'
    attempts = json.loads(attempts_path.read_text(encoding='utf-8')) if attempts_path.exists() else []
    if args.fetch:
        for sid in EXTRA_IDS:
            path = ROOT / 'sources' / f'{sid}.html'
            if path.exists():
                continue
            url = f'https://ngdc.cncb.ac.cn/biosample/browse/{sid}'
            record = {'id': sid, 'url': url, 'attempted_utc': datetime.now(timezone.utc).isoformat()}
            try:
                request = urllib.request.Request(url, headers={'User-Agent': 'Research metadata verification/1.0'})
                with urllib.request.urlopen(request, timeout=35) as response:
                    data = response.read(250_001)
                    if len(data) > 250_000:
                        raise ValueError('Stopped at 250 KB metadata limit')
                    text = data.decode('utf-8')
                    if 'BioSample Detail' not in text or sid not in text:
                        raise ValueError('Response is not the expected public BioSample page; no challenge handled')
                    path.write_bytes(data)
                    record.update(status='downloaded', http_status=response.status,
                                  final_url=response.url, bytes=len(data), sha256=digest(path))
            except Exception as exc:
                record.update(status='failed', error=f'{type(exc).__name__}: {exc}')
            attempts.append(record)
            write('fetch_log.json', attempts)
            print(json.dumps(record, ensure_ascii=True), flush=True)

    prior_records = json.loads((PRIOR / 'source_records.json').read_text(encoding='utf-8'))
    archived = []
    for filename, sid in [('sra_sesamin_records.xml', 'SRA_XML'), ('PRJCA021826.html', 'NGDC_PROJECT')]:
        path = PRIOR / 'sources' / filename
        record = next(r for r in prior_records if r.get('filename') == filename and r.get('status') == 'downloaded')
        archived.append({'id': sid, 'url': record['url'], 'local_path': str(path.relative_to(ROOT.parent.parent)),
                         'path_base': 'next_steps', 'bytes': path.stat().st_size, 'sha256': digest(path),
                         'mode': 'read_existing_archive', 'original_fetch_utc': record['attempted_utc']})
    manifest_path = PRIOR / 'raw_read_manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    tree = ET.parse(PRIOR / 'sources' / 'sra_sesamin_records.xml')
    packages = {p.find('EXPERIMENT').get('accession'): p for p in tree.findall('.//EXPERIMENT_PACKAGE')}
    records = []
    checks = []
    def check(name, passed):
        checks.append({'check': name, 'passed': bool(passed)})
    abstracts = set()
    for row in manifest:
        sid = row['source_sample']
        html_path = PRIOR / 'sources' / f'{sid}.html' if sid in OLD_IDS else ROOT / 'sources' / f'{sid}.html'
        if not html_path.exists():
            raise FileNotFoundError(f'Need archived page {sid}; run --fetch first')
        attrs_html = parse_html(html_path)
        archived.append({'id': sid, 'url': f'https://ngdc.cncb.ac.cn/biosample/browse/{sid}',
                         'local_path': str(html_path.relative_to(ROOT.parent.parent)), 'path_base': 'next_steps',
                         'bytes': html_path.stat().st_size, 'sha256': digest(html_path),
                         'mode': 'read_existing_archive' if sid in OLD_IDS else 'new_download'})
        p = packages[row['mirror_experiment']]
        sample = p.find('SAMPLE')
        sample_path = f".//EXPERIMENT_PACKAGE[SAMPLE/@accession='{row['mirror_sample']}']/SAMPLE"
        exp_path = f".//EXPERIMENT[@accession='{row['mirror_experiment']}']"
        attrs = {a.findtext('TAG'): a.findtext('VALUE') for a in sample.findall('SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')}
        sample_evidence = evidence('SRA_XML', sample_path + '/SAMPLE_ATTRIBUTES')
        condition_evidence = [evidence(sid, "table/tr[th='样品名称']/td"),
                              evidence('SRA_XML', f".//STUDY[@accession='DRP018055']/DESCRIPTOR/STUDY_ABSTRACT"),
                              evidence(sid, "table[@id='attribute_table']/tbody/tr[th='处理方法']/td (empty)")]
        abstract = p.findtext('STUDY/DESCRIPTOR/STUDY_ABSTRACT')
        abstracts.add(abstract)
        f = {}
        for key in ['source_project', 'source_archive', 'mirror_project', 'mirror_study', 'source_sample', 'mirror_sample',
                    'biosample_mirror', 'source_experiment', 'mirror_experiment', 'source_run', 'mirror_run']:
            f[key] = field(row[key], 'known', 'public_record_identifier', [evidence('SRA_XML', exp_path + ' and enclosing EXPERIMENT_PACKAGE')],
                           'Repository identifiers are not animal or culture donor identifiers.')
        f['sample_label'] = field(row['sample_label'], 'known', 'sample_level', [evidence(sid, "table/tr[th='样品名称']/td"), evidence('SRA_XML', sample_path + '/TITLE')])
        f['label_prefix'] = field(row['group_label'], 'known', 'sample_label_only', [evidence(sid, "table/tr[th='样品名称']/td")], 'A label prefix; not explicit treatment metadata.')
        for key, tag in [('organism','organism'), ('age','age'), ('sex','sex'), ('tissue','tissue'), ('biomaterial_provider','biomaterial_provider')]:
            value = sample.findtext('SAMPLE_NAME/SCIENTIFIC_NAME') if tag == 'organism' else attrs[tag]
            f[key] = field(value, 'known', 'sample_level_reported', [sample_evidence, evidence(sid, 'sample properties table')])
        for key, xpath in [('library_strategy','DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY'),
                           ('library_source','DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE'),
                           ('library_selection','DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SELECTION'),
                           ('library_construction_protocol','DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_CONSTRUCTION_PROTOCOL'),
                           ('platform','PLATFORM/ILLUMINA/INSTRUMENT_MODEL')]:
            f[key] = field(p.find('EXPERIMENT').findtext(xpath), 'known', 'library_level_reported', [evidence('SRA_XML', exp_path + '/' + xpath)])
        f['layout'] = field('PAIRED', 'known', 'library_level_reported', [evidence('SRA_XML', exp_path + '/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED')], 'Paired-end reads do not imply a paired biological design.')
        f['treatment_group'] = field(row['condition_inference'], 'inferred', 'sample_level_mapping_from_study', condition_evidence,
                                     'Frozen label interpretation; do not relabel from expression results. HFD is not proof of an in vivo high-fat diet experiment.')
        f['cell_model'] = field('mouse primary hepatocytes', 'known', 'study_level_only', [evidence('SRA_XML', './/STUDY/DESCRIPTOR/STUDY_ABSTRACT')],
                                'Study description is explicit; sample-level tissue is hepatocyte.')
        f['strain'] = field(None, 'unknown', 'sample_level', [sample_evidence, evidence(sid, "table[@id='attribute_table']/tbody/tr[th='品种']/td")], raw_value=attrs['strain'])
        f['isolate_or_donor_id'] = field(None, 'unknown', 'sample_level', [sample_evidence], 'isolate=missing; BioSample accession is not a donor identity.', attrs['isolate'])
        f['collection_date'] = field(None, 'unknown', 'sample_level', [sample_evidence], 'A missing-data reason is not a collection date.', attrs['collection_date'])
        for key in ['independent_biological_replicate', 'technical_replicate_structure', 'pooled_animals', 'donor_count',
                    'paired_donors_across_conditions', 'culture_batch', 'sequencing_batch',
                    'sesamin_dose', 'palmitate_dose', 'oleate_dose', 'pa_oa_ratio',
                    'sesamin_exposure_duration', 'fatty_acid_exposure_duration', 'treatment_order',
                    'vehicle', 'vehicle_matching', 'synchronized_harvest']:
            f[key] = field(None, 'unknown', 'not_reported_in_reviewed_metadata', [sample_evidence, evidence('SRA_XML', exp_path + '/DESIGN/DESIGN_DESCRIPTION (empty)'), evidence(sid, "table[@id='attribute_table']/tbody/tr[th='处理方法']/td (empty)")])
        records.append({'record_key': sid, 'fields': f, 'official_sample_attributes_verbatim': attrs, 'ngdc_table_fields': attrs_html})
        check(sid + ': label crossmatch', attrs_html.get('样品名称') == row['sample_label'] and sample.findtext('TITLE') == f"{row['sample_label']} ({sid})")
        check(sid + ': accession crossmatch', sample.get('accession') == row['mirror_sample'] and attrs['ngdc_sample_id'] == sid)
        check(sid + ': empty treatment', '处理方法' in attrs_html and attrs_html['处理方法'] == '')
        check(sid + ': demographic crossmatch', attrs_html.get('年龄') == '8 Week(s)' and attrs_html.get('性别') == attrs['sex'] and attrs_html.get('组织器官') == attrs['tissue'])
        check(sid + ': missing donor identity', attrs['isolate'] == 'missing')
        check(sid + ': no treatment attribute in SRA', not any(k.lower() in ('treatment', 'condition', 'dose', 'donor') for k in attrs))
        check(sid + ': design description empty', not p.findtext('EXPERIMENT/DESIGN/DESIGN_DESCRIPTION'))
        check(sid + ': paired reads', p.find('EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED') is not None)
    check('nine distinct samples', len(records) == len({r['record_key'] for r in records}) == 9)
    check('three labels per prefix', all(sum(r['fields']['label_prefix']['value'] == g for r in records) == 3 for g in ['ctrl','HFD','sess']))
    check('consistent explicit study abstract', len(abstracts) == 1 and 'mouse primary hepatocytes' in next(iter(abstracts)))
    check('all treatment mappings remain inferred', all(r['fields']['treatment_group']['status'] == 'inferred' for r in records))
    check('no donor independence inferred', all(r['fields']['independent_biological_replicate']['status'] == 'unknown' for r in records))
    result = {'audit_date_local': '2026-10-08', 'timezone': 'Asia/Taipei', 'status_definitions': {
        'known': 'Explicitly reported in the source at the stated scope; not independently experimentally verified.',
        'inferred': 'Reasonable interpretation from labels plus study description; not explicit sample treatment metadata.',
        'unknown': 'Unavailable or a missing-value placeholder in the limited reviewed sources.'},
        'study': {'source_project': 'PRJCA021826', 'source_archive': 'CRA013797', 'mirror_project': 'PRJDB40955', 'mirror_study': 'DRP018055',
            'design': field(['normal hepatocytes','PA/OA treated hepatocytes','sesamin + PA/OA hepatocytes'], 'known', 'study_level', [evidence('SRA_XML','.//STUDY/DESCRIPTOR/STUDY_ABSTRACT'), evidence('NGDC_PROJECT',"table/tr[th='描述信息']/td")]),
            'public_libraries': 9, 'known_independent_animal_count': None,
            'publication_link': field(None, 'unknown', 'study_identity', [evidence('SRA_XML', './/STUDY/STUDY_LINKS'), evidence('NGDC_PROJECT', 'project detail and linked resources')], 'No direct accession-to-PMID/DOI link in these reviewed project/sample/SRA records. Prior stage06 considered PMID 38679210 a plausible match; no publication methods are transferred into these samples here.')},
        'samples': records,
        'allowed_next_operations': ['label-agnostic download/checksum verification under project storage plan', 'read-level QC', 'label-agnostic expression quantification', 'exploratory plots retaining inferred labels after quantification'],
        'not_supported_now': ['Treat nine libraries as nine independent mice.', 'Claim three biological replicates per group.', 'Use read-pair status as evidence of paired donors.', 'Assert causal sesamin rescue or independent validation from these metadata.', 'Use sample suffix 1/2/3 as matched donor identities across conditions.', 'Fit confirmatory biological significance tests without resolving the experimental units.'],
        'relabeling_rule': 'Keep all nine records and the prior inferred mappings unchanged; expression results cannot be used to choose a favorable group assignment.',
        'raw_inputs': {'manifest_path': str(manifest_path.relative_to(ROOT.parent.parent)), 'manifest_sha256': digest(manifest_path)},
        'execution': {'new_sample_pages': len(EXTRA_IDS), 'raw_sequence_payload_downloaded': False, 'expression_quantification_executed': False, 'installs': 0, 'Claude_calls': 0}}
    write('sample_design.json', result)
    write('source_records.json', {'sources': archived, 'new_download_log': 'fetch_log.json', 'web_tool_attempts': [
        {'url':'https://ngdc.cncb.ac.cn/biosample/browse/SAMC3214403','result':'not accessible via web tool; archived through ordinary HTTP'},
        {'url':'https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA021826','result':'not accessible via web tool; reused prior archived official page'}]})
    write('validation.json', {'passed': all(c['passed'] for c in checks), 'check_count': len(checks), 'checks': checks,
                              'outputs': {n: digest(ROOT/n) for n in ['sample_design.json','source_records.json']}})
    print(json.dumps({'passed': all(c['passed'] for c in checks), 'checks':len(checks), 'samples':len(records), 'new_pages':len(EXTRA_IDS)}, ensure_ascii=True))
    if not all(c['passed'] for c in checks):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
