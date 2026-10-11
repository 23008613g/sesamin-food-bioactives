"""Install checksum-verified replacements and preserve failed transfers in D project."""
from common08 import *

def main():
    stage = HERE/'download_recovery'
    resolution_path = stage/'resolution.json'
    assert not resolution_path.exists(), 'Recovery already finalized; inspect its record'
    original = json.loads((HERE/'execution/remaining_download.json').read_text(encoding='utf-8'))
    assert original['status'] == 'failed'
    manifest = json.loads((PREV/'fastq_manifest.json').read_text(encoding='utf-8'))
    rows = {e['filename']: e for e in manifest['entries']}
    replacements = []
    for name, recordpath in [
        ('DRR1010539_2.fastq.gz', stage/'replacement_integrity.json'),
        ('DRR1010540_2.fastq.gz', stage/'second_file/replacement_integrity.json')]:
        rec = json.loads(recordpath.read_text(encoding='utf-8'))
        new = Path(rec['new_file'])
        old = DATA/'fastq'/name
        quarantine = DATA/'failed_downloads'/(name.removesuffix('.fastq.gz')+'_attempt1')
        for path in (new, old, quarantine):
            assert path.resolve().is_relative_to(DATA.resolve())
        assert new.resolve() == (DATA/'download_recovery'/name).resolve()
        assert rec['passed'] and new.is_file() and not Path(str(new)+'.aria2').exists()
        assert new.stat().st_size == rec['actual_bytes'] == rows[name]['bytes']
        assert digest(new, 'md5') == rec['actual_md5'] == rec['official_md5'] == rows[name]['expected_md5']
        assert digest(new) == rec['sha256']
        assert old.is_file() and not quarantine.exists()
        bad_md5 = digest(old, 'md5')
        assert bad_md5 != rows[name]['expected_md5'], 'Do not quarantine a valid original'
        replacements.append({'name':name, 'source':str(new), 'final_path':str(old),
            'quarantine':str(quarantine), 'replacement_record':str(recordpath),
            'replacement_sha256':rec['sha256'], 'failed_original_md5':bad_md5,
            'failed_original_sha256':digest(old), 'failed_original_bytes':old.stat().st_size})
    journal = {'started_utc':now(), 'status':'validated_before_moves', 'replacements':replacements, 'moves':[]}
    write_json(stage/'finalization_journal.json', journal)
    for item in replacements:
        old, new, quarantine = map(Path, (item['final_path'],item['source'],item['quarantine']))
        quarantine.mkdir(parents=True, exist_ok=False)
        for source in (old, Path(str(old)+'.aria2')):
            if not source.exists(): continue
            target = quarantine/source.name
            assert source.resolve().is_relative_to(DATA.resolve()) and target.resolve().is_relative_to(DATA.resolve())
            assert not target.exists()
            source.replace(target)
            journal['moves'].append({'from':str(source),'to':str(target)})
            write_json(stage/'finalization_journal.json', journal)
        assert not old.exists() and new.resolve().is_relative_to(DATA.resolve()) and old.resolve().is_relative_to(DATA.resolve())
        new.replace(old)
        journal['moves'].append({'from':str(new),'to':str(old)})
        write_json(stage/'finalization_journal.json', journal)
    results=[]
    for name, e in rows.items():
        path=DATA/'fastq'/name
        actual_md5=digest(path,'md5')
        assert path.stat().st_size == e['bytes'] and actual_md5 == e['expected_md5'] and not Path(str(path)+'.aria2').exists()
        results.append({'file':str(path),'url':e['url'],'expected_bytes':e['bytes'],
            'actual_bytes':path.stat().st_size,'expected_md5':e['expected_md5'],
            'actual_md5':actual_md5,'sha256':digest(path),'passed':True})
    assert len(results)==18
    combined={'passed':True,'verified_utc':now(),'files':results,
        'recovery_applied_to':[r['name'] for r in replacements],
        'original_download_command_status':'failed',
        'recovery_records':[r['replacement_record'] for r in replacements],
        'all_final_inputs_independently_rehashed':True}
    write_json(stage/'all18_final_integrity.json',combined)
    write_json(HERE/'remaining_download_integrity.json',dict(combined,files=[r for r in results if not Path(r['file']).name.startswith('DRR1010532_')]))
    journal.update(status='passed',completed_utc=now(),passed=True,failed_bytes_deleted=False,
        final_FASTQ_files_verified=18,scientific_parameters_changed=False,
        cause='Original checksum failures unresolved; user reported subsequent network loss, and interrupted replacement was resumed. No causal attribution made.')
    write_json(resolution_path,journal)
    write_json(stage/'finalization_journal.json',journal)
    adjustments=json.loads((HERE/'execution_adjustments.json').read_text(encoding='utf-8'))
    adjustments['recorded_changes'].append({'phase':'remaining transfer recovery',
        'failure':'Two R2 transfers failed official ENA MD5; isolated replacement transfer was also interrupted during reported network loss.',
        'repair':'Re-download in isolated paths using four HTTPS connections; resume interrupted transfer; confirm official metadata and independent MD5/SHA256; preserve failed originals; rehash all18 final files.',
        'resolution':'download_recovery/resolution.json','scientific_parameters_changed':False,
        'failed_inputs_used_for_QC_or_quant':False,'root_cause_established':False})
    write_json(HERE/'execution_adjustments.json',adjustments)
    print(json.dumps({'passed':True,'final_FASTQ_files_verified':18,'failed_transfers_preserved':2}))

if __name__=='__main__': main()
