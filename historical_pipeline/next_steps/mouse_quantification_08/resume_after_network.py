"""Resume the isolated replacement after network loss; keep previous attempts."""
from common08 import *
import shutil

name = 'DRR1010539_2.fastq.gz'
destination = DATA / 'download_recovery'
stage = HERE / 'download_recovery'
row = next(e for e in json.loads((PREV/'fastq_manifest.json').read_text(encoding='utf-8'))['entries'] if e['filename'] == name)
prior = json.loads((HERE/'execution/recover_DRR1010539_R2_resume4.json').read_text(encoding='utf-8'))
assert prior['status'] == 'failed'
attempt = 1
while (HERE/'execution'/f'recover_DRR1010539_R2_network{attempt}.json').exists():
    previous = json.loads((HERE/'execution'/f'recover_DRR1010539_R2_network{attempt}.json').read_text(encoding='utf-8'))
    assert previous['status'] == 'failed', 'Do not duplicate running/passed transfer'
    attempt += 1
key = f'recover_DRR1010539_R2_network{attempt}'
run_logged(key, [shutil.which('aria2c'), '--dir='+str(destination), '--out='+name,
    '--max-connection-per-server=4', '--split=4', '--min-split-size=20M',
    '--continue=true', '--check-integrity=true', '--checksum=md5='+row['expected_md5'],
    '--auto-file-renaming=false', '--allow-overwrite=false', '--file-allocation=none',
    '--summary-interval=30', '--console-log-level=warn', '--max-tries=8', '--retry-wait=10', row['url']])
path = destination/name
actual = digest(path, 'md5')
assert path.stat().st_size == row['bytes'] and actual == row['expected_md5']
assert not Path(str(path)+'.aria2').exists()
write_json(stage/'replacement_integrity.json', {'passed':True, 'verified_utc':now(),
    'new_file':str(path), 'actual_bytes':path.stat().st_size, 'actual_md5':actual,
    'sha256':digest(path), 'official_md5':row['expected_md5'], 'execution_key':key})
print(json.dumps({'replacement_integrity_passed':True, 'file':name}))
