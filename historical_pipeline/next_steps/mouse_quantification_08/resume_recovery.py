"""Resume only the isolated replacement transfer; all failed originals retained."""
from common08 import *
import shutil

name='DRR1010539_2.fastq.gz'
stage=HERE/'download_recovery'
destination=DATA/'download_recovery'
row=next(e for e in json.loads((PREV/'fastq_manifest.json').read_text(encoding='utf-8'))['entries'] if e['filename']==name)
prior=json.loads((HERE/'execution/recover_DRR1010539_R2.json').read_text(encoding='utf-8'))
assert prior['status']=='failed'
assert (stage/'failed_transfer.json').exists() and (destination/name).exists()
write_json(stage/'resume_reason.json',{'recorded_utc':now(),'reason':'Isolated single-connection replacement averaged about0.5MiB/s, so root intentionally stopped only its verified aria2 process and resumed the same isolated transfer using4connections. This is not evidence of the original checksum failure cause.','original_failed_raw_file_untouched':True,'parameters_changed':'Transfer connections1→4 only; same official URL, MD5, bytes and analysis plan','first_helper_edit_failure':'A syntax typo was caught before execution, corrected before any recovery file actions.'})
run_logged('recover_DRR1010539_R2_resume4',[shutil.which('aria2c'),'--dir='+str(destination),'--out='+name,'--max-connection-per-server=4','--split=4','--min-split-size=20M','--continue=true','--check-integrity=true','--checksum=md5='+row['expected_md5'],'--auto-file-renaming=false','--allow-overwrite=false','--file-allocation=none','--summary-interval=30','--console-log-level=warn','--max-tries=8','--retry-wait=10',row['url']])
path=destination/name
actual=digest(path,'md5')
assert path.stat().st_size==row['bytes'] and actual==row['expected_md5']
write_json(stage/'replacement_integrity.json',{'passed':True,'verified_utc':now(),'new_file':str(path),'actual_bytes':path.stat().st_size,'actual_md5':actual,'sha256':digest(path),'official_md5':row['expected_md5'],'execution_key':'recover_DRR1010539_R2_resume4'})
print(json.dumps({'replacement_integrity_passed':True,'file':name}))
