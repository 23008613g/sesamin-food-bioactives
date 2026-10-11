"""Recover one checksum-failed transfer without overwriting or hiding the failed bytes."""
from common08 import *
import argparse, csv, io, shutil, urllib.request

RUN='DRR1010539'
NAME=RUN+'_2.fastq.gz'
RECOVERY=DATA/'download_recovery'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--finalize',action='store_true')
    args=p.parse_args()
    manifest=json.loads((PREV/'fastq_manifest.json').read_text(encoding='utf-8'))
    row=next(e for e in manifest['entries'] if e['filename']==NAME)
    RECOVERY.mkdir(exist_ok=True)
    stage=HERE/'download_recovery'
    stage.mkdir(exist_ok=True)
    new=RECOVERY/NAME
    old=DATA/'fastq'/NAME
    if not args.finalize:
        url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession='+RUN+'&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv'
        req=urllib.request.Request(url,headers={'User-Agent':'sesamin-research-integrity/1.0'})
        with urllib.request.urlopen(req,timeout=60) as resp:
            payload=resp.read()
        (stage/'ena_refreshed.tsv').write_bytes(payload)
        ena=list(csv.DictReader(io.StringIO(payload.decode()),delimiter='\t'))
        assert len(ena)==1 and ena[0]['run_accession']==RUN
        official=dict(zip([u.rsplit('/',1)[-1] for u in ena[0]['fastq_ftp'].split(';')],zip(ena[0]['fastq_md5'].split(';'),map(int,ena[0]['fastq_bytes'].split(';')))))
        assert official[NAME]==(row['expected_md5'],row['bytes']), 'Official metadata changed; do not assume same source'
        failure={'checked_utc':now(),'file':str(old),'actual_bytes':old.stat().st_size,'actual_md5':digest(old,'md5'),'actual_sha256':digest(old),'expected_md5':row['expected_md5'],'expected_bytes':row['bytes'],'fresh_ENA_metadata_url':url,'fresh_ENA_sha256':digest(stage/'ena_refreshed.tsv'),'official_metadata_unchanged':True,'cause':'Transfer checksum failure; exact root cause not established','failed_input_used_for_QC_or_quant':False}
        assert failure['actual_md5']!=failure['expected_md5']
        write_json(stage/'failed_transfer.json',failure)
        # Fresh path and one connection: keep the prior failure untouched while testing recovery.
        assert not new.exists()
        run_logged('recover_DRR1010539_R2',[shutil.which('aria2c'),'--dir='+str(RECOVERY),'--out='+NAME,'--max-connection-per-server=1','--split=1','--max-concurrent-downloads=1','--check-integrity=true','--checksum=md5='+row['expected_md5'],'--auto-file-renaming=false','--allow-overwrite=false','--file-allocation=none','--summary-interval=30','--console-log-level=warn','--max-tries=8','--retry-wait=10',row['url']])
        actual=digest(new,'md5')
        assert new.stat().st_size==row['bytes'] and actual==row['expected_md5']
        write_json(stage/'replacement_integrity.json',{'passed':True,'verified_utc':now(),'new_file':str(new),'actual_bytes':new.stat().st_size,'actual_md5':actual,'sha256':digest(new),'official_md5':row['expected_md5']})
    else:
        original_download=json.loads((HERE/'execution/remaining_download.json').read_text(encoding='utf-8'))
        assert original_download['status'] in ('failed','passed'), 'Do not move files while original aria2 is active'
        replacement=json.loads((stage/'replacement_integrity.json').read_text(encoding='utf-8'))
        assert replacement['passed'] and digest(new)==replacement['sha256']
        failed=DATA/'failed_downloads/DRR1010539_R2_attempt1'
        failed.mkdir(parents=True,exist_ok=False)
        moved=[]
        for source in (old,Path(str(old)+'.aria2')):
            destination=failed/source.name
            assert source.resolve().is_relative_to(DATA.resolve()) and destination.resolve().is_relative_to(DATA.resolve())
            assert not destination.exists()
            if source.exists():
                source.replace(destination)
                moved.append({'from':str(source),'to':str(destination)})
        assert new.resolve().is_relative_to(DATA.resolve()) and old.resolve().is_relative_to(DATA.resolve()) and not old.exists()
        new.replace(old)
        results=[]
        for e in manifest['entries']:
            if e['run']==sorted({x['run'] for x in manifest['entries']})[0]: continue
            path=DATA/'fastq'/e['filename']
            actual_md5=digest(path,'md5')
            assert path.stat().st_size==e['bytes'] and actual_md5==e['expected_md5'] and not Path(str(path)+'.aria2').exists()
            results.append({'file':str(path),'url':e['url'],'expected_bytes':e['bytes'],'actual_bytes':path.stat().st_size,'expected_md5':e['expected_md5'],'actual_md5':actual_md5,'sha256':digest(path),'passed':True})
        write_json(HERE/'remaining_download_integrity.json',{'passed':True,'verified_utc':now(),'files':results,'recovery_applied_to':[NAME],'original_download_command_status':original_download['status'],'recovery_record':'download_recovery/replacement_integrity.json','all_final_inputs_independently_rehashed':True})
        write_json(stage/'resolution.json',{'resolved_utc':now(),'passed':True,'original_failed_transfer_preserved':moved,'replacement_final_path':str(old),'replacement_sha256':digest(old),'original_transfer_record':'execution/remaining_download.json','failed_bytes_deleted':False})
    print(json.dumps({'phase':'finalize' if args.finalize else 'replacement_download','passed':True}))

if __name__=='__main__': main()
