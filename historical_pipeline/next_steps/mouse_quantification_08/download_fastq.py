"""Download frozen study inputs; aria2 validates ENA's MD5; retain sources."""
from common08 import *
import argparse, shutil

def download_runs(runs, key):
    manifest=json.loads((PREV/'fastq_manifest.json').read_text(encoding='utf-8'))
    entries=[e for e in manifest['entries'] if e['run'] in runs]
    assert len(entries)==2*len(runs)
    assert shutil.disk_usage('D:/').free>100_000_000_000
    directory=DATA/'fastq'
    directory.mkdir(parents=True,exist_ok=True)
    inputfile=HERE/(key+'.aria2')
    inputfile.write_text('\n'.join(f"{e['url']}\n  out={e['filename']}\n  checksum=md5={e['expected_md5']}" for e in entries)+'\n',encoding='utf-8')
    aria=shutil.which('aria2c')
    assert aria
    connections='8' if len(runs)>1 else '4'
    run_logged(key,[aria,'--input-file='+str(inputfile),'--dir='+str(directory),'--max-concurrent-downloads=2','--max-connection-per-server='+connections,'--split='+connections,'--min-split-size=20M','--continue=true','--check-integrity=true','--auto-file-renaming=false','--allow-overwrite=false','--file-allocation=none','--summary-interval=30','--console-log-level=warn','--download-result=full','--max-tries=8','--retry-wait=10','--connect-timeout=30','--timeout=60'])
    results=[]
    for e in entries:
        path=directory/e['filename']
        measured={'file':str(path),'url':e['url'],'expected_bytes':e['bytes'],'actual_bytes':path.stat().st_size,'expected_md5':e['expected_md5'],'actual_md5':digest(path,'md5'),'sha256':digest(path)}
        measured['passed']=measured['expected_bytes']==measured['actual_bytes'] and measured['actual_md5']==measured['expected_md5']
        results.append(measured)
    write_json(HERE/(key+'_integrity.json'),{'verified_utc':now(),'passed':all(r['passed'] for r in results),'files':results})
    assert all(r['passed'] for r in results)
    print(json.dumps({'downloaded_runs':runs,'verified_files':len(results)},ensure_ascii=False),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--remaining',action='store_true')
    args=parser.parse_args()
    manifest=json.loads((PREV/'fastq_manifest.json').read_text(encoding='utf-8'))
    allruns=sorted({e['run'] for e in manifest['entries']})
    if args.remaining:
        gate=json.loads((HERE/'pilot_gate.json').read_text(encoding='utf-8'))
        assert gate['passed']
        download_runs(allruns[1:],'remaining_download')
    else:
        download_runs(allruns[:1],'pilot_download')
