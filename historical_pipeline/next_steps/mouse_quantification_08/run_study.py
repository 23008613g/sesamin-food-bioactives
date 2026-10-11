"""Full-file QC then fixed-parameter paired-end quantification, one run at a time."""
from common08 import *
import argparse, csv, math, sys
from concurrent.futures import ThreadPoolExecutor

def wait_and_verify_download(run):
    """Permit download of later pairs to overlap the one-at-a-time quant job."""
    manifest=json.loads((PREV/'fastq_manifest.json').read_text(encoding='utf-8'))
    entries=[e for e in manifest['entries'] if e['run']==run]
    assert len(entries)==2
    deadline=time.monotonic()+4*3600
    while True:
        ready=all((DATA/'fastq'/e['filename']).is_file() and (DATA/'fastq'/e['filename']).stat().st_size==e['bytes'] and not (DATA/'fastq'/(e['filename']+'.aria2')).exists() for e in entries)
        if ready:
            break
        download=HERE/'execution/remaining_download.json'
        if download.exists():
            state=json.loads(download.read_text(encoding='utf-8'))
            assert state['status']!='failed', 'Remaining download failed; inspect original log'
        if time.monotonic()>deadline:
            raise TimeoutError('Four-hour bounded wait expired for '+run)
        time.sleep(10)
    results=[]
    for e in entries:
        path=DATA/'fastq'/e['filename']
        actual_md5=digest(path,'md5')
        assert actual_md5==e['expected_md5'], 'ENA checksum mismatch'
        results.append({'file':str(path),'url':e['url'],'expected_bytes':e['bytes'],'actual_bytes':path.stat().st_size,'expected_md5':e['expected_md5'],'actual_md5':actual_md5,'sha256':digest(path),'passed':True})
    write_json(HERE/'per_run_download_integrity'/(run+'.json'),{'passed':True,'verified_utc':now(),'files':results,'download_complete_basis':'Exact size and no aria2 control file, followed by independently verified full-file MD5'})

def prepare_run(run):
    index=json.loads((HERE/'index_validation.json').read_text(encoding='utf-8'))
    assert index['passed']
    manifests=[HERE/'pilot_download_integrity.json',HERE/'remaining_download_integrity.json',HERE/'per_run_download_integrity'/(run+'.json')]
    checked={f['file']:f for p in manifests if p.exists() for f in json.loads(p.read_text(encoding='utf-8'))['files'] if f['passed']}
    mates=[DATA/'fastq'/f'{run}_{mate}.fastq.gz' for mate in (1,2)]
    for m in mates:
        assert str(m) in checked and m.stat().st_size==checked[str(m)]['actual_bytes']
        assert digest(m,'md5')==checked[str(m)]['expected_md5'], 'FASTQ no longer matches ENA MD5'
    input_hashes={str(m):digest(m) for m in mates}
    for m in mates:
        assert input_hashes[str(m)]==checked[str(m)]['sha256']
    input_hashes[str(HERE/'fastq_qc/fastq_qc.py')]=digest(HERE/'fastq_qc/fastq_qc.py')
    qcfile=HERE/'qc'/(run+'.json')
    qcrec=run_once('qc_'+run,[sys.executable,'-X','utf8',HERE/'fastq_qc/fastq_qc.py',*mates,'--output',qcfile],input_fingerprints=input_hashes,output_paths=[qcfile])
    qc=json.loads(qcfile.read_text(encoding='utf-8'))
    assert qc['passed'] and qc['full_scan_complete']
    return index,mates,input_hashes,qcfile,qcrec,qc

def quant_run(run,prepared=None):
    index,mates,input_hashes,qcfile,qcrec,qc=prepared or prepare_run(run)
    out=DATA/'quant'/run
    out.parent.mkdir(parents=True,exist_ok=True)
    quant_hashes=dict(input_hashes)
    quant_hashes[str(EXE)]=digest(EXE)
    assert quant_hashes[str(EXE)]==json.loads((HERE/'tool_validation.json').read_text(encoding='utf-8'))['executable_sha256']
    quant_hashes[index['index']]=digest(index['index'])
    assert quant_hashes[index['index']]==index['index_sha256'], 'Index changed'
    quant_hashes[str(qcfile)]=digest(qcfile)
    key='quant_'+run
    prior=HERE/'execution'/(key+'.json')
    if prior.exists() and json.loads(prior.read_text(encoding='utf-8'))['status']=='failed':
        key+='_retry1' # Original failed log retained; any second failure still stops.
    rec=run_once(key,[EXE,'quant','-i',index['index'],'-o',out,'-t','8','--plaintext',*mates],input_fingerprints=quant_hashes,output_paths=[out/'abundance.tsv',out/'run_info.json'])
    info=read_kallisto_info(out/'run_info.json')
    assert info['n_processed']==qc['total_pairs']>0
    assert 0<info['n_pseudoaligned']<=info['n_processed']
    assert 0<=info['n_unique']<=info['n_pseudoaligned']
    assert info['kallisto_version']=='0.51.1' and info['n_bootstraps']==0
    with (out/'abundance.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==info['n_targets'] and len({r['target_id'] for r in rows})==len(rows)
    assert all(math.isfinite(float(r[k])) and float(r[k])>=0 for r in rows for k in ('est_counts','tpm','length','eff_length'))
    total_tpm=math.fsum(float(r['tpm']) for r in rows)
    total_counts=math.fsum(float(r['est_counts']) for r in rows)
    assert abs(total_tpm-1e6)<10
    assert abs(total_counts-info['n_pseudoaligned'])<max(10,info['n_pseudoaligned']*0.001)
    rate=info['n_pseudoaligned']/info['n_processed']
    flags=[]
    # Engineering diagnostic threshold frozen before any study expression result.
    if rate<0.5:
        flags.append('pseudoalignment_below_50_percent_requires_diagnostic_review')
    for name,stat in qc['mates'].items():
        sample=stat['sample_first_pairs_only']
        if sample['q30_base_fraction_assuming_33']<0.8:
            flags.append(name+'_first100k_q30_below_80pct')
        if sample['fraction_containing_motif']>0.05:
            flags.append(name+'_first100k_exact_adapter_motif_over_5pct')
    result={'run':run,'completed_utc':now(),'passed':not flags,'numerical_and_structure_checks_passed':True,'flags_requiring_review':flags,'run_info':info,'pseudoalignment_fraction':rate,'qc_file':str(qcfile),'qc_sha256':digest(qcfile),'abundance_file':str(out/'abundance.tsv'),'abundance_sha256':digest(out/'abundance.tsv'),'raw_run_info_sha256':digest(out/'run_info.json'),'total_tpm':total_tpm,'total_estimated_counts':total_counts,'transcripts':len(rows),'index_sha256':index['index_sha256'],'quant_execution_key':key,'quant_wall_seconds':rec['wall_seconds'],'quant_peak_working_set_bytes_observed':rec['peak_working_set_bytes_observed'],'full_qc_wall_seconds':qcrec['wall_seconds'],'sample_design_status':'Treatment mapping inferred; biological independence unknown; no validation/causal claim'}
    write_json(HERE/'run_results'/(run+'.json'),result)
    print(json.dumps({'run':run,'pseudoalignment_fraction':rate,'pairs':qc['total_pairs'],'flags':flags,'transcripts':len(rows)}),flush=True)
    if flags:
        raise RuntimeError(f'Diagnostic review required: {flags}; no automatic deletion or relabelling')
    return result

def wait_and_prepare(run):
    wait_and_verify_download(run)
    return prepare_run(run)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--remaining',action='store_true')
    args=p.parse_args()
    manifest=json.loads((PREV/'fastq_manifest.json').read_text(encoding='utf-8'))
    runs=sorted({e['run'] for e in manifest['entries']})
    if args.remaining:
        assert json.loads((HERE/'pilot_gate.json').read_text(encoding='utf-8'))['passed']
        # At most2 independent QC jobs overlap downloads; quant remains sequential8threads.
        with ThreadPoolExecutor(max_workers=2) as pool:
            pending={run:pool.submit(wait_and_prepare,run) for run in runs[1:]}
            for run in runs[1:]:
                quant_run(run,pending[run].result())
    else:
        result=quant_run(runs[0])
        write_json(HERE/'pilot_gate.json',{'passed':result['passed'],'pilot_run':runs[0],'checked_utc':now(),'run_result_sha256':digest(HERE/'run_results'/(runs[0]+'.json')),'basis':'Full MD5/gzip/paired-read QC, actual processed counts, abundance integrity and technical diagnostic flags; no expression direction criterion','thresholds_predeclared':'pseudoalignment<50%, sampledQ30<80%, sampledexactadapter>5% require review, not sample deletion; first100k quality estimates are nonrandom indicators','permission':'User authorized continuation through all9 after technical gate','biological_validation':False})

if __name__=='__main__':
    main()
