"""Authorized project-local extraction and actual bundled kallisto test."""
from common08 import *
import csv, math, shutil, tarfile

def main():
    plan=json.loads((PREV/'processing_plan.json').read_text(encoding='utf-8'))
    archive=PREV/'sources/kallisto_windows-v0.51.1.tar.gz'
    assert digest(archive)==plan['tool']['archive_sha256']
    assert archive.stat().st_size==plan['tool']['archive_bytes']
    assert shutil.disk_usage('D:/').free >= plan['proposed_free_space_gate_bytes']
    write_json(HERE/'authorization.json', {'recorded_utc':now(),'user_statement':'是否允許我在D槽專案目錄配置kallisto 0.51.1便攜版，完成試跑，技術檢查通過後繼續全部9份文庫的定量？（允许）另外这次任务可以对壓縮包，进行配置或執行','authorized_scope':'Project-local tool configuration/execution; pinned references; pilot then all9 after technical QC; local data processing','system_changes':False,'new_external_fees':False})
    target=DATA/'tools/kallisto-0.51.1'
    already_extracted=target.exists()
    target.mkdir(parents=True,exist_ok=True)
    with tarfile.open(archive,'r:gz') as tf:
        members=tf.getmembers()
        for m in members:
            destination=(target/m.name).resolve()
            if not destination.is_relative_to(target.resolve()) or not (m.isdir() or m.isfile()):
                raise RuntimeError(f'Unsafe member: {m.name}')
        if already_extracted:
            for m in members:
                if m.isfile():
                    installed=target/m.name
                    assert installed.is_file() and installed.stat().st_size==m.size
                    assert hashlib.sha256(tf.extractfile(m).read()).hexdigest()==digest(installed)
        else:
            tf.extractall(target, filter='data')
    assert EXE.is_file()
    if not (HERE/'execution/tool_version.json').exists():
        run_logged('tool_version', [EXE,'version'])
    # This Windows binary rejects --help. Bare subcommand emits documented usage.
    for sub in ('quant','index'):
        run_once('tool_'+sub+'_usage', [EXE,sub], allowed_codes=(0,1))
        usage=(DATA/('logs/tool_'+sub+'_usage.stdout.log')).read_text()+(DATA/('logs/tool_'+sub+'_usage.stderr.log')).read_text()
        assert 'Usage: kallisto '+sub in usage
        if sub=='quant':
            assert '--plaintext' in usage and '--threads' in usage
    version=(DATA/'logs/tool_version.stdout.log').read_text()+(DATA/'logs/tool_version.stderr.log').read_text()
    assert '0.51.1' in version
    toy=DATA/'software_test'
    toy.mkdir(exist_ok=True)
    test=EXE.parent/'test'
    run_once('toy_index',[EXE,'index','-i',toy/'test.idx','-k','31',test/'transcripts.fasta.gz'])
    run_once('toy_quant',[EXE,'quant','-i',toy/'test.idx','-o',toy/'quant','-t','8','--plaintext',test/'reads_1.fastq.gz',test/'reads_2.fastq.gz'])
    info=read_kallisto_info(toy/'quant/run_info.json')
    assert info['n_processed']>0 and info['n_pseudoaligned']>0
    with (toy/'quant/abundance.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert rows and all(math.isfinite(float(r[k])) and float(r[k])>=0 for r in rows for k in ('est_counts','tpm'))
    report={'passed':True,'completed_utc':now(),'archive_sha256':digest(archive),'executable_sha256':digest(EXE),'executable':str(EXE),'version_output':version.strip(),'extracted':True,'executable_run':True,'help_compatibility_note':'Initial quant --help rejected with exit1; original failure retained. Bare quant/index emit usage; planned options checked from actual binary.','toy_test_scope':'Bundled synthetic/test software data, not this research study','toy_run_info':info,'toy_transcript_rows':len(rows),'global_environment_changes':False}
    write_json(HERE/'tool_validation.json',report)
    print(json.dumps(report,ensure_ascii=False),flush=True)

if __name__=='__main__':
    main()
