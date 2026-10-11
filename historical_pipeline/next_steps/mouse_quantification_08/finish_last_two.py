"""Continue only the two unfinished runs after independently verified recovery."""
from run_study import *

def finish():
    recovery=json.loads((HERE/'download_recovery/resolution.json').read_text(encoding='utf-8'))
    assert recovery['passed'] and recovery['final_FASTQ_files_verified']==18
    assert json.loads((HERE/'pilot_gate.json').read_text(encoding='utf-8'))['passed']
    for number in range(1010532,1010539):
        result=json.loads((HERE/'run_results'/f'DRR{number}.json').read_text(encoding='utf-8'))
        assert result['passed']
        assert digest(result['abundance_file'])==result['abundance_sha256']
        assert digest(result['qc_file'])==result['qc_sha256']
    runs=['DRR1010539','DRR1010540']
    with ThreadPoolExecutor(max_workers=2) as pool:
        prepared={run:pool.submit(wait_and_prepare,run) for run in runs}
        for run in runs:
            quant_run(run,prepared[run].result())

if __name__=='__main__': finish()
