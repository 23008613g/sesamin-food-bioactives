"""Build the actual pinned reference index after reference/tool validation."""
from common08 import *

def main():
    ref=json.loads((HERE/'reference_build/reference_validation.json').read_text(encoding='utf-8'))
    tool=json.loads((HERE/'tool_validation.json').read_text(encoding='utf-8'))
    assert ref['passed'] and tool['passed'] and ref['release']==116
    fasta=Path(ref['combined_fasta'])
    assert fasta.is_file()
    index=DATA/'reference/mouse_ensembl116.idx'
    if index.exists():
        raise RuntimeError('Index already exists; inspect before rerun')
    rec=run_logged('mouse_index',[EXE,'index','-i',index,'-k','31',fasta])
    assert index.stat().st_size>0
    write_json(HERE/'index_validation.json',{'passed':True,'completed_utc':now(),'index':str(index),'index_bytes':index.stat().st_size,'index_sha256':digest(index),'fasta':str(fasta),'fasta_sha256':digest(fasta),'reference_validation_sha256_at_build':digest(HERE/'reference_build/reference_validation.json'),'version':tool['version_output'],'k':31,'execution':rec})

if __name__=='__main__':
    main()
