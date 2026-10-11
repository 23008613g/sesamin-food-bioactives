"""Read-only GSE248903 count-table and all-nine-sample metadata audit.

No differential expression, enrichment, causal inference or model calls.
"""
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
import urllib.request

from verify_targets import fetch

HERE=Path(__file__).resolve().parent
OUT=HERE/'sesamin_inputs'
OUT.mkdir(exist_ok=True)
URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248903/suppl/GSE248903_total_count.anno.txt.gz'
FILE=OUT/'GSE248903_total_count.anno.txt.gz'


def main():
    if not FILE.exists():
        with urllib.request.urlopen(URL,timeout=35) as response:
            data=response.read(20_000_001)
        if len(data)>20_000_000:
            raise ValueError('File over explicit 20 MB single-count-file cap')
        FILE.write_bytes(data)
    assert FILE.stat().st_size==15515659
    csv.field_size_limit(10_000_000)
    summaries={}
    genes=set()
    duplicate_ids=0
    annotation_missing=0
    rows=0
    width_errors=0
    all_zero=0
    with gzip.open(FILE,'rt',encoding='utf-8-sig',newline='') as handle:
        reader=csv.reader(handle,delimiter='\t')
        header=next(reader)
        columns=header[1:10]
        expected=[f'{g}_{i}' for g in ['Model','Control','Treat'] for i in [1,2,3]]
        assert columns==expected
        summaries={c:{'nonzero_genes':0,'sum_counts':0,'max_count':0} for c in columns}
        for row in reader:
            rows+=1
            if len(row)!=len(header):
                width_errors+=1
                raise ValueError(f'row width mismatch {rows}')
            duplicate_ids+=row[0] in genes
            genes.add(row[0])
            annotation_missing+=row[10] in ['', '-', 'NA', 'N/A']
            values=[]
            for c,v in zip(columns,row[1:10],strict=True):
                value=float(v)
                if not math.isfinite(value) or value<0 or not value.is_integer():
                    raise ValueError(f'Invalid integer count {rows} {c} {v}')
                value=int(value)
                values.append(value)
                s=summaries[c]
                s['sum_counts']+=value
                s['nonzero_genes']+=value>0
                s['max_count']=max(s['max_count'],value)
            all_zero+=not any(values)
    assert duplicate_ids==0
    metadata=[]
    manifest=[]
    for gsm in range(7922593,7922602):
        r=fetch('GSM'+str(gsm))
        metadata.append(r)
        f=r['fields']
        column=f['!Sample_description'][0]
        assert column in columns
        group={'Model':'OA_model','Control':'control','Treat':'OA_plus_sesamin'}[column.rsplit('_',1)[0]]
        manifest.append({'gsm':r['accession'],'matrix_column':column,'group':group,
                         'geo_title':f['!Sample_title'][0], 'geo_source':f['!Sample_source_name_ch1'][0],
                         'replicate_label':column.rsplit('_',1)[1],
                         'independent_culture_batch':'not_verified', 'count_file':FILE.name})
    assert {r['matrix_column'] for r in manifest}==set(columns)
    manifest.sort(key=lambda x:columns.index(x['matrix_column']))
    with (HERE/'sesamin_sample_manifest.csv').open('w',encoding='utf-8-sig',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    (HERE/'sesamin_sample_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    record={'executed_utc':datetime.now(timezone.utc).isoformat(), 'interpreter':sys.executable,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'input':{'url':URL,'file':str(FILE.relative_to(HERE)),'bytes':FILE.stat().st_size,
                     'sha256':hashlib.sha256(FILE.read_bytes()).hexdigest()},
            'shape':{'gene_rows':rows,'sample_columns':len(columns),'annotation_columns':len(header)-10},
            'duplicate_gene_ids':duplicate_ids,'all_zero_rows':all_zero,'missing_gene_symbol_rows':annotation_missing,
            'integer_nonnegative_finite_counts':True,'row_width_errors':width_errors,
            'sample_summaries':summaries, 'mapped_GSMs':len(manifest),
            'scope':'File/mapping readiness only; no DE, no independent validation, no new biology or AI score',
            'metadata_limits':['GEO treatment text loses unit symbols: literal 400  M OA and 16   M sesamin; do not silently repair units',
                               'OA 24 h followed by sesamin another 36 h stated; model/control corresponding waiting time and continued OA status unresolved',
                               '3 submitted samples per group do not establish independent culture batches',
                               'No sesamin-only RNA group; OA steatosis is not by itself verified lipotoxic injury'],
            'download_note':'Initial 10 MB one-off read was stopped before saving; HEAD confirmed 15515659 bytes, then full file retrieved under a revised 20 MB cap. No partial raw file retained.',
            'actual_model_calls':0}
    (HERE/'sesamin_input_audit.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record,ensure_ascii=True))


if __name__=='__main__':
    main()
