"""Download three small public processed tables and audit structure only.

No effect estimates, differential expression, pathway ranking or AI scoring.
Existing downloads are read as-is and are never overwritten by this script.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import re
import sys
import urllib.request

FILES = [
    ('GSE200765', 'series_matrix', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200765/matrix/GSE200765_series_matrix.txt.gz', 14),
    ('GSE270357', 'TPM', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/GSE270357/suppl/GSE270357_TPMs_allsamples.txt.gz', 9),
    ('GSE247407', 'filename_gene_count', 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE247nnn/GSE247407/suppl/GSE247407_gene_count.txt.gz', 30),
]

SAMPLE_PATTERNS = {
    'GSE200765': r'GSM[0-9]+',
    'GSE270357': r'PA(?:_Cur(?:_CuCl2)?)?[1-3]',
    'GSE247407': r'(?:F95P|F95|M95|H961|H20)[CF][1-3]',
}


def inspect(accession, kind, url, expected, folder):
    path = folder / url.rsplit('/', 1)[-1]
    r = {'accession':accession, 'kind':kind, 'url':url, 'local_file':'downloads/'+path.name,
         'expected_sample_columns':expected, 'status':'not_attempted'}
    try:
        if path.exists():
            data = path.read_bytes()
            r['retrieval'] = 'existing_file_read_without_overwrite'
        else:
            req = urllib.request.Request(url, headers={'User-Agent':'PublicDataFeasibilityAudit/1.0'})
            with urllib.request.urlopen(req, timeout=40) as response:
                data = response.read(25*1024*1024+1)
                r['http_status'] = response.status
                r['last_modified'] = response.headers.get('Last-Modified')
            if len(data) > 25*1024*1024:
                raise ValueError('Compressed download exceeds 25 MiB limit')
            with path.open('xb') as out:
                out.write(data)
            r['retrieval'] = 'downloaded_new_file'
        r['bytes'] = len(data)
        r['sha256'] = hashlib.sha256(data).hexdigest()
        ids, header, rows, bad_width = set(), None, 0, 0
        sample_indices = []
        nonnumeric, missing, nonfinite, negative, fractional, duplicates = 0, 0, 0, 0, 0, 0
        in_table = kind != 'series_matrix'
        total = 0
        with gzip.GzipFile(fileobj=io.BytesIO(data)) as gz:
            with io.TextIOWrapper(gz, encoding='utf-8-sig') as stream:
                for line in stream:
                    total += len(line.encode('utf-8'))
                    if total > 256*1024*1024:
                        raise ValueError('Uncompressed text exceeds 256 MiB limit')
                    if line.startswith('!series_matrix_table_begin'):
                        in_table = True
                        continue
                    if line.startswith('!series_matrix_table_end'):
                        break
                    if not in_table or not line.strip():
                        continue
                    row = next(csv.reader([line], delimiter='\t'))
                    if header is None:
                        header = row
                        sample_indices = [i for i,name in enumerate(header) if re.fullmatch(SAMPLE_PATTERNS[accession], name)]
                        r['sample_columns'] = [header[i] for i in sample_indices]
                        r['annotation_columns'] = [name for i,name in enumerate(header) if i not in sample_indices]
                        continue
                    rows += 1
                    if len(row) != len(header):
                        bad_width += 1
                    if row[0] in ids:
                        duplicates += 1
                    ids.add(row[0])
                    for i in sample_indices:
                        if i >= len(row):
                            missing += 1
                            continue
                        value = row[i]
                        if value.strip().lower() in ('','na','nan','null'):
                            missing += 1
                            continue
                        try:
                            x = float(value)
                        except ValueError:
                            nonnumeric += 1
                            continue
                        if not math.isfinite(x):
                            nonfinite += 1
                        elif x < 0:
                            negative += 1
                        elif x != math.floor(x):
                            fractional += 1
        r.update({'header':header, 'feature_rows':rows, 'duplicate_feature_ids':duplicates,
                  'observed_sample_columns':len(sample_indices),
                  'bad_width_rows':bad_width, 'missing_values':missing,
                  'nonnumeric_values':nonnumeric, 'nonfinite_values':nonfinite,
                  'negative_values':negative, 'fractional_nonnegative_values':fractional})
        good = header and len(sample_indices) == expected and rows > 0 and bad_width == 0 and nonnumeric == 0 and nonfinite == 0
        r['status'] = 'structure_read' if good else 'structure_requires_review'
        r['interpretation'] = 'Structural counts only. No normalization, outlier, donor independence, phenotype, biological effect or replication validation.'
    except Exception as exc:
        r['status'] = 'failed'
        r['error'] = f'{type(exc).__name__}: {exc}'
    return r


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    folder = args.output.resolve().parent / 'downloads'
    folder.mkdir(parents=True, exist_ok=True)
    records = [inspect(*item, folder) for item in FILES]
    result = {'retrieved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'decision_date_taipei':'2026-10-07', 'command':sys.argv,
              'interpreter':sys.executable, 'python_version':platform.python_version(),
              'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope':'Processed-table download and structural/numeric-format inspection only; all exploration',
              'records':records}
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([{k:r.get(k) for k in ['accession','status','bytes','feature_rows','observed_sample_columns','error']} for r in records],ensure_ascii=True,indent=2))


if __name__ == '__main__':
    main()
