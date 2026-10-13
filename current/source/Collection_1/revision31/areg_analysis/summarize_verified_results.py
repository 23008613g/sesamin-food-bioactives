"""Summarize the freshly recomputed cell table; never substitutes for raw rerun."""
from pathlib import Path
import hashlib
import json
import pandas as pd

HERE = Path(__file__).resolve().parent
cells = pd.read_csv(HERE / 'organoid_areg_cells.csv.gz')
filters = {
    'Deposited broad Cholangiocytes': cells.source_broad_label.eq('Cholangiocytes'),
    'Cholangiocytes + Ductal cells excluding SMC': cells.source_subtype_label.isin(['Cholangiocytes', 'Ductal cells']),
    'Hepatocytes': cells.source_broad_label.eq('Hepatocytes'),
}
rows = []
for name, mask in filters.items():
    selected = cells[mask]
    rows.append({'population': name, 'n_cells': len(selected),
                 'n_raw_nonzero': int(selected.AREG_raw_UMI.gt(0).sum()),
                 'raw_detection_pct': 100 * selected.AREG_raw_UMI.gt(0).mean(),
                 'mean_log1p_CP10k': selected.AREG_log1p_CP10k.mean()})
pd.DataFrame(rows).to_csv(HERE / 'areg_broad_exclusion_pooled.csv', index=False, float_format='%.17g')
lib = pd.read_csv(HERE / 'areg_population_by_library.csv')
paired = []
for accession, df in lib.groupby('library_accession'):
    hep = df[df.population.eq('Hepatocytes')].iloc[0]
    for name in list(filters)[:2]:
        target = df[df.population.eq(name)].iloc[0]
        paired.append({'library_accession': accession, 'population': name,
                       'n_cells': target.n_cells, 'n_hepatocyte_cells': hep.n_cells,
                       'raw_detection_pct': target.raw_detection_pct,
                       'hepatocyte_raw_detection_pct': hep.raw_detection_pct,
                       'delta_detection_percentage_points': target.raw_detection_pct - hep.raw_detection_pct,
                       'mean_log1p_CP10k': target.mean_log1p_CP10k,
                       'hepatocyte_mean_log1p_CP10k': hep.mean_log1p_CP10k,
                       'delta_mean_log1p_CP10k': target.mean_log1p_CP10k - hep.mean_log1p_CP10k})
pd.DataFrame(paired).to_csv(HERE / 'areg_broad_vs_hepatocyte_by_library.csv', index=False, float_format='%.17g')
freeze = json.loads((HERE / 'input_hashes_before_run.json').read_text('utf-8'))
original = json.loads((HERE.parents[1] / 'organoid_expression/decompression_log.json').read_text('utf-8'))
source = next(f for f in freeze['files'] if f['path'].endswith('GSE207889_os.h5ad'))
assert source['sha256'] == original['sha256']
audit = {'original_decompression_log_sha256_matches_current_input': True,
         'object_sha256': source['sha256'],
         'original_download_deposited_md5': original['deposited_compressed_md5'],
         'original_download_actual_md5': original['actual_compressed_md5'],
         'qualification': 'Compressed file MD5 copied from the original download/decompression log; not freshly re-downloaded.',
         'all_8_broad_vs_hepatocyte_detection_and_mean_differences_positive': all(r['delta_detection_percentage_points'] > 0 and r['delta_mean_log1p_CP10k'] > 0 for r in paired)}
(HERE / 'additional_source_validation.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
outputs = []
for p in sorted(HERE.iterdir()):
    if p.is_file() and p.name != 'output_hashes.json':
        outputs.append({'name': p.name, 'bytes': p.stat().st_size,
                        'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()})
(HERE / 'output_hashes.json').write_text(json.dumps(outputs, indent=2), encoding='utf-8')
print(json.dumps({'pooled': rows, 'additional_checks': audit}, indent=2))
