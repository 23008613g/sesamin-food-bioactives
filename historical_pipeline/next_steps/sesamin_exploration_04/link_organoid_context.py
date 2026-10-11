"""Descriptive disease-only reference for every sesamin exploratory candidate."""
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'organoid_expression'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = datetime.now(timezone.utc).isoformat()
    source_audit = json.loads((SOURCE / 'aggregation_audit.json').read_text())
    files = ['all_gene_library_subtype_counts.tsv.gz', 'library_subtype_summary.csv']
    for name in files:
        expected = next(x['sha256'] for x in source_audit['outputs'] if x['file'] == name)
        assert sha(SOURCE / name) == expected
    data = pd.read_csv(SOURCE / files[0], sep='\t', keep_default_na=False).set_index('gene_id')
    meta = pd.read_csv(SOURCE / files[1])
    candidates = pd.read_csv(HERE / 'priority_all.csv', keep_default_na=False)
    assert data.index.is_unique and candidates['gene_id'].is_unique and len(candidates) == 25
    for row in meta.itertuples():
        assert int(data[row.count_column].sum()) == row.total_counts
    records, comparisons = [], []
    for candidate in candidates.itertuples():
        present = candidate.gene_id in data.index
        libs = {}
        for row in meta.itertuples():
            value = int(data.at[candidate.gene_id, row.count_column]) if present else None
            cpm = value / row.total_counts * 1e6 if present and row.total_counts > 0 else None
            records.append({'gene_id': candidate.gene_id, 'bulk_symbol': candidate.symbol,
                            'organoid_symbol': data.at[candidate.gene_id, 'gene_symbol'] if present else None,
                            'gsm': row.gsm, 'condition': row.group, 'subtype': row.subtype,
                            'cells_aggregated': row.cells_aggregated, 'total_UMI': row.total_counts,
                            'candidate_UMI': value, 'CPM': cpm,
                            'status': 'observed_positive' if present and value > 0 else 'observed_zero' if present else 'feature_not_covered'})
            if row.subtype == 'ALL_HEPATOCYTE_LIKE':
                libs[row.gsm] = cpm
        item = {'gene_id': candidate.gene_id, 'bulk_symbol': candidate.symbol,
                'feature_covered': present, 'bulk_OA_log2_effect': candidate.log2_OA_vs_control,
                'bulk_sesamin_log2_effect': candidate.log2_sesamin_vs_OA,
                'same_gene_symbol': bool(data.at[candidate.gene_id, 'gene_symbol'] == candidate.symbol) if present else None}
        item.update({gsm + '_CPM': value for gsm, value in libs.items()})
        if present:
            control = np.mean([libs['GSM6322577'], libs['GSM6322578']])
            pa = np.mean([libs['GSM6322579'], libs['GSM6322580']])
            delta = float(np.log2((pa + 1) / (control + 1)))
            item.update({'CTRL_technical_library_mean_CPM': control, 'PA_technical_library_mean_CPM': pa,
                         'organoid_PA_log2_CPM_ratio': delta,
                         'libraries_with_detected_UMI': sum(v > 0 for v in libs.values()),
                         'disease_direction_same_sign': bool(np.sign(delta) == np.sign(candidate.log2_OA_vs_control))})
        comparisons.append(item)
    pd.DataFrame(records).to_csv(HERE / 'organoid_candidate_context.csv', index=False)
    combined = pd.DataFrame(comparisons)
    combined.to_csv(HERE / 'organoid_context_summary.csv', index=False)
    selected = combined.loc[combined['bulk_symbol'].isin(['CXCL8', 'VNN1'])].to_dict('records')
    audit = {'started_utc': started, 'completed_utc': datetime.now(timezone.utc).isoformat(),
             'interpreter': sys.executable, 'script_sha256': sha(Path(__file__)),
             'context_plan_sha256': sha(HERE / 'context_plan.json'),
             'candidate_sha256': sha(HERE / 'priority_all.csv'),
             'source_hashes': {name: sha(SOURCE / name) for name in files},
             'requested_candidates': len(candidates), 'matched_exact_ids': int(combined['feature_covered'].sum()),
             'unmatched_ids': combined.loc[~combined['feature_covered'], 'gene_id'].tolist(),
             'context_rows': len(records), 'library_subtype_columns': len(meta),
             'spotlight_after_literature_review': selected,
             'limitations': json.loads((HERE / 'context_plan.json').read_text())['limitations'],
             'artifact_hashes': {name: sha(HERE / name) for name in ['organoid_candidate_context.csv', 'organoid_context_summary.csv']}}
    (HERE / 'organoid_context_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({'matched': audit['matched_exact_ids'], 'missing': audit['unmatched_ids'],
                      'context_rows': len(records), 'spotlight': selected}))


if __name__ == '__main__':
    main()
