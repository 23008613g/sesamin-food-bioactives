"""Read-only independent numeric cross-check of R output against main analysis.

Run using the existing project Python environment. Writes audit/example outputs
only below results/limma; does not refit or select biological hypotheses.
"""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats

stage = Path(sys.argv[1])
results = stage / 'results'
out = results / 'limma'
annotation = pd.read_csv(results / 'rat_gene_annotation.csv', dtype={'rat_entrez': str})
assert annotation.rat_entrez.is_unique
checks = {}
examples = []
for variant in ('deposited', 'qn10', 'qn20'):
    main = pd.read_csv(results / f'rat_gene_results_{variant}.csv.gz', dtype={'rat_entrez': str})
    limma = pd.read_csv(out / f'rat_limma_{variant}.csv.gz', dtype={'gene_id': str})
    assert len(main) == len(limma) == 14604
    assert main.rat_entrez.is_unique and limma.gene_id.is_unique
    assert main.rat_entrez.tolist() == limma.gene_id.tolist()
    errors = {
        'logFC_vs_main_difference': float(np.max(np.abs(limma.logFC - main.difference))),
        'mean_control_vs_main': float(np.max(np.abs(limma.mean_control - main.mean_control))),
        'mean_sesamin_vs_main': float(np.max(np.abs(limma.mean_sesamin - main.mean_sesamin))),
        'p_from_scipy_t': float(np.max(np.abs(limma['P.Value'] - 2 * stats.t.sf(np.abs(limma.t), limma['df.total'])))),
        'BH_from_scipy': float(np.max(np.abs(limma['adj.P.Val'] - stats.false_discovery_control(limma['P.Value'], method='bh')))),
        'BY_from_scipy': float(np.max(np.abs(limma['adj.P.Val.BY'] - stats.false_discovery_control(limma['P.Value'], method='by')))),
    }
    halfwidth = stats.t.ppf(.975, limma['df.total']) * limma.SE_moderated
    errors['CI_lower_from_scipy'] = float(np.max(np.abs(limma.CI95_moderated_low - (limma.logFC - halfwidth))))
    errors['CI_upper_from_scipy'] = float(np.max(np.abs(limma.CI95_moderated_high - (limma.logFC + halfwidth))))
    assert all(error < 1e-10 for error in errors.values()), errors
    checks[variant] = {
        'status': 'passed', 'n_genes': len(limma), 'all_gene_ids_in_identical_order': True,
        'max_absolute_errors': errors,
        'limma_BH_lt_005': int((limma['adj.P.Val'] < .05).sum()),
        'limma_BY_lt_005': int((limma['adj.P.Val.BY'] < .05).sum()),
    }
    merged = limma.merge(annotation[['rat_entrez', 'gene_symbol']], left_on='gene_id', right_on='rat_entrez', validate='one_to_one')
    merged = merged.merge(main[['rat_entrez', 'welch_bh', 'welch_by', 'exact_p', 'exact_bh', 'exact_by']], on='rat_entrez', validate='one_to_one')
    selected = merged[merged.gene_symbol.str.upper().isin(['VNN1', 'AREG'])].copy()
    assert len(selected) == 2
    selected['variant'] = variant
    selected['example_role'] = 'Pre-existing manuscript example, not a new discovery or fitting filter'
    examples.append(selected[['variant', 'rat_entrez', 'gene_symbol', 'logFC', 'P.Value', 'adj.P.Val', 'adj.P.Val.BY', 'CI95_moderated_low', 'CI95_moderated_high', 'welch_bh', 'welch_by', 'exact_p', 'exact_bh', 'exact_by', 'example_role']])

report = {
    'status': 'passed',
    'python': sys.executable,
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'purpose': 'Independent output verification; no additional statistical method or model fit',
    'variants': checks,
}
(out / 'independent_crosscheck.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
example_table = pd.concat(examples, ignore_index=True)
example_table.to_csv(out / 'preexisting_VNN1_AREG_examples.csv', index=False)
print(json.dumps(report, indent=2))
print(example_table.to_string(index=False))
