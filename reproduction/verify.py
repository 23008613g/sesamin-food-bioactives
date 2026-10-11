"""Numerically compare recomputed tables to stage 15 references, retaining missingness."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args()
    results = args.results.resolve()
    records = []
    total_numeric = 0
    for expected in sorted((ROOT/'expected').glob('*/*')):
        module = expected.parent.name
        actual = results/module/('outputs' if module=='rat' else '')/expected.name
        assert actual.exists(), str(actual)
        a = pd.read_csv(actual)
        b = pd.read_csv(expected)
        assert list(a.columns)==list(b.columns), (expected.name,'column names or order')
        assert a.shape==b.shape, (expected.name,'shape',a.shape,b.shape)
        columns = []
        for col in a.columns:
            x, y = a[col], b[col]
            assert np.array_equal(x.isna(), y.isna()), (expected.name,col,'missingness')
            numeric = is_numeric_dtype(x) and is_numeric_dtype(y)
            if numeric:
                xx, yy = x.to_numpy(dtype=float), y.to_numpy(dtype=float)
                valid = np.isfinite(xx)&np.isfinite(yy)
                assert np.array_equal(np.isfinite(xx),np.isfinite(yy)), (expected.name,col,'finite values')
                xx, yy = xx[valid], yy[valid]
                # The arbitrary sign of each rat PCA component is equivalent.
                sign_flip = expected.name=='all20_array_pca.csv' and col in ['PC1','PC2'] and not np.allclose(xx,yy,rtol=1e-9,atol=1e-10) and np.allclose(xx,-yy,rtol=1e-9,atol=1e-10)
                if sign_flip:
                    xx = -xx
                integer_column = np.equal(xx,np.floor(xx)).all() and np.equal(yy,np.floor(yy)).all()
                if integer_column:
                    assert np.array_equal(xx,yy), (expected.name,col,'integer/boolean mismatch')
                else:
                    assert np.allclose(xx,yy,rtol=1e-9,atol=1e-10), (expected.name,col,float(np.max(abs(xx-yy))))
                total_numeric += len(xx)
                columns.append({'column':col,'finite_values':len(xx),'NA_values':int(x.isna().sum()),
                                'max_absolute_error':float(np.max(abs(xx-yy))) if len(xx) else 0.,
                                'integer_exact':bool(integer_column),'PCA_sign_aligned':bool(sign_flip)})
            else:
                assert np.array_equal(x.fillna('<NA>').astype(str),y.fillna('<NA>').astype(str)), (expected.name,col,'labels')
        records.append({'file':f'{module}/{expected.name}','rows':len(a),'columns':len(a.columns),
                        'status':'PASS','numeric_comparisons':columns})
    cross=pd.read_csv(results/'robustness/cross_species_configuration_consensus.csv')
    means=cross[cross.aggregation=='mean']
    rat=pd.read_csv(results/'rat/outputs/candidate_results_all25_two_normalizations.csv')
    exact=pd.read_csv(results/'rat/outputs/posthoc_exact_sensitivity_all25.csv')
    organoid=pd.read_csv(results/'organoid/weighting_sensitivity_all25.csv')
    covered=organoid[organoid.covered]
    key_results={
        'candidate_count_human_mouse':15,'configurations_per_species':64,'pairs':4096,
        'both_concordant_baseline':means.loc[means.both_baseline_same_nonzero,'human_symbol'].tolist(),
        'both_concordant_all4096':means.loc[means.both_model_and_treatment_concordant_n.eq(4096),'human_symbol'].tolist(),
        'rat_welch_fdr05_by_method':{m:d.loc[d.q_bh_fixed25.lt(.05),'human_symbol'].tolist() for m,d in rat.groupby('method')},
        'rat_exact_fdr05_by_method':{m:d.loc[d.permutation_q_bh25.lt(.05),'human_symbol'].tolist() for m,d in exact.groupby('method')},
        'rat_exact_allocations':sorted(exact.allocations.unique().tolist()),
        'organoid_weighting_sign_agreement':int(covered.total_direction_agrees.sum()),
        'organoid_weighting_denominator':len(covered),
        'organoid_weighting_disagreements':covered.loc[~covered.total_direction_agrees.astype(bool),'symbol'].tolist(),
    }
    external=json.loads((results/'external/analysis_summary.json').read_text(encoding='utf-8'))
    external_checks=json.loads((results/'external/calculation_checks.json').read_text(encoding='utf-8'))
    focused=json.loads((results/'external/focused_self_audit.json').read_text(encoding='utf-8'))
    assert external_checks['all_passed'] and all(external_checks['checks'].values())
    assert focused['all_passed'] and all(focused['checks'].values())
    paired=pd.read_csv(results/'external/paired_descriptive_summary_all25.csv')
    key_results['external']={'donors':external['n_donors'],'cells':external['n_cells'],
        'matched_candidates':external['n_candidates_matched'],'bounded_checks':len(external_checks['checks']),
        'focused_checks':len(focused['checks']),
        'AREG_descriptive_rows':json.loads(paired.loc[paired.symbol.eq('AREG')].to_json(orient='records'))}
    report={'status':'PASS','tolerance':{'float_rtol':1e-9,'float_atol':1e-10,'integer_boolean':'exact',
             'missingness_labels_row_and_column_order':'exact','rat_PCA_component_global_sign':'may align'},
            'reference_tables':len(records),'reference_rows':sum(r['rows'] for r in records),
            'numeric_values_compared':total_numeric,'key_results':key_results,'tables':records,
            'verification_type':'Recomputed values compared cell by cell; not merely copied-file SHA verification.'}
    (results/'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='tables'},indent=2))


if __name__ == '__main__':
    main()
