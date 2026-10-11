"""Read-only public raw-count/author-annotation reanalysis, four GSE207889 libraries.

Run with the existing project .venv. No downloads, installation, model inference,
UMAP recomputation, reclustering, new QC exclusions, or significance testing.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy
from scipy.io import mmread

BASE = Path(__file__).resolve().parent
N = BASE.parent.parent
ORIG = N / 'organoid_expression'
SMALL = N / 'organoid_inputs'
CANDIDATES = N / 'publication_package_13/data/human_candidates_all25.csv'
HEP_SUBTYPES = ['AH', 'FH1', 'FH2', 'HB1', 'HB2']
SPOTLIGHT = ['VNN1', 'AREG', 'JAG2', 'CXCL8']
INPUT_RECORDS = []
OUTPUTS = []


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def inp(path, expected=None):
    path = Path(path)
    digest = sha(path)
    if expected is not None and digest != expected:
        raise ValueError(f'Input hash mismatch: {path.name}')
    stat = path.stat()
    INPUT_RECORDS.append({'path': str(path), 'sha256': digest, 'bytes': stat.st_size,
                          'mtime_ns': stat.st_mtime_ns, 'prior_audit_matches': expected is not None})
    return path


def decode(node):
    if isinstance(node, h5py.Dataset):
        return node.asstr()[...] if h5py.check_string_dtype(node.dtype) is not None else node[...]
    if 'codes' in node and 'categories' in node:
        cats = decode(node['categories'])
        return np.array([cats[int(i)] if i >= 0 else None for i in node['codes'][...]], dtype=object)
    raise ValueError(f'Unsupported encoding {node.name}')


def table(name, frame):
    path = BASE / name
    frame.to_csv(path, index=False, encoding='utf-8', na_rep='NA',
                 compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None,
                 float_format='%.12g')
    OUTPUTS.append({'file': name, 'rows': len(frame), 'columns': list(frame.columns),
                    'bytes': path.stat().st_size, 'sha256': sha(path)})


def feature_table(path):
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        return [line.rstrip('\n').split('\t') for line in f]


def summarize_genes(cells, raw, cp10k, lognorm, candidates, grouping, scope):
    out = []
    grouper = grouping[0] if len(grouping) == 1 else grouping
    for keys, indices in cells.groupby(grouper, sort=False, dropna=False).groups.items():
        keys = keys if isinstance(keys, tuple) else (keys,)
        idx = np.asarray(list(indices), dtype=int)
        assert len(idx) > 0
        for j, cand in candidates.iterrows():
            values = raw[idx, j]
            covered = bool(cand['covered'])
            present = int(np.sum(values > 0)) if covered else None
            out.append({**dict(zip(grouping, keys)), 'scope': scope,
                        'gene_id': cand['gene_id'], 'symbol': cand['symbol'],
                        'raw_feature_symbol': cand['raw_feature_symbol'], 'covered': covered,
                        'n_cells': len(idx), 'n_detected_cells': present,
                        'detected_pct': present / len(idx) * 100 if covered else None,
                        'mean_raw_counts': float(values.mean()) if covered else None,
                        'mean_cp10k': float(cp10k[idx, j].mean()) if covered else None,
                        'mean_log1p_cp10k': float(lognorm[idx, j].mean()) if covered else None,
                        'sum_raw_counts': int(values.sum()) if covered else None,
                        'pooling': 'All retained cells in this stratum, cell-count weighted; not equal-library weighted'})
    return pd.DataFrame(out)


def composition(cells, grouping, denominator_scope):
    rows = []
    levels = sorted(cells[grouping].unique())
    for gsm, sub in cells.groupby('gsm', sort=False):
        n = len(sub)
        for level in levels:
            count = int((sub[grouping] == level).sum())
            rows.append({'gsm': gsm, 'condition': sub.condition.iloc[0], grouping: level,
                         'n_cells': count, 'denominator_n_cells': n,
                         'cell_fraction_pct': count/n*100, 'denominator_scope': denominator_scope})
    return pd.DataFrame(rows)


def main():
    started = time.monotonic()
    BASE.mkdir(parents=True, exist_ok=True)
    join_audit = json.loads(inp(ORIG/'annotation_join_audit.json').read_text(encoding='utf-8'))
    agg_audit = json.loads(inp(ORIG/'aggregation_audit.json').read_text(encoding='utf-8'))
    dec_audit = json.loads(inp(ORIG/'decompression_log.json').read_text(encoding='utf-8'))
    small_audit = json.loads(inp(N/'organoid_metadata_audit.json', join_audit['small_input_audit_sha256']).read_text(encoding='utf-8'))
    manifest = pd.read_csv(inp(N/'organoid_sample_manifest.csv', join_audit['manifest_sha256']), dtype=str)
    cells = pd.read_csv(inp(ORIG/'four_library_annotations.csv.gz', join_audit['annotation_file_sha256']), dtype=str)
    candidates = pd.read_csv(inp(CANDIDATES))[['gene_id', 'symbol']].copy()
    assert len(candidates) == 25 and candidates.gene_id.is_unique
    assert len(cells) == 24034 and not cells.duplicated(['gsm', 'barcode']).any()
    assert cells.original_obs_index.is_unique
    cells.insert(0, 'cell_key', cells.gsm+'|'+cells.barcode)
    h5path = inp(ORIG/'sources/GSE207889_os.h5ad', dec_audit['sha256'])
    embedding = {'source': str(h5path), 'key': 'obsm/X_umap',
                 'scope': 'Four-library subset of author-provided full-atlas UMAP coordinates',
                 'recomputed': False, 'new_clustering': False, 'random_seed': None,
                 'seed_note': 'Not applicable: deposited coordinates reused, no stochastic embedding run'}
    with h5py.File(h5path, 'r') as h5:
        obs = h5['obs']
        ids = decode(obs[obs.attrs['_index']])
        if len(set(ids)) != len(ids):
            raise ValueError('Full-atlas observation IDs are not unique')
        id_to_pos = dict(zip(ids, range(len(ids))))
        positions = np.array([id_to_pos[x] for x in cells.original_obs_index], dtype=int)
        cells['author_obs_row'] = positions
        umap = h5['obsm/X_umap'][...]
        assert umap.shape == (82467, 2) and np.isfinite(umap).all()
        cells['umap_1'] = umap[positions, 0]
        cells['umap_2'] = umap[positions, 1]
        for label in ['condition', 'detailed_condition', 'celltype_category', 'celltype_scTypedb']:
            assert (decode(obs[label])[positions].astype(str) == cells[label].to_numpy()).all(), label
        for field in ['n_genes', 'n_genes_by_counts', 'total_counts', 'pct_mito', 'doublet_score']:
            if field in obs:
                cells['author_'+field] = decode(obs[field])[positions]
        embedding.update({'full_atlas_n_cells': len(ids), 'selected_n_cells': len(cells),
                          'all_atlas_coordinates_finite': bool(np.isfinite(umap).all()),
                          'selected_axis_min': umap[positions].min(axis=0).tolist(),
                          'selected_axis_max': umap[positions].max(axis=0).tolist(),
                          'saved_umap_parameters': {}, 'author_colors': {}})
        if 'umap' in h5['uns']:
            def collect(name, node):
                if isinstance(node, h5py.Dataset):
                    value = node[()]
                    embedding['saved_umap_parameters'][name] = value.item() if hasattr(value, 'item') else str(value)
            h5['uns/umap'].visititems(collect)
        for key in ['celltype_category', 'celltype_scTypedb']:
            catnode = obs[key]
            color_key = key+'_colors'
            if isinstance(catnode, h5py.Group) and color_key in h5['uns']:
                embedding['author_colors'][key] = dict(zip(decode(catnode['categories']).tolist(), decode(h5['uns'][color_key]).tolist()))
    small_by_name = {Path(x['local_file']).name: x for x in small_audit['downloads']}
    feature_rows = None
    gene_index = None
    raw_selected = np.full((len(cells), len(candidates)), np.nan)
    matrix_audits = []
    old_count_path = ORIG/'all_gene_library_subtype_counts.tsv.gz'
    old_count_expected = next(x['sha256'] for x in agg_audit['outputs'] if x['file'] == old_count_path.name)
    old_counts = pd.read_csv(inp(old_count_path, old_count_expected), sep='\t', index_col=0)
    for _, lib in manifest.iterrows():
        gsm = lib.gsm
        t0 = time.monotonic()
        fpath = SMALL/lib.features_url.rsplit('/',1)[-1]
        bpath = SMALL/lib.barcode_url.rsplit('/',1)[-1]
        feat = feature_table(inp(fpath, small_by_name[fpath.name]['sha256']))
        inp(bpath, small_by_name[bpath.name]['sha256'])
        with gzip.open(bpath, 'rt', encoding='utf-8') as f:
            barcodes = [s.strip() for s in f]
        assert len(barcodes) == len(set(barcodes))
        assert all(len(x)==3 and x[2]=='Gene Expression' for x in feat)
        assert len(set(x[0] for x in feat)) == len(feat)
        if feature_rows is None:
            feature_rows = feat
            gene_index = {x[0]: i for i,x in enumerate(feat)}
            candidates['covered'] = candidates.gene_id.isin(gene_index)
            candidates['raw_feature_symbol'] = [feat[gene_index[g]][1] if g in gene_index else None for g in candidates.gene_id]
            candidates['feature_row_0based'] = [gene_index.get(g, np.nan) for g in candidates.gene_id]
        else:
            assert feat == feature_rows
        idx = cells.index[cells.gsm==gsm].to_numpy()
        lookup = {b:i for i,b in enumerate(barcodes)}
        raw_cols = np.array([lookup[b] for b in cells.loc[idx, 'barcode']], dtype=int)
        assert len(set(raw_cols)) == len(idx)
        cells.loc[idx, 'raw_barcode_column_0based'] = raw_cols
        expected = next(x for x in agg_audit['libraries'] if x['gsm']==gsm)
        mpath = inp(ORIG/'sources'/lib.matrix_url.rsplit('/',1)[-1], expected['matrix_sha256'])
        with gzip.open(mpath, 'rb') as f:
            full = mmread(f)
        assert full.shape == (len(feat), len(barcodes))
        assert full.nnz == expected['declared_nnz']
        assert np.issubdtype(full.data.dtype, np.integer) and (full.data >= 0).all()
        n_entries = full.nnz
        full = full.tocsc()
        full.sum_duplicates()
        duplicates = n_entries-full.nnz
        full.eliminate_zeros()
        assert int(full.sum()) == expected['all_raw_matrix_counts']
        selected = full[:, raw_cols]
        totals = np.asarray(selected.sum(axis=0)).ravel().astype(np.int64)
        detected = selected.getnnz(axis=0).astype(np.int64)
        mito_rows = [j for j,x in enumerate(feat) if x[1].startswith('MT-')]
        mito_counts = np.asarray(selected[mito_rows,:].sum(axis=0)).ravel().astype(np.int64)
        assert (totals>0).all()
        cells.loc[idx, 'raw_total_counts'] = totals
        cells.loc[idx, 'raw_n_genes_by_counts'] = detected
        cells.loc[idx, 'raw_mito_counts'] = mito_counts
        cells.loc[idx, 'raw_pct_mito_MT_prefix'] = mito_counts/totals*100
        covered_j = np.where(candidates.covered.to_numpy())[0]
        feature_j = candidates.loc[covered_j, 'feature_row_0based'].astype(int).to_numpy()
        values = selected[feature_j,:].toarray().T
        raw_selected[np.ix_(idx,covered_j)] = values
        n_checks = 0
        for subtype in HEP_SUBTYPES:
            mask = cells.loc[idx,'celltype_scTypedb'].to_numpy()==subtype
            actual = np.asarray(selected[:,mask].sum(axis=1)).ravel()
            previous = old_counts.loc[[x[0] for x in feat], gsm+'|'+subtype].to_numpy()
            assert np.array_equal(actual, previous), (gsm, subtype)
            n_checks += len(feat)
        hep_mask = cells.loc[idx,'celltype_category'].to_numpy()=='Hepatocytes'
        assert int(totals[hep_mask].sum()) == expected['selected_hepatocyte_counts']
        matrix_audits.append({'gsm': gsm, 'raw_matrix_shape': list(full.shape),
                              'raw_nnz_entries': n_entries, 'duplicate_coordinate_entries_merged': duplicates,
                              'annotated_cells_joined': len(idx), 'annotation_barcode_missing': 0,
                              'all_features': len(feat), 'mitochondrial_MT_prefix_features': len(mito_rows),
                              'selected_total_counts': int(totals.sum()), 'whole_raw_matrix_counts': int(full.sum()),
                              'prior_all_gene_subtype_count_values_exactly_matched': n_checks,
                              'seconds': round(time.monotonic()-t0,2)})
        print(json.dumps(matrix_audits[-1]), flush=True)
        del full, selected
    for col in ['raw_barcode_column_0based','raw_total_counts','raw_n_genes_by_counts','raw_mito_counts']:
        cells[col] = cells[col].astype(np.int64)
    # The deposited field named pct_mito is numerically a 0--1 fraction.
    # Its identity is checked against the same MT- numerator / author total.
    mito_fraction_residual = np.abs(cells.raw_mito_counts/cells.author_total_counts-cells.author_pct_mito)
    assert mito_fraction_residual.max() < 1e-7
    cells['author_pct_mito_times100'] = cells.author_pct_mito*100
    cells['raw_vs_author_total_count_difference_pct'] = (cells.raw_total_counts-cells.author_total_counts)/cells.raw_total_counts*100
    cp10k = raw_selected/cells.raw_total_counts.to_numpy()[:,None]*1e4
    lognorm = np.log1p(cp10k)
    assert np.isfinite(lognorm[:, candidates.covered]).all()
    assert np.isnan(lognorm[:,~candidates.covered]).all()
    for j,cand in candidates.iterrows():
        cells[cand.symbol+'__counts'] = raw_selected[:,j]
        cells[cand.symbol+'__log1p_cp10k'] = lognorm[:,j]
    table('cell_level_umap_qc_expression.csv.gz', cells)
    table('candidate_gene_coverage.csv', candidates)
    table('candidate25_library_subtype_all10.csv', summarize_genes(cells, raw_selected, cp10k, lognorm, candidates,
                                                                 ['gsm','condition','celltype_category','celltype_scTypedb'], 'All author subtypes, each library separately'))
    pooled10 = summarize_genes(cells, raw_selected, cp10k, lognorm, candidates, ['celltype_category','celltype_scTypedb'], 'All author subtypes, four-library cell pool')
    table('candidate25_subtype_all10_pooled.csv', pooled10)
    table('candidate25_condition_subtype_all10.csv', summarize_genes(cells, raw_selected, cp10k, lognorm, candidates,
                                                                    ['condition','celltype_category','celltype_scTypedb'], 'All author subtypes, condition cell pool'))
    hep = cells[cells.celltype_category=='Hepatocytes']
    table('candidate25_hepatocyte_subtype_pooled.csv', summarize_genes(hep, raw_selected, cp10k, lognorm, candidates,
                                                                     ['celltype_scTypedb'], 'Five author hepatocyte-like subtypes, four-library cell pool'))
    table('candidate25_library_hepatocyte_subtype.csv', summarize_genes(hep, raw_selected, cp10k, lognorm, candidates,
                                                                       ['gsm','condition','celltype_scTypedb'], 'Five author hepatocyte-like subtypes, each library separately'))
    cat4 = summarize_genes(cells, raw_selected, cp10k, lognorm, candidates, ['celltype_category'], 'Four author broad categories, four-library cell pool')
    table('spotlight4_category4_pooled.csv', cat4[cat4.symbol.isin(SPOTLIGHT)])
    table('spotlight4_subtype10_pooled.csv', pooled10[pooled10.symbol.isin(SPOTLIGHT)])
    table('cell_composition_library_category4.csv', composition(cells, 'celltype_category', 'All author-retained cells in this library'))
    table('cell_composition_library_subtype10.csv', composition(cells, 'celltype_scTypedb', 'All author-retained cells in this library'))
    table('cell_composition_library_hepatocyte_subtype5.csv', composition(hep, 'celltype_scTypedb', 'Author Hepatocytes category cells in this library'))
    qc=[]
    for gsm, sub in cells.groupby('gsm',sort=False):
        row={'gsm':gsm, 'condition':sub.condition.iloc[0], 'n_cells':len(sub),
             'n_hepatocyte_like_cells':int((sub.celltype_category=='Hepatocytes').sum())}
        for col in ['raw_total_counts','raw_n_genes_by_counts','raw_pct_mito_MT_prefix']:
            for q in [0,.1,.25,.5,.75,.9,1]:
                row[f'{col}_q{int(q*100):02d}']=float(sub[col].quantile(q))
        for raw, author in [('raw_total_counts','author_total_counts'), ('raw_n_genes_by_counts','author_n_genes_by_counts'), ('raw_pct_mito_MT_prefix','author_pct_mito_times100')]:
            a=sub[author].to_numpy(dtype=float)
            b=sub[raw].to_numpy(dtype=float)
            row[raw+'_author_exact_matches']=int(np.sum(a==b))
            row[raw+'_author_max_abs_difference']=float(np.max(np.abs(a-b)))
        row['author_pct_mito_original_unit']='fraction; multiplied by 100 for percentage-point comparison'
        row['author_mito_fraction_numerator_identity_max_abs_difference']=float(mito_fraction_residual.loc[sub.index].max())
        qc.append(row)
    table('qc_library_summary.csv', pd.DataFrame(qc))
    invariants = {
        'unique_library_barcode': bool(cells.cell_key.is_unique),
        'all_umap_finite': bool(np.isfinite(cells[['umap_1','umap_2']]).all().all()),
        'unique_original_obs_ids': bool(cells.original_obs_index.is_unique),
        'selected_cell_counts': cells.groupby('gsm',sort=False).size().to_dict(),
        'selected_total_cells': len(cells), 'selected_hepatocyte_like_cells': len(hep),
        'n_categories': cells.celltype_category.nunique(), 'n_subtypes': cells.celltype_scTypedb.nunique(),
        'covered_candidate_count': int(candidates.covered.sum()),
        'missing_candidates': candidates.loc[~candidates.covered,['gene_id','symbol']].to_dict('records'),
        'no_missing_candidate_filled_zero': True,
        'all250_pooled_gene_subtype_rows': len(pooled10),
        'all_matrix_aggregate_values_match_prior': sum(x['prior_all_gene_subtype_count_values_exactly_matched'] for x in matrix_audits),
        'all_detected_pct_in_range': bool(pooled10.detected_pct.dropna().between(0,100).all()),
        'no_new_quality_exclusions': True,
        'author_pct_mito_unit': '0-to-1 fraction, verified as raw MT-prefix count divided by author total count',
        'author_mito_fraction_identity_max_abs_difference': float(mito_fraction_residual.max()),
        'max_raw_vs_author_mito_percentage_point_difference': float((cells.raw_pct_mito_MT_prefix-cells.author_pct_mito_times100).abs().max()),
        'max_raw_vs_author_total_count_difference_pct': float(cells.raw_vs_author_total_count_difference_pct.max()),
    }
    for rec in INPUT_RECORDS:
        stat=Path(rec['path']).stat()
        rec['size_mtime_unchanged_after_read'] = stat.st_size==rec['bytes'] and stat.st_mtime_ns==rec['mtime_ns']
        assert rec['size_mtime_unchanged_after_read']
    summary={'status':'PASS_DESCRIPTIVE_ANALYSIS', 'client_date':'2026-10-09',
             'executed_utc':datetime.now(timezone.utc).isoformat(), 'command':[sys.executable,*sys.argv],
             'script_sha256':sha(__file__), 'environment':{'python':sys.version, 'executable':sys.executable,
              'platform':platform.platform(), 'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'h5py':h5py.__version__},
             'embedding':embedding, 'checks':invariants,'matrix_audits':matrix_audits,
             'methods':{
              'counts':'Original 33,538-gene MTX integer values; no author processed X/layer used as counts',
              'selection':'All author-retained/annotated cells from GSM6322577-80; exact GSM+10x-barcode matching',
              'normalization':'For each cell, log1p(raw gene UMI count / total UMI across all 33,538 measured genes in that cell * 10,000); natural logarithm',
              'detection':'Raw count > 0; fraction uses all retained cells in the specified library/condition/type, including measured zeros',
              'dot_mean':'Arithmetic mean of per-cell log1p(CP10k), including measured zeros; not log1p of a mean and not z-scored',
              'missing':'Gene IDs absent from raw features remain NA, including count/mean/detection; never converted to measured zeros',
              'composition':'Cell counts divided by retained cells in the explicitly named denominator scope; never UMI shares',
              'mitochondrial':'MT- prefix in original feature gene symbols; percentage of full raw gene UMI total',
              'author_QC_comparison':'Deposited pct_mito is a fraction, verified as MT-prefix counts divided by author total_counts (maximum residual <1e-7); multiply by 100 before comparing percentages. Deposited total_counts/n_genes differ slightly from the full 33,538-feature raw matrices. The precise filtering stage responsible is not established here.',
              'pooling':'Pooled summaries are weighted by observed cell numbers, not by equal library weights; library summaries retained',
             },
             'limits':['One differentiation, four technical libraries; cells/libraries are not independent biological replicates',
              'Control/PA background only, no sesamin in these four libraries',
              'Author-provided cell labels and UMAP reused, no independent cell annotation or embedding reconstruction',
              'Subset occupies original full-atlas embedding; UMAP distances/visual segregation are not effect-size or significance measures',
              'Sampling/QC/library effects can alter observed cell composition and detection',
              'HB1 counts are 4/69/9/128 across the four libraries; sparse subgroup descriptions require care',
              'Feature NA is not absence of expression; per-cell zero is an observed count zero, not proof of biological absence',
              'No differential-expression significance, causal, therapeutic or new-mechanism inference'],
             'exploratory_failures':['Initial console inspection stopped when Windows GBK could not encode the micro symbol from a manifest; no analysis outputs affected. Subsequent inspection used ASCII-safe JSON.', 'First draft QC comparison subtracted deposited pct_mito fraction from raw percentage. Self-review and independent technical review caught the unit mismatch; final tables convert deposited fraction to percent first and verify its numerator/denominator identity. Raw counts, normalization and figure data did not change.'],
             'inputs':INPUT_RECORDS,'outputs':OUTPUTS,'elapsed_seconds':round(time.monotonic()-started,2)}
    (BASE/'analysis_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    lines=['# Single-cell extension: executed analysis','',
           'Status: PASS for descriptive source-table construction and computational checks.',
           '',f"Four author-annotated libraries contribute {len(cells):,} cells: {len(hep):,} in the author Hepatocytes category, four broad categories and ten author subtypes.",
           'The full author UMAP contains 82,467 cells; this analysis reuses its original coordinates for the four-library subset. No UMAP was recomputed and no clustering or reannotation was performed.',
           '', '## Main verification', '',
           f"All 670,760 prior all-gene library/subtype count sums were exactly reproduced from the original matrices. {int(candidates.covered.sum())}/25 candidate gene IDs occur in all four raw feature tables; missing IDs remain NA.",
           'Raw counts, total UMIs, detected genes and MT-prefix mitochondrial percentages were calculated from the 33,538-gene matrices. Each normalized value is log1p(count / all-gene cell UMI total × 10,000). Means include measured zeros.',
           'Per-cell metrics, author QC values, 1,000 library × subtype × gene records, 250 pooled subtype × gene records, five-hepatocyte-subtype records and explicitly denominated cell-composition records are saved as CSV.',
           '', '## Spotlight genes: pooled descriptive summaries', '']
    for gene in SPOTLIGHT:
        rows=pooled10[pooled10.symbol==gene].sort_values('mean_log1p_cp10k',ascending=False)
        top=rows.iloc[0]
        lines.append(f"- {gene}: largest pooled mean log1p(CP10k) in author subtype {top.celltype_scTypedb}: {top.mean_log1p_cp10k:.4f}; detected in {top.detected_pct:.2f}% of {int(top.n_cells):,} cells. This descriptive ranking is not a subtype-specificity or differential-expression test.")
    lines+=['','## QC comparison with deposited author fields','']
    for r in qc:
        lines.append(f"- {r['gsm']}: raw total count max absolute difference {r['raw_total_counts_author_max_abs_difference']:.8g}; detected genes difference {r['raw_n_genes_by_counts_author_max_abs_difference']:.8g}; MT-prefix percentage-point difference after converting author fraction to percent {r['raw_pct_mito_MT_prefix_author_max_abs_difference']:.8g}.")
    lines += ['', 'The deposited field pct_mito is a fraction: across all selected cells, raw MT-prefix counts divided by deposited author total_counts match this fraction within 7.45e-9. Multiplication by 100 aligns percentage units. The remaining mitochondrial percentage differences follow the distinct total-count denominators. Deposited total-count and detected-gene values are slightly lower than full raw-matrix values; the precise source filtering stage is not established here. All new per-cell normalization and QC plots consistently use the full raw-matrix denominator.']
    lines+=['','## Interpretation limits','']+['- '+x for x in summary['limits']]
    lines+=['','## Re-run','', '`<project>/.venv/Scripts/python.exe <project>/next_steps/single_cell_extension_14/analysis/derive_single_cell.py`',
            '', 'Input hashes are checked against prior audits when available; input sizes and modification times are checked again after analysis. Existing raw sources and Stage 13 artifacts are read only. Dependency versions and output hashes are recorded in analysis_summary.json.']
    (BASE/'analysis_summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({'status':summary['status'],'checks':invariants,'elapsed_seconds':summary['elapsed_seconds']},ensure_ascii=True),flush=True)


if __name__=='__main__':
    main()
