"""Library-aware descriptive localization and symmetric mixture decomposition.

No embedding, statistical significance test, cell exclusion or batch correction.
See the frozen protocol for the exact estimands and mathematical definitions.
"""
from __future__ import annotations
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import traceback

import numpy as np
import pandas as pd
import scipy

ROOT=Path(__file__).resolve().parents[1]
HERE=Path(sys.argv[1]).resolve()
INPUT=ROOT/'inputs/organoid'
LIBRARIES=['GSM6322577','GSM6322578','GSM6322579','GSM6322580']
CONDITIONS=['CTRL','PA']
CATEGORIES=['Hepatocytes','Cholangiocytes','Fibroblasts','Hepatic stellate cells']
SUBTYPES=['AH','FH1','FH2','HB1','HB2','Cholangiocytes','Ductal cells','Fibroblasts','Hepatic stellate cells','Smooth muscle cells']
SCHEMES=['equal_library_actual','pooled_cell_actual','separate_average_reference_diagnostic']
SPOTLIGHT=['VNN1','AREG','JAG2','CXCL8']
FROZEN_PROTOCOL='4ab6ea1a388b1c5fc66e1e3c172acbb6baeb7d6391d630357a25995e68034edd'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def safe(value):
    if isinstance(value,dict): return {str(k):safe(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [safe(v) for v in value]
    if isinstance(value,np.generic): value=value.item()
    if isinstance(value,float) and not math.isfinite(value): return None
    return value


def save_json(name,value):
    (HERE/name).write_text(json.dumps(safe(value),indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def csv(name,rows):
    frame=pd.DataFrame(rows)
    frame.to_csv(HERE/name,index=False,na_rep='NA',float_format='%.17g')
    return frame


def direction(value):
    if pd.isna(value): return 'NA'
    return 'positive' if value>1e-12 else 'negative' if value < -1e-12 else 'zero'


def pattern(values):
    signs=[direction(x) for x in values]
    return 'NA' if 'NA' in signs else 'all_positive' if all(x=='positive' for x in signs) else 'all_negative' if all(x=='negative' for x in signs) else 'all_zero' if all(x=='zero' for x in signs) else 'mixed_or_contains_zero'


def main():
    started=time.perf_counter()
    started_utc=datetime.now(timezone.utc).isoformat()
    assert sha(HERE/'protocol.md')==FROZEN_PROTOCOL, 'Frozen protocol has changed'
    filenames=['cell_level_umap_qc_expression.csv.gz','candidate_gene_coverage.csv',
       'candidate25_library_subtype_all10.csv','cell_composition_library_category4.csv']
    source=[{'file':x,'path':str(INPUT/x),'sha256':sha(INPUT/x),'bytes':(INPUT/x).stat().st_size} for x in filenames]
    data=pd.read_csv(INPUT/filenames[0])
    coverage=pd.read_csv(INPUT/'candidate_gene_coverage.csv')
    genes=coverage.symbol.tolist()
    cover=dict(zip(coverage.symbol,coverage.covered))
    gene_ids=dict(zip(coverage.symbol,coverage.gene_id))
    assert len(genes)==25 and len(set(genes))==25
    assert sum(cover.values())==24 and [g for g in genes if not cover[g]]==['AOC4P']
    assert len(data)==24034 and data.cell_key.is_unique
    assert set(data.gsm)==set(LIBRARIES) and set(data.condition)==set(CONDITIONS)
    assert set(data.celltype_category)==set(CATEGORIES) and set(data.celltype_scTypedb)==set(SUBTYPES)
    assert data.groupby('gsm').condition.nunique().eq(1).all()
    assert data.groupby('gsm').size().to_dict()==dict(zip(LIBRARIES,[3093,5510,7399,8032]))
    assert (data.raw_total_counts>0).all() and np.isfinite(data.raw_total_counts).all()
    assert np.isfinite(data[['raw_n_genes_by_counts','raw_pct_mito_MT_prefix']]).all().all()
    assert data.raw_pct_mito_MT_prefix.between(0,100).all()
    lib_condition={lib:str(data.loc[data.gsm==lib,'condition'].iloc[0]) for lib in LIBRARIES}
    assert [lib_condition[x] for x in LIBRARIES]==['CTRL','CTRL','PA','PA']
    masks_lib={lib:(data.gsm==lib).to_numpy() for lib in LIBRARIES}
    masks_type={k:(data.celltype_category==k).to_numpy() for k in CATEGORIES}
    norm={};raw={};normalization_residual={}
    for gene in genes:
        counts=data[f'{gene}__counts'].to_numpy()
        stored=data[f'{gene}__log1p_cp10k'].to_numpy()
        if cover[gene]:
            assert np.isfinite(counts).all() and (counts>=0).all() and (counts==np.floor(counts)).all()
            values=np.log1p(counts/data.raw_total_counts.to_numpy()*10000)
            normalization_residual[gene]=float(np.max(np.abs(values-stored)))
            assert normalization_residual[gene]<1e-10
            norm[gene]=values;raw[gene]=counts
        else:
            assert np.isnan(counts).all() and np.isnan(stored).all()
            norm[gene]=np.full(len(data),np.nan);raw[gene]=counts

    def summarize(mask,gene):
        n=int(np.sum(mask));assert n>0
        if not cover[gene]:
            return {'gene_id':gene_ids[gene],'symbol':gene,'covered':False,'n_cells':n,
                'n_detected_cells':np.nan,'detected_pct':np.nan,'mean_log1p_cp10k':np.nan,'sum_raw_counts':np.nan}
        y=norm[gene][mask];x=raw[gene][mask]
        detected=int((x>0).sum())
        return {'gene_id':gene_ids[gene],'symbol':gene,'covered':True,'n_cells':n,'n_detected_cells':detected,
             'detected_pct':detected/n*100,'mean_log1p_cp10k':float(y.mean()),'sum_raw_counts':int(x.sum())}

    subtype_rows=[];category_rows=[];lib_rows=[];sampling=[];qc=[]
    for lib in LIBRARIES:
        ml=masks_lib[lib];condition=lib_condition[lib];nf=int(ml.sum())
        f=data.loc[ml]
        q={'gsm':lib,'condition':condition,'n_cells':nf}
        for col in ['raw_total_counts','raw_n_genes_by_counts','raw_pct_mito_MT_prefix']:
            values=f[col].to_numpy()
            q.update({f'{col}_mean':float(values.mean()),f'{col}_median':float(np.median(values)),
                f'{col}_min':float(values.min()),f'{col}_max':float(values.max())})
        for category in CATEGORIES:
            mask=ml&masks_type[category];n=int(mask.sum());assert n>0
            sampling.append({'gsm':lib,'condition':condition,'celltype_category':category,'n_cells':n,
                'denominator_n_cells':nf,'fraction':n/nf,'cell_fraction_pct':n/nf*100})
            q[f'fraction_{category}']=n/nf
            for gene in genes:
                category_rows.append({'gsm':lib,'condition':condition,'celltype_category':category,**summarize(mask,gene)})
        for subtype in SUBTYPES:
            mask=ml&(data.celltype_scTypedb==subtype).to_numpy()
            assert mask.sum()>0
            for gene in genes:
                subtype_rows.append({'gsm':lib,'condition':condition,'celltype_scTypedb':subtype,**summarize(mask,gene)})
        for gene in genes: lib_rows.append({'gsm':lib,'condition':condition,**summarize(ml,gene)})
        qc.append(q)
    sub=csv('library_subtype_candidates_all25.csv',subtype_rows)
    cats=csv('library_category_candidates_all25.csv',category_rows)
    libs=csv('library_candidates_all25.csv',lib_rows)
    sample=csv('category_sampling_counts.csv',sampling)
    csv('library_qc_context.csv',qc)
    assert len(sub)==1000 and len(cats)==400 and len(libs)==100
    assert np.allclose(sample.groupby('gsm').fraction.sum(),1,rtol=0,atol=1e-14)

    old=pd.read_csv(INPUT/'candidate25_library_subtype_all10.csv')
    joined=sub.merge(old,on=['gsm','celltype_scTypedb','symbol'],validate='one_to_one',suffixes=('_new','_old'))
    assert len(joined)==1000
    comparisons={}
    for col in ['n_cells','n_detected_cells','detected_pct','mean_log1p_cp10k','sum_raw_counts']:
        a=joined[f'{col}_new'].to_numpy();b=joined[f'{col}_old'].to_numpy()
        assert np.array_equal(np.isnan(a),np.isnan(b))
        valid=np.isfinite(a)
        residual=float(np.max(np.abs(a[valid]-b[valid])))
        assert np.allclose(a[valid],b[valid],rtol=1e-9,atol=1e-10),col
        comparisons[col]={'max_abs_difference':residual,'finite_values_compared':int(valid.sum()),'NA_entries_retained':int((~valid).sum())}
    original_composition=pd.read_csv(INPUT/'cell_composition_library_category4.csv')
    for row in sample.itertuples():
        ref=original_composition[(original_composition.gsm==row.gsm)&(original_composition.celltype_category==row.celltype_category)].iloc[0]
        assert row.n_cells==ref.n_cells and row.denominator_n_cells==ref.denominator_n_cells
        assert np.isclose(row.cell_fraction_pct,ref.cell_fraction_pct,rtol=1e-9,atol=1e-10)

    # Within-library comparisons are kept for every fixed candidate, not just AREG.
    contrast_rows=[];robustness=[]
    cat_lookup=cats.set_index(['gsm','celltype_category','symbol'])
    for gene in genes:
        gened=[];genem=[]
        for lib in LIBRARIES:
            ch=cat_lookup.loc[(lib,'Cholangiocytes',gene)];he=cat_lookup.loc[(lib,'Hepatocytes',gene)]
            d=float(ch.detected_pct-he.detected_pct);m=float(ch.mean_log1p_cp10k-he.mean_log1p_cp10k)
            gened.append(d);genem.append(m)
            contrast_rows.append({'symbol':gene,'gene_id':gene_ids[gene],'covered':bool(cover[gene]),'gsm':lib,'condition':lib_condition[lib],
              'cholangiocyte_cells':int(ch.n_cells),'hepatocyte_cells':int(he.n_cells),
              'cholangiocyte_detection_pct':ch.detected_pct,'hepatocyte_detection_pct':he.detected_pct,
              'detection_difference_percentage_points':d,'detection_direction':direction(d),
              'cholangiocyte_mean_log1p_cp10k':ch.mean_log1p_cp10k,'hepatocyte_mean_log1p_cp10k':he.mean_log1p_cp10k,
              'mean_expression_difference':m,'mean_expression_direction':direction(m)})
        ch_pool=summarize(masks_type['Cholangiocytes'],gene);he_pool=summarize(masks_type['Hepatocytes'],gene)
        pool_d=ch_pool['detected_pct']-he_pool['detected_pct'];pool_m=ch_pool['mean_log1p_cp10k']-he_pool['mean_log1p_cp10k']
        pattern_d=pattern(gened);pattern_m=pattern(genem)
        robustness.append({'symbol':gene,'gene_id':gene_ids[gene],'covered':bool(cover[gene]),
            'detection_within_library_pattern':pattern_d,'mean_expression_within_library_pattern':pattern_m,
            'pooled_cell_detection_difference_pp':pool_d,'equal_library_detection_difference_pp':float(np.mean(gened)),
            'min_library_detection_difference_pp':float(np.min(gened)),'max_library_detection_difference_pp':float(np.max(gened)),
            'pooled_cell_mean_expression_difference':pool_m,'equal_library_mean_expression_difference':float(np.mean(genem)),
            'min_library_mean_expression_difference':float(np.min(genem)),'max_library_mean_expression_difference':float(np.max(genem)),
            'complete_Simpson_reversal_detection':bool((pattern_d=='all_positive' and pool_d<0)or(pattern_d=='all_negative' and pool_d>0)) if cover[gene] else None,
            'complete_Simpson_reversal_mean_expression':bool((pattern_m=='all_positive' and pool_m<0)or(pattern_m=='all_negative' and pool_m>0)) if cover[gene] else None})
    contrast=csv('cholangiocyte_vs_hepatocyte_per_library_all25.csv',contrast_rows)
    robust=csv('cholangiocyte_vs_hepatocyte_robustness_all25.csv',robustness)

    profiles=[];terms=[];totals=[];identity_max=0.;direct_max=0.;reverse_max=0.;cov_max=0.
    for gene in genes:
        condition_profile={}
        for condition in CONDITIONS:
            chosen=[l for l in LIBRARIES if lib_condition[l]==condition]
            assert len(chosen)==2
            nlib=np.array([int(masks_lib[l].sum()) for l in chosen],dtype=float)
            ncat=np.array([[int((masks_lib[l]&masks_type[k]).sum()) for k in CATEGORIES] for l in chosen],dtype=float)
            p=ncat/nlib[:,None]
            assert np.isfinite(p).all() and np.allclose(p.sum(axis=1),1,rtol=0,atol=1e-14)
            if cover[gene]:
                mu=np.array([[float(cat_lookup.loc[(l,k,gene)].mean_log1p_cp10k) for k in CATEGORIES] for l in chosen])
                m=np.array([float(norm[gene][masks_lib[l]].mean()) for l in chosen])
                assert np.allclose(np.sum(p*mu,axis=1),m,rtol=0,atol=1e-12)
                direct_el=float(m.mean())
                direct_pool=float(norm[gene][(data.condition==condition).to_numpy()].mean())
                p_el=p.mean(axis=0)
                mu_el=np.sum(p*mu,axis=0)/p.sum(axis=0)
                p_pool=ncat.sum(axis=0)/nlib.sum()
                mu_pool=np.sum(ncat*mu,axis=0)/ncat.sum(axis=0)
                mu_ref=mu.mean(axis=0)
                covariance=float(np.mean(np.sum((p-p_el)*(mu-mu_ref),axis=1)))
                values={'equal_library_actual':(p_el,mu_el,direct_el),
                        'pooled_cell_actual':(p_pool,mu_pool,direct_pool),
                        'separate_average_reference_diagnostic':(p_el,mu_ref,float(np.sum(p_el*mu_ref)))}
                assert np.isclose(direct_el-values[SCHEMES[2]][2],covariance,rtol=0,atol=1e-12)
            else:
                covariance=np.nan;direct_el=np.nan;direct_pool=np.nan
                values={scheme:(p.mean(axis=0) if scheme!=SCHEMES[1] else ncat.sum(axis=0)/nlib.sum(),np.full(4,np.nan),np.nan) for scheme in SCHEMES}
            condition_profile[condition]={'values':values,'covariance':covariance,'direct_el':direct_el,'direct_pool':direct_pool}
            for scheme,(ps,mus,direct) in values.items():
                for idx,category in enumerate(CATEGORIES):
                    profiles.append({'symbol':gene,'gene_id':gene_ids[gene],'covered':bool(cover[gene]),'condition':condition,
                        'scheme':scheme,'celltype_category':category,'category_fraction':ps[idx],
                        'category_mean_log1p_cp10k':mus[idx],'category_contribution_to_profile_mean':ps[idx]*mus[idx],
                        'profile_mean':direct,'observed_cells_in_category':int(ncat[:,idx].sum()),'technical_libraries':2})
        for scheme in SCHEMES:
            p0,mu0,direct0=condition_profile['CTRL']['values'][scheme]
            p1,mu1,direct1=condition_profile['PA']['values'][scheme]
            composition=(p1-p0)*(mu1+mu0)/2
            within=(mu1-mu0)*(p1+p0)/2
            actual_el_delta=condition_profile['PA']['direct_el']-condition_profile['CTRL']['direct_el']
            cov_delta=condition_profile['PA']['covariance']-condition_profile['CTRL']['covariance']
            if cover[gene]:
                mean0=float(np.dot(p0,mu0));mean1=float(np.dot(p1,mu1));delta=mean1-mean0
                ctotal=float(composition.sum());wtotal=float(within.sum());residual=delta-ctotal-wtotal
                scalar=sum((float(p1[i])-float(p0[i]))*(float(mu1[i])+float(mu0[i]))/2+
                           (float(mu1[i])-float(mu0[i]))*(float(p1[i])+float(p0[i]))/2 for i in range(4))
                reverse_c=float(np.sum((p0-p1)*(mu0+mu1)/2));reverse_w=float(np.sum((mu0-mu1)*(p0+p1)/2))
                assert abs(residual)<1e-12 and abs(delta-scalar)<1e-12
                assert abs(mean0-direct0)<1e-12 and abs(mean1-direct1)<1e-12
                assert abs(reverse_c+ctotal)<1e-12 and abs(reverse_w+wtotal)<1e-12
                identity_max=max(identity_max,abs(residual),abs(delta-scalar))
                direct_max=max(direct_max,abs(mean0-direct0),abs(mean1-direct1))
                reverse_max=max(reverse_max,abs(reverse_c+ctotal),abs(reverse_w+wtotal))
                if scheme==SCHEMES[2]:
                    diagnostic_residual=actual_el_delta-(ctotal+wtotal+cov_delta)
                    assert abs(diagnostic_residual)<1e-12
                    cov_max=max(cov_max,abs(diagnostic_residual))
                else: diagnostic_residual=np.nan
            else:
                mean0=mean1=delta=ctotal=wtotal=residual=diagnostic_residual=np.nan
            for idx,k in enumerate(CATEGORIES):
                terms.append({'symbol':gene,'gene_id':gene_ids[gene],'covered':bool(cover[gene]),'scheme':scheme,'celltype_category':k,
                   'control_fraction':p0[idx],'PA_fraction':p1[idx],'control_mean_log1p_cp10k':mu0[idx],'PA_mean_log1p_cp10k':mu1[idx],
                   'composition_term':composition[idx],'within_category_term':within[idx],'category_total_difference':composition[idx]+within[idx]})
            totals.append({'symbol':gene,'gene_id':gene_ids[gene],'covered':bool(cover[gene]),'scheme':scheme,
               'control_profile_mean':mean0,'PA_profile_mean':mean1,'total_PA_minus_control':delta,
               'composition_total':ctotal,'within_category_total':wtotal,'decomposition_identity_residual':residual,
               'actual_equal_library_delta':actual_el_delta,
               'diagnostic_covariance_delta':cov_delta if scheme==SCHEMES[2] else np.nan,
               'diagnostic_actual_EL_reconstruction_residual':diagnostic_residual,
               'interpretation':'Reference mixture only; add covariance delta to reconstruct actual EL contrast' if scheme==SCHEMES[2] else 'Actual condition-mean contrast under stated library weighting'})
    profile_df=csv('condition_category_profiles_all25.csv',profiles)
    term_df=csv('kitagawa_terms_all25.csv',terms)
    total_df=csv('kitagawa_summary_all25.csv',totals)
    assert len(profile_df)==600 and len(term_df)==300 and len(total_df)==75
    assert total_df.loc[total_df.symbol=='AOC4P','total_PA_minus_control'].isna().all()
    assert total_df.loc[total_df.covered,'total_PA_minus_control'].notna().all()
    assert np.allclose(profile_df.groupby(['symbol','condition','scheme']).category_fraction.sum(),1,rtol=0,atol=1e-14)
    comparison=[]
    for gene in genes:
        el=total_df[(total_df.symbol==gene)&(total_df.scheme==SCHEMES[0])].iloc[0]
        pool=total_df[(total_df.symbol==gene)&(total_df.scheme==SCHEMES[1])].iloc[0]
        comparison.append({'symbol':gene,'covered':bool(cover[gene]),
          'equal_library_total':el.total_PA_minus_control,'pooled_cell_total':pool.total_PA_minus_control,
          'equal_library_composition':el.composition_total,'pooled_cell_composition':pool.composition_total,
          'equal_library_within':el.within_category_total,'pooled_cell_within':pool.within_category_total,
          'total_direction_equal_library':direction(el.total_PA_minus_control),'total_direction_pooled':direction(pool.total_PA_minus_control),
          'total_direction_agrees':direction(el.total_PA_minus_control)==direction(pool.total_PA_minus_control) if cover[gene] else None})
    comparison_df=csv('weighting_sensitivity_all25.csv',comparison)
    unchanged=[{'file':r['file'],'sha256_before':r['sha256'],'sha256_after':sha(Path(r['path'])),'unchanged':r['sha256']==sha(Path(r['path']))} for r in source]
    assert all(r['unchanged'] for r in unchanged)
    validation={'status':'PASS','cell_count':24034,'fixed_candidate_count':25,'covered_candidate_count':24,'missing_candidate':'AOC4P',
        'source_read_only_checks':unchanged,'normalized_value_max_abs_difference_from_stage14':max(normalization_residual.values()),
        'stage14_subtype_summary_comparison':comparisons,'category_fraction_sum_checks':'PASS',
        'decomposition_max_abs_identity_residual':identity_max,'profile_vs_direct_mean_max_abs_residual':direct_max,
        'condition_swap_sign_reversal_max_abs_residual':reverse_max,'diagnostic_covariance_identity_max_abs_residual':cov_max,
        'library_subtype_rows':1000,'library_category_rows':400,'library_overall_rows':100,
        'condition_category_profile_rows':600,'decomposition_term_rows':300,'decomposition_summary_rows':75,
        'hypothesis_tests_run':False,'confidence_intervals_computed':False,'correlations_computed':False,
        'UMAP_recomputed':False,'new_exclusions':0,'raw_h5ad_or_MTX_read':False}
    save_json('numerical_validation.json',validation)
    highlights={gene:{
        'localization_per_library':contrast.loc[contrast.symbol==gene].to_dict(orient='records'),
        'localization_robustness':robust.loc[robust.symbol==gene].iloc[0].to_dict(),
        'decomposition':total_df.loc[total_df.symbol==gene].to_dict(orient='records')
        } for gene in SPOTLIGHT}
    sensitivity= comparison_df.loc[comparison_df.covered]
    summary={'status':'PASS_DESCRIPTIVE_ANALYSIS','date':'2026-10-09','primary_scheme':'equal_library_actual',
       'sensitivity_scheme':'pooled_cell_actual','diagnostic_scheme':'separate_average_reference_diagnostic',
       'fixed_candidates':genes,'missing':'AOC4P','highlights':highlights,
       'weighted_total_sign_agreement':int(sensitivity.total_direction_agrees.sum()),'weighted_total_sign_denominator':24,
       'weighted_total_sign_disagreements':sensitivity.loc[~sensitivity.total_direction_agrees.astype(bool),'symbol'].tolist(),
       'library_qc':qc,'validation':validation,
       'limits':['Four technical libraries from one differentiation; no independent biological replication',
          'Control/PA background only; no sesamin treatment',
          'Mean log1p(CP10k) decomposition is descriptive algebra, not causal adjustment or a protective mechanism',
          'Detection differences may depend on sequencing depth, recovery and cell mixtures; no batch effect correction or significance test',
          'Author-provided broad and subtype labels retained; no new annotation or clustering'],
       'skills_used':['validate-data: exact estimand definitions, completeness and numerical crosschecks'],
       'protocol_sha256':FROZEN_PROTOCOL,'inputs':source}
    save_json('analysis_summary.json',summary)
    lines=['# Library-aware single-cell extension','',
      'This executed descriptive analysis retains all 25 fixed candidates and every one of the 24,034 source cells. Twenty-four candidates are measured; AOC4P remains NA. Four technical libraries from one differentiation contain control/PA conditions without sesamin.','',
      '## Per-library cellular localization','',
      '| Candidate | Detection difference: Cholangiocytes minus Hepatocytes, C1 / C2 / P1 / P2 (percentage points) | Mean-expression sign in all four libraries |',
      '|---|---|---|']
    for gene in SPOTLIGHT:
        rows=contrast[contrast.symbol==gene];r=robust[robust.symbol==gene].iloc[0]
        numbers=' / '.join(f'{v:+.4f}' for v in rows.detection_difference_percentage_points)
        lines.append(f'| {gene} | {numbers} | {r.mean_expression_within_library_pattern} |')
    lines+=['','Full four-library results are retained for all 25 candidates, including contrary and zero differences. The broad Cholangiocytes category and its Ductal cells subtype are not interchangeable labels.','',
      '## Actual condition-mean decomposition','',
      'Primary equal-library means give each library weight 1/2. Composition and within-category terms use the symmetric Kitagawa formula and sum to the actual PA-minus-control mean difference. Pooled-cell weighting is a sensitivity analysis. Terms are in mean log1p(CP10k) units, not fold changes or causal effect proportions.','',
      '| Candidate | Weighting | Total difference | Composition | Within-category |',
      '|---|---|---:|---:|---:|']
    for gene in SPOTLIGHT:
        for scheme in SCHEMES[:2]:
            r=total_df[(total_df.symbol==gene)&(total_df.scheme==scheme)].iloc[0]
            lines.append(f'| {gene} | {scheme} | {r.total_PA_minus_control:+.8f} | {r.composition_total:+.8f} | {r.within_category_total:+.8f} |')
    lines += ['',f'Across all 24 measured candidates, total-difference direction agrees between the two actual-mean weightings for {summary["weighted_total_sign_agreement"]}/24. Discordant candidates: '+(', '.join(summary['weighted_total_sign_disagreements']) or 'none')+'.','',
      'The diagnostic separate-average reference mixture is not presented as the observed equal-library total. Its library composition-expression covariance residual is reported explicitly and reconstructs the actual equal-library contrast.','',
      '## Technical-library context and checks','',
      'All four libraries retain separate QC and expression summaries. No four-point correlation coefficient, P value, confidence interval, new QC cutoff, embedding or batch correction was computed. Per-cell normalization and all 1,000 library-by-subtype candidate summaries were independently reconstructed and matched the stage 14 tables.','',
      f'Maximum decomposition identity residual: {identity_max:.3g}; maximum direct-mean reconstruction residual: {direct_max:.3g}; diagnostic covariance reconstruction residual: {cov_max:.3g}. Control/PA swapping reverses both terms exactly within numerical precision. Inputs remained unchanged.','',
      'These results strengthen the transparency of cellular localization and weighting sensitivity. They do not supply independent biological validation, establish a sesamin response in single cells, or separate technical-library, cell-recovery and biological effects.']
    (HERE/'analysis_report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    elapsed=time.perf_counter()-started
    outputs=[{'file':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(HERE.iterdir()) if p.is_file() and p.name not in ['execution_log.json','artifact_manifest.json','execution_failure.json']]
    log={'status':'completed','started_utc':started_utc,'finished_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':elapsed,
      'client_date':'2026-10-09','command':[sys.executable,str(Path(__file__).resolve())],
      'environment':{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__},
      'script_sha256':sha(Path(__file__)),'protocol_sha256':FROZEN_PROTOCOL,'inputs':source,'outputs':outputs,
      'dependencies_installed':False,'external_transmissions':False,'failure_record_present':(HERE/'execution_failure.json').exists()}
    save_json('execution_log.json',log)
    manifest=[{'file':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='artifact_manifest.json']
    save_json('artifact_manifest.json',manifest)
    print(json.dumps({'status':summary['status'],'elapsed_seconds':round(elapsed,2),
       'AREG':safe(highlights['AREG']['localization_robustness']),
       'spotlight_decomposition':safe(total_df[(total_df.symbol.isin(SPOTLIGHT))&(total_df.scheme.isin(SCHEMES[:2]))].to_dict(orient='records')),
       'all24_weighting_sign_agreement':summary['weighted_total_sign_agreement'],'identity_max_residual':identity_max}))


if __name__=='__main__':
    try:
        main()
    except Exception:
        save_json('execution_failure.json',{'time_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(Path(__file__)),
            'command':[sys.executable,str(Path(__file__).resolve())],'traceback':traceback.format_exc(),
            'interpretation':'Failed attempt; subsequent successful results, if present, are separately logged.'})
        raise
