"""Independent-formula checks and bounded report for stage 15 sensitivity."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import statistics

import numpy as np
import pandas as pd

P = Path(__file__).resolve().parent
N = P.parents[1]
M = Path('D:/ResearchData/sesamin_mouse_2026-10-08/matrices')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    s = json.loads((P / 'analysis_summary.json').read_text(encoding='utf-8'))
    primary = pd.read_csv(P / 'frozen_primary_candidates.csv')
    submitted = pd.read_csv(M / 'candidate_mouse_descriptive.tsv', sep='\t').set_index('original_human_id')
    configs = pd.read_csv(P / 'enumerated_configurations.csv')
    result = pd.read_csv(P / 'all_configuration_effects.csv.gz')
    c = pd.read_csv(P / 'cross_species_configuration_consensus.csv')
    checks = []
    def check(name, condition):
        checks.append({'check': name, 'pass': bool(condition)})
        assert condition, name
    # Independent mouse path: stdlib mean/median + math.log2 from archived sample table.
    contrasts = [('model_vs_control', 1, 0), ('treatment_vs_model', 2, 1), ('treatment_vs_control', 2, 0)]
    checked_mouse = 0
    for r in primary.itertuples(index=False):
        values = [float(submitted.loc[r.original_human_id, f'DRR{i}_TPM']) for i in range(1010532, 1010541)]
        for cfg in configs.itertuples(index=False):
            groups = [[int(i) for i in g.split(',')] for g in cfg.retained_indices.split('|')]
            for agg, fn in [('mean', statistics.mean), ('median', statistics.median)]:
                centers = [fn([values[i] for i in group]) for group in groups]
                rows = result[(result.species == 'mouse') & (result.metric == 'TPM') &
                              (result.human_symbol == r.human_symbol) & (result.configuration == cfg.configuration) &
                              (result.aggregation == agg)].set_index('contrast')
                for contrast, a, b in contrasts:
                    v = math.log2((centers[a] + .1) / (centers[b] + .1))
                    assert math.isclose(v, float(rows.loc[contrast, 'log2_ratio']), abs_tol=1e-12)
                    checked_mouse += 1
    check('5760 mouse configuration contrast values independently recomputed', checked_mouse == 5760)
    # Recompute consensus by straightforward nested loops, independent of array broadcasting.
    checked_joint = 0
    for agg in ['mean', 'median']:
        for gene in primary.human_symbol:
            rr = result[(result.aggregation == agg) & (result.human_symbol == gene)]
            h = rr[(rr.species == 'human') & (rr.metric == 'median_ratio_reestimated')].pivot(index='configuration', columns='contrast', values='direction')
            m = rr[(rr.species == 'mouse') & (rr.metric == 'TPM')].pivot(index='configuration', columns='contrast', values='direction')
            joint = 0
            for hi, mi in itertools.product(range(64), repeat=2):
                if all(h.at[hi, name] != 0 and h.at[hi, name] == m.at[mi, name]
                       for name in ['model_vs_control', 'treatment_vs_model']):
                    joint += 1
            expected = int(c[(c.aggregation == agg) & (c.human_symbol == gene)].iloc[0].both_model_and_treatment_concordant_n)
            assert joint == expected
            checked_joint += 1
    check('30 candidate-summary joint fractions independently enumerated', checked_joint == 30)
    # Verify selected human combinations from raw counts using log-space normalization.
    hn = [f'{g}_{i}' for g in ['Control', 'Model', 'Treat'] for i in [1, 2, 3]]
    raw = pd.read_csv(N / 'candidate_screen_03/sesamin_inputs/GSE248903_total_count.anno.txt.gz', sep='\t', usecols=['name'] + hn).set_index('name')[hn]
    checked_human = 0
    for cid in [0, 1, 4, 16, 21, 63]:
        cfg = configs.loc[configs.configuration == cid].iloc[0]
        groups = [[int(i) for i in g.split(',')] for g in cfg.retained_indices.split('|')]
        retained = [i for g in groups for i in g]
        active = raw.iloc[:, retained]
        positive = active.gt(0).all(axis=1)
        logref = np.log(active.loc[positive])
        ratios = np.exp(logref.sub(logref.mean(axis=1), axis=0))
        factors = ratios.median(axis=0)
        factors = factors / math.exp(float(np.log(factors).mean()))
        normalized = active / factors
        for r in primary.itertuples(index=False):
            centers = [statistics.mean(normalized.loc[r.original_human_id, [hn[i] for i in g]].tolist()) for g in groups]
            rows = result[(result.species == 'human') & (result.metric == 'median_ratio_reestimated') &
                          (result.human_symbol == r.human_symbol) & (result.configuration == cid) &
                          (result.aggregation == 'mean')].set_index('contrast')
            for contrast, a, b in contrasts:
                value = math.log2((centers[a]+1)/(centers[b]+1))
                assert math.isclose(value, float(rows.loc[contrast, 'log2_ratio']), abs_tol=1e-12)
                checked_human += 1
    check('270 human configuration contrasts independently recomputed', checked_human == 270)
    vals = pd.read_csv(P / 'candidate_library_values.csv')
    low = []
    for gene in primary.human_symbol:
        h = vals[(vals.human_symbol == gene) & (vals.species == 'human')]
        m = vals[(vals.human_symbol == gene) & (vals.species == 'mouse')]
        low.append({'human_symbol': gene, 'human_raw_positive_libraries': int(h.raw_or_estimated_count.gt(0).sum()),
                    'human_CPM_ge1_libraries': int(h.CPM.ge(1).sum()),
                    'mouse_TPM_positive_libraries': int(m.TPM.gt(0).sum()),
                    'mouse_TPM_ge1_libraries': int(m.TPM.ge(1).sum()),
                    'mouse_TPM_min': m.TPM.min(), 'mouse_TPM_max': m.TPM.max(),
                    'mouse_all9_TPM_below1': bool(m.TPM.lt(1).all())})
    low = pd.DataFrame(low)
    low.to_csv(P / 'candidate_abundance_flags.csv', index=False)
    check('5 mouse low-expression flags reproduce archived candidate labels',
          set(low.loc[low.mouse_all9_TPM_below1, 'human_symbol']) == {'NEGR1', 'BCAS1', 'NGFR', 'JAG2', 'RAB3B'})
    sums = pd.read_csv(P / 'configuration_effect_summary.csv')
    hmain = sums[(sums.species == 'human') & (sums.metric == 'median_ratio_reestimated') &
                 (sums.aggregation == 'mean') & sums.contrast.isin(['model_vs_control', 'treatment_vs_model'])]
    check('all 15 human candidates retain both primary directions in all64 configurations',
          len(hmain) == 30 and hmain.baseline_direction_retained_n.eq(64).all())
    for output in s['outputs']:
        check('output_hash:' + output['file'], sha(P / output['file']) == output['sha256'])
    audit = {'status': 'PASS_BOUNDED_SELF_CHECK', 'utc': datetime.now(timezone.utc).isoformat(),
             'check_script_sha256': sha(__file__), 'analysis_script_sha256': s['script_sha256'],
             'checks': checks, 'numeric_values_recomputed': checked_mouse + checked_human,
             'cross_species_summary_rows_recomputed': checked_joint,
             'limits': 'Same Codex worker using a separate calculation path; not an independent expert review. PCA was numerically reconstructed; no figure generated in this branch.'}
    (P / 'calculation_checks.json').write_text(json.dumps(audit, indent=2) + '\n', encoding='utf-8')
    mean = c[c.aggregation == 'mean']
    lines = []
    for r in mean.itertuples(index=False):
        lines.append(f'| {r.human_symbol} | {r.treatment_vs_model_concordant_n}/4096 | {r.both_model_and_treatment_concordant_n}/4096 |')
    text = '''# Stage 15 candidate robustness and full-filtered-transcriptome QC

Status: executed and formula-checked; share with stated caveats. The frozen 15 primary one-to-one orthologues were retained without reselection. All 18 libraries remain. Original inputs and earlier outputs are unchanged. No biological hypothesis tests, confidence intervals, new experiments or AI performance scores were produced.

## What was added

The earlier human analysis already tested 27 simultaneous one-library-per-group omissions; this extension does not call that analysis new. It enumerates every configuration retaining two or three libraries per group (64 per dataset), adds mouse-library sensitivity, and compares the full 64 × 64 Cartesian product across species without inventing paired libraries. Arithmetic means are primary; medians are a sensitivity. Human normalization factors are re-estimated in each configuration, with frozen factors and CPM as method checks. Mouse TPM is primary; estimated-count median-ratio normalization is a secondary measurement sensitivity, without tximport length correction or a differential-expression model.

The 4096 fractions are **fractions of the enumerated configurations**, not probabilities of replication, confidence levels, independent biological observations or evidence from new data. Configurations overlap and share libraries. The human candidates were selected on these data; stability does not establish independent validation. The model group is shared by two contrasts.

## Main result

VNN1 was the only one of the 15 primary candidates whose model and treatment contrast directions remained concordant between the human and mouse datasets across every configuration pair (4096/4096). AREG retained both directions in 3072/4096 (75.0%); JAG2 did so in 2560/4096 (62.5%). These are complete mean-based enumeration results. All 15 human candidates retained their two primary signs in all 64 human configurations, which is expected to be favored by their previous selection and stability filtering. The instability of AREG/JAG2 in this extension arose in the mouse treatment contrast.

Mouse sess versus HFD percentage differences, calculated directly from TPM means without a pseudocount, ranged from −48.4% to −27.2% for VNN1; −23.2% to +9.8% for AREG; and −76.6% to +78.1% for JAG2. Thus a consistent reduction is supported for VNN1 within the observed libraries, whereas AREG and JAG2 can change direction with library composition. These min–max ranges are sensitivity ranges, not confidence intervals.

The all-library mean-to-median change preserves the primary treatment direction of all three spotlight genes on TPM. Across the full configurations, using medians yields 4096/4096 for VNN1, 2816/4096 for AREG and 3072/4096 for JAG2 for both directions. Median-ratio-normalized mouse estimated counts preserve the VNN1 reduction in all configurations; the all-library median JAG2 treatment direction becomes positive. This reinforces caution for JAG2, which is below 1 TPM in all nine mouse libraries. RAB3B is also below 1 TPM throughout; stability alone does not overcome low abundance.

Four candidates—VNN1, PALLD, RAB3B and SOCS2—retain intervention-direction agreement in all 4096 configurations; only VNN1 also retains the model-direction agreement. These are overlapping summaries, not newly discovered genes.

| Human symbol | Treatment-direction concordance | Both-direction concordance |
|---|---:|---:|
''' + '\n'.join(lines) + '''

## Transcriptome QC

PCA uses all genes passing an expression filter independent of submitted group labels, separately within each species. Human: 13,014 genes with CPM ≥1 in at least three libraries, log2(median-ratio normalized count+1), gene-centered without variance scaling; PC1 and PC2 explain 45.0% and 16.7%. Mouse: 15,340 genes with TPM ≥1 in at least three libraries, log2(TPM+1), the same centering; PC1 and PC2 explain 40.6% and 14.2%. The previous human top-500-variable-gene PCA is not overwritten, and its variance fractions must not be mixed with these.

Human pairwise Pearson correlations range from 0.9622 to 0.9915; mouse correlations range from 0.8759 to 0.9718. No exact duplicate profiles were found. The mouse ctrl libraries have greater spread in the first two principal components, including DRR1010534 separated along PC1. This is an observation, not a diagnosis of technical failure or biological independence. No library was removed. PCA separation cannot by itself establish treatment attribution, recovery, or a causal sesamin effect. The human treatment libraries do not globally return to the control position in this PCA; candidate-specific directionality must not be generalized to transcriptome-wide normalization.

## Traceability and limits

- `protocol.json` and `protocol_freeze.json`: pre-calculation rules and SHA.
- `run_robustness.py`: deterministic analysis; `analysis_summary.json`: version, invocation, eight input hashes, output hashes, original-result reproductions and QC.
- `check_and_report.py` / `calculation_checks.json`: 6,030 numeric effect checks using a separate calculation path, plus 30 complete cross-species summary checks. Same worker self-check, not independent human or expert review.
- `candidate_library_values.csv` / `candidate_abundance_flags.csv`: all 270 species–gene–library observations and low-expression flags.
- `all_configuration_effects.csv.gz` / `configuration_effect_summary.csv`: every result, including zero, unstable and discordant candidates.
- `cross_species_configuration_consensus.csv`: mean and median results for all 15, with baseline signs and counts.
- `human_*` / `mouse_*`: PCA scores, all eight nontrivial component loadings, variance fractions, exact gene filter, sample QC and full Pearson/Spearman matrices.

The exact independent culture/animal units and mouse RNA treatment mapping remain unresolved. The mouse HFD string is a submitted group label, not evidence that these are liver tissues from high-fat-fed animals. These calculations cannot remedy that provenance gap. They strengthen transparency and prioritization within existing data; they do not establish a new mechanism. VNN1 was previously reported from the same human dataset, so its strongest stability is not biological novelty. Any later external intervention study must remain a distinct analysis, with its own design and mapping restrictions.

## Suggested manuscript wording

**Methods.** We fixed the 15 previously eligible one-to-one orthologues and enumerated every configuration retaining either two or three of the three submitted libraries per group (64 configurations in each dataset). Human median-ratio factors were re-estimated per configuration, while mouse gene TPM remained on its archived scale. Direction concordance was evaluated across all 4096 human–mouse configuration pairs. Group medians, frozen human normalization factors, CPM, and normalized mouse estimated counts were examined as descriptive sensitivity analyses. Enumeration fractions and min–max effect ranges were not interpreted as probabilities, independent replication or confidence intervals.

**Results.** VNN1 retained concordant model- and treatment-associated directions across all 4096 configuration pairs, whereas AREG and JAG2 retained both directions in 75.0% and 62.5%, respectively. The mouse VNN1 sess-group mean remained 27.2–48.4% below the HFD-label mean across the enumerated configurations. AREG and JAG2 treatment effects crossed zero under library omission; JAG2 also remained below 1 TPM in every library. These analyses distinguish a stable observed VNN1 pattern from weaker, library-dependent candidates while preserving all 15 candidates and all libraries.

**Discussion.** Stability within the selection and comparison datasets improves the precision of follow-up priorities but is not independent validation. The strongest signal, VNN1, was already reported from the source human dataset; AREG remains a context-sensitive lead and JAG2 warrants particular caution because of low abundance and measurement sensitivity.
'''
    (P / 'report.md').write_text(text, encoding='utf-8')
    print(json.dumps({'status': audit['status'], 'numeric_values': audit['numeric_values_recomputed'], 'checks': len(checks)}))


if __name__ == '__main__':
    main()
