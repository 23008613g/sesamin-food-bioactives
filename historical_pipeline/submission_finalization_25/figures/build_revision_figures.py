"""Revise the study diagram and stability display from frozen numerical inputs."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'scripts'))
from audit_panel_alignment import require_matplotlib_panel_alignment

WIDTH = 183 / 25.4
INK = '#203644'
BLUE = '#42758D'
ORANGE = '#B87345'
GREY = '#97A2A9'
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 7.3, 'axes.labelsize': 7.3, 'axes.titlesize': 8,
    'xtick.labelsize': 6.8, 'ytick.labelsize': 7, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.linewidth': .65, 'pdf.fonttype': 42,
    'svg.fonttype': 'none', 'legend.frameon': False, 'text.color': INK,
    'axes.labelcolor': INK, 'figure.facecolor': 'white', 'savefig.facecolor': 'white'})

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def label(ax, name):
    ax.text(0, 1, name, transform=ax.transAxes + ScaledTranslation(-18/72, 12/72, ax.figure.dpi_scale_trans),
            fontsize=10, weight='bold', va='bottom', ha='left')

def export(fig, number, axes, exemptions=()):
    qa = HERE / 'qa'; qa.mkdir(exist_ok=True)
    out = HERE / 'main' / f'Figure_{number:02d}'
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, axes=axes, exemptions=exemptions,
        json_out=qa / (out.name + '.alignment.json'), overlay_svg=qa / (out.name + '.alignment.svg'),
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=len(axes) > 1, strict=True)
    fig.savefig(out.with_suffix('.pdf'), dpi=600)
    fig.savefig(out.with_suffix('.svg'), dpi=600)
    fig.savefig(out.with_suffix('.png'), dpi=600)
    fig.savefig(out.with_suffix('.tiff'), dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)

def study_design(external):
    fig, ax = plt.subplots(figsize=(WIDTH, 165 / 25.4))
    fig.subplots_adjust(left=.035, right=.965, bottom=.035, top=.965)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    ax.text(.015, .985, 'Study design and evidence roles', fontsize=11, weight='bold', va='top')
    ax.text(.015, .934, 'Fixed candidates assessed across experimental and cellular contexts', fontsize=8, va='top')
    ax.text(.025, .855, 'Sesamin-related comparisons', fontsize=9, color=BLUE, weight='bold')
    ax.text(.54, .855, 'Context without sesamin', fontsize=9, color=ORANGE, weight='bold')
    left = [
        ('Human HepG2 cells', 'GSE248903 | 9 libraries', '25 fixed candidates; 15 primary orthologues', 'Culture independence and timing unresolved'),
        ('Mouse hepatocytes', 'PRJCA021826 / CRA013797 | 9 libraries', '64 library configurations per species', 'Treatment labels remain provisional'),
        ('Independent rat liver', 'GSE42693 | 5 control + 5 sesamin rats', 'Effect direction and small-sample sensitivity', 'Healthy liver; original body-weight selection'),
    ]
    right = [
        ('Primary human hepatocytes', 'GSE247407 | 30 libraries; 5 donor labels', 'Fatty-acid-associated expression background', 'Donor-label relationships partly unresolved'),
        ('Human liver organoids', 'GSE207889 | 24,034 retained cells', 'AREG localization and cell-mixture accounting', '4 technical libraries; one differentiation'),
    ]
    if external:
        right.append(('Independent adult human liver', 'GSE115469 | 8,444 cells; 5 donors',
                      'Donor-level cell-source comparison', 'Author-processed scale; sparse subsets flagged'))
    for x, cards, accent in [(.015, left, BLUE), (.525, right, ORANGE)]:
        for row, card in enumerate(cards):
            top = .811 - row * .202
            ax.add_patch(FancyBboxPatch((x, top - .169), .46, .169,
                boxstyle='round,pad=0.004,rounding_size=0.009', edgecolor='#CFD8DD', facecolor='#F7F9FA', lw=.7))
            ax.text(x + .014, top - .015, card[0], fontsize=8.3, weight='bold', color=accent, va='top')
            ax.text(x + .014, top - .053, card[1], fontsize=7.1, va='top')
            ax.text(x + .014, top - .090, card[2], fontsize=7.1, va='top')
            ax.text(x + .014, top - .128, card[3], fontsize=6.9, va='top', color='#53616A')
    ax.text(.015, .157, 'Interpretation', fontsize=8.6, weight='bold')
    ax.text(.015, .119, 'Within-study stability, cross-model agreement and cellular source answer different questions.', fontsize=7.5)
    ax.text(.015, .080, 'All fixed candidates and contrary results retained; expression units not pooled across species.', fontsize=7.1)
    ax.text(.015, .041, 'Study links indicate analytical roles, not causal mechanisms or demonstrated protection.', fontsize=7.1)
    export(fig, 1, [ax])
    (HERE / 'source_data/study_design.json').write_text(json.dumps({
        'intervention_related': left, 'contextual': right, 'external_atlas_included': external,
        'scientific_role': 'study schematic, not experimental measurements'}, indent=2), encoding='utf-8')

def stability():
    paths = {name: HERE / 'source_data' / (name + '.csv') for name in
             ['configuration_effect_summary', 'cross_species_configuration_consensus', 'frozen_primary_candidates']}
    inputs = {name: pd.read_csv(path) for name, path in paths.items()}
    summary = inputs['configuration_effect_summary']
    candidates = inputs['frozen_primary_candidates']
    genes = candidates.human_symbol.tolist()
    ranges = summary[(summary.species == 'mouse') & (summary.metric == 'TPM') &
                     (summary.aggregation == 'mean') & (summary.contrast == 'treatment_vs_model')].set_index('human_symbol').loc[genes]
    consensus = inputs['cross_species_configuration_consensus']
    consensus = consensus[consensus.aggregation == 'mean'].set_index('human_symbol').loc[genes]
    assert len(genes) == 15 and ranges.configurations.eq(64).all() and consensus.configuration_pairs.eq(4096).all()
    counts = []
    for species, metric in [('human', 'median_ratio_reestimated'), ('mouse', 'TPM')]:
        for contrast in ['model_vs_control', 'treatment_vs_model']:
            frame = summary[(summary.species == species) & (summary.metric == metric) &
                            (summary.aggregation == 'mean') & (summary.contrast == contrast)].set_index('human_symbol').loc[genes]
            assert len(frame) == 15 and frame.configurations.eq(64).all()
            counts.append(frame.baseline_direction_retained_n.to_numpy(dtype=int))
    counts = np.asarray(counts).T
    assert counts.shape == (15, 4) and (counts[:, :2] == 64).all()
    assert ((counts >= 0) & (counts <= 64)).all()
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 168 / 25.4),
        gridspec_kw={'width_ratios': [1.18, 1.20, .82]})
    fig.subplots_adjust(left=.14, right=.978, top=.827, bottom=.265, wspace=.40)
    fig.text(.035, .975, 'Within-study stability and cross-species agreement', fontsize=10, weight='bold', va='top')
    fig.text(.035, .934, 'All 15 primary orthologues | 64 overlapping configurations per species', fontsize=7.5, va='top')
    a, b, c = axes
    a.axvline(0, color='#BBC6CC', lw=.6, zorder=0)
    low = dict(zip(candidates.human_symbol, candidates.all_libraries_TPM_lt1))
    for i, (gene, row) in enumerate(ranges.iterrows()):
        color = BLUE if gene == 'VNN1' else ORANGE if gene in ['AREG', 'JAG2'] else GREY
        assert row.min_log2_ratio <= row.full_data_log2_ratio + 1e-12 <= row.max_log2_ratio + 1e-12
        a.hlines(i, row.min_log2_ratio, row.max_log2_ratio, color=color, lw=1.3)
        a.scatter(row.full_data_log2_ratio, i, s=14, color=color, zorder=3)
        n = int(consensus.loc[gene, 'both_model_and_treatment_concordant_n'])
        c.barh(i, n / 4096 * 100, height=.52, color=color)
        c.text(104, i, f'{n / 4096 * 100:.1f}'.rstrip('0').rstrip('.'), fontsize=6.4, va='center', ha='left')
    a.set_yticks(range(15), [g + (' *' if low[g] else '') for g in genes])
    a.set_xlim(float(ranges.min_log2_ratio.min()) - .13, float(ranges.max_log2_ratio.max()) + .13)
    a.set_title('Mouse effect range', pad=12)
    a.set_xlabel('sess / HFD log2 ratio\n0.1 TPM pseudocount', labelpad=7)
    b.imshow(counts, vmin=0, vmax=64, cmap='Blues', aspect='auto', interpolation='nearest', extent=(-.5, 3.5, 14.5, -.5))
    for row in range(15):
        for col in range(4):
            b.text(col, row, str(counts[row, col]), ha='center', va='center', fontsize=6.7,
                   color='white' if counts[row, col] >= 40 else INK)
    b.set_xticks(range(4), ['Human\nModel', 'Human\nTreat', 'Mouse\nModel', 'Mouse\nTreat'])
    b.tick_params(axis='x', labelsize=6.1, length=0, pad=5)
    b.set_title('Baseline sign retained', pad=12)
    b.set_xlabel('Configurations out of 64\nEach contrast within its species', labelpad=9)
    c.set_title('Both directions agree', pad=12)
    c.set_xlim(0, 133); c.set_xticks([0, 50, 100])
    c.set_xlabel('Human-mouse pairs (%)\nOut of 4096 pairs', labelpad=7)
    for ax, panel in zip(axes, 'abc'):
        ax.set_ylim(14.6, -.6); ax.tick_params(axis='y', length=0)
        label(ax, panel)
    for ax in [b, c]:
        ax.set_yticks([])
    for spine in b.spines.values():
        spine.set_visible(False)
    fig.text(.035, .143, 'Ranges are min-max sensitivity summaries, not confidence intervals. Model = model/control; Treat = treatment/model.', fontsize=6.9)
    fig.text(.035, .108, 'All human signs remain stable after prior selection. Overlapping configurations are not independent biological samples.', fontsize=6.9)
    fig.text(.035, .073, 'Cross-species agreement requires matching signs, whereas baseline retention only asks whether a sign changes within a species.', fontsize=6.9)
    fig.text(.035, .038, '* All nine mouse libraries below 1 TPM. Source treatment labels remain provisional. No hypothesis test is added.', fontsize=6.9)
    exemptions = [{'panels': list('abc'), 'checks': ['panel-width'],
                   'reason': 'Intentional range / four-column count matrix / agreement bar widths; all plot heights and shared row edges retained.'}]
    export(fig, 3, list(axes), exemptions)
    out = pd.DataFrame(counts, columns=['human_model_sign_n_of64', 'human_treatment_sign_n_of64',
                                     'mouse_model_sign_n_of64', 'mouse_treatment_sign_n_of64'])
    out.insert(0, 'human_symbol', genes)
    out['both_directions_agree_pairs_n_of4096'] = consensus.both_model_and_treatment_concordant_n.to_numpy()
    out.to_csv(HERE / 'source_data/figure3_species_sign_counts.csv', index=False)
    (HERE / 'qa/figure3_numerical_checks.json').write_text(json.dumps({
        'status': 'PASS', 'genes': genes, 'human_all15_all64_primary_signs_retained': True,
        'configurations_per_species': 64, 'pairs': 4096,
        'inputs': {name: {'path': str(path.relative_to(HERE)), 'sha256': sha(path)} for name, path in paths.items()},
        'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
        'matplotlib': matplotlib.__version__, 'analysis_scope': 'Reshaped existing results, no new statistical inference'}, indent=2), encoding='utf-8')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--external', action='store_true')
    parser.add_argument('--only', choices=['design', 'stability', 'all'], default='all')
    args = parser.parse_args()
    if args.only in ['design', 'all']:
        study_design(args.external)
    if args.only in ['stability', 'all']:
        stability()
    print('Rendered requested revised figures from frozen input tables.')
