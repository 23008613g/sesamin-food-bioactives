"""Plot actual submitted-library data; no significance or mechanism schematic."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
QA = Path('<LOCAL_PATH>')
sys.path.insert(0, str(QA))
from audit_panel_alignment import require_matplotlib_panel_alignment

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                     'font.size': 8, 'axes.titlesize': 10, 'axes.labelsize': 9,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
                     'svg.fonttype': 'none', 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.8, 'legend.frameon': False})
GROUPS = ['Control', 'OA', 'OA + sesamin']
COLORS = ['#666666', '#B56A44', '#387FA8']
MARKERS = ['o', '^', 's']


def export(fig, name):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=str(HERE / (name + '.alignment.json')),
                                       tolerance_pt=1.5, strict=True)
    fig.savefig(HERE / (name + '.pdf'))
    fig.savefig(HERE / (name + '.svg'))
    fig.savefig(HERE / (name + '.png'), dpi=300)
    plt.close(fig)


def main():
    summary = json.loads((HERE / 'analysis_summary.json').read_text())
    points = pd.read_csv(HERE / 'pca_scores.csv')
    fig, ax = plt.subplots(figsize=(7.2047244, 4.3307087))
    fig.subplots_adjust(left=0.13, right=0.96, bottom=0.17, top=0.76)
    for group, color, marker in zip(GROUPS, COLORS, MARKERS, strict=True):
        subset = points.loc[points['group'] == group]
        assert len(subset) == 3
        ax.scatter(subset['PC1'], subset['PC2'], color=color, marker=marker, s=48,
                   label=group, linewidths=0.6, edgecolors='white')
    ev = summary['pca_explained_variance']
    ax.set_xlabel(f'PC1 ({ev[0] * 100:.1f}% variance)')
    ax.set_ylabel(f'PC2 ({ev[1] * 100:.1f}% variance)')
    ax.margins(x=0.16, y=0.20)
    fig.suptitle('GSE248903: all nine deposited libraries', y=0.96, fontsize=10)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), ncol=3)
    export(fig, 'sample_pca')

    samples = pd.read_csv(HERE / 'known_panel_samples.csv')
    assert (samples['normalized_count'] >= 0).all()
    assert np.allclose(np.log2(samples['normalized_count'] + 1), samples['log2_normalized_plus1'])
    spec = json.loads((HERE / 'protocol.json').read_text(encoding='utf-8'))
    genes = spec['known_claim_panel']
    fig, ax = plt.subplots(figsize=(7.2047244, 4.3307087))
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.17, top=0.76)
    for j, (group, color, marker) in enumerate(zip(GROUPS, COLORS, MARKERS, strict=True)):
        x, y = [], []
        for i, gene in enumerate(genes):
            subset = samples.loc[(samples['symbol'] == gene) & (samples['group'] == group)].sort_values('sample')
            assert len(subset) == 3, 'Gene-symbol duplicates require explicit separate plotting'
            x.extend(i + (j - 1) * 0.23 + np.array([-0.045, 0.0, 0.045]))
            y.extend(subset['log2_normalized_plus1'])
        ax.scatter(x, y, color=color, marker=marker, s=18, alpha=0.85,
                   label=group, linewidths=0.35, edgecolors='white')
    ax.set_xticks(range(len(genes)), genes)
    ax.set_ylabel('log2(normalized count + 1)')
    ax.set_ylim(bottom=-0.6, top=float(samples['log2_normalized_plus1'].max()) + 1)
    ax.set_xlim(-0.7, len(genes) - 0.3)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05), ncol=3)
    fig.suptitle('Prior-claim genes: library-level expression', y=0.96, fontsize=10)
    export(fig, 'known_gene_panel')
    print('Exported two real-data figures: PDF, SVG, PNG; no inferential statistics.')


if __name__ == '__main__':
    main()
