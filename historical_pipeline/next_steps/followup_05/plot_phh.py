"""Show both spotlight genes and every PHH donor-label, including discordance."""
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, '<LOCAL_PATH>')
from audit_panel_alignment import require_matplotlib_panel_alignment

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
                     'font.size': 8, 'axes.labelsize': 9, 'xtick.labelsize': 8,
                     'ytick.labelsize': 8, 'svg.fonttype': 'none', 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.8})


def main():
    frame = pd.read_csv(HERE / 'phh_donor_effects.csv')
    assert (frame[['Control_mean_norm', 'FFA_mean_norm']] >= 0).all().all()
    expected = np.log2((frame['FFA_mean_norm'] + 1) / (frame['Control_mean_norm'] + 1))
    assert np.allclose(expected, frame['log2_FFA_Control_mean_norm_plus1'])
    donors = ['F95', 'M95', 'H20', 'H961', 'F95P']
    fig, ax = plt.subplots(figsize=(7.2047244, 5.1181102))
    fig.subplots_adjust(left=0.23, right=0.96, bottom=0.23, top=0.82)
    positions, labels = [], []
    for i, (gene, color, marker) in enumerate([('CXCL8', '#387FA8', 'o'), ('VNN1', '#94628D', 's')]):
        subset = frame.loc[frame['symbol'] == gene].set_index('donor_label').loc[donors]
        y = np.arange(5) + 6 * i
        ax.scatter(subset['log2_FFA_Control_mean_norm_plus1'], y, s=38, color=color, marker=marker,
                   zorder=3, linewidths=0.5, edgecolors='white')
        positions.extend(y)
        labels.extend([gene + ' / ' + donor for donor in donors])
    ax.axvline(0, color='#999999', linewidth=0.8, zorder=0)
    ax.set_yticks(positions, labels)
    ax.set_ylim(10.7, -0.7)
    ax.set_xlim(-1.0, 2.2)
    ax.set_xlabel('Within-label FFA / Control log2 ratio (normalized means + 1)')
    fig.suptitle('PHH fatty-acid responses differ across donor labels', y=0.96, fontsize=10)
    fig.text(0.5, 0.89, 'Every label shown; this study contains no sesamin treatment', ha='center', fontsize=8)
    fig.text(0.5, 0.06, 'Three libraries per condition averaged within each label; no significance testing.', ha='center', fontsize=7)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=str(HERE / 'phh_spotlight.alignment.json'),
                                       tolerance_pt=1.5, strict=True)
    fig.savefig(HERE / 'phh_spotlight.pdf')
    fig.savefig(HERE / 'phh_spotlight.svg')
    fig.savefig(HERE / 'phh_spotlight.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
