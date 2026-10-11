"""Display the independently analyzed adult-liver atlas without adding inference."""
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from build_revision_figures import HERE, WIDTH, INK, BLUE, ORANGE, GREY, label, export, sha

ROOT = HERE.parent
SOURCE = ROOT / 'external_single_cell/results'

def main():
    names = ['figure_AREG_paired_donors.csv', 'figure_AREG_all_cell_context.csv', 'donor_broad_results_all25.csv']
    sources = []
    for name in names:
        target = HERE / 'source_data' / name
        shutil.copy2(SOURCE / name, target)
        sources.append({'path': 'source_data/' + name, 'sha256': sha(target)})
    paired = pd.read_csv(HERE / 'source_data/figure_AREG_paired_donors.csv')
    context = pd.read_csv(HERE / 'source_data/figure_AREG_all_cell_context.csv')
    details = pd.read_csv(HERE / 'source_data/donor_broad_results_all25.csv')
    # Zero-cell donor/type combinations remain NA in the full source table;
    # only donor/type combinations with observed cells define plotted means.
    details = details[(details.symbol == 'AREG') & details.available & (details.n_cells > 0)]
    assert len(paired) == 5 and paired.donor.nunique() == 5 and len(context) == 15
    assert paired.paired_available.all()
    assert (paired.delta_nonzero_percentage_points > 0).sum() == 4
    assert (paired.delta_mean_source_log2CPM > 0).sum() == 4
    assert np.allclose(paired.chol_nonzero_pct - paired.hep_nonzero_pct,
                       paired.delta_nonzero_percentage_points, rtol=0, atol=1e-9)
    assert np.allclose(paired.chol_mean_source_log2CPM - paired.hep_mean_source_log2CPM,
                       paired.delta_mean_source_log2CPM, rtol=0, atol=1e-9)
    # Independently reconstruct displayed equal-donor summaries from donor tables.
    for row in context.itertuples():
        values = details[details.cell_type == row.cell_type]
        assert len(values) == row.n_donors_with_available_expression
        assert np.isclose(values.pct_nonzero_source.mean(), row.equal_donor_nonzero_pct, atol=1e-8)
        assert np.isclose(values.mean_source_log2CPM.mean(), row.equal_donor_mean_source_log2CPM, atol=1e-10)
    fig = plt.figure(figsize=(WIDTH, 203 / 25.4))
    grid = fig.add_gridspec(2, 2, left=.325, right=.975, top=.844, bottom=.165,
                            height_ratios=[1, 2.35], hspace=.53, wspace=.43)
    a = fig.add_subplot(grid[0, 0]); b = fig.add_subplot(grid[0, 1]); c = fig.add_subplot(grid[1, :])
    fig.text(.035, .98, 'AREG cellular context in independent adult human liver', fontsize=10, weight='bold', va='top')
    fig.text(.035, .945, 'GSE115469 | Five donors | 8,444 author-retained cells | No sesamin exposure', fontsize=7.4, va='top')
    handles = [Line2D([0], [0], marker='o', color='none', markerfacecolor=color,
                markeredgecolor='none', markersize=5, label=text)
               for text, color in [('Hepatocytes', BLUE), ('Cholangiocytes', ORANGE)]]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.035, .919), ncol=2,
               fontsize=7.4, handletextpad=.4)
    for i, row in enumerate(paired.itertuples()):
        for ax, h_value, c_value in [(a, row.hep_nonzero_pct, row.chol_nonzero_pct),
                                    (b, row.hep_mean_source_log2CPM, row.chol_mean_source_log2CPM)]:
            ax.plot([h_value, c_value], [i, i], color='#B4BFC6', lw=1)
            ax.scatter(h_value, i, s=23, color=BLUE, edgecolor='none', zorder=3)
            ax.scatter(c_value, i, s=23, color=ORANGE, edgecolor='none', zorder=3)
    donor_labels = [f"{r.donor.replace('TLH', '')}{' *' if not r.both_types_ge20 else ''}  ({r.n_cholangiocytes}/{r.n_hepatocytes})"
                    for r in paired.itertuples()]
    a.set_yticks(range(5), donor_labels); a.set_ylabel('Donor (cholangiocytes/hepatocytes)', fontsize=7, labelpad=5)
    b.set_yticks(range(5), [r.donor.replace('TLH', '') for r in paired.itertuples()])
    for ax, title, xlabel in [(a, 'Nonzero source entries', 'Cells with source value > 0 (%)'),
                               (b, 'Mean source expression', 'Author-normalized log2CPM')]:
        ax.set_title(title, pad=12); ax.set_xlabel(xlabel, fontsize=7, labelpad=6)
        ax.set_ylim(4.65, -.65); ax.tick_params(axis='y', length=0, pad=4)
    a.set_xlim(-.8, 17); a.set_xticks([0, 5, 10, 15])
    b.set_xlim(-.008, .17); b.set_xticks([0, .05, .10, .15])
    type_labels = {
        'Central_venous_LSECs': 'Central venous LSECs',
        'Cholangiocytes': 'Cholangiocytes', 'Erythroid_Cells': 'Erythroid cells',
        'Hepatic_Stellate_Cells': 'Hepatic stellate cells', 'Hepatocytes': 'Hepatocytes',
        'Inflammatory_Macrophage': 'Inflammatory macrophages', 'Mature_B_Cells': 'Mature B cells',
        'NK-like_Cells': 'NK-like cells', 'Non-inflammatory_Macrophage': 'Non-inflammatory macrophages',
        'Periportal_LSECs': 'Periportal LSECs', 'Plasma_Cells': 'Plasma cells',
        'Portal_endothelial_Cells': 'Portal endothelial cells', 'alpha-beta_T_Cells': 'Alpha-beta T cells',
        'gamma-delta_T_Cells_1': 'Gamma-delta T cells 1', 'gamma-delta_T_Cells_2': 'Gamma-delta T cells 2'}
    for i, row in enumerate(context.itertuples()):
        donor_values = details[details.cell_type == row.cell_type].sort_values('donor')
        offsets = np.linspace(-.17, .17, len(donor_values))
        c.scatter(donor_values.pct_nonzero_source, i + offsets, s=12, color='#BBC3C8',
                  linewidths=0, zorder=2)
        color = BLUE if row.cell_type == 'Hepatocytes' else ORANGE if row.cell_type == 'Cholangiocytes' else INK
        c.scatter(row.equal_donor_nonzero_pct, i, marker='D', s=18, color=color, linewidths=0, zorder=3)
        c.text(94, i, f'{row.n_donors_with_cells}/5', fontsize=6.8, va='center')
    c.set_yticks(range(15), [type_labels[t] for t in context.cell_type])
    c.tick_params(axis='y', length=0, pad=5); c.set_ylim(14.6, -.6); c.set_xlim(-2, 103)
    c.set_xticks([0, 20, 40, 60, 80]); c.set_xlabel('Cells with nonzero source value (%)', labelpad=6)
    c.set_title('All source cell types: donor values and equal-donor mean', pad=12)
    c.text(94, -1.15, 'Donors', fontsize=6.8, ha='left')
    for ax, name in zip([a, b, c], 'abc'):
        label(ax, name)
    fig.text(.035, .101, 'Paired differences favor cholangiocytes in 4/5 donors; P1 is retained despite the opposite observed direction.', fontsize=7)
    fig.text(.035, .068, '* Fewer than 20 cholangiocytes. Grey points are donors; diamonds are equally weighted available-donor means.', fontsize=7)
    fig.text(.035, .035, 'Nonzero refers to the deposited normalized matrix, not verified raw-UMI detection. No inferential tests are shown.', fontsize=7)
    export(fig, 8, [a, b, c])
    (HERE / 'qa/figure8_numerical_checks.json').write_text(json.dumps({
        'status': 'PASS', 'donors': 5, 'cell_types': 15, 'paired_direction_positive': 4,
        'paired_direction_negative': 1, 'source_scale': 'deposited author-normalized log2CPM',
        'source_files': sources, 'all15_celltype_equal_donor_summaries_reconstructed': True,
        'hypothesis_tests_added': False, 'donors_excluded': []}, indent=2), encoding='utf-8')
    print('Figure 8 rendered and plotted donor/mean values numerically verified.')

if __name__ == '__main__':
    main()
