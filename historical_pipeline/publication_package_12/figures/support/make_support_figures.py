"""Draw scientific figures from previously executed result tables only."""
from pathlib import Path
import hashlib
import json
import platform
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
NEXT = HERE.parents[2]
PROJECT = NEXT.parent
SKILL = Path.home() / '.agents/skills/nature-figure/scripts'
sys.path.insert(0, str(SKILL))
from audit_panel_alignment import require_matplotlib_panel_alignment

DATA = HERE / 'source_data'
DATA.mkdir(exist_ok=True)
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'], 'font.size': 7.5,
    'axes.titlesize': 8, 'axes.labelsize': 7.5,
    'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.linewidth': .65,
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'legend.frameon': False, 'savefig.facecolor': 'white',
})
BLUE, ORANGE, NEUTRAL = '#3E6D9C', '#C47738', '#6C7580'
GROUP = {'Control': NEUTRAL, 'OA': ORANGE, 'OA + sesamin': BLUE}
MARK = {'Control': 'o', 'OA': 's', 'OA + sesamin': '^'}
GENES = ['VNN1', 'AREG', 'JAG2', 'CXCL8']
INPUTS = []
WIDTH_IN = 180 / 25.4

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def register(path):
    path = Path(path)
    INPUTS.append({'path': str(path.relative_to(PROJECT)), 'sha256': sha(path), 'bytes': path.stat().st_size})
    return path

def read_json(path):
    return json.loads(register(path).read_text(encoding='utf-8'))

def read_csv(path):
    return pd.read_csv(register(path))

def panel_label(ax, letter):
    ax.annotate(letter, xy=(0, 1), xycoords='axes fraction',
                xytext=(-23, 10), textcoords='offset points',
                fontsize=9, fontweight='bold', va='bottom', ha='left')

def export(fig, stem, **alignment):
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig, json_out=HERE / f'{stem}.alignment.json',
        overlay_svg=HERE / f'{stem}.alignment.svg',
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True,
        **alignment)
    # Fixed page geometry preserves exactly 180 mm width.
    fig.savefig(HERE / f'{stem}.pdf')
    fig.savefig(HERE / f'{stem}.svg')
    fig.savefig(HERE / f'{stem}.png', dpi=600)
    fig.savefig(HERE / f'{stem}.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)

def figure1(summary):
    # One workflow canvas; box arrows are analytic flow, not molecular mechanisms.
    fig, ax = plt.subplots(figsize=(WIDTH_IN, 150/25.4))
    fig.subplots_adjust(left=.035, right=.965, bottom=.02, top=.975)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis('off')

    def box(x, y, w, h, title, text, edge=BLUE, fill='#F3F6FA'):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.006,rounding_size=0.012',
                                   facecolor=fill, edgecolor=edge, linewidth=.8))
        ax.text(x+.018, y+h-.028, title, ha='left', va='top', fontsize=8, fontweight='bold', color=edge)
        ax.text(x+.018, y+h-.066, text, ha='left', va='top', fontsize=7.4, linespacing=1.45)

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1,y1), (x2,y2), arrowstyle='-|>', mutation_scale=9,
                                     linewidth=.9, color=NEUTRAL))

    ax.text(.02, .985, 'Public transcriptomes: selection, comparison and context', fontsize=10,
            fontweight='bold', va='top')
    ax.text(.02, .93, 'SESAMIN-RELATED DATA', fontsize=7.5, color=NEUTRAL, fontweight='bold')
    box(.02, .731, .445, .172, 'HepG2 cells | GSE248903',
        '9 submitted libraries: 3 per group\nControl / OA / OA + sesamin\nCandidate selection', BLUE)
    box(.535, .731, .445, .172, 'Mouse hepatocytes | PRJCA021826',
        'CRA013797: 9 libraries; 3 per label\nctrl / HFD / sess; treatment mapping inferred\nCross-species comparison', ORANGE, '#FBF5EF')
    arrow(.24, .724, .24, .679)
    arrow(.755, .724, .755, .679)
    box(.02, .462, .445, .204, 'Descriptive candidate selection',
        f"{summary['gene_rows']:,} records → {summary['filtered_gene_rows']:,} filtered\n"
        f"{summary['base_descriptive_rule_n']} pass initial heuristic + CPM check\n"
        f"{summary['priority_rule_n']} pass leave-library-out stability rule\nNo differential-expression significance test", BLUE)
    box(.535, .462, .445, .204, 'Human–mouse orthology of 25 candidates',
        '15 one-to-one primary comparisons\nACAT2: 2 mappings, secondary only\n9 candidates not comparable (NA)\nNo substitute mouse gene for CXCL8', ORANGE, '#FBF5EF')
    arrow(.471, .563, .529, .563)
    ax.text(.02, .409, 'DISEASE-ONLY CONTEXT: NO SESAMIN INTERVENTION', fontsize=7.5,
            fontweight='bold', color=NEUTRAL)
    arrow(.24, .455, .24, .431)
    # Context boxes are linked by their shared role, not by a causal arrow.
    box(.02, .171, .445, .203, 'Primary human hepatocytes | GSE247407',
        '30 libraries / 5 donor labels\nFFA versus control within each label\nDonor-label heterogeneity\nLabel independence remains unresolved', NEUTRAL, '#F5F6F7')
    box(.535, .171, .445, .203, 'Liver organoids | GSE207889',
        '4 technical libraries / 1 differentiation\nPA versus control; single-cell data\nAuthor-defined hepatocyte-like subtypes\nSubtype and composition context', NEUTRAL, '#F5F6F7')
    ax.text(.02, .111, 'Output: traceable candidate evidence and testable follow-up questions',
            fontsize=8.4, fontweight='bold', color=BLUE)
    ax.text(.02, .063, 'Libraries are not assumed to be independent biological replicates.\n'
            'Contextual agreement does not establish sesamin efficacy or a protective mechanism.',
            fontsize=7.1, color=NEUTRAL, linespacing=1.5)
    export(fig, 'fig1_design')

def figure_s1(summary, pca, qc, mouse):
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, 87/25.4))
    fig.subplots_adjust(left=.09, right=.98, bottom=.20, top=.82, wspace=.40)
    a, b = axes
    for group in ['Control', 'OA', 'OA + sesamin']:
        rows = pca.loc[pca.group == group]
        a.scatter(rows.PC1, rows.PC2, c=GROUP[group], marker=MARK[group], s=30,
                  edgecolors='white', linewidths=.35, label=group, zorder=3)
    a.set_xlim(-23, 16)
    a.set_ylim(-10.5, 11)
    a.set_xlabel(f"PC1 ({100*summary['pca_explained_variance'][0]:.1f}%)")
    a.set_ylabel(f"PC2 ({100*summary['pca_explained_variance'][1]:.1f}%)")
    a.set_title('HepG2 library structure', pad=14)
    a.legend(loc='upper left', ncol=1, handletextpad=.4, borderaxespad=.2, fontsize=7)
    a.text(.01, -.29, '9 libraries; cached PCA coordinates', transform=a.transAxes, fontsize=7, color=NEUTRAL)
    y = np.arange(len(mouse))
    for group, color, marker in [('ctrl', NEUTRAL, 'o'), ('HFD', ORANGE, 's'), ('sess', BLUE, '^')]:
        idx = mouse['original_group'] == group
        b.scatter(mouse.loc[idx, 'pseudoalignment_percent'], y[idx], c=color, marker=marker,
                  s=25, edgecolors='white', linewidths=.35, zorder=3)
    b.set_yticks(y, mouse.original_label)
    b.invert_yaxis()
    b.set_ylim(8.6, -.6)
    b.set_xlim(80, 91)
    b.set_xticks([80, 85, 90])
    b.set_xlabel('Pseudoaligned read pairs (%)')
    b.set_title('Mouse library pseudoalignment', pad=14)
    for i, row in mouse.iterrows():
        b.annotate(f'{row.pseudoalignment_percent:.2f}', (row.pseudoalignment_percent, i),
                   xytext=(5,0), textcoords='offset points', va='center', fontsize=7)
    b.text(.01, -.29, '9 libraries; original sample labels', transform=b.transAxes, fontsize=7, color=NEUTRAL)
    for letter, ax in zip('ab', axes):
        panel_label(ax, letter)
    export(fig, 'figS1_quality', panel_ids=['a','b'], require_panel_labels=True)

def figure_s2(pseudo, omission):
    fig, axes = plt.subplots(2, 4, figsize=(WIDTH_IN, 150/25.4), sharey='row')
    fig.subplots_adjust(left=.09, right=.985, bottom=.16, top=.865, wspace=.33, hspace=.93)
    labels = [('OA_vs_control','OA / Control', ORANGE, 'o'),
              ('sesamin_vs_OA','OA + Sesamin / OA', BLUE, 's'),
              ('sesamin_vs_control','OA + Sesamin / Control', NEUTRAL, '^')]
    for j, gene in enumerate(GENES):
        a, b = axes[:,j]
        rows = pseudo.loc[pseudo.symbol == gene].sort_values('pseudocount')
        assert len(rows) == 4 and set(rows.pseudocount) == {0, .5, 1, 5}
        for field, label, color, marker in labels:
            a.plot(range(4), rows[field], color=color, marker=marker, ms=3.7, lw=.9, label=label)
        a.axhline(0, lw=.5, ls='--', color='#A5ADB5', zorder=0)
        a.set_xticks(range(4), ['0', '0.5', '1', '5'])
        a.set_ylim(-1.5, 2.1)
        a.set_xlabel('Pseudocount p')
        a.set_title(gene, fontstyle='italic', pad=12)
        donor_order = ['F95','M95','H20','H961','F95P']
        rows = omission.loc[omission.symbol == gene].set_index('omitted_label').loc[donor_order]
        assert len(rows) == 5
        b.scatter(range(5), rows.median_effect_other4, c=BLUE, s=22, edgecolors='white', linewidths=.35, zorder=3)
        b.axhline(0, lw=.5, ls='--', color='#A5ADB5', zorder=0)
        b.set_xticks(range(5), donor_order, rotation=45, ha='right', rotation_mode='anchor')
        b.set_xlim(-.5, 4.5)
        b.set_ylim(-.65, .55)
        b.set_xlabel('Omitted donor label')
        b.set_title(gene, fontstyle='italic', pad=12)
    axes[0,0].set_ylabel('Descriptive log2 ratio')
    axes[1,0].set_ylabel('Median log2 ratio, other 4 labels')
    for letter, ax in zip('abcdefgh', axes.flat):
        panel_label(ax, letter)
    fig.text(.09, .965, 'HepG2: sensitivity to the additive pseudocount', fontsize=8.5, fontweight='bold')
    handles, texts = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, texts, loc='upper center', bbox_to_anchor=(.54,.944), ncol=3,
               handlelength=1.8, columnspacing=1.7, borderaxespad=0, fontsize=7)
    fig.text(.09, .49, 'PHH: leave-one-label-out FFA / Control summaries (no sesamin)', fontsize=8.5, fontweight='bold')
    fig.text(.09, .035, 'All settings shown; omission subsets overlap and are not independent replicates or confidence intervals.',
             fontsize=7, color=NEUTRAL)
    export(fig, 'figS2_sensitivity', panel_ids=list('abcdefgh'), require_panel_labels=True)

def main():
    summary = read_json(NEXT / 'sesamin_exploration_04/analysis_summary.json')
    read_json(NEXT / 'sesamin_exploration_04/protocol.json')
    pca = read_csv(NEXT / 'sesamin_exploration_04/pca_scores.csv')
    qc = read_csv(NEXT / 'sesamin_exploration_04/sample_qc.csv')
    assert len(pca) == len(qc) == 9
    assert list(pca['sample']) == list(qc['sample'])
    register(NEXT / 'manuscript_draft_11/manuscript_zh.md')
    register(NEXT / 'mouse_quantification_08/results_summary.json')
    rows = []
    for i, run in enumerate(range(1010532,1010541)):
        r = read_json(NEXT / f'mouse_quantification_08/run_results/DRR{run}.json')
        frac = r['pseudoalignment_fraction']
        assert abs(frac - r['run_info']['n_pseudoaligned']/r['run_info']['n_processed']) < 1e-12
        rows.append({'run':r['run'], 'original_label': ['ctrl','HFD','sess'][i//3]+str(i%3+1),
                     'original_group':['ctrl','HFD','sess'][i//3],
                     'n_processed_pairs':r['run_info']['n_processed'],
                     'n_pseudoaligned_pairs':r['run_info']['n_pseudoaligned'],
                     'pseudoalignment_percent': frac*100})
    mouse = pd.DataFrame(rows)
    pca.merge(qc, on=['sample','group'], validate='one_to_one').to_csv(DATA/'figS1_human_library_qc.csv',index=False)
    mouse.to_csv(DATA/'figS1_mouse_library_qc.csv',index=False)
    pseudo_all = read_csv(NEXT/'followup_05/bulk_pseudocount_sensitivity.csv')
    omission_all = read_csv(NEXT/'followup_05/phh_label_omission.csv')
    pseudo = pseudo_all.loc[pseudo_all.symbol.isin(GENES)].copy()
    omission = omission_all.loc[omission_all.symbol.isin(GENES)].copy()
    assert len(pseudo) == 16 and len(omission) == 20
    pseudo.to_csv(DATA/'figS2_pseudocounts.csv',index=False)
    omission.to_csv(DATA/'figS2_label_omission.csv',index=False)
    pd.DataFrame([
        ['GSE248903','HepG2','9 submitted libraries','Candidate selection','Sesamin-related treatment labels; independent cultures unresolved'],
        ['PRJCA021826 / CRA013797','Mouse hepatocytes','9 submitted libraries','Cross-species comparison','Treatment mapping inferred; independent biological units unresolved'],
        ['GSE247407','Primary human hepatocytes','30 libraries / 5 donor labels','Disease-only context','No sesamin; donor-label independence unresolved'],
        ['GSE207889','Liver organoid single-cell','4 technical libraries / 1 differentiation','Disease-only context','No sesamin; author-defined hepatocyte-like subtypes']
    ], columns=['accession','model','unit','role','boundary']).to_csv(DATA/'fig1_dataset_roles.csv',index=False)
    pd.DataFrame([['All gene records',summary['gene_rows']],['Expression-filtered',summary['filtered_gene_rows']],
                  ['Initial heuristic and CPM agreement',summary['base_and_CPM_agreement_n']],
                  ['Stability-rule candidates',summary['priority_rule_n']],['One-to-one primary',15],
                  ['One-to-many candidate genes (ACAT2)',1],['Secondary mouse mappings for ACAT2',2],
                  ['Not comparable',9]],columns=['stage','count']).to_csv(DATA/'fig1_counts.csv',index=False)
    figure1(summary)
    figure_s1(summary,pca,qc,mouse)
    figure_s2(pseudo,omission)
    record = {
        'date':'2026-10-09', 'interpreter':sys.executable, 'python':platform.python_version(),
        'matplotlib':matplotlib.__version__, 'numpy':np.__version__, 'pandas':pd.__version__,
        'input_fingerprints':INPUTS, 'plot_script_sha256':sha(Path(__file__)),
        'randomness':'None', 'biological_analysis_recomputed':False, 'libraries_excluded':0,
        'illustrative_genes':GENES, 'source_pseudocount_rows':len(pseudo_all),'display_pseudocount_rows':len(pseudo),
        'source_PHh_omission_rows':len(omission_all),'display_PHh_omission_rows':len(omission),
        'subset_reason':'The four illustrative genes were assigned before drawing; all their recorded settings and all three HepG2 contrasts are shown.',
        'pca_explained_variance':summary['pca_explained_variance'][:2],
        'mouse_pseudoalignment_percent_range':[float(mouse.pseudoalignment_percent.min()),float(mouse.pseudoalignment_percent.max())],
        'source_data_fingerprints':[{ 'file':str(p.relative_to(HERE)), 'sha256':sha(p)} for p in sorted(DATA.glob('*.csv'))]
    }
    (HERE/'execution_record.json').write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'figures':3,'libraries':{'human':9,'mouse':9},'sensitivity_rows':[len(pseudo),len(omission)]}))

if __name__ == '__main__':
    main()
