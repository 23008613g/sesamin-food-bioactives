"""Programmatic Fig1 workflow and source-data Fig3 revision, Python only.

First run: pass --archive and --genome-stage to populate verified source copies.
Later runs reuse review/figure_updates/source_data without the old project.
No pixel editing, image-generation service, new hypothesis test or data exclusion.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from zipfile import ZipFile

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures_new'
QA = ROOT / 'review/figure_updates'
DATA = QA / 'source_data'
WIDTH_MM = 174
INK, BLUE, ORANGE, GREY = '#223744', '#39748F', '#B96F3E', '#7D8B92'
CONTEXT_TEXT, MOUSE_TEXT = '#995325', '#536773'
AUDITS = []
CHECKS = []
sys.dont_write_bytecode = True
mpl.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 9,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'axes.linewidth': .65, 'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'legend.frameon': False, 'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': GREY,
})


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def check(name, condition):
    CHECKS.append({'name': name, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(name)


def copy_verified_bytes(data, destination):
    digest = hashlib.sha256(data).hexdigest()
    if destination.exists():
        if sha(destination) != digest:
            raise ValueError(f'Existing frozen source differs: {destination.name}')
    else:
        destination.write_bytes(data)
    return digest


def prepare_sources(archive, genome_stage):
    inputs = {}
    manifest_path = QA / 'source_manifest.json'
    if archive is not None:
        archive_before = sha(archive)
        names = ['human_candidate_samples.csv', 'mouse_candidate_samples.csv',
                 'rat_samples.csv', 'Fig3_all_observations.csv']
        with ZipFile(archive) as z:
            for name in names:
                member = f'figures/source_data/{name}'
                data = z.read(member)
                inputs[name] = {'source_kind': 'archived submitted figure input', 'archive_member': member,
                                'sha256': copy_verified_bytes(data, DATA / name), 'bytes': len(data)}
            data = z.read('figures/build_figures.py')
            inputs['original_figure_code.txt'] = {'source_kind': 'read-only plotting provenance',
                'archive_member': 'figures/build_figures.py',
                'sha256': copy_verified_bytes(data, DATA / 'original_figure_code.txt'), 'bytes': len(data)}
        check('source_archive_unchanged', sha(archive) == archive_before)
        inputs['archive'] = {'sha256': archive_before, 'bytes': archive.stat().st_size}
    if genome_stage is not None:
        for src, name in [(genome_stage / 'analysis_summary.json', 'genome_analysis_summary.json'),
                          (genome_stage / 'review/cross_species_execution.json', 'cross_species_execution.json'),
                          (genome_stage / 'results/rat_gene_expression_deposited.csv.gz', 'rat_gene_expression_deposited.csv.gz')]:
            before = sha(src)
            data = src.read_bytes()
            inputs[name] = {'source_kind': 'executed genome-wide results',
                            'source_relative_path': str(src.relative_to(genome_stage)).replace('\\', '/'),
                            'sha256': copy_verified_bytes(data, DATA / name), 'bytes': len(data)}
            check(f'{name}_source_unchanged', sha(src) == before)
    if inputs:
        if manifest_path.exists():
            old = json.loads(manifest_path.read_text(encoding='utf-8'))
            old['inputs'].update(inputs)
            inputs = old['inputs']
        dump(manifest_path, {'policy': 'Source originals are read-only; copies preserve exact bytes.', 'inputs': inputs})
    if not manifest_path.is_file():
        raise FileNotFoundError('First run requires --archive and --genome-stage')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    for name, item in manifest['inputs'].items():
        if name != 'archive':
            check(f'{name}_frozen_hash', sha(DATA / name) == item['sha256'])
    return manifest


def export(fig, stem, axes, qa_tools, audit_python, require_labels):
    sys.path.insert(0, str(qa_tools))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    fig.canvas.draw()
    alignment = require_matplotlib_panel_alignment(
        fig, axes=axes, panel_ids=list('abcdefghijklmnopqrstuvwxyz')[:len(axes)],
        json_out=QA / f'{stem}.alignment.json', overlay_svg=QA / f'{stem}.alignment.svg',
        tolerance_pt=1.5, gutter_tolerance_pt=1.5,
        require_panel_labels=require_labels, strict=True)
    # Preserve explicit physical dimensions; do not crop the canvas.
    fig.savefig(OUT / f'{stem}.pdf')
    fig.savefig(OUT / f'{stem}.svg')
    fig.savefig(OUT / f'{stem}.eps')
    fig.savefig(OUT / f'{stem}.png', dpi=300)
    with Image.open(OUT / f'{stem}.png') as im:
        check(stem + '_300dpi_preview', min(im.info['dpi']) > 299)
        check(stem + '_physical_width', abs(im.width / 300 * 25.4 - WIDTH_MM) < .1)
    width, height = fig.get_size_inches() * 25.4
    plt.close(fig)
    text = subprocess.run([str(audit_python), str(qa_tools / 'audit_pdf_text.py'),
                           str(OUT / f'{stem}.pdf'), '--min-pt', '5', '--json'],
                          capture_output=True, text=True, encoding='utf-8')
    (QA / f'{stem}.text_audit.json').write_text(text.stdout, encoding='utf-8')
    collision = subprocess.run([str(audit_python), str(qa_tools / 'audit_figure_collisions.py'),
                                str(OUT / f'{stem}.pdf'), '--json-out', str(QA / f'{stem}.collision.json'),
                                '--overlay-pdf', str(QA / f'{stem}.collision_overlay.pdf')],
                               capture_output=True, text=True, encoding='utf-8')
    (QA / f'{stem}.collision.log').write_text(collision.stdout + collision.stderr, encoding='utf-8')
    AUDITS.append({'figure': stem, 'width_mm': float(width), 'height_mm': float(height),
                   'alignment': alignment.get('verdict'), 'pdf_text_exit': text.returncode,
                   'collision_exit': collision.returncode})
    print(json.dumps(AUDITS[-1]), flush=True)


def fig1(qa_tools, audit_python):
    genome = json.loads((DATA / 'genome_analysis_summary.json').read_text(encoding='utf-8'))
    cross = json.loads((DATA / 'cross_species_execution.json').read_text(encoding='utf-8'))
    h = pd.read_csv(DATA / 'human_candidate_samples.csv')
    check('workflow_candidates_25', h.symbol.nunique() == 25)
    check('workflow_rat_genes_14604', genome['n_genes'] == 14604)
    check('workflow_GO_sets_1311', genome['go_results']['deposited']['n_gene_sets'] == 1311)
    check('workflow_shared_genes_7912', cross['counts']['shared_genes'] == 7912)
    def luminance(colour):
        c = [int(colour[i:i+2], 16) / 255 for i in (1, 3, 5)]
        c = [v / 12.92 if v <= .04045 else ((v + .055) / 1.055) ** 2.4 for v in c]
        return .2126*c[0] + .7152*c[1] + .0722*c[2]
    contrasts = {colour: (luminance('#F7F9FA')+.05)/(luminance(colour)+.05)
                 for colour in (INK, BLUE, CONTEXT_TEXT, MOUSE_TEXT, '#56666F')}
    check('Fig1_text_contrast_at_least_4_5', min(contrasts.values()) >= 4.5)
    dump(QA / 'Fig1.text_contrast.json', {'background': '#F7F9FA', 'ratios': contrasts})
    fig, ax = plt.subplots(figsize=(WIDTH_MM / 25.4, 119 / 25.4))
    fig.subplots_adjust(left=.025, right=.975, bottom=.025, top=.975)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis('off')
    geometry = []

    def card(x, y, w, hgt, title, lines, tone=BLUE, group=''):
        ax.add_patch(FancyBboxPatch((x, y), w, hgt, boxstyle='round,pad=.003,rounding_size=.007',
                                    facecolor='#F7F9FA', edgecolor='#C7D2D8', linewidth=.65))
        # Text resides wholly within its vector card; no white masks or image edits.
        ax.text(x+.014, y+hgt-.023, title, fontsize=8.6, color=tone, weight='bold', va='top')
        for i, line in enumerate(lines):
            ax.text(x+.014, y+hgt-.074-i*.038, line, fontsize=7.7, color=INK, va='top')
        geometry.append({'group': group, 'title': title, 'x': x, 'y': y, 'width': w, 'height': hgt})

    ax.text(.008, .974, 'Candidate-focused comparisons', fontsize=9.2, weight='bold', color=BLUE, va='top')
    card(.008, .736, .469, .186, 'HepG2 cells',
         ['GSE248903 | 9 libraries', '25 fixed candidates', 'Tumour-derived; culture units unresolved'], group='candidate')
    card(.523, .736, .469, .186, 'Mouse hepatocyte libraries',
         ['CRA013797 | 9 libraries', 'ctrl / HFD / sess source labels', 'Exposure mapping unconfirmed'], tone=MOUSE_TEXT, group='candidate')

    ax.text(.008, .677, 'Candidate-independent assessment', fontsize=9.2, weight='bold', color=BLUE, va='top')
    card(.008, .439, .469, .186, 'Rat liver after dietary sesamin',
         ['GSE42693 | 5 control + 5 sesamin', '14,604 genes | 1,311 GO sets', 'All eligible genes; 15-day exposure'], group='whole_genome')
    card(.523, .439, .469, .186, 'Human-rat comparison',
         ['Full retained HepG2 transcriptome + rat liver', '7,912 one-to-one shared genes', 'No candidate filter'], group='whole_genome')

    ax.text(.008, .380, 'Cellular context without sesamin', fontsize=9.2, weight='bold', color=CONTEXT_TEXT, va='top')
    card(.008, .094, .309, .234, 'Primary human hepatocytes',
         ['GSE247407 | 30 libraries', '5 donor labels', 'Fatty-acid context', 'Donor-label relation unresolved'], tone=CONTEXT_TEXT, group='context')
    card(.346, .094, .309, .234, 'Human liver organoids',
         ['GSE207889 | 24,034 cells', '4 technical libraries', '1 differentiation', 'AREG cellular context'], tone=CONTEXT_TEXT, group='context')
    card(.683, .094, .309, .234, 'Adult human liver atlas',
         ['GSE115469 | 8,444 cells', '5 donors', 'Donor and immune-cell context', 'AREG cellular context'], tone=CONTEXT_TEXT, group='context')
    ax.text(.008, .032, 'Study-specific expression evidence; no normal-human sesamin validation',
            fontsize=7.7, color='#56666F', va='top')
    # Fig1 is a single diagram axis, not comparable quantitative panels.
    # Explicitly verify peer cards within each band have common sizes and edges.
    for group in ('candidate', 'whole_genome', 'context'):
        same = [item for item in geometry if item['group'] == group]
        check(f'Fig1_{group}_card_alignment', all(same[0][key] == item[key]
              for item in same for key in ('y', 'height', 'width')))
    dump(QA / 'Fig1.card_geometry.json', {'role': 'single-axis vector schematic', 'cards': geometry,
          'no_causal_arrows': True, 'source_dataset_count': 6,
          'seventh_box_is_analysis_not_dataset': True})
    dump(QA / 'Fig1.contract.json', {
        'core_conclusion': 'Six public datasets contribute distinct candidate, candidate-independent and cellular-context evidence; none supplies normal-human sesamin validation.',
        'archetype': 'schematic-led workflow', 'backend': 'Python/matplotlib vector primitives',
        'source_counts': {'HepG2_candidates': 25, 'rat_genes': 14604, 'rat_GO_sets': 1311, 'shared_genes': 7912},
        'no_sesamin_context': ['GSE247407', 'GSE207889', 'GSE115469'],
        'new_dataset_or_mechanism': False, 'generative_image_service': False,
        'human_author_approval': 'Not established by automated rendering or internal agent QA'})
    export(fig, 'Fig1', [ax], qa_tools, audit_python, require_labels=False)


def letter(ax, char):
    ax.text(0, 1, char, transform=ax.transAxes + ScaledTranslation(-22/72, 12/72, ax.figure.dpi_scale_trans),
            fontsize=10, fontweight='bold', ha='left', va='bottom')


def fig3(qa_tools, audit_python):
    h = pd.read_csv(DATA / 'human_candidate_samples.csv').rename(columns={'library': 'sample', 'median_ratio_count': 'normalized_count'}).query("symbol == 'VNN1'")
    m = pd.read_csv(DATA / 'mouse_candidate_samples.csv').rename(columns={'symbol': 'human_symbol', 'library': 'run'}).query("human_symbol == 'VNN1'").copy()
    r = pd.read_csv(DATA / 'rat_samples.csv').query("human_symbol == 'VNN1' and method == 'deposited_log2' and primary_comparison")
    check('Fig3_all_28_observations', len(h) == 9 and len(m) == 9 and len(r) == 10)
    genome = pd.read_csv(DATA / 'rat_gene_expression_deposited.csv.gz', dtype={'rat_entrez': str}).set_index('rat_entrez')
    sample_order = [f'GSM{i}' for i in range(1048036, 1048046)]
    rat_current = genome.loc['29142', sample_order].to_numpy(dtype=float)
    rat_old = r.set_index('sample').loc[sample_order, 'aggregate_log2_signal'].to_numpy(dtype=float)
    check('Fig3_rat_values_match_full_genome', np.allclose(rat_current, rat_old, atol=1e-12, rtol=0))
    specs = [(h, 'normalized_count', ['Control', 'Model', 'Treat'], ['Control', 'Model', 'Treat'], 'Human HepG2', 'Normalized counts', 'sample'),
             (m, 'TPM', ['ctrl', 'HFD', 'sess'], ['ctrl', 'HFD', 'sess'], 'Mouse source labels', 'TPM', 'run'),
             (r, 'aggregate_log2_signal', ['control', 'sesamin'], ['Control', 'Sesamin'], 'Dietary rat liver', 'log2(signal + 1)', 'sample')]
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH_MM/25.4, 89/25.4))
    fig.subplots_adjust(left=.105, right=.985, bottom=.32, top=.86, wspace=.68)
    display = []
    for ax, (data, field, order, ticks, head, unit, idcol), panel in zip(axes, specs, 'abc'):
        for j, group in enumerate(order):
            subset = data.loc[data.group == group].sort_values(idcol)
            check(f'Fig3_{panel}_{group}_n', len(subset) == (5 if panel == 'c' else 3))
            check(f'Fig3_{panel}_{group}_finite', np.isfinite(subset[field]).all())
            colour = [GREY, ORANGE, BLUE][j] if len(order) == 3 else [GREY, BLUE][j]
            ax.scatter(j+np.linspace(-.12, .12, len(subset)), subset[field], s=25,
                       c=colour, edgecolors='white', linewidths=.35, zorder=3)
            ax.hlines(subset[field].mean(), j-.24, j+.24, color=INK, linewidth=1.1, zorder=4)
            display.extend({'panel': panel, 'observation_id': row[idcol], 'group': group,
                            'value': row[field], 'unit': unit} for row in subset.to_dict('records'))
        ax.set_xlim(-.5, len(order)-.5)
        ax.set_xticks(range(len(order)), ticks)
        ax.set_title(head, pad=13)
        ax.set_ylabel(unit, labelpad=5)
        letter(ax, panel)
        ax.tick_params(axis='x', length=0, pad=5)
    axes[0].set_ylim(0, max(h.normalized_count)*1.15)
    axes[1].set_ylim(0, max(m.TPM)*1.15)
    # Rat autoscaling is identical to the source plotting code; no axis transform.
    fig.text(.105, .153, 'All 28 observations retained; short lines show group means', fontsize=8)
    fig.text(.105, .097, 'Rat: five animals per group; cell datasets: three libraries per label', fontsize=8)
    fig.text(.105, .041, 'Cell-library independence and mouse exposure mapping remain unresolved', fontsize=8)
    shown = pd.DataFrame(display)
    original = pd.read_csv(DATA / 'Fig3_all_observations.csv')
    check('Fig3_identifiers_groups_units_unchanged', shown.drop(columns='value').equals(original.drop(columns='value')))
    check('Fig3_values_unchanged', np.allclose(shown.value, original.value, rtol=0, atol=1e-12))
    shown.to_csv(QA / 'Fig3_display_values.csv', index=False)
    dump(QA / 'Fig3.contract.json', {
        'core_conclusion': 'The pre-existing VNN1 example shows study-specific expression patterns across unlike systems.',
        'archetype': 'quantitative grid', 'reuse_level': 'source-data redraw retaining original scales and geometry',
        'panels': [{'panel': ch, 'n_observations': len(spec[0]), 'unit': spec[5],
                    'ylim': list(ax.get_ylim()), 'center': 'arithmetic group mean',
                    'uncertainty': 'All individual observations shown; no CI or SD bar added',
                    'biological_unit': 'individual rat' if ch == 'c' else 'library; culture-unit independence unresolved'}
                   for ch, spec, ax in zip('abc', specs, axes)],
        'removed': 'Old 25-candidate-family BH annotation only', 'new_tests': False,
        'new_rat_comparison_max_error': float(np.max(np.abs(rat_current-rat_old))),
        'display_max_error_vs_archive': float(np.max(np.abs(shown.value-original.value))),
        'HepG2_scope_for_caption': 'Tumour-derived model, not normal hepatocytes',
        'mouse_scope_for_caption': 'Source labels only; treatment mapping unconfirmed',
        'rat_inference_for_caption': 'Report full-genome results elsewhere, with method dependence; these plotted values do not change.'})
    export(fig, 'Fig3', list(axes), qa_tools, audit_python, require_labels=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--genome-stage', type=Path)
    parser.add_argument('--qa-tools', type=Path, required=True)
    parser.add_argument('--audit-python', type=Path, required=True)
    args = parser.parse_args()
    for path in (OUT, QA, DATA):
        path.mkdir(parents=True, exist_ok=True)
    for name in ('audit_panel_alignment.py', 'audit_pdf_text.py', 'audit_figure_collisions.py', 'validate_figure.py'):
        if not (args.qa_tools / name).is_file():
            raise FileNotFoundError(f'Missing existing QA tool: {name}')
    runtime = subprocess.run([str(args.audit_python), '-c', 'import fitz,sys,json; print(json.dumps({"python":sys.executable,"pymupdf":fitz.VersionBind}))'],
                             capture_output=True, text=True, encoding='utf-8')
    if runtime.returncode:
        raise RuntimeError('Existing PDF audit interpreter failed: ' + runtime.stderr)
    manifest = prepare_sources(args.archive, args.genome_stage)
    static = subprocess.run([sys.executable, str(args.qa_tools/'validate_figure.py'), str(Path(__file__)), '--json'],
                            capture_output=True, text=True, encoding='utf-8')
    (QA / 'source_static_audit.json').write_text(static.stdout, encoding='utf-8')
    fig1(args.qa_tools, args.audit_python)
    fig3(args.qa_tools, args.audit_python)
    result = {'status': 'rendered; inspect actual PNG and review any PDF warnings before delivery',
              'completed_utc': datetime.now(timezone.utc).isoformat(), 'command': sys.argv,
              'python': sys.executable, 'python_version': platform.python_version(),
              'matplotlib': mpl.__version__, 'numpy': np.__version__, 'pandas': pd.__version__,
              'audit_runtime': json.loads(runtime.stdout), 'script_sha256': sha(Path(__file__)),
              'static_audit_exit': static.returncode, 'checks': CHECKS, 'audits': AUDITS,
              'outputs': [{'file': path.name, 'sha256': sha(path), 'bytes': path.stat().st_size}
                          for path in sorted(OUT.iterdir()) if path.name.startswith(('Fig1.', 'Fig3.'))]}
    dump(QA / 'execution.json', result)
    if static.returncode or any(a['pdf_text_exit'] or a['collision_exit'] for a in AUDITS):
        raise RuntimeError('Figure audit requires correction or documented review; see review/figure_updates')


if __name__ == '__main__':
    main()
