"""Copy verified plotting inputs and official QA helpers without source changes."""
from pathlib import Path
import hashlib
import json
import shutil

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
DATA = HERE / 'source_data'
DATA.mkdir(exist_ok=True)
for name in ['main', 'qa', 'scripts']:
    (HERE / name).mkdir(exist_ok=True)
specs = {
 'human_candidate_samples.csv': 'candidate_analysis/human_candidate_samples.csv',
 'mouse_candidate_samples.csv': 'candidate_analysis/mouse_candidate_samples.csv',
 'candidate_direction_order.csv': 'candidate_analysis/candidate_direction_order.csv',
 'rat_samples.csv': 'rat_analysis/results/candidate_sample_values_all25.csv',
 'rat_results.csv': 'rat_analysis/results/candidate_results_all25_welch_exact_two_normalizations.csv',
 'areg_subtype_by_library.csv': 'areg_analysis/areg_subtype_by_library.csv',
 'areg_excluding_smc.csv': 'areg_analysis/areg_excluding_smc.csv',
 'areg_matched_scale.csv': 'areg_analysis/areg_matched_scale.csv',
 'areg_population_by_library.csv': 'areg_analysis/areg_population_by_library.csv',
 'adult_donor_context.csv': 'areg_analysis/areg_adult_donor_context.csv',
 'adult_celltype_donor.csv': 'areg_analysis/areg_adult_all_celltypes_by_donor.csv',
 'adult_celltype_summary.csv': 'areg_analysis/areg_adult_immune_context.csv',
}
records = []
for output, source in specs.items():
    src = STAGE / source
    dst = DATA / output
    shutil.copyfile(src, dst)
    digest = hashlib.sha256(src.read_bytes()).hexdigest()
    assert digest == hashlib.sha256(dst.read_bytes()).hexdigest()
    records.append({'source_relative_to_stage': source, 'plot_input': output, 'sha256': digest, 'bytes': dst.stat().st_size})
helper_root = Path('C:/Users/AI/.agents/skills/nature-figure/scripts')
for name in ['validate_figure.py', 'audit_panel_alignment.py', 'audit_pdf_text.py', 'audit_figure_collisions.py']:
    shutil.copyfile(helper_root / name, HERE / 'scripts' / name)
(HERE / 'input_manifest.json').write_text(json.dumps(records, indent=2) + '\n', encoding='utf-8')
print('Copied', len(records), 'verified source tables; source files unchanged')
