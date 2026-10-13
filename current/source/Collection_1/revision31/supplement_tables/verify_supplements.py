"""Bounded content/source checks; Word/PDF visual QA belongs to document build."""
from pathlib import Path
import hashlib
import json
import re
import pandas as pd

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
OLD = STAGE.parent / 'pfhn_submission_28'
checks = {}


def verify(name, outcome):
    checks[name] = bool(outcome)
    if not outcome:
        raise AssertionError(name)


def rows(text, number, lang):
    prefix = '## Supplementary Table S' if lang == 'en' else '## 补充表S'
    heading = re.search('^' + re.escape(prefix + str(number) + '.') + r'.*$', text, re.M).group(0)
    after = text.split(heading, 1)[1]
    block = re.split(r'^## ', after, maxsplit=1, flags=re.M)[0]
    return [line for line in block.splitlines() if line.startswith('|')]


report = {}
for lang in ['en', 'zh']:
    p = STAGE / 'manuscripts' / f'supplement_{lang}.md'
    text = p.read_text('utf-8')
    old = (OLD / 'manuscripts' / p.name).read_text('utf-8')
    prefix = '## Supplementary Table S' if lang == 'en' else '## 补充表S'
    nums = re.findall('^' + re.escape(prefix) + r'(\d+)\.', text, re.M)
    verify(lang + '_tables_S1_through_S15', nums == [str(n) for n in range(1, 16)])
    links = re.findall(r'!\[[^]]+\]\(([^)]+)\)', text)
    verify(lang + '_17_figure_links', len(links) == 17)
    verify(lang + '_figure_files_exist', all((p.parent / x).exists() for x in links))
    verify(lang + '_figures_order_S01_to_S17', [Path(x).stem for x in links] == [f'Figure_S{n:02}' for n in range(1, 18)])
    verify(lang + '_no_old_title', 'Food bioactives in fatty liver research' not in text and '食源性活性成分在脂肪肝研究中的分子证据评估' not in text)
    verify(lang + '_one_author_placeholder', text.count('AUTHOR_INPUT_NEEDED') == 1)
    for num in list(range(1, 9)) + [10]:
        verify(f'{lang}_archived_table_S{num}_rows_unchanged', rows(text, num, lang) == rows(old, num, lang))
    verify(lang + '_S11_eight_raw_rerun_rows', len(rows(text, 11, lang)) == 10)
    verify(lang + '_S12_all15_both_species', len(rows(text, 12, lang)) == 34)
    verify(lang + '_S13_three_thresholds', len(rows(text, 13, lang)) == 5)
    verify(lang + '_S14_twelve_subtype_library_rows', len(rows(text, 14, lang)) == 14)
    verify(lang + '_S15_eight_group_library_rows', len(rows(text, 15, lang)) == 10)
    # New table display precision must match executed sources.
    matched = pd.read_csv(STAGE / 'areg_analysis/areg_matched_scale.csv')
    for i, r in matched.iterrows():
        verify(f'{lang}_S11_delta_{i}_matches_raw_recalculation', f'{r.PA_minus_control:+.6f}' in '\n'.join(rows(text, 11, lang)))
    verify(lang + '_negative_SMC_excluded_direction_visible', '-0.124007' in text and '-0.048424' in text)
    verify(lang + '_adult_negative_and_sparse_P1_kept', '| P1TLH | 8 | 96 | -1.0417 | -0.00159 |' in text)
    verify(lang + '_adult_sparse_P4_kept', '| P4TLH | 13 | 43 |' in text)
    verify(lang + '_Welch_and_exact_main_q_retained', '0.00257' in text and '0.1984' in text)
    verify(lang + '_new_BH15_is_separate', '0.1190' in text and 'BH15' in text)
    report[lang] = {'characters': len(text), 'words_whitespace': len(text.split()),
                    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}

for n in range(1, 15):
    src = OLD / f'figures/supplement/Figure_S{n:02}.png'
    dest = STAGE / f'figures/supplement/Figure_S{n:02}.png'
    verify(f'old_Figure_S{n:02}_image_unmodified', hashlib.file_digest(src.open('rb'), 'sha256').hexdigest() == hashlib.file_digest(dest.open('rb'), 'sha256').hexdigest())

(HERE / 'supplement_checks.json').write_text(json.dumps({'all_passed': all(checks.values()), 'n_checks': len(checks),
    'checks': checks, 'files': report, 'scope': 'Text, table transcriptions, source image integrity and links. DOCX/PDF layout not checked here.'}, indent=2), encoding='utf-8')
print(json.dumps({'all_passed': all(checks.values()), 'n_checks': len(checks), 'files': report}, indent=2))
