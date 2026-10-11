"""Build a deterministic orthology ledger from saved official responses."""
from __future__ import annotations
import csv
import hashlib
import io
import json
import re
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

HERE = Path(__file__).resolve().parent
RAW = HERE / 'raw'
collected = json.loads((HERE / 'collected.json').read_text(encoding='utf-8'))
records = json.loads((HERE / 'source_records.json').read_text(encoding='utf-8'))
checks = []

def check(name, ok):
    checks.append({'name': name, 'passed': bool(ok)})
    if not ok:
        raise ValueError(name)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

for r in records:
    check('hash_' + r['key'], digest(HERE / r['file']) == r['sha256'])
check('candidate_count', len(collected['candidates']) == 25)
check('unique_original_ids', len({r['original_human_id'] for r in collected['candidates']}) == 25)
check('rest_release_116', collected['metadata']['info_data']['releases'] == [116])
check('rest_software_116', collected['metadata']['info_software']['release'] == 116)
registry = ET.fromstring((RAW / 'biomart_registry.txt').read_text(encoding='utf-8'))
mart_db = [r.attrib for r in registry if r.attrib.get('name') == 'ENSEMBL_MART_ENSEMBL']
check('same_biomart_release', len(mart_db) == 1 and mart_db[0]['database'] == 'ensembl_mart_116')
bm_lines = list(csv.reader(io.StringIO((RAW / 'biomart_25.txt').read_text(encoding='utf-8')), delimiter='\t'))
check('biomart_10_fields', all(len(r) == 10 for r in bm_lines if r))
bm_rows = [r for r in bm_lines[1:] if r]
bm_by_id = {}
for row in bm_rows:
    bm_by_id.setdefault(row[0], []).append(row)
check('biomart_covers_all_25', set(bm_by_id) == {r['original_human_id'] for r in collected['candidates']})

entries = []
for c in collected['candidates']:
    gid = c['original_human_id']
    lookup = c.get('human_lookup') or {}
    hom_resp = c.get('homology_response') or {}
    e = {'original_human_id': gid, 'input_symbol': c['input_symbol'],
         'official_current_human_symbol': lookup.get('display_name'),
         'human_gene_version': lookup.get('version'), 'human_biotype': lookup.get('biotype'),
         'human_assembly': lookup.get('assembly_name'),
         'human_id_resolves': c['human_lookup_status'] == 200 and lookup.get('id') == gid and lookup.get('species') == 'homo_sapiens' and lookup.get('object_type') == 'Gene',
         'mouse_orthologues': [], 'primary_comparison_eligible': False,
         'original_id_changed': False, 'original_symbol_changed': bool(lookup.get('display_name') and lookup['display_name'] != c['input_symbol'])}
    if not e['human_id_resolves']:
        e['status'] = 'human_id_unresolved'
        e['archive_response'] = c.get('archive_response')
    elif c['homology_http_status'] != 200 or not isinstance(hom_resp.get('data'), list):
        e['status'] = 'homology_api_failure'
    else:
        blocks = hom_resp['data']
        check('exact_homology_block_' + gid, len(blocks) == 1 and blocks[0]['id'] == gid)
        hs = blocks[0]['homologies']
        check('same_rest_biomart_targets_' + gid, {h['target']['id'] for h in hs} == {b[2] for b in bm_by_id[gid] if b[2]})
        for h in hs:
            mid = h['target']['id']
            ml = c['mouse_lookups'][mid]
            bms = [b for b in bm_by_id[gid] if b[2] == mid]
            check('single_biomart_pair_' + gid + '_' + mid, len(bms) == 1)
            b = bms[0]
            check('same_orthology_type_' + gid + '_' + mid, h['type'] == b[4])
            check('exact_species_' + gid + '_' + mid, h['source']['id'] == gid and h['source']['species'] == 'homo_sapiens' and h['source']['taxon_id'] == 9606 and h['target']['species'] == 'mus_musculus' and h['target']['taxon_id'] == 10090)
            mouse_ok = ml.get('id') == mid and ml.get('species') == 'mus_musculus' and ml.get('object_type') == 'Gene' and ml.get('assembly_name') == 'GRCm39'
            check('mouse_lookup_valid_' + mid, mouse_ok)
            check('symbol_matches_' + mid, ml.get('display_name') == b[3])
            confidence = int(b[5]) if b[5] else None
            raw_rest_identity = [h['source']['perc_id'], h['target']['perc_id']]
            raw_bm_identity = [float(b[6]), float(b[7])]
            check('identity_pair_values_agree_' + gid + '_' + mid, sorted(raw_rest_identity) == sorted(raw_bm_identity))
            e['mouse_orthologues'].append({
                'mouse_id': mid, 'mouse_symbol': ml.get('display_name'), 'mouse_gene_version': ml.get('version'),
                'mouse_biotype': ml.get('biotype'), 'mouse_assembly': ml.get('assembly_name'),
                'orthology_type': h['type'], 'official_orthology_confidence': confidence,
                'rest_source_perc_id': h['source']['perc_id'], 'rest_target_perc_id': h['target']['perc_id'],
                'biomart_identity_raw': {bm_lines[0][6]: b[6], bm_lines[0][7]: b[7]},
                'goc_score': float(b[8]) if b[8] else None, 'wga_coverage': float(b[9]) if b[9] else None,
                'taxonomy_level': h.get('taxonomy_level'), 'method_link_type': h.get('method_link_type'),
                'source_protein_id': h['source'].get('protein_id'), 'target_protein_id': h['target'].get('protein_id'),
                'mouse_id_verified': mouse_ok})
        if not hs:
            e['status'] = 'current_human_id_no_mouse_orthologue_in_ensembl116'
            e['exclusion_reason'] = 'Successful exact-ID lookup and successful orthology query; empty mouse orthologue set. This does not assert absence of analogous mouse biology.'
        elif len(hs) == 1 and hs[0]['type'] == 'ortholog_one2one':
            e['primary_comparison_eligible'] = all(o['official_orthology_confidence'] == 1 and o['mouse_id_verified'] for o in e['mouse_orthologues'])
            e['status'] = 'one2one_high_confidence' if e['primary_comparison_eligible'] else 'one2one_confidence_pending_or_low'
        else:
            e['status'] = 'non_one2one_mapping'
            e['exclusion_reason'] = 'Not one-to-one; retain all targets without choosing a matching expression direction.'
    entries.append(e)

primary = [e for e in entries if e['primary_comparison_eligible']]
check('all_primary_mouse_ids_unique', len({e['mouse_orthologues'][0]['mouse_id'] for e in primary}) == len(primary))
check('cxcl8_no_proxy', next(e for e in entries if e['original_human_id'] == 'ENSG00000169429')['mouse_orthologues'] == [])
check('acat2_two_targets_not_primary', next(e for e in entries if e['original_human_id'] == 'ENSG00000120437')['status'] == 'non_one2one_mapping')
summary = {'total_candidates': len(entries), 'primary_eligible': len(primary), 'status_counts': dict(Counter(e['status'] for e in entries)),
           'primary_symbols': [e['input_symbol'] for e in primary], 'expression_analysis_performed': False,
           'sequence_files_downloaded': False, 'functional_analog_substitutions': 0, 'id_replacements': 0,
           'mapping_only_not_validation': True}
release = {'ensembl_rest_data_release': 116, 'ensembl_rest_api': collected['metadata']['info_rest'],
           'ensembl_biomart_database': mart_db[0]['database'], 'biomart_actual_host': next(r['final_url'] for r in records if r['key'] == 'biomart_registry').split('/')[2],
           'human_assembly': collected['metadata']['assembly_human']['assembly_name'],
           'mouse_assembly': collected['metadata']['assembly_mouse']['assembly_name'],
           'mouse_assembly_accession': collected['metadata']['assembly_mouse']['assembly_accession'],
           'mouse_genebuild_last_geneset_update': collected['metadata']['assembly_mouse']['genebuild_last_geneset_update'],
           'reference_freeze_for_later_quantification': 'Ensembl release 116, Mus musculus GRCm39; sequence/annotation files not downloaded or hashed in this task'}
output = {'summary': summary, 'release': release, 'input_sha256': collected['candidate_sha256'],
          'mapping_rule_sha256': collected['mapping_rule_sha256'], 'candidates': entries,
          'identity_field_note': 'REST source/target perc_id and BioMart headers are retained verbatim as separate fields. Values match as a pair; no direction-label conversion is inferred.',
          'scope': 'Gene-identity/homology readiness only. No conserved treatment response, causal effect or efficacy is inferred.'}
(HERE / 'orthology.json').write_text(json.dumps(output, indent=2, ensure_ascii=False), encoding='utf-8')
lines = ['# 人—鼠候选基因映射（stage07）', '',
         '已执行官方元数据核对，未读取小鼠表达结果。25 个原始人基因全部保留；15 个符合冻结的主要比较规则，1 个是一对多，9 个在本版没有返回鼠直系同源。', '',
         '## 版本与规则', '',
         '- Ensembl REST 数据与软件 release 116；BioMart 实际重定向到 `jun2026.archive.ensembl.org`，registry 为 `ensembl_mart_116`。',
         '- 人参考为 GRCh38；鼠为 GRCm39 / GCA_000001635.9，API 返回最近 geneset 更新为 2025-11。后续应使用 release 116 的注释及序列，不用滑动 latest。此步未下载 FASTA/GTF。',
         '- 主要映射规则在拉取候选映射前写入 [mapping_rule.json](mapping_rule.json)：原始 ENSG 正常解析、恰好一个官方 one-to-one 鼠同源、BioMart orthology confidence=1，鼠 ENSMUSG 与符号经独立 lookup 核对。',
         '- 未按鼠表达方向挑选基因；未替换旧 ID；未用功能类似基因冒充直系同源。confidence 是 Ensembl 的注释可信标记，不是治疗效应可信度。', '',
         '## 全部 25 个候选', '',
         '| 原始人 ID | 输入符号（当前官方符号有差异时另列） | 鼠同源 | 类型/状态 | 官方 confidence | 主要比较 |',
         '|---|---|---|---|---|---|']
for e in entries:
    name = e['input_symbol']
    if e['original_symbol_changed']:
        name += ' → ' + e['official_current_human_symbol'] + '（同一 ENSG 的当前命名）'
    targets = '<br>'.join(o['mouse_symbol'] + ' / ' + o['mouse_id'] for o in e['mouse_orthologues']) or '本版未返回'
    types = ', '.join(sorted({o['orthology_type'] for o in e['mouse_orthologues']})) or '当前人 ID 有效；无鼠直系同源返回'
    conf = ', '.join(str(o['official_orthology_confidence']) for o in e['mouse_orthologues']) or '不适用'
    lines.append(f"| {e['original_human_id']} | {name} | {targets} | {types} | {conf} | {'是' if e['primary_comparison_eligible'] else '否'} |")
lines += ['', '## 对研究的影响', '',
          '- **CXCL8：不能用这套小鼠数据做该基因的一对一验证。** 原始 ENSG 查询成功，但 REST 与 BioMart 均未给鼠直系同源；不采用 Cxcl1/Cxcl2 替代。没有对应并不表示小鼠没有类似炎症过程。',
          '- **ACAT2：同时对应鼠 Acat2 与 Acat3，均被标为 one-to-many / high confidence。** 主要一对一比较排除它，两个鼠基因均单列保留，不取表达更好看的那个。',
          '- **VNN1 → Vnn1** 为官方 high-confidence one-to-one，可进入预定映射集；这不预先说明芝麻素反应是否一致。',
          '- 其余未映射项是有效人 ID 的空鼠直系同源返回；本轮没有弃用人 ID，也没有候选查询 API 失败。非编码转录本的缺少直系同源不能被解释为表达量为零。',
          '- 15 个主要可映射基因：' + '、'.join(e['input_symbol'] for e in primary) + '。', '',
          '## 复现与追溯', '',
          '- 使用项目 `.venv/Scripts/python.exe`；两个脚本均只依赖 Python 标准库。',
          '- `fetch_orthology.py` 保存完整响应及每个 URL、最终 URL、UTC 时间、HTTP 状态、SHA256；成功与失败缓存均保留。再次运行默认复用原缓存。需要未来重新查询时可另指定 `--refresh-suffix`，不可将新版当作此冻结版。',
          '- `build_orthology.py` 仅从保存响应生成 [orthology.json](orthology.json)、此表和 [validation.json](validation.json)。复现本版不要更改映射规则。',
          '- [source_records.json](source_records.json) 对应 75 个元数据响应，共 74,530 bytes；[collected.json](collected.json) 为合并原始元数据。',
          '- [documentation_records.json](documentation_records.json) 保存文档抓取结果。homology API 和 archive API 文档获取成功；归档站质量控制与 BioMart 文档 HTTP 403，已记录，不伪称拿到全文。官方质量控制页面和 BioMart 官方文档另经 web 搜索核读。',
          '- 采集最初为顺序小请求，完成部分缓存后主动中止并改为最多 3 并发；已存在原始响应经哈希核对后复用，没有把中止过程作为完成或丢弃失败证据。', '',
          '## 官方依据', '',
          '- [Ensembl homology API](https://rest.ensembl.org/documentation/info/homology_species_gene_id)：按原始 ENSG 与目标物种查询。',
          '- [Ensembl orthology quality controls](https://static.ensembl.org/info/genome/compara/Ortholog_qc_manual.html)：confidence 综合序列 identity、基因顺序/全基因组比对等标准；本项目直接使用官方 flag，不自行杜撰阈值。',
          '- [BioMart 官方 REST 文档](https://jun2026-bacteria.ensembl.org/info/data/biomart/biomart_restful.html)：此文档镜像说明 XML 接口，实际查询服务器和 release 以本地源记录为准。', '',
          '本结果仅完成身份与同源关系核查，不构成跨物种表达复现、芝麻素功效或新机制验证。', '']
(HERE / 'README.md').write_text('\n'.join(lines), encoding='utf-8')
validation = {'passed': all(c['passed'] for c in checks), 'check_count': len(checks), 'checks': checks,
              'summary': summary, 'output_sha256': digest(HERE / 'orthology.json'),
              'rule_sha256': digest(HERE / 'mapping_rule.json'), 'builder_sha256': digest(Path(__file__))}
(HERE / 'validation.json').write_text(json.dumps(validation, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps({'passed': validation['passed'], 'checks': len(checks), 'summary': summary}, ensure_ascii=False))
