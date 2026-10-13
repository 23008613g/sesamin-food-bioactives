"""Build PFHN supplement revisions from stage28, preserving old tables/figures.

New numerical display tables are transcribed from executed stage31 CSV outputs.
Does not edit main manuscripts, images, source data, or earlier stages.
"""
from pathlib import Path
import hashlib
import json
import re
import shutil
import pandas as pd

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
OLD = STAGE.parent / 'pfhn_submission_28'
OUT = STAGE / 'manuscripts'
OUT.mkdir(exist_ok=True)


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(str(x) for x in row) + ' |' for row in rows])


def section(text, heading, replacement):
    pat = r'^' + re.escape(heading) + r'\n.*?(?=^#{1,3} |\Z)'
    result, n = re.subn(pat, lambda m: replacement.rstrip() + '\n\n', text,
                        count=1, flags=re.M | re.S)
    assert n == 1, heading
    return result


def sign(v):
    return {1: '↑', -1: '↓', 0: '0'}[int(v)]


def yes(v, lang):
    return ('Yes' if v else 'No') if lang == 'en' else ('是' if v else '否')


title_en = 'Hepatic transcriptomic evidence for the sesame lignan sesamin: cross-study comparison and cellular context'
title_zh = '芝麻木脂素芝麻素的肝脏转录组证据：跨研究比较与细胞背景分析'
sources = {
    'candidate_direction_order.csv': STAGE / 'candidate_analysis/candidate_direction_order.csv',
    'adjacent_threshold_summary.csv': STAGE / 'candidate_analysis/adjacent_threshold_summary.csv',
    'areg_subtype_by_library.csv': STAGE / 'areg_analysis/areg_subtype_by_library.csv',
    'areg_excluding_smc.csv': STAGE / 'areg_analysis/areg_excluding_smc.csv',
    'areg_matched_scale.csv': STAGE / 'areg_analysis/areg_matched_scale.csv',
    'areg_adult_donor_context.csv': STAGE / 'areg_analysis/areg_adult_donor_context.csv',
    'rat_all25_reanalysis.csv': STAGE / 'rat_analysis/results/candidate_results_all25_welch_exact_two_normalizations.csv',
}
for name, p in sources.items():
    shutil.copyfile(p, HERE / name)
cand = pd.read_csv(sources['candidate_direction_order.csv'])
threshold = pd.read_csv(sources['adjacent_threshold_summary.csv'])
subtypes = pd.read_csv(sources['areg_subtype_by_library.csv'])
exclude = pd.read_csv(sources['areg_excluding_smc.csv'])
matched = pd.read_csv(sources['areg_matched_scale.csv'])
lib_names = {'GSM6322577': 'C1', 'GSM6322578': 'C2', 'GSM6322579': 'P1', 'GSM6322580': 'P2'}

en = (OLD / 'manuscripts/supplement_en.md').read_text('utf-8')
zh = (OLD / 'manuscripts/supplement_zh.md').read_text('utf-8')
en = en.replace('Food bioactives in fatty liver research: assessing molecular evidence through a sesamin case study', title_en)
zh = zh.replace('食源性活性成分在脂肪肝研究中的分子证据评估：以芝麻素为例', title_zh)
en = en.replace('Tables S1–S11 and Figures S1–S14', 'Tables S1–S15 and Figures S1–S17')
zh = zh.replace('表S1–S11与图S1–S14', '表S1–S15与图S1–S17')

en = section(en, '### Study design and data roles', '''### Study design and data roles

The study reanalyzed public data, with no newly performed biological experiments. GSE248903 supplied three HepG2 libraries each labelled Control, Model and Treat for exploratory candidate selection. Its deposited protocol describes OA exposure followed by sesamin in Treat. GSE42693 supplied a separately documented dietary sesamin comparison in healthy rat liver [1]. Neither dataset alone establishes a treatment mechanism in human fatty liver.

Mouse PRJCA021826/CRA013797 and its DRP018055/PRJDB40955 mirror represent the same nine libraries, not independent studies. The project descriptions identify normal, PA/OA and Ses-treated PA/OA hepatocytes, but the retrieved sample records do not explicitly map ctrl1–3, HFD1–3 and sess1–3 to those conditions. Numerical comparisons therefore retain the deposited labels as supplementary exploration and are not evidence of a confirmed sesamin-response replication. HFD is not expanded into a demonstrated in-vivo high-fat-diet group.

GSE247407 (30 primary human hepatocyte libraries), GSE207889 (four control/PA technical organoid libraries from one differentiation) and GSE115469 (five adult-liver donors) supply cellular or fatty-acid context without sesamin intervention. Figure 1 and main Table 1 distinguish these roles. The full source audit is provided as Online Resource 2, revision31/metadata/sample_metadata_audit.csv, with 77 records: nine human, nine mouse, 20 rat, 30 primary hepatocyte, four organoid and five adult-liver records.''')

en = section(en, '### Sample-design and source checks', '''### Sample-design and source checks

The bounded metadata audit was repeated on 11 October 2026 using GEO series/sample records, official NGDC and ENA project/sample/experiment records, the original rat paper and the organoid correction. Cached original records are distinguished from fresh retrievals in revision31/metadata/metadata_evidence_log.md and its source manifest. Human, mouse and organoid library numbers were not equated to independently verified cultures or animals.

All nine GSE248903 records explicitly retain Control, Model and Treat labels. The common protocol describes 24 h OA exposure followed by a further 36 h sesamin for Treat. Deposited OA and sesamin strings lack unit prefixes; they are not silently repaired to micromolar. Matched final harvest times, continued OA, washing/medium changes, vehicle matching and independent culture identities remain unresolved. Fresh GEO records supply no linked publication; related papers [2,3] do not establish an accession-to-library link and their doses are not imported.

All nine mouse mirror records link back to the original SAMC/CRX identifiers and PRJCA021826. They identify male, eight-week-old source animals and hepatocyte material, but lack animal/isolation identities, pooling, independent culture units, doses, duration, vehicle and explicit label-to-condition mapping. Project-level mention of sesamin does not resolve those library-level gaps. This source-labelled comparison remains descriptive and supplementary; no biological-replicate DESeq2 model was added.

GSE247407 describes five patient labels, three repeats per condition and seven days of Control/FFA culture. F95 and F95P remain separate labels whose donor relationship and within-label replication are unresolved. All 30 libraries are retained.

The corrected organoid methods [8,9] specify fatty-acid preparation and carrier controls: PA stock in ethanol is conjugated with BSA before dilution; OA stock preparation uses NaOH and BSA before dilution. PA treatment is described from day 21 to day 25, whereas the four GEO records state five treatment days. This discrepancy remains recorded. The source Statistical analysis defines scRNA-seq replicates as suspensions from one differentiation; a general library-preparation phrase does not establish independent differentiations. The correction identifies original Fig. 4B labelling errors and omitted preparation details, not an AREG annotation correction. Correction text was read; a visual verification of individual corrected figure labels was not achieved and is not claimed.''')

threshold_en = '''

An additional post-hoc sensitivity repeated the candidate rule at neighboring fold thresholds of 1.4 and 1.6, holding the expression filter, sign reversal, strictly intermediate Treat ordering, CPM agreement and at-least-22-of-27 stability requirement fixed. The plan was recorded before this rerun with earlier results already known. Final set sizes were 39, 25 and 15 at thresholds 1.4, 1.5 and 1.6, respectively; VNN1, AREG and JAG2 remained in all three sets (Table S13). The original 25 candidates remain the reporting family; the sensitivity does not replace them or constitute held-out validation. Raw counts were reread for this calculation; a fitted differential-expression model was not introduced.
'''
en = en.replace('### Mouse quantification and orthology', threshold_en + '\n### Mouse quantification and orthology', 1)
en = en.replace('### Independent dietary rat comparison', '''### Independent dietary rat comparison

The revision reread the deposited series matrix and GPL1355 annotation. GEO documents GCOS 1.4 global scaling, and matrix values (minimum 0.0221933, median 27.04825, maximum 18,018.4) are unlogged signal. The original article additionally describes downstream GeneSpring normalization and filtering; those values are not substituted for the deposited matrix. We retain log2(signal + 1), with no double-log correction, raw CEL reprocessing or new detection-call filtering.
''', 1)
en = en.replace('### Primary hepatocyte and organoid context', '''The fixed candidate-to-probe audit retained 15 measurable candidates among the 25: nine lack a returned rat orthologue and one lacks eligible one-to-one orthology. Each measurable candidate has one eligible GPL1355 probe under the single-Entrez rule, so the per-animal median-probe rule reduces to that probe. VNN1 maps to rat Vnn1 (Entrez 29142), probe 1389253_at. No probe was selected by effect size, and no sequence-level remapping was performed. A new explicitly post-hoc BH15 sensitivity used only measurable candidates while retaining BH25 as primary. VNN1 exact-label q was 0.1984 under BH25 and 0.1190 under BH15; no candidate passed exact-label q <0.05 under either family or normalization. Full probe, animal, NA and dual-inference records are supplied in revision31/rat_analysis/results/.

### Primary hepatocyte and organoid context''', 1)

en = section(en, '### Added matched-scale analysis and original label crosswalk', '''### Raw-source AREG verification, label provenance and exclusion sensitivity

The earlier matched-scale analysis reweighted archived category summaries. In this post-hoc revision, with previous results already known, we reread the original GSE207889 H5AD obs fields and all four MTX matrices, features and barcode files. The protocol and source hashes were recorded before calculation. The broad label is deposited in obs/celltype_category and the subtype in obs/celltype_scTypedb; neither was generated by relabelling cells in this reanalysis. The H5AD hash also matched the original download/decompression record.

Library accession plus barcode gave unique matches for all 3,093, 5,510, 7,399 and 8,032 author-retained cells, respectively. AREG counts and per-cell totals were extracted from all 33,538 measured features, and log1p(CP10k) was recomputed with measured zeros included. No cell count was forced to match the manuscript and no additional QC threshold, downsampling or annotation was applied. There were no unmatched retained cells. Absent features or empty groups are represented as NA, not zero; every AREG subtype used here was represented in all four libraries.

The deposited broad Cholangiocytes category contains source subtypes Cholangiocytes (213), Ductal cells (2,965) and Smooth muscle cells (241), total 3,419. Hepatocytes contains AH (10,672), FH1 (4,125), FH2 (482), HB1 (210) and HB2 (2,608), total 18,097. Labels identify deposited groups, not independently validated mature identities. Table S14 separates the three source subtypes. The additional exclusion sensitivity combines Cholangiocytes and Ductal cells (3,178 cells) without SMC; both the original category and exclusion result remain visible (Table S15).

Raw AREG counts were nonzero in 188/213 Cholangiocytes, 2,453/2,965 Ductal cells and 180/241 SMC. The original broad category had 2,821/3,419 nonzero cells (82.51%); excluding SMC gave 2,641/3,178 (83.10%), compared with 7,290/18,097 (40.28%) in Hepatocytes. Both the original broad category and the SMC-excluded group retained higher raw detection and mean log1p(CP10k) than Hepatocytes within each library. These are descriptive comparisons, not cell-level significance tests or evidence of specificity.

Let m_lk and n_lk be the zero-inclusive mean log1p(CP10k) and retained-cell count in library l and group k. Whole-library means are sum_k(n_lk m_lk)/sum_k(n_lk); the Hepatocytes mean uses its corresponding cells. Equal-library condition summaries average the two technical-library means; pooled-cell sensitivity instead weights them by relevant cell count. On this fixed scale, raw recomputation reproduced PA-minus-control differences of +0.054898 in whole organoid and +0.039797 in Hepatocytes. The original broad Cholangiocytes and SMC-excluded groups had negative differences, -0.129474 and -0.124007 (Table S11). Thus the positive whole/hepatocyte comparison is not a claim of uniform change across cell populations.

The earlier hepatocyte pooled-UMI CPM contrast was also recovered (-0.111612 on the log2 ratio scale). It changes weighting and transformation, so its contrast with the positive mean log1p(CP10k) is not a population-only sign reversal or a causal composition effect. Prior all-25-candidate summaries remain available; the new raw-source extraction and its 95 passing arithmetic/identity checks are specifically an AREG audit, not a complete new 25-gene analysis or biological validation. The raw matrix, author-label and output hashes, cell maps and old/new reconciliations are provided in revision31/areg_analysis/.

For the adult atlas, the deposited normalized CSV and original annotations were reread for AREG, preserving all 8,444 cell IDs and five donors. The values were not normalized again. All five donor differences agree with the earlier stored results to recorded precision, including the negative P1 difference and sparse P1/P4 groups (Table S8). Primary means remain +0.071352 source log2CPM and +8.279730 percentage points for Cholangiocytes minus Hepatocytes. The >=20-cell sensitivity retains three donors, without deleting either sparse donor from the primary analysis. Immune-cell expression remains in the complete context tables; source nonzero percentages are not relabelled raw-UMI detection.''')

en = en.replace('Main Fig. 3a and Supplementary Figs. S9, S3 and S4', 'The candidate displays and Supplementary Figs. S9, S3 and S4')
en = en.replace('the full-filtered-gene PCA in main Fig. 2', 'the full-filtered-gene PCA in Fig. S15')
en = en.replace('mouse intervention −37.7%', 'mouse sess/HFD label contrast −37.7%')
en = en.replace('Mouse intervention −37.7%', 'Mouse sess/HFD label contrast −37.7%')
en = en.replace('| Intervention agreement | Both contrasts agree |', '| sess/HFD label-sign agreement | Both label contrasts agree |')
en = en.replace('Stable in cell data; opposite rat direction;', 'Stable human and mouse label directions; mouse intervention assignment unresolved; opposite rat direction;')
en = en.replace('Remaining full-data intervention-concordant candidates;', 'Remaining full-data sess/HFD-label-concordant candidates;')
en = en.replace('Stable intervention directions also deserve retention;', 'Stable sess/HFD label directions also deserve retention;')
en = en.replace('Intervention direction is stable,', 'The sess/HFD label direction is stable,')
en = en.replace('Full-data intervention signs oppose human;', 'Full-data sess/HFD label signs oppose human;')
en = en.replace('Full-data model and intervention signs oppose human', 'Full-data HFD/ctrl and sess/HFD label signs oppose human')
en = en.replace('Mouse intervention −19.7%; −14.6%', 'Mouse sess/HFD label contrast −19.7%; −14.6%')
en = en.replace('Full-data mouse intervention −2.9%; −58.3%', 'Full-data mouse sess/HFD label contrast −2.9%; −58.3%')
en = en.replace('Full-data mouse intervention −13.6%; −0.086%', 'Full-data mouse sess/HFD label contrast −13.6%; −0.086%')
en = en.replace('intervention range −76.6% to +78.1%', 'sess/HFD range −76.6% to +78.1%')
en = en.replace('Library-dependent intervention signal; cholangiocyte-category detection exceeds hepatocytes in all four organoid libraries.', 'Library-dependent mouse label signal; the deposited broad category and its SMC-excluded sensitivity exceed Hepatocytes in raw detection in all four organoid libraries.')
en = en.replace('“Both” means nonzero concordance for model/control and intervention/model.', '“Both” means nonzero concordance of the two within-study label contrasts. Mouse labels do not establish confirmed intervention assignments.')
en = en.replace('Accompanying tables are display transcriptions, not rerun analyses;', 'Tables S1–S10 retain archived values except where explicitly verified against this revision; Tables S11–S15 transcribe the new raw-source or candidate-rule calculations;')
en = en.replace('This verifies the recorded processed-matrix/candidate-count scope; it does not rerun FASTQ alignment or quantification, raw CEL preprocessing, all-gene organoid UMI extraction, UMAP fitting or installation in a newly provisioned environment.', 'That earlier run verifies its recorded processed-input scope. The present revision additionally reread the human count matrix, rat deposited signal/probe annotation and the AREG raw organoid matrices/author annotations. Neither stage reran FASTQ quantification, raw CEL preprocessing, UMAP fitting or installation in a newly provisioned environment.')
en = en.replace('We compared the model and intervention directions across the Cartesian product', 'We compared the two within-study numerical label directions across the Cartesian product')
en = en.replace('Mouse intervention −62.5%', 'Mouse sess/HFD label contrast −62.5%')
en = en.replace('A locally frozen protocol specified an external descriptive check', 'For the original external-atlas analysis, a locally frozen protocol specified an external descriptive check')
en = en.replace('The prespecified ≥20-cells sensitivity includes P2, P3 and P5', 'The ≥20-cells sensitivity specified for the original atlas analysis includes P2, P3 and P5')
en = en.replace('VNN1 increases, opposing the original human/mouse intervention directions.', 'VNN1 increases, opposing human Treat/Model and the exploratory mouse sess/HFD label direction; the mouse exposure assignment remains unresolved.')
en = en.replace('Grey, orange and teal identify the respective control, model and treatment groups.', 'Grey, orange and teal identify Control/Model/Treat or ctrl/HFD/sess source labels; mouse exposure assignments remain unresolved.')
en = en.replace('The public v24 archive predates this addition and is not claimed to contain it.', 'The public v24 archive predates this addition and all revision31 modules and is not claimed to contain them.')

# Chinese audit methods: targeted replacement, retaining correct original formulas.
zh = section(zh, '### 设计、单位与人源归一化', '''### 设计、单位与人源归一化

正文表1区分六套资料的角色。2026年10月11日的限定范围元数据复核形成77条逐样本记录：人9、小鼠9、大鼠20、原代肝细胞30、类器官4及成人肝脏5。Online Resource 2的revision31/metadata/sample_metadata_audit.csv及metadata_evidence_log.md逐项区分原始字段、已确认信息与缺口，并区分新获取与复用的官方记录。

GSE248903九库明确标为Control、Model和Treat。共同协议描述OA 24小时后Treat另加36小时芝麻素，但OA“400 M”和芝麻素“16 M”缺浓度前缀，不能自动改成μM。对照和模型的同步最终收样、持续OA、洗脱换液、溶剂匹配及独立培养身份仍未确认。新GEO记录未给出关联论文；相关论文的条件不能据主题相似直接移入测序描述。因此仅作探索性文库层面候选筛选，不将Treat差值称为已确立的因果获益。

小鼠PRJCA021826/CRA013797与DRP018055/PRJDB40955镜像经原SAMC/CRX映射确认为同九库。项目级文字描述正常、PA/OA及芝麻素相关处理，但逐样本记录只有ctrl、HFD、sess标签，未提供明确对应表。材料为雄性、8周龄来源肝细胞；动物、分离、混样、培养单位、浓度、时序、溶剂均未解决。其数值保留为补充标签层面探索，不作为已确认的芝麻素反应重复验证，HFD亦不展开为确定的体内高脂饮食组。

GSE247407描述五名患者标签、每条件三次重复及七天培养，但F95与F95P的身份关系及标签内重复性质仍未独立确认。全部30库保留作无芝麻素的脂肪酸背景。

记非负计数K_ij，i为基因、j为文库。CPM分母使用各库全部分配到基因的计数，CPM≥1且至少3库达标给出13,014条表达记录。原九库中15,412基因均为正，以g_i=(∏_j K_ij)^(1/9)计算几何均值，初因子为median_i(K_ij/g_i)，再将因子几何均值缩为1；N_ij=K_ij/s_j。该调整依赖参考分布整体稳定或近似平衡，不估计每细胞绝对RNA；未拟合DESeq2检验模型。

两个主对比为log2[(Model均值+1)/(Control均值+1)]和log2[(Treat均值+1)/(Model均值+1)]，共用Model组。初始要求反号、绝对值均≥log2(1.5)、Treat严格在另外两均值之间，且CPM方向与排序一致；98条满足。各组省略1库构成27种重叠组合，在原表达过滤基因集上重估因子；至少22种同时保留原方向和完整规则才纳入，最终25候选。组合不独立，后续同数据使用不是保留测试集验证。

本轮在已知旧结果后冻结事后方案，从人源计数重新执行相邻1.4、1.5、1.6倍阈值，其余规则固定，最终分别39、25、15候选，VNN1、AREG、JAG2三者均保留（表S13）。原25候选仍是报告集合，不改成更有利子集，不以此声称预注册或独立验证。''')

zh = zh.replace('### 独立大鼠比较与统计敏感性', '''### 独立大鼠比较与统计敏感性

本轮重新读取存档GSE42693矩阵和GPL1355注释。GEO明确GCOS 1.4整体缩放；矩阵最小、中位、最大为0.0221933、27.04825、18,018.4，属于未取对数signal。原论文另有GeneSpring后续归一化和过滤，不能混称为本矩阵处理。保留log2(signal+1)，未发现双重log错误；没有重跑原CEL或新增detection-call过滤。
''', 1)
zh = zh.replace('### 脂肪酸背景与细胞汇总', '''固定25候选中15可测、9无返回大鼠同源、1不符合一对一规则。15个可测者各有一个符合单一Entrez规则的探针，因此逐动物多探针中位数在本数据退化为单探针，没有按效应选择。VNN1对应Vnn1、Entrez 29142、探针1389253_at；未作序列重映射。本轮另列事后BH15可测家族敏感性，但不替代主要BH25：VNN1精确标签q分别0.1984与0.1190，两种归一化及两种家族均没有精确标签q<0.05候选。完整探针、逐动物、NA及两套统计见revision31/rat_analysis/results/。

### 脂肪酸背景与细胞汇总''', 1)
zh = zh.replace('### 独立成人肝脏细胞背景', '''类器官原文及更正[8,9]补充脂肪酸配制和载体对照：PA乙醇储液与BSA结合后稀释，OA以NaOH配制并与BSA结合。原文明确第21天至第25天四天处理，而GEO写五天，差异保留不强行协调。Statistical analysis将scRNA-seq重复定义为同一次分化的细胞悬液，不能用另一处笼统“biological replicates”文字替代独立分化依据。更正针对原图4B标签及配制遗漏，不是AREG细胞注释更正；已读文字，未成功逐一视觉比较更正前后图标签。

### 独立成人肝脏细胞背景''', 1)
zh = section(zh, '### 新增统一尺度分析与原始标签对应', '''### 原始来源AREG重算、标签核对与去SMC敏感性

旧统一尺度分析由存档类别均值重新加权。本轮在已知旧结果后先记录事后方案及输入哈希，再重新读取GSE207889原H5AD的obs字段与四份MTX、features及barcode文件。大类来自原文件obs/celltype_category，亚型来自obs/celltype_scTypedb，不是我们重新合并标注；对象哈希亦与原下载解压记录一致。

以文库号加barcode唯一连接3,093、5,510、7,399、8,032个作者保留细胞，共24,034个，无重复或未匹配。逐细胞总UMI取全部33,538测量特征，重新提取AREG raw UMI并计算ln(1+10,000x/T)，含实测零。没有强凑旧细胞数、新质控门槛、下采样或新注释；空组或缺失特征应为NA而非零，本次AREG各亚型在四库均有细胞。

原宽类别Cholangiocytes包含Cholangiocytes亚型213、Ductal cells 2,965、Smooth muscle cells（SMC）241，共3,419。Hepatocytes含AH 10,672、FH1 4,125、FH2 482、HB1 210及HB2 2,608，共18,097。它们是原存储标签，不是独立验证的成熟身份。表S14拆分三原亚型；新敏感性只组合Cholangiocytes与Ductal cells、排除SMC，共3,178细胞，同时保留原宽类别（表S15）。

AREG raw UMI>0的细胞分别188/213、2,453/2,965及180/241。原宽类检出2,821/3,419（82.51%），去SMC为2,641/3,178（83.10%），Hepatocytes为7,290/18,097（40.28%）。原宽类与去SMC组在全部四库的检出率及含零均值都高于Hepatocytes；这是描述性观察，不是细胞层面显著性或胆管特异性证据。

记m_lk为逐库逐组含零log1p(CP10k)均值，n_lk为细胞数，全库均值sum_k(n_lk m_lk)/sum_k(n_lk)。主分析每条件两技术库等权；合并细胞敏感性按相应群体细胞数加权。原矩阵重算复现whole PA减对照+0.054898及Hepatocytes +0.039797；原宽Cholangiocytes及去SMC组则为−0.129474、−0.124007（表S11）。因此不能把whole/hepatocyte正方向写成全部群体方向一致。

旧肝细胞样pooled-UMI CPM对比亦复现为−0.111612（log2比值尺度）。它改变RNA权重和变换，不是仅换细胞群体导致的方向反转，也不证明构成因果机制。旧全25候选表继续保留；本轮raw提取和95项算术/身份检查限定于AREG，不声称全部25基因新重算或生物验证。逐细胞映射、原始及分析标签、哈希、旧新比较均在revision31/areg_analysis/。

成人图谱从原存档归一化CSV及注释重新读取AREG，保留8,444细胞与全部五供者，不再次归一化。逐供者结果在旧记录精度内一致，保留P1反向和P1/P4稀疏组（表S8）。主要胆管减肝细胞供者等权差为+0.071352 source log2CPM及+8.279730个原值非零比例百分点；两类≥20细胞敏感性只用三供者，但不从主结果删去稀疏供者。免疫细胞更高表达仍完整保留；成人原归一化值非零比例不能称raw UMI检出率。''')
zh = zh.replace('与正文图2的全部13,014过滤基因PCA不同', '与图S15的全部13,014过滤基因PCA不同')
zh = zh.replace('本次文本改写未重跑FASTQ', '本轮矩阵重算未重跑FASTQ')
zh = zh.replace('在查看候选表达前，本地冻结方案规定对GSE115469', '原成人图谱分析在查看候选表达前本地冻结方案，规定对GSE115469')
zh = zh.replace('配置比例基于均值汇总。“两项”指模型/对照及处理/模型均为非零且人鼠一致。', '配置比例基于均值汇总。“两项”指各研究两项标签对比均为非零且人鼠数值符号一致；鼠标签的暴露赋值尚未确认，不能解释成已确认的处理反应重复。')
zh = zh.replace('| 候选 | 处理方向一致 | 两项方向一致 |', '| 候选 | sess/HFD标签方向一致 | 两项标签方向一致 |')
zh = zh.replace('鼠处理−', '鼠sess/HFD标签对比−')
zh = zh.replace('细胞数据内稳定，大鼠相反；', '人和鼠标签数值内稳定，鼠暴露对应未确认，大鼠方向相反；')
zh = zh.replace('处理信号依赖文库；类器官胆管大类检出在四库均较高，应明确细胞来源。', '鼠标签信号依赖文库；原宽类别和去SMC组在四库raw检出均高于Hepatocytes，应明确细胞来源。')
zh = zh.replace('处理范围−76.6%', 'sess/HFD标签范围−76.6%')
zh = zh.replace('处理方向也稳定且保留；模型方向在全部配置均与人相反', 'sess/HFD标签方向也稳定且保留；HFD/ctrl标签方向在全部配置均与人相反')
zh = zh.replace('处理方向稳定，但全数据模型方向相反', 'sess/HFD标签方向稳定，但全数据HFD/ctrl标签方向相反')
zh = zh.replace('其余全数据处理方向一致者；模型方向及／或文库敏感性', '其余全数据sess/HFD标签方向一致者；HFD/ctrl标签方向及／或文库敏感性')
zh = zh.replace('全数据处理方向与人相反', '全数据sess/HFD标签方向与人相反')
zh = zh.replace('全数据模型及处理两项均与人相反', '全数据HFD/ctrl及sess/HFD标签方向均与人相反')
zh = zh.replace('配套表格为已有结果的展示转录，不是重新分析；', '表S1–S10保留存档值并注明本轮核验，表S11–S15转录新原始来源或候选规则计算；')
zh = zh.replace('该检查验证记录的处理矩阵／候选计数范围，没有重跑FASTQ比对或定量、原始CEL预处理、类器官全基因UMI提取、UMAP拟合或全新环境安装。', '该较早检查只验证其记录的处理后输入范围；本轮另重新读取人源计数、大鼠signal/探针及AREG类器官原矩阵和原注释。两阶段均未重跑FASTQ定量、原CEL预处理、UMAP拟合或全新环境安装。')
zh = zh.replace('VNN1上升与原人鼠细胞下降相反；', 'VNN1上升与人Treat/Model及探索性鼠sess/HFD标签方向相反，鼠暴露对应仍未明确；')

pop_en = {'Whole organoid': 'Whole organoid', 'Hepatocytes': 'Hepatocytes',
          'Deposited broad Cholangiocytes': 'Deposited broad Cholangiocytes',
          'Cholangiocytes + Ductal cells excluding SMC': 'Cholangiocytes + Ductal cells; no SMC'}
pop_zh = {'Whole organoid': '全类器官', 'Hepatocytes': 'Hepatocytes',
          'Deposited broad Cholangiocytes': '原宽Cholangiocytes',
          'Cholangiocytes + Ductal cells excluding SMC': 'Cholangiocytes+Ductal；去SMC'}
for lang in ['en', 'zh']:
    text = en if lang == 'en' else zh
    names = pop_en if lang == 'en' else pop_zh
    rows = [[names[r.population], ('Equal library' if r.weighting == 'equal_library' else 'Pooled cell') if lang == 'en' else ('文库等权' if r.weighting == 'equal_library' else '合并细胞'),
             f'{r.control_mean:.6f}', f'{r.PA_mean:.6f}', f'{r.PA_minus_control:+.6f}'] for r in matched.itertuples()]
    head = '## Supplementary Table S11. Matched-scale AREG means' if lang == 'en' else '## 补充表S11. 统一尺度的AREG均值'
    note = ('All values were recomputed from raw counts in this revision. Means include measured zeros and use ln(1 + 10000x/T), with T covering all 33,538 features. Two technical libraries per condition receive equal weight in the primary analysis; pooled-cell weighting is a sensitivity. Whole-organoid cells number 8,603 control and 15,431 PA; Hepatocytes 6,283 and 11,814; original broad Cholangiocytes 1,213 and 2,206; SMC-excluded 1,109 and 2,069. These are descriptive technical-library summaries, without biological confidence intervals. Earlier all-25 summaries remain in analysis/matched_scale_areg/outputs/; the new AREG raw calculations are in revision31/areg_analysis/.' if lang == 'en' else
            '本轮均由原始计数重算，均值含实测零，ln(1+10000x/T)的T为全部33,538特征。主要分析每条件两技术库等权，合并细胞为敏感性。对照/PA细胞数依次为全类器官8,603/15,431，Hepatocytes 6,283/11,814，原宽Cholangiocytes 1,213/2,206，去SMC 1,109/2,069。不计算生物学置信区间。旧全25表保留于analysis/matched_scale_areg/outputs/，新AREG raw计算见revision31/areg_analysis/。')
    text = section(text, head, head + '\n\n' + note + '\n\n' + table(['Population', 'Weighting', 'Control mean', 'PA mean', 'PA minus control'] if lang == 'en' else ['群体', '权重', '对照均值', 'PA均值', 'PA减对照'], rows))

    head12 = '## Supplementary Table S12. Contrast directions and three-group ordering of all 15 primary orthologues' if lang == 'en' else '## 补充表S12. 全部15主要同源候选的方向与三组排序'
    intro12 = ('The fixed 15 orthologues are displayed without reselection. Human units are median-ratio normalized counts; mouse units are TPM and are not interchangeable. Directions use human pseudocount 1 and mouse pseudocount 0.1. “Between” requires the third-group mean to fall strictly between the other two. Mouse labels remain unconfirmed exposure assignments; arrows are numerical signs, not significance or treatment validation. All human candidates have intermediate Treat ordering by construction. VNN1 does not reproduce that ordering in mouse labels: sess is below both ctrl and HFD.' if lang == 'en' else
               '固定15个同源候选全部列出，不重新筛选。人源为中位数比值归一化计数，鼠为TPM，不可互换。方向分别加人1、鼠0.1伪计数；“中间”要求第三组均值严格位于另两组之间。鼠标签的暴露对应仍未确认，箭头不是显著性或处理验证。人源Treat中间排序来自筛选规则；VNN1鼠sess低于ctrl及HFD，没有重现中间排序。')
    hrows = [[r.human_symbol, f'{r.control_mean:.3f}', f'{r.model_mean:.3f}', f'{r.treat_mean:.3f}', sign(r.model_control_sign) + '/' + sign(r.treat_model_sign), yes(r.treat_between, lang)] for r in cand.itertuples()]
    mrows = [[r.human_symbol, f'{r.ctrl_mean_TPM:.6f}', f'{r.HFD_mean_TPM:.6f}', f'{r.sess_mean_TPM:.6f}', sign(r.HFD_ctrl_sign) + '/' + sign(r.sess_HFD_sign), yes(r.sess_between, lang), yes(r.all_libraries_TPM_lt1, lang)] for r in cand.itertuples()]
    s12 = head12 + '\n\n' + intro12 + '\n\n### S12a. Human HepG2\n\n' + table(['Gene', 'Control', 'Model', 'Treat', 'M/C; T/M', 'Treat between'] if lang == 'en' else ['基因', 'Control', 'Model', 'Treat', 'M/C；T/M', 'Treat中间'], hrows)
    s12 += '\n\n### S12b. Mouse source labels\n\n' + table(['Human gene', 'ctrl TPM', 'HFD TPM', 'sess TPM', 'HFD/ctrl; sess/HFD', 'sess between', 'All nine <1 TPM'] if lang == 'en' else ['人基因', 'ctrl TPM', 'HFD TPM', 'sess TPM', 'HFD/ctrl；sess/HFD', 'sess中间', '九库均<1 TPM'], mrows)
    head13 = '## Supplementary Table S13. Neighboring human candidate thresholds' if lang == 'en' else '## 补充表S13. 人源候选相邻阈值敏感性'
    s13 = head13 + '\n\n' + ('Post-hoc sensitivity with the other selection conditions and 22/27 stability requirement held fixed. The original 25-candidate reporting set is unchanged.' if lang == 'en' else '事后敏感性固定其他筛选条件及22/27稳定性要求，不改变原25候选报告集合。') + '\n\n'
    s13 += table(['Fold threshold', 'Initial + CPM agreement', 'Final candidates', 'Original 25 retained', 'VNN1 / AREG / JAG2'] if lang == 'en' else ['倍数阈值', '初选且CPM一致', '最终数', '保留原25数', 'VNN1 / AREG / JAG2'],
                 [[f'{r.fold_threshold:.1f}', r.cpm_agree_n, r.selected_n, r.frozen25_retained, '/'.join(yes(x, lang) for x in [r.VNN1, r.AREG, r.JAG2])] for r in threshold.itertuples()])
    head14 = '## Supplementary Table S14. AREG in the three source subtypes of the deposited broad category' if lang == 'en' else '## 补充表S14. 原宽类别三个来源亚型的逐库AREG'
    s14 = head14 + '\n\n' + ('C1/C2 are GSM6322577/78 control libraries; P1/P2 are GSM6322579/80 PA libraries. Counts are raw UMI >0; means include zeros on log1p(CP10k). SMC retains the deposited Smooth muscle cells label. All ten subtypes and 40 library/subtype rows are retained in the CSV; the three source subtypes relevant to the broad-category audit are printed here.' if lang == 'en' else 'C1/C2为GSM6322577/78对照，P1/P2为GSM6322579/80 PA库。非零为raw UMI>0，均值为含零log1p(CP10k)。SMC保留原Smooth muscle cells标签。CSV保留全部10亚型×4库共40行，此表打印与宽类别核对有关的三原亚型。') + '\n\n'
    selected = subtypes[subtypes.source_subtype_label.isin(['Cholangiocytes', 'Ductal cells', 'Smooth muscle cells'])]
    s14 += table(['Library', 'Source subtype', 'Cells', 'Raw nonzero', 'Detection %', 'Mean log1p(CP10k)'] if lang == 'en' else ['文库', '原亚型', '细胞数', 'raw非零数', '检出%', '含零均值'],
                 [[lib_names[r.library_accession], r.source_subtype_label.replace('Smooth muscle cells', 'SMC'), r.n_cells, int(r.n_raw_nonzero), f'{r.raw_detection_pct:.2f}', f'{r.mean_log1p_CP10k:.6f}'] for r in selected.itertuples()])
    head15 = '## Supplementary Table S15. Original broad category and exclusion of Smooth muscle cells' if lang == 'en' else '## 补充表S15. 原宽类别与排除SMC的敏感性'
    s15 = head15 + '\n\n' + ('The post-hoc exclusion keeps source-labelled Cholangiocytes + Ductal cells, without reannotating identities. Both versions and all four technical libraries remain. No P values or biological confidence intervals are calculated; Table S11 retains both negative within-population PA-minus-control differences.' if lang == 'en' else '事后排除只保留原Cholangiocytes+Ductal cells，不重新注释身份。原宽类别及排除版的四库全部保留，不计算P值或生物置信区间；表S11同时保留这两群PA减对照的负向差值。') + '\n\n'
    s15 += table(['Library', 'Group', 'Cells', 'Raw nonzero', 'Detection %', 'Mean log1p(CP10k)'] if lang == 'en' else ['文库', '组别', '细胞数', 'raw非零数', '检出%', '含零均值'],
                 [[lib_names[r.library_accession], ('Original broad' if r.population.startswith('Deposited') else 'Excluding SMC') if lang == 'en' else ('原宽类别' if r.population.startswith('Deposited') else '去SMC'), r.n_cells, int(r.n_raw_nonzero), f'{r.raw_detection_pct:.2f}', f'{r.mean_log1p_CP10k:.6f}'] for r in exclude.itertuples()])
    marker = '## Source data and reproducibility scope' if lang == 'en' else '## 源数据与复现范围'
    text = text.replace(marker, '\n\n'.join([s12, s13, s14, s15]) + '\n\n' + marker, 1)
    if lang == 'en':
        en = text
    else:
        zh = text

en = en.replace('## Source data and reproducibility scope', '''## Source data and reproducibility scope

The revision31 folder of Online Resource 2 contains the new 77-record metadata audit, count-based candidate/ordering reanalysis, deposited-signal rat rerun, and raw-source AREG audit. Protocols explicitly disclose prior knowledge of earlier results. Human/rat/AREG reruns and retained legacy displays are distinguished in each module; none is a new wet-laboratory experiment. Source files were read-only. The original independent-directory reproduction record remains historical evidence of its narrower executed scope, rather than a claim that all subsequent modules were included in that earlier test.

The original dietary study's Table 2 [1] provides physiological context: control versus sesamin liver triacylglycerol was 49.1 ± 3.1 versus 25.8 ± 2.0 micromol/g, and serum triacylglycerol 396 ± 74 versus 198 ± 30 micromol/dL (means ± SE, seven animals per group; distinct original comparison superscripts, P <0.05). These are source-paper observations, not new phenotype measurements. Only five selected animals per group entered the arrays, and no individual phenotype-to-array mapping was available. No gene–phenotype correlation was calculated from group means, and a dietary percentage was not converted to human sesame intake. The checked source and measurement distinctions are retained in revision31/metadata/nutrition_context_evidence.md.
''', 1)
zh = zh.replace('## 源数据与复现范围', '''## 源数据与复现范围

Online Resource 2的revision31目录包含新77条元数据审核、计数层面候选/排序复算、大鼠存档signal重算及AREG原始来源核验。所有方案明确此前结果已知，区分本轮已执行、人鼠旧记录及保留图；不是新增湿实验。原始来源只读，旧独立目录复现记录仍只证明其原已执行范围，不能扩大成已验证后续所有模块。

原膳食论文表2[1]提供生理背景：对照与芝麻素组肝三酰甘油49.1±3.1与25.8±2.0 μmol/g，血清三酰甘油396±74与198±30 μmol/dL，均值±SE、每组7动物，原文不同上标表示P<0.05。这是原作者观察，不是本研究新测表型；转录组每组只含选取的5动物，缺逐动物表型—阵列对应。未用组均值拼接基因—表型相关性，也未换算成人芝麻摄入量。核验来源与尺度见revision31/metadata/nutrition_context_evidence.md。
''', 1)

# Preserve S1-S14 image content; add relocated old main displays without changing them.
cap15en = ('Whole-transcriptome sample structure. a, Nine human libraries and 13,014 filtered genes, log2(median-ratio count + 1). b, Nine mouse libraries and 15,340 filtered genes, log2(TPM + 1). c, Twenty rat arrays and all 31,099 probes, log2(signal + 1). Each point is one submitted library or array. Archived coordinates and explained variance are reused without refitting; features were centered without variance scaling. Panels have separate axes and are not a pooled cross-species embedding. Mouse source labels do not confirm intervention assignments. Rat episesamin and sesamolin are not additional sesamin replicates. The older top-500-gene human PCA in Fig. S1 uses a different feature basis.')
cap16en = ('Candidate directions under overlapping library configurations. a, Full-data mouse label log2 ratios and minimum–maximum ranges across 64 configurations, using a 0.1 TPM pseudocount; ranges are not confidence intervals. b, Configurations retaining the full-data numerical sign for each species and contrast, each out of 64. Labels 64 (0) mark constant-zero states, not nonzero responses. c, Nonzero agreement of both within-study label contrasts across 4096 overlapping human–mouse configuration pairs. Each group retains two or three libraries; human normalization is re-estimated per configuration and mouse TPM is unchanged. Asterisks denote <1 TPM in all mouse libraries. Source-label exposure assignments remain unresolved. These descriptive fractions are neither validation rates nor independent replication and do not reapply the full original candidate-selection rule.')
cap17en = ('Organoid cellular context on the deposited UMAP. a, Ten source subtype labels. b, AREG log1p(CP10k) on the same coordinates and retained cells. All 24,034 cells from four control/PA technical libraries are shown. The embedding is a subset of the author\'s full atlas, not a newly fitted UMAP. Per-cell normalization uses all measured genes; zeros are retained. HSC, hepatic stellate cells; SMC, smooth muscle cells. AH, FH1, FH2, HB1 and HB2 retain author labels without new identity assignments. The image is retained from the previous manuscript; the raw-source AREG recomputation independently reproduced its supporting expression summaries. This dataset contains no sesamin intervention.')
cap15zh = '全转录组样本结构。a，人源9库、13,014过滤基因，log2(中位数比值计数+1)；b，小鼠9库、15,340过滤基因，log2(TPM+1)；c，大鼠20阵列、31,099探针，log2(signal+1)。每点一文库或阵列，沿用存档坐标及解释率，不重拟合；特征中心化但不标准化方差。三面板分别建模，不是跨物种合并嵌入。鼠原标签不等于已确认的暴露赋值，表芝麻素和芝麻林素不算芝麻素重复。图S1旧人源最高变异500基因PCA使用不同特征基底。'
cap16zh = '重叠文库配置下的候选方向。a，鼠原标签全数据log2比值及64配置最小—最大范围，伪计数0.1 TPM，范围不是置信区间。b，各物种各对比保留全数据数值符号的配置数，每项共64；64(0)代表恒零，不是非零反应。c，4096个重叠人鼠配置配对中，两项标签对比同时非零同向的计数。每组保留2或3库，人源因子逐配置重估，鼠TPM不变。星号表示鼠九库均<1 TPM。鼠标签的暴露对应仍未解决；这些比例不是验证率、独立复制，也不重新执行完整候选筛选规则。'
cap17zh = '原作者UMAP上的类器官细胞背景。a，10个原亚型标签；b，同坐标同细胞的AREG log1p(CP10k)。展示四份Control/PA技术库全部24,034细胞，是原完整图谱坐标子集，不是新拟合UMAP。逐细胞分母使用全部测量基因并保留零。HSC为肝星状细胞，SMC为平滑肌细胞；AH、FH1、FH2、HB1、HB2不扩展为新身份。本图沿用旧稿，本轮从原计数复算支持表达摘要。没有芝麻素干预。'
for num, ce, cz in [(15, cap15en, cap15zh), (16, cap16en, cap16zh), (17, cap17en, cap17zh)]:
    en += f'\n\n<!-- pagebreak -->\n\n### Figure S{num}\n\n![Figure S{num}](../figures/supplement/Figure_S{num:02}.png)\n\nFigure S{num}. **' + ce.split('. ', 1)[0] + '.** ' + ce.split('. ', 1)[1] + '\n'
    zh += f'\n\n<!-- pagebreak -->\n\n### 图S{num}\n\n![图S{num}](../figures/supplement/Figure_S{num:02}.png)\n\n图S{num}。' + cz + '\n'

for lang, text in [('en', en), ('zh', zh)]:
    (OUT / f'supplement_{lang}.md').write_text(text, encoding='utf-8')

manifest = {'source_stage': str(OLD), 'new_table_sources': {}, 'output_files': {},
            'retained_old_SI_figures': list(range(1, 15)), 'relocated_main_figures': {'2': 'S15', '5': 'S16', '6': 'S17'},
            'scope': 'Supplement text and source table copies only; main text and images untouched.'}
for name, p in sources.items():
    manifest['new_table_sources'][name] = {'path': str(p), 'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
for lang in ['en', 'zh']:
    p = OUT / f'supplement_{lang}.md'
    manifest['output_files'][p.name] = {'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest(), 'characters': len(p.read_text('utf-8'))}
(HERE / 'supplement_revision_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps(manifest['output_files'], indent=2))
