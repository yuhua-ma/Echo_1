# SHANK3 ASO 结果阅读指南

面向：Echo Ma · 读 `shank3-aso-mechanism/results/` + 最终报告  
范围：ASO_4_5_2_4、ASO_4_5_2_6 · GRCh38 · SHANK3  
来源：PR [#1](https://github.com/yuhua-ma/Echo_1/pull/1) · `cursor/shank3-aso-mechanism-7245`

---

## 推荐阅读顺序（先看这个）

1. **`results/PIPELINE_SUMMARY.txt`** — 一页总览（匹配数、化学假设、top RBP、top 机制）
2. **`01_sequence_qc.tsv`** — 序列是否干净 + 化学标签（`likely_uniform_MOE` / RNase-H）
3. **`02_shank3_mapping.tsv`** — ASO 是否真的对上 SHANK3 transcript
4. **`04_rbp_motif_overlap.tsv` + `04_rbp_clip_overlap.tsv`** — RBP 证据强弱（motif ≠ 蛋白结合）
5. **`06_ranked_rbps.tsv`** — top-15 候选（启发式排序）
6. **`07_mechanism_hypotheses.tsv`** — 机制假说排名（不是结论）
7. **`03_rna_structure.tsv` / `03_accessibility.tsv`** — 结构上下文（近似值）
8. **`07_validation_recommendations.tsv`** — 下一步实验清单

报告入口：`report/FINAL_SHANK3_ASO_MECHANISM_REPORT.html`（或 `.pdf`）

---

## 先记住的三件事

- **表里的分数 / 排名 = 假说优先级**，不是实验证明。
- **Motif hit ≠ ASO 作为配体结合了那个 RBP 蛋白**；只表示靶 RNA 窗口里有该 RBP 偏好的序列模式。
- **当前 CLIP 表几乎是空的**：只有 ENCODE 元数据查询，**没有**位点解析的 eCLIP peak 重叠（未编造 peak）。

---

## 细胞类型怎么理解（neuron + neural_precursor_cell）

配置（`_config.yml`）只限定两种细胞语境：

- `neuron`
- `neural_precursor_cell`（NPC）

**它实际管什么**

- 表达先验（`04_rbp_expression.tsv`）：只在这两种细胞类型上标 `expressed=TRUE/FALSE`（二元 prior，不是定量 RNAseq）。
- 功能/排名语境：优先考虑“在神经元/NPC 里讲得通”的 RBP。

**它不证明什么**

- **不是** neuron-specific CLIP 证据。
- CLIP 表里的 `cell_type` 只是配置字符串（`neuron,neural_precursor_cell`），**不等于**已在这些细胞测到 peak。
- 没有你自己实验体系（原代神经元 / iPSC-NPC 等）的定量表达。

一句话：细胞类型缩小了“先验候选池”，**没有**把证据升级成神经元位点结合。

---

## 化学标签速查（贯穿全表）

| 字段 / 概念 | 当前取值 | 读法 |
|---|---|---|
| `chemistry_architecture` | `likely_uniform_MOE` | 更像整链 MOE 占位/位阻，不像经典 gapmer |
| `chemistry_gapmer` | `not_supported_by_current_supplier_notation` | 现有供应商记号**不支持** gapmer 假设 |
| `rnase_h_direct_SHANK3_mRNA` / `rnase_h_compatible` | `low_support` | 直接 RNase-H 切 SHANK3 mRNA **支持度低** |
| backbone | phosphorothioate (PS) | 蛋白/RBP 接触更说得通；**PS 位点图谱尚未从供应商字符串解析** |
| 5m-dC | `reported_but_needs_confirmation` | 有报告，未实验确认 |
| FAM | 有，但 `chemistry_fam_in_unlabeled_mechanism=FALSE` | 解释无标记机制时**排除 FAM** |

表型假设（排名用）：**increased SHANK3 expression**（表达升高）→ 与“RNase-H 敲低”方向不一致。

---

## Step 01 — 序列 QC

**在做什么**  
清洗 ASO 碱基、算 GC/发夹风险，并把化学注释写进表。

**输入**  
`_config.yml` 里的 ASO 序列 + chemistry 块。

**打开这些文件**

- `results/01_sequence_qc.tsv`
- `results/01_aso_sequences.fasta`

**怎么读主要列**

- `cleaned_sequence` / `reverse_complement` / `rna_compatible_sequence`：后续比对与结构用的规范序列
- `gc_percent`、`homopolymers`、`hairpin_potential`、`self_complementarity_score`、`qc_warnings`：序列风险旗标
- `chemistry_*`：见上表；`status=ok` 只表示碱基合法

**不证明**

- 不证明递送、稳定性、体内活性
- 不证明 5m-dC / PS map 已确认（看 `qc_warnings`）

---

## Step 02 — Transcript / 基因组定位

**在做什么**  
用 ASO 的 **reverse complement** 在 SHANK3 transcript / 基因座上找匹配；顺带做 **chr22** 轻量 off-target 扫描。

**输入**  
Step 01 清洗序列；GRCh38（TxDb / BSgenome / GENCODE 缓存）。

**打开这些文件**

- `results/02_shank3_mapping.tsv`（主表）
- `results/02_target_windows.fasta` / `.bed`（±200 nt 窗口，供结构与 motif）
- `results/02_offtargets_chr22.tsv`

**怎么读主要列**

- `exact_shank3_match` / `match_class`：是否精确 transcript RC 命中（当前两支 ASO 都有 exact 行）
- `transcript_id`、`feature_type`、`chrom/start/end`：位点与 isoform
- `rank` / `rank_score`：同 ASO 多命中时的排序
- off-target 表：仅 **chr22 genomic exact**，不是全转录组

**不证明**

- 精确匹配 ≠ 已确认细胞内 on-target engagement
- chr22 无命中 ≠ 全基因组无 off-target
- 窗口坐标是映射产物，不是实验 peak

---

## Step 03 — ViennaRNA 结构（重点）

**在做什么**  
用 ViennaRNA 估计：ASO 自身折叠、靶窗口 MFE 结构、位点附近是否易单链、近似杂交能量。

**输入**  
`02_target_windows.fasta` + ASO 规范 RNA 序列；需本机有 `RNAfold` / `RNAplfold`（缺失则记录命令、不编造）。

**打开这些文件**

- `results/03_rna_structure.tsv`
- `results/03_accessibility.tsv`
- `figures/03_structure/`（若已生成）
- `logs/vienna_status.txt`

**工具各干什么**

- **RNAfold**：最小自由能（MFE）二级结构 → `structure_dotbracket`、`mfe_kcal_mol`
- **RNAplfold**：局部窗口内的 **unpaired probability**（平均进 `mean_unpaired_probability`）
- **RNAcofold**（若有）：ASO↔靶的近似杂交能量 → `approx_hybridization_energy`

**怎么读主要列**

- `context`：`aso_self` vs `target_window`
- `site_structure_class`：`loop_or_unpaired` / `stem` / `bulge_or_junction` / `mixed` — 位点处 stem-loop 语境的粗分类
- `mean_unpaired_probability`：越高 → 计算上更“开放/可及”（仍是模型）
- `chemistry_caveat`：表里已写明 — **规范 RNA 近似，不是 MOE/PS/5m-dC 的生化实测**

**为什么不是精确生化测量**

- 折叠按 **canonical RNA** 热力学参数
- **MOE、PS、5m-dC 未建模**（能量、杂交都会偏）
- 无 DMS/SHAPE ± ASO 实验约束

**不证明**

- 不证明细胞里真实结构
- 不证明 ASO 结合后一定 remodel
- 不证明 accessibility 高就一定好结合

---

## Step 04 — RBP 证据（CLIP / motif / 表达 / 功能）

### 4A · CLIP 重叠（重点）

**在做什么**  
试图把 ASO 靶位点与公开 eCLIP / POSTAR 等 peak 做基因组重叠；**做不到就不编造**。

**输入**  
映射位点坐标；ENCODE eCLIP 搜索 API；POSTAR3 等（多数仅 provenance）。

**打开**  
`results/04_rbp_clip_overlap.tsv`

**`overlap_0nt / 10nt / 30nt` 是什么意思**

| 列 | 含义（若为 TRUE） |
|---|---|
| `overlap_0nt` | peak 与位点区间直接相交（最强 CLIP 层） |
| `overlap_10nt` | 在位点外扩 10 nt 内碰到 peak |
| `overlap_30nt` | 在位点外扩 30 nt 内碰到 peak |
| `nearest_distance_nt` | 最近 peak 距离（无重叠时） |

证据强度规则（代码注释）：strong≈0 nt 直接重叠；moderate≈10–30 nt 或 motif+表达；weak≈仅距离/仅 motif/仅表达。

**当前运行的真实状态（请按表读）**

- `status` = `metadata_only_no_peak_overlap`
- `peak_start/peak_end` = NA；`rbp` = NA
- 三个 overlap 列均为 **FALSE**；`evidence_strength` = **none**
- 含义：**查过 ENCODE 元数据，未自动下载并比对 SHANK3 位点的 peak BED**；POSTAR3 无批量 peak dump → **重叠留空，不是“无 RBP”，而是“本流水线未解析到 peak”**

**数据库 vs 真实 peak**

- 有：ENCODE 搜索 JSON 缓存（元数据）
- 无：位点解析的 eCLIP peak 重叠表
- 因此：**没有 neuron-specific CLIP 命中可引用**

**不证明**

- 不证明位点上没有 RBP（只是没测到/没下到 peak）
- 不证明某个 RBP 在该位点结合

---

### 4B · Motif 重叠（关键 — 易误读）

**在做什么**  
在 **靶 RNA 窗口（含 ASO footprint 区域）** 里扫描公开共识 motif（ATtRACT/RBPDB 风格、经策展的正则）。

**打开**  
`results/04_rbp_motif_overlap.tsv`

**怎么读**

- `window_id` / `match_start` / `match_end`：命中在 **RNA 窗口序列**上的坐标
- `motif_sequence_or_regex` / `motif_id`：用的是哪条共识模式
- `evidence_strength`：当前多为 **`weak`**
- `warning`：共识近似；**该位点未经 CLIP 验证**

**必须分清的两句话**

1. Motif hit = 靶 RNA / ASO footprint 区域的序列，**符合**某 RBP 已知序列偏好。  
2. Motif hit **≠** “ASO 2.4 这条化学分子作为蛋白配体结合了该 RBP”。

**不证明**

- 不证明 RBP 蛋白真的占着这个位点
- 不证明 ASO–RBP 直接蛋白结合
- **仅有 motif = 弱证据**，需 CLIP / RIP / 功能遗传学升级

---

### 4C · 表达先验 & 功能注释

**打开**

- `results/04_rbp_expression.tsv`
- `results/04_rbp_function.tsv`

**怎么读**

- 表达：`expression_source=curated_neuronal_prior_not_quant_RNAseq`；`expression_units=binary_prior`
- 功能：`regulatory_direction_hint`（repressive / activating_stabilizing / …）、`rna_process`
- 细胞类型列：仅 neuron / NPC 两行 prior

**不证明**

- 不是 ENCODE/GTEx 定量
- 功能标签不是 SHANK3 位点特异性证明

---

## Step 05 — TF / 染色质（条件步骤）

**在做什么**  
仅当靶点被判为核内调控 / enhancer / promoter / antisense / intron 调控语境时，才扫 DNA 上的 TF motif。

**当前结果**

- `results/05_tf_chromatin.tsv`：`status=skipped_not_regulatory_context`
- 主命中是 **transcript exon spliced** → **跳过**（不做 RNA 上的 TF 推断）

**不证明**  
跳过 ≠ 转录调控不可能；只是本流水线按 locus class **未启用**该层。

---

## Step 06 — Top-15 RBP 排名

**在做什么**  
把 motif / 表达 prior / 功能方向 /（若有）CLIP 合成可解释启发式分数；标记 `in_top_report`（默认 top **15**）。

**输入**  
Step 04 四张表（+ 结构可选用）。

**打开**  
`results/06_ranked_rbps.tsv`

**怎么读主要列**

- `score` + `score_rationale`：加分明细（透明启发式）
- `has_clip_0nt` / `has_motif` / `has_expression_prior`：证据组成（当前多为 **无 CLIP**）
- `rbp_class`、`regulatory_direction_hint`：和“表达升高”表型怎么对齐的假设
- `rank`、`in_top_report`：报告里重点看的 15 个
- `warning`：明确写着 **不是实验验证**

**计分直觉（简化）**

- CLIP 0 nt 或 motif+表达 prior → 更高
- 仅表达 prior / 结构相关 → 较低
- 与表型方向矛盾时有软惩罚

**当前两支 ASO 的 top 模式（摘要）**  
PTBP2、RBFOX1、RBFOX3、SRSF1、AGO2…（完整以表为准；两边高度相似）

**不证明**

- 高分 ≠ 机制确认
- top-15 ≠ 这 15 个都在位点上结合
- 很多高分主要靠 **表达 prior ± 功能标签 ± motif**，不是 CLIP

---

## Step 07 — 机制假说 + 验证建议

**在做什么**  
结合化学（`likely_uniform_MOE`、RNase-H `low_support`、PS）与 top RBP，给机制标签打假说分；并列出验证实验。

**打开**

- `results/07_mechanism_hypotheses.tsv`
- `results/07_validation_recommendations.tsv`
- `results/PIPELINE_SUMMARY.txt`
- `report/FINAL_SHANK3_ASO_MECHANISM_REPORT.*`

**怎么读机制表**

- `mechanism_label` + `hypothesis_score` + `rank` + `is_top`
- 当前 top：**RBP occlusion**（位阻/占据 repressor 位点 → 与表达升高相容的一种假说）
- `rnase_h_status=low_support`；`rnase_h_supported=FALSE`
- `rationale` / `caveat`：为什么加分 + 明确 “Hypothesis only”

**化学如何拉动排名（直觉）**

- uniform MOE → 更偏向 steric occupancy（occlusion / displacement），**压低** RNase-H cleavage
- 表型 = 表达升高 → 与 on-target RNase-H 敲低不一致
- PS → 蛋白接触/RBP 干扰更说得通（但仍缺 PS map）

**验证表怎么用**  
按 top RBP 展开实验菜单（CLIP/RIP、KD/OE、EMSA、SHAPE/DMS、RNase-H assay 等）。  
把它当 **待做清单**，不是已完成结果。

**不证明**

- 任何 `mechanism_label` 都不是已证实 MoA
- 数据库 + 启发式 ≠ 药理结论

---

## 一页诚实差距清单

当前流水线**有**：

- SHANK3 transcript 精确 RC 匹配与 ±200 nt 窗口
- ViennaRNA 规范 RNA 结构近似
- 策展 motif 扫描 + 神经元/NPC 表达先验 + 功能标签
- 化学感知的机制假说排序与验证建议

当前流水线**没有 / 未完成**：

- 位点解析的 neuron/NPC eCLIP peak 重叠（仅 metadata）
- 定量 RNAseq 表达、细胞内 on-target 确认
- 带真实 MOE/PS/5m-dC 的结构/杂交测量
- 供应商 PS linkage 字符串解析出的位点图谱
- 5m-dC 实验确认
- 全基因组 / 全转录组 off-target

---

## 文件速查

| Step | 主要打开 |
|---|---|
| 总览 | `results/PIPELINE_SUMMARY.txt` |
| 01 | `01_sequence_qc.tsv` |
| 02 | `02_shank3_mapping.tsv`、`02_offtargets_chr22.tsv` |
| 03 | `03_rna_structure.tsv`、`03_accessibility.tsv` |
| 04 | `04_rbp_clip_overlap.tsv`、`04_rbp_motif_overlap.tsv`、`04_rbp_expression.tsv`、`04_rbp_function.tsv` |
| 05 | `05_tf_chromatin.tsv`（当前 skipped） |
| 06 | `06_ranked_rbps.tsv` |
| 07 | `07_mechanism_hypotheses.tsv`、`07_validation_recommendations.tsv` |
| 报告 | `report/FINAL_SHANK3_ASO_MECHANISM_REPORT.html` |

路径前缀（仓库内）：`shank3-aso-mechanism/`
