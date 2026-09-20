# PI 预测表 vs 我们的 SHANK3 ASO 结果 — 对照说明

面向：Echo Ma  
状态：**完成**（已读 PI xlsx）— 2026-09-20  
PI 表：`docs/inputs/Shank3_prediction_UTR.xlsx`（原名 `Shank3_prediction_UTR (2).xlsx`）  
我方：[PR #1](https://github.com/yuhua-ma/Echo_1/pull/1) · `cursor/shank3-aso-mechanism-7245` · `shank3-aso-mechanism/results/`

---

## 0. 30 秒结论

- PI 表 **不是 ASO 设计表**，而是 **SHANK3 UTR 上的 RBP motif 目录**（ATtRACT 风格：motif + Off Set + 实验来源 + Exon250/CDS/Intron 分）。
- **没有** ASO ID、没有化学、没有机制假说排名 → 和我们流水线是 **不同问题**，只能比「UTR RBP 先验」这一层。
- **诚实重合（gene 名）**：我方 top-15 里 **7/15** 出现在 PI（`RBFOX1, SRSF1, AGO2, ELAVL2, ELAVL4, FMR1, CELF4`）。
- **关键分歧**：我方 #1 `PTBP2`、#3 `RBFOX3` **不在** PI；PI 顶频是 `PTBP1`（63 行）+ 一堆 hnRNP/SR。
- **位点**：在假定 PI Off Set ≈ 某条 ~1.9 kb 3′UTR（与 `ENST00000262795` 3′UTR 长度接近）时，两支 ASO 的 UTR footprint **附近有大量 PI motif**，但 footprint 上的 PI 基因多为 `PTBP1`/hnRNP/`SRSF*`，**不是**我方 motif 三件套的精确同位点命中。
- 两边都是预测；**都不是实验金标准**。PI = 全 UTR motif 扫描先验；我方 = ASO 窗 + 化学感知机制排序（CLIP 仍空）。

---

## 1. PI 表里有什么？

### Sheet

| Sheet | 规模 | 内容 |
|---|---|---|
| `SHANK3_UTR` | **912** 行 × 13 列 | 主表：每个 RBP×motif×文献/实验 一行 |
| `Frequency_gene-SHANK3_UTR` | 131 个 gene + 计数 | 透视：每个 Gene Name 出现次数 |

### 列（主表）

`Gene Name` · `Gene Id` · `Organism` · `Motif` · `Len` · `Pubmed` · `Experiment` · `Domain` · `Off Set` · `GO terms` · `Exon250` · `CDS` · `Intron`

### 怎么读

- **Gene Name** = 预测的 RBP（约 **120** 个 unique，含 `HNRNPH1*` 这类带 `*` 标签 → 归一化后 ~120）。
- **Motif / Len / Off Set**：在 SHANK3 UTR 序列上的 motif 与位置；`Off Set` 可多值（`;` 分隔）；范围约 **3–1916**（≈ 一条 ~1.9 kb UTR）。
- **Experiment / Pubmed / Domain**：motif 证据来自 RNAcompete / SELEX / X-ray / NMR / RIP-chip 等（**几乎不是**位点 CLIP；CLIP-like 行极少）。
- **Exon250 / CDS / Intron**：数值分（可正可负）— 像区域富集/偏好分数，**不是** ASO 亲和力。
- **没有**：ASO 序列、isoform ENST、基因组 chr 坐标、机制标签、化学。

### 频次 Top（PI 自己的“热度”）

`PTBP1`(63) › `SRSF1`(51–52) › `PPIE`(44) › `HNRNPH1/H2` › `HNRNPF` › `MBNL1` › `HNRNPA1` › `CELF2` › `ELAVL1/2` › `RBFOX1`(11) › `FMR1`(9) …

→ 这是 **UTR motif 目录的密度排序**，不是 ASO MoA 排名。

---

## 2. 我们这边（对照轴）

| 项 | 取值 |
|---|---|
| ASO | `ASO_4_5_2_4` / `ASO_4_5_2_6` |
| 化学 | `likely_uniform_MOE`；RNase-H `low_support` |
| 匹配 | exact SHANK3 transcript RC **13** 行；多为 **3′UTR**（Ensembl 附加粗分） |
| Motif 命中 | **`FMR1`, `RBFOX1`, `RBFOX3`**（weak；窗内共识） |
| Top-15 | `PTBP2, RBFOX1, RBFOX3, SRSF1, AGO2, ELAVL2, ELAVL4, FMR1, IGF2BP1, PUM2, UPF1, CELF4, MOV10, ATXN2, CSTF2` |
| 机制 top | **RBP occlusion** |
| CLIP | `metadata_only_no_peak_overlap` |

---

## 3. Overlap（同意 / 不同意）

机器表：`internal/pi-vs-ours-overlap.tsv`（协调员可选用；用户主看本文）。

### 3.1 靶点 / isoform / UTR

| | PI | 我们 |
|---|---|---|
| 对象 | 整段 SHANK3 **UTR** motif 扫描 | **两条 ASO** 的 transcript RC 位点 ±200 nt |
| ASO | **无** | 有 |
| UTR vs CDS | 表名 + Off Set 暗示 **UTR** | 流水线标 `transcript_exon_spliced`；附加 Ensembl：**多为 3′UTR** |
| 坐标能否直接叠 | 仅 Off Set（相对 UTR） | transcript 坐标；换算成 3′UTR offset 后，**仅在 ~1.9 kb UTR 假设下**可与 PI 比 |

**同意（弱）**：两边都指向 SHANK3 **UTR / 3′UTR 语境**，不是 intron/enhancer 主叙事。  
**不同意 / 不可比**：PI **没有** ASO footprint；不能说“PI 预测了 ASO_4_5_2_4 的靶点”。

在 `ENST00000262795` 上粗算 3′UTR offset（CDS 止于 ~5247 → UTR 长 ~1900，与 PI max 1916 接近）：

- `ASO_4_5_2_4` → UTR **218–235**
- `ASO_4_5_2_6` → UTR **182–199**

该假设下 footprint 上出现的 PI 基因偏 `PTBP1`、hnRNP、`NOVA1`、`SRSF1/9` 等；**`RBFOX1`/`FMR1` 不在 footprint 正中**（`RBFOX1` 的 PI offsets 主要在 **114 / 1043** 一带；`FMR1` 更散、偏远端）。

> 警告：若 PI 用的不是同一条 UTR 异构体，数字重合只是巧合。应用前问清 PI 的 UTR 序列 / ENST。

### 3.2 RBP 名单

**两边都有（top-15 ∩ PI）— 7 个**

`RBFOX1` · `SRSF1` · `AGO2` · `ELAVL2` · `ELAVL4` · `FMR1` · `CELF4`

**我方 top-15 有、PI 无 — 8 个**

`PTBP2` · `RBFOX3` · `IGF2BP1` · `PUM2` · `UPF1` · `MOV10` · `ATXN2` · `CSTF2`

**PI 很热、我方没进 top（或未收录）**

`PTBP1` · `HNRNPH1/H2/H3` · `HNRNPF` · `PPIE` · `MBNL1` · `HNRNPA1` · `CELF2` · `NOVA1` · `RBFOX2`(PI 有，我方未排) …

**家族错位（重要）**

| 我们 | PI | 读法 |
|---|---|---|
| `PTBP2` #1 | `PTBP1` 顶频 | 同家族不同 paralog；神经元常更关心 PTBP2 |
| `RBFOX3` motif | 仅 `RBFOX1`/`RBFOX2` | 我方多报了一个 isoform-family 成员 |
| `PUM2` / `IGF2BP1` | 有 `PUM1` / `IGF2BP3` 等 | 近亲在 PI，精确 ID 不一致 |

**Motif 层**：我方 `FMR1`+`RBFOX1` 在 PI；`RBFOX3` 不在。且 **motif 字符串并不相同**（我方 RBFOX=`GCATG`；PI=`UGCAUG`/`AGCAUG`… — 同家族；FMR1 两边 motif 差得更远）。

### 3.3 机制类别

| | PI | 我们 |
|---|---|---|
| 机制 | **无**（只到 “UTR 上可能有这些 RBP motif”） | **RBP occlusion** top；RNase-H 垫底 |
| 化学 | **无** | `likely_uniform_MOE` 拉动 steric |

→ **谈不上机制同意/不同意**；PI 停在结合位点先验，我们多走了化学 + 表型兼容排序。

---

## 4. 差异怎么解释（别过度解读）

1. **问题不同**：PI = 全 UTR RBP motif 目录；我们 = **两条 ASO** 的 mapping + 局部窗 + 机制假说。
2. **化学**：我们压低 RNase-H、抬高 occlusion；PI 完全不编码化学 → 机制层只能听我们（仍是假说）。
3. **细胞类型**：我们用 neuron/NPC **表达 prior** 抬 `PTBP2` 等；PI 不按细胞过滤 → 更像泛组织 motif 库。
4. **motif ≠ CLIP**：两边主要都是 motif/生化 SELEX 类证据；我们 CLIP 表仍空；PI 也几乎无位点 CLIP。
5. **坐标**：只有在「同一条 ~1.9 kb 3′UTR」假设下，Off Set 才能和 ASO footprint 比；否则只比 gene 名。

---

## 5. 实操建议（信谁 / 重查 / 问 PI）

### 更该信什么（下一步实验）

- **信我方**：ASO → SHANK3 transcript 的 **exact RC 坐标**（序列事实）；化学方向（别按 gapmer 敲低排主线）。
- **信 PI（作为 prior）**：SHANK3 UTR 上 **高密度** 的 `PTBP1`/SR/hnRNP/CELF/ELAVL/RBFOX1 等 motif 生态 — 说明 UTR 很“RBP 忙”，和 occlusion 故事 **相容但不证明**。
- **两边交叉更优先跟的 RBP**：`RBFOX1`、`SRSF1`、`ELAVL2/4`、`FMR1`、`CELF4`、`AGO2`（名称层重合）。
- **单独追问的家族**：`PTBP2`(我们) vs `PTBP1`(PI)；补做 **两者** 或查神经元 isoform。
- **不要信**：任一分数 = MoA；PI 频次 ≠ ASO 靶点重要性。

### 建议重查

1. PI 的 UTR 序列 / ENST / 基因组版本是否 = GRCh38 · 我们用的 isoform。  
2. 把 PI Off Set 映射到基因组，与我方 `02_shank3_mapping.tsv` 做 0/10/30 nt 重叠。  
3. 补我方 **locus eCLIP peak**（关掉 metadata-only）。  
4. 验证实验：KD/OE 或 RIP 优先交叉名单；并行 PTBP1/2。

### 建议直接问 PI

- 这条 UTR 对应哪个 **ENST**？Off Set=1 是 3′UTR 第一个碱基吗？  
- 导出工具是否就是 **ATtRACT**（或同类）？有无过滤阈值？  
- `Exon250/CDS/Intron` 分数定义？  
- 是否做过 / 见过 **SHANK3 位点 CLIP**？  
- 是否知道我们这两条 ASO 的序列与 **uniform MOE** 化学？

---

## 6. 文件指针

| 文件 | 路径 |
|---|---|
| 本对照 | `docs/pi-vs-ours-comparison.md` |
| PI 表（已存） | `docs/inputs/Shank3_prediction_UTR.xlsx` |
| Overlap TSV | `internal/pi-vs-ours-overlap.tsv` |
| 我方阅读指南 | `docs/results-reading-guide.md` |
| 我方结果 | `shank3-aso-mechanism/results/`（PR 分支） |

---

*PI = expert prior（UTR motif 目录），不是 ground truth。我方 = ASO 中心假说生成器。交叉名单可开实验；全表“一致”不成立。*
