# SHANK3 ASO Mechanism-of-Action Bioinformatics Workflow
#
# Reproducible R / targets pipeline for chemistry-aware ASO target mapping,
# RNA structure, RBP evidence, TF/chromatin (conditional), and mechanism ranking.
#
# Species: Homo sapiens | Genome: GRCh38 | Gene: SHANK3
# Chemistry: PS backbone, 2'-MOE, all C = 5m-dC; architecture unknown;
#            RNase-H compatibility unknown. Do NOT assume gapmer / RNase-H.

## Project layout

```
shank3-aso-mechanism/
├── _targets.R
├── _config.yml
├── README.md
├── renv.lock
├── R/                         # analysis modules
├── report/                    # RMarkdown reports
├── data/raw/                  # cached downloads + provenance
├── data/processed/
├── results/
├── figures/
└── logs/
```

## Requirements

- R >= 4.3
- Optional: ViennaRNA (`RNAfold`, `RNAplfold`) — structure step records commands and continues if missing
- Bioconductor packages: Biostrings, GenomicRanges, IRanges, rtracklayer, GenomicFeatures,
  Rsamtools, BiocFileCache, AnnotationHub, biomaRt (and optionally BSgenome.Hsapiens.UCSC.hg38)

## Install

```bash
cd shank3-aso-mechanism
Rscript -e 'install.packages("renv"); renv::restore()'
# or, without renv snapshot:
Rscript -e 'source("R/00_bootstrap.R"); bootstrap_packages()'
```

ViennaRNA (Debian/Ubuntu):

```bash
sudo apt-get install -y vienna-rna
# or build from https://www.tbi.univie.ac.at/RNA/
```

## Run

```bash
cd shank3-aso-mechanism
Rscript -e 'targets::tar_make()'
```

Render final report after pipeline completion:

```bash
Rscript -e 'rmarkdown::render("report/FINAL_SHANK3_ASO_MECHANISM_REPORT.Rmd",
  output_format = c("html_document", "pdf_document"),
  output_dir = "report")'
```

## Ranking knobs

In `_config.yml`:

```yaml
ranking:
  top_n_rbps: 15   # RBPs flagged in_top_report / shown in summaries (default 15)
```

## Hard constraints (do not relax)

- Analysis / statistics / plotting / reporting in **R only**
- ViennaRNA only via `system2()`
- Do **not** assume ASOs are gapmers or that RNase-H is active
- Do **not** infer a confirmed mechanism from sequence/database evidence alone
- GRCh38/hg38 only — never mix genome builds
- If no exact SHANK3 transcript match: use the ambiguity statement in `_config.yml`
- If ViennaRNA is missing: warn, record command, continue — never fabricate structure results
- Unreachable databases → empty/partial tables with `status` / `warning` columns

## Outputs (key)

| Step | Primary outputs |
|------|-----------------|
| 01 Sequence QC | `results/01_sequence_qc.tsv`, `results/01_aso_sequences.fasta` |
| 02 Mapping | `results/02_shank3_mapping.tsv`, windows fasta/bed |
| 03 Structure | `results/03_*.tsv`, `figures/03_structure/` |
| 04 RBP | `results/04_rbp_*.tsv` |
| 05 TF/chromatin | `results/05_*.tsv` / `.bed` (conditional) |
| 06 RBP rank | `results/06_ranked_rbps.tsv` |
| 07 Mechanisms | `results/07_mechanism_hypotheses.tsv` + final report |

## Chemistry note

Revised chemistry: core **MOE** (all bases in MOE brackets), phosphorothioate backbone,
architecture **likely_uniform_MOE**, gapmer **not_supported_by_current_supplier_notation**,
RNase-H direct SHANK3 mRNA **low_support**, 5m-dC **reported_but_needs_confirmation**,
FAM present in supplier construct but excluded from unlabeled mechanism interpretation.
PS linkage pattern awaits supplier notation (map not invented).
Phenotype for both ASOs is **increased SHANK3 expression**.

## Citation / provenance

Every downloaded resource is recorded under `data/raw/PROVENANCE.tsv` with annotation release,
database version, URL, and download date.
