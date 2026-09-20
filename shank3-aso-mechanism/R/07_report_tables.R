# R/07_report_tables.R — validation recommendations + end-of-run summary

.validation_experiments <- function() {
  c(
    "ASO + RBP CLIP / eCLIP (or RIP-qPCR) at the mapped site",
    "RBP knockdown / CRISPRi and measure SHANK3 RNA/protein (directionality)",
    "RBP overexpression rescue of ASO phenotype",
    "RNA immunoprecipitation after ASO to test occlusion vs recruitment",
    "Competitive EMSA / binding assays with PS/MOE/5m-dC ASO chemistry",
    "DMS-MaPseq or SHAPE-MaP ± ASO for structure remodeling",
    "RNase-H cleavage assay in vitro with cognate chemistry (architecture test)",
    "Nascent transcription (GRO-seq / TT-seq) if transcriptional mechanism hypothesized",
    "Isoform-resolved RNA-seq / splicing reporter if splicing RBP hypothesized",
    "Polysome profiling / ribosome footprinting if translation regulation hypothesized"
  )
}

build_validation_table <- function(cfg, rbp_rank) {
  suppressPackageStartupMessages(library(dplyr))
  root <- cfg$project_root
  top <- rbp_rank$ranked %>% filter(in_top_report)
  exps <- .validation_experiments()
  rows <- list()
  for (i in seq_len(nrow(top))) {
    rbp <- top$rbp[[i]]
    aso_id <- top$aso_id[[i]]
    hyp <- paste0("ASO engages ", rbp, " pathway (", top$rbp_class[[i]], "; ",
                  top$regulatory_direction_hint[[i]], ")")
    for (e in exps) {
      tests <- dplyr::case_when(
        grepl("CLIP|RIP-qPCR", e) ~ "Tests whether RBP occupies the ASO site and whether ASO alters occupancy (occlusion/displacement)",
        grepl("knockdown|CRISPRi", e) ~ "Tests whether RBP is necessary for basal SHANK3 repression/activation matching ASO phenotype",
        grepl("overexpression", e) ~ "Tests sufficiency / epistasis of RBP relative to ASO",
        grepl("immunoprecipitation after ASO", e) ~ "Distinguishes occlusion vs recruitment of the RBP",
        grepl("EMSA", e) ~ "Tests direct ASO–RBP or ASO–RNA–RBP competition with matched chemistry",
        grepl("DMS|SHAPE", e) ~ "Tests RNA-structure remodeling hypothesis at the site",
        grepl("RNase-H cleavage", e) ~ "Tests whether chemistry/architecture supports RNase-H (expected low_support under likely_uniform_MOE; gapmer not supported by current supplier notation)",
        grepl("Nascent transcription", e) ~ "Tests transcriptional/enhancer regulation hypothesis",
        grepl("Isoform-resolved|splicing", e) ~ "Tests splicing or processing hypothesis",
        grepl("Polysome|ribosome", e) ~ "Tests translation regulation hypothesis",
        TRUE ~ "Maps to mechanism hypothesis testing"
      )
      rows[[length(rows) + 1]] <- data.frame(
        aso_id = aso_id,
        rbp = rbp,
        rbp_rank_score = top$score[[i]],
        hypothesis_context = hyp,
        experiment = e,
        hypothesis_tested = tests,
        stringsAsFactors = FALSE
      )
    }
  }
  df <- bind_rows(rows)
  out <- file.path(root, "results/07_validation_recommendations.tsv")
  write_tsv(df, out)
  list(validation = df, tsv = out)
}

print_pipeline_summary <- function(cfg, qc_result, map_result, rbp_rank, mech_result, struct_result) {
  mapping <- map_result$mapping
  ranked <- rbp_rank$ranked
  hyp <- mech_result$hypotheses

  exact_n <- sum(mapping$exact_shank3_match, na.rm = TRUE)
  approx_n <- sum(mapping$approximate_shank3_match, na.rm = TRUE)
  ot_path <- file.path(cfg$project_root, "results/02_offtargets_chr22.tsv")
  ot_n <- if (file.exists(ot_path)) {
    ot <- tryCatch(readr::read_tsv(ot_path, show_col_types = FALSE), error = function(e) data.frame())
    nrow(ot)
  } else 0L

  cat("\n========== SHANK3 ASO MoA PIPELINE SUMMARY ==========\n")
  cat("Genome build: GRCh38 only | Gene: SHANK3\n")
  chem <- chemistry_annotation(cfg)
  cat("Chemistry: core=", chem$core_chemistry,
      " | architecture=", chem$architecture,
      " | gapmer=", chem$gapmer,
      " | RNase-H(SHANK3)=", chem$rnase_h_direct_SHANK3_mRNA,
      " | 5m-dC=", chem$five_methyl_dC,
      " | PS_linkage=", chem$ps_linkage_status,
      " | FAM_in_mech=", chem$fam_label_in_unlabeled_mechanism, "\n", sep = "")
  cat("Exact SHANK3 matches (rows): ", exact_n, "\n", sep = "")
  cat("Approximate SHANK3 matches (rows): ", approx_n, "\n", sep = "")
  cat("Transcriptome/genomic off-target hits recorded (chr22 screen): ", ot_n, "\n", sep = "")
  if (exact_n == 0) {
    cat("AMBIGUITY: ", ambiguity_statement(cfg), "\n", sep = "")
  }
  for (aso in cfg$asos) {
    id <- aso$id
    top <- ranked %>% dplyr::filter(aso_id == id, in_top_report) %>% dplyr::arrange(dplyr::desc(score))
    top_mech <- hyp %>% dplyr::filter(aso_id == id, is_top)
    cat("\n--- ", id, " ---\n", sep = "")
    cat("Top RBPs: ", paste(top$rbp, collapse = ", "), "\n", sep = "")
    cat("Top mechanism hypothesis: ", top_mech$mechanism_label[1] %||% "unresolved", "\n", sep = "")
    cat("RNase-H supported/unsupported/unknown: ",
        chem$rnase_h_direct_SHANK3_mRNA,
        " (gapmer=", chem$gapmer, "; likely_uniform_MOE favors steric/occupancy)\n", sep = "")
  }
  cat("\nMissing info for definitive interpretation:\n")
  cat("  - Supplier PS linkage notation string (ps_linkage_pattern awaits parse; map not invented)\n")
  cat("  - Experimental confirmation of 5m-dC (currently reported_but_needs_confirmation)\n")
  cat("  - Experimental confirmation of on-target engagement\n")
  cat("  - Neuron/NPC eCLIP peak overlaps at the precise site\n")
  cat("  - Quantitative expression of candidate RBPs in the treated cell system\n")
  cat("  - Structure probing ± ASO with matched chemistry\n")
  cat("  - Functional genetics (RBP KD/OE epistasis with ASO)\n")
  cat("  - FAM is present in supplier construct but excluded from unlabeled mechanism interpretation\n")
  cat("====================================================\n\n")

  summary_path <- file.path(cfg$project_root, "results/PIPELINE_SUMMARY.txt")
  sink(summary_path)
  cat("Exact SHANK3 match rows: ", exact_n, "\n", sep = "")
  cat("Approximate SHANK3 match rows: ", approx_n, "\n", sep = "")
  cat("Off-target rows (chr22): ", ot_n, "\n", sep = "")
  cat("Chemistry architecture: ", chem$architecture, "\n", sep = "")
  cat("Gapmer: ", chem$gapmer, "\n", sep = "")
  cat("RNase-H direct SHANK3 mRNA: ", chem$rnase_h_direct_SHANK3_mRNA, "\n", sep = "")
  cat("5m-dC: ", chem$five_methyl_dC, "\n", sep = "")
  cat("PS linkage status: ", chem$ps_linkage_status, "\n", sep = "")
  cat("FAM in unlabeled mechanism: ", chem$fam_label_in_unlabeled_mechanism, "\n", sep = "")
  for (aso in cfg$asos) {
    id <- aso$id
    top <- ranked %>% dplyr::filter(aso_id == id, in_top_report) %>% dplyr::arrange(dplyr::desc(score))
    top_mech <- hyp %>% dplyr::filter(aso_id == id, is_top)
    cat(id, " top RBPs: ", paste(top$rbp, collapse = ", "), "\n", sep = "")
    cat(id, " top mechanism: ", top_mech$mechanism_label[1] %||% "unresolved", "\n", sep = "")
  }
  cat("RNase-H: ", chem$rnase_h_direct_SHANK3_mRNA, "\n", sep = "")
  sink()

  invisible(list(exact = exact_n, approx = approx_n, offtargets = ot_n, chemistry = chem))
}

write_report_index <- function(cfg) {
  root <- cfg$project_root
  files <- c(
    "report/final_report.Rmd",
    "report/target_mapping_report.Rmd",
    "report/structure_report.Rmd",
    "report/rbp_evidence_report.Rmd",
    "report/tf_chromatin_report.Rmd",
    "report/FINAL_SHANK3_ASO_MECHANISM_REPORT.Rmd"
  )
  writeLines(files, file.path(root, "results/07_report_index.txt"))
  files
}
