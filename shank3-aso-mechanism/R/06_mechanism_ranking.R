# R/06_mechanism_ranking.R — RBP prioritization + chemistry-aware mechanism ranking
#
# RBP scoring rules (transparent):
#   +3  motif hit overlapping target window AND neuronal/NPC expression prior
#   +2  motif hit in window OR strong functional match to increased-expression phenotype
#   +1  neuronal/NPC expression prior only OR structure-related function
#   -2  contradictory direction vs increased SHANK3 expression phenotype
#
# Mechanism labels:
#   RBP occlusion, RBP displacement, RNA-structure remodeling, RNA stabilization,
#   RNase-H-mediated degradation, transcriptional/enhancer regulation,
#   splicing or processing, translation regulation, off-target mechanism, unresolved
#
# Chemistry/mechanism rules (revised supplier-informed chemistry):
#   1. Architecture likely_uniform_MOE; gapmer not_supported_by_current_supplier_notation
#   2. rnase_h_direct_SHANK3_mRNA = low_support → RNase-H cleavage strongly down-weighted
#   3. Phenotype = increased SHANK3 → on-target RNase-H knockdown directionally inconsistent
#   4. PS backbone → protein binding / RBP occlusion–displacement plausible
#      (PS map not invented when ps_linkage_pattern awaits supplier notation)
#   5. Uniform MOE / all bases in MOE brackets → steric/occupancy >> RNase-H
#   6. 5m-dC reported_but_needs_confirmation → caveat only; do not treat as confirmed
#   7. Structure remodeling remains hypothetical without experimental DMS/SHAPE + ASO
#   8. Sequence/database evidence alone ≠ confirmed mechanism
#   9. FAM present in supplier construct but excluded from unlabeled mechanism interpretation

rank_rbps <- function(cfg, rbp_result, struct_result = NULL) {
  suppressPackageStartupMessages(library(dplyr))
  root <- cfg$project_root
  motif <- rbp_result$motif
  expr <- rbp_result$expression
  func <- rbp_result$function_annot
  clip <- rbp_result$clip

  aso_ids <- unique(c(expr$aso_id, motif$aso_id, func$aso_id))
  rows <- list()

  for (aso_id in aso_ids) {
    rbps <- unique(c(func$rbp[func$aso_id == aso_id], motif$rbp[motif$aso_id == aso_id]))
    rbps <- rbps[!is.na(rbps)]
    for (rbp in rbps) {
      score <- 0
      reasons <- c()
      has_motif <- nrow(motif %>% filter(aso_id == !!aso_id, rbp == !!rbp)) > 0
      has_expr <- nrow(expr %>% filter(aso_id == !!aso_id, rbp == !!rbp, expressed)) > 0
      frow <- func %>% filter(aso_id == !!aso_id, rbp == !!rbp)
      direction <- if (nrow(frow)) frow$regulatory_direction_hint[[1]] else "unknown"
      process <- if (nrow(frow)) frow$rna_process[[1]] else "unknown"
      has_clip0 <- any(clip$aso_id == aso_id & isTRUE(clip$overlap_0nt) &
                         !is.na(clip$rbp) & clip$rbp == rbp)

      if (has_clip0) {
        score <- score + 3
        reasons <- c(reasons, "+3 CLIP 0nt overlap")
      }
      if (has_motif && has_expr) {
        score <- score + 3
        reasons <- c(reasons, "+3 motif+expression")
      } else if (has_motif) {
        score <- score + 2
        reasons <- c(reasons, "+2 motif")
      } else if (has_expr) {
        score <- score + 1
        reasons <- c(reasons, "+1 expression prior")
      }
      if (grepl("stabil", direction) || grepl("activat", direction)) {
        score <- score + 2
        reasons <- c(reasons, "+2 activating/stabilizing function (phenotype-compatible if occlusion relieves repression OR displacement of repressor — context needed)")
      }
      if (grepl("repress", direction)) {
        score <- score + 2
        reasons <- c(reasons, "+2 repressive function (occlusion could increase expression)")
      }
      if (grepl("structure", direction) || grepl("remodel", process)) {
        score <- score + 1
        reasons <- c(reasons, "+1 structure-related")
      }
      if (grepl("splic", process)) {
        score <- score + 1
        reasons <- c(reasons, "+1 splicing/processing")
      }
      # Soft penalty: activating stabilizer + assuming occlusion without CLIP is ambiguous
      if (!has_clip0 && grepl("activating_stabilizing", direction) && has_motif) {
        score <- score - 2
        reasons <- c(reasons, "-2 ambiguous: stabilizing RBP occlusion would predict decrease without additional context")
      }

      rbp_class_val <- case_when(
        grepl("repress", direction) ~ "repressive",
        grepl("activat|stabil", direction) ~ "activating_stabilizing",
        grepl("structure", direction) ~ "structure_related",
        grepl("splic", process) ~ "splicing",
        TRUE ~ "other"
      )

      rows[[length(rows) + 1]] <- data.frame(
        aso_id = aso_id,
        rbp = rbp,
        score = score,
        rbp_class = rbp_class_val,
        regulatory_direction_hint = direction,
        rna_process = process,
        has_clip_0nt = has_clip0,
        has_motif = has_motif,
        has_expression_prior = has_expr,
        score_rationale = paste(reasons, collapse = " | "),
        phenotype = "increased_SHANK3_expression",
        status = "scored",
        warning = "Scores are transparent heuristics; not experimental validation",
        stringsAsFactors = FALSE
      )
    }
  }

  ranked <- bind_rows(rows) %>%
    group_by(aso_id) %>%
    arrange(desc(score), rbp) %>%
    mutate(rank = row_number()) %>%
    ungroup()

  # Ensure diversity in top set: try include repressive, activating, structure-related
  top <- ranked %>%
    group_by(aso_id) %>%
    group_modify(~ {
      d <- .x
      pick <- d %>% filter(rank <= 5)
      # force-include best of each class if not present
      for (cl in c("repressive", "activating_stabilizing", "structure_related")) {
        if (!any(pick$rbp_class == cl)) {
          extra <- d %>% filter(rbp_class == cl) %>% slice_head(n = 1)
          pick <- bind_rows(pick, extra) %>% distinct(rbp, .keep_all = TRUE)
        }
      }
      pick %>% slice_head(n = 5) %>% mutate(in_top_report = TRUE)
    }) %>%
    ungroup()

  ranked <- ranked %>%
    left_join(top %>% dplyr::select(aso_id, rbp, in_top_report), by = c("aso_id", "rbp")) %>%
    mutate(in_top_report = ifelse(is.na(in_top_report), FALSE, in_top_report))

  out <- file.path(root, "results/06_ranked_rbps.tsv")
  write_tsv(ranked, out)
  log_message("Step 06 RBP ranking complete")
  list(ranked = ranked, top = top, tsv = out)
}

rank_mechanisms <- function(cfg, qc_result, map_result, struct_result, rbp_result, rbp_rank, tf_result) {
  suppressPackageStartupMessages(library(dplyr))
  root <- cfg$project_root
  chem <- chemistry_annotation(cfg)
  mapping <- map_result$mapping
  ranked <- rbp_rank$ranked

  labels <- c(
    "RBP occlusion", "RBP displacement", "RNA-structure remodeling", "RNA stabilization",
    "RNase-H-mediated degradation", "transcriptional/enhancer regulation",
    "splicing or processing", "translation regulation", "off-target mechanism", "unresolved"
  )

  rows <- list()
  for (aso in cfg$asos) {
    aso_id <- aso$id
    msub <- mapping %>% filter(aso_id == !!aso_id)
    exact <- any(msub$exact_shank3_match, na.rm = TRUE)
    approx <- any(msub$approximate_shank3_match, na.rm = TRUE)
    top_rbps <- ranked %>% filter(aso_id == !!aso_id, in_top_report) %>% arrange(desc(score))
    has_repressor <- any(top_rbps$rbp_class == "repressive")
    has_struct <- any(grepl("structure", top_rbps$rbp_class)) ||
      any(struct_result$structure$aso_id == aso_id &
            !is.na(struct_result$structure$site_structure_class))
    regulatory_tf <- isTRUE(tf_result$ran)

    # Apply 8 rules to score each label (higher = more plausible as hypothesis only)
    scores <- setNames(rep(0, length(labels)), labels)
    rationale <- setNames(vector("list", length(labels)), labels)

    add <- function(lab, pts, why) {
      scores[[lab]] <<- scores[[lab]] + pts
      rationale[[lab]] <<- c(rationale[[lab]], paste0(if (pts >= 0) "+" else "", pts, " ", why))
    }

    # Rule 1–2: gapmer not supported; RNase-H direct SHANK3 mRNA low support
    add("RNase-H-mediated degradation", -6,
        paste0("Rule1/2: gapmer=", chem$gapmer,
               "; rnase_h_direct_SHANK3_mRNA=", chem$rnase_h_direct_SHANK3_mRNA,
               "; architecture=", chem$architecture))
    # Rule 3: phenotype increased expression
    add("RNase-H-mediated degradation", -3,
        "Rule3: phenotype is increased SHANK3 — on-target RNase-H cleavage typically decreases RNA")
    add("RNA stabilization", 2, "Rule3: phenotype-compatible (stabilization / occupancy)")
    add("RBP occlusion", 3, "Rule3/4: occlusion of repressor can increase expression; PS favors protein engagement")
    add("RBP displacement", 2, "Rule4: PS ASO may displace RBP complexes")
    if (identical(chem$ps_linkage_status, "awaiting_supplier_notation")) {
      add("unresolved", 1,
          "Rule4 caveat: PS linkage pattern awaits supplier notation — PS map not invented")
    }
    # Rule 5: likely uniform MOE / all bases in MOE brackets → steric/occupancy
    if (isTRUE(chem$moe) || grepl("uniform_MOE|MOE", chem$architecture, ignore.case = TRUE)) {
      add("RBP occlusion", 3,
          "Rule5: likely_uniform_MOE / MOE brackets favor steric occupancy over RNase-H cleavage")
      add("RBP displacement", 2, "Rule5: uniform MOE occupancy can displace RBPs")
      add("RNA-structure remodeling", 1, "Rule5: steric MOE binding may remodel local structure (hypothesis)")
      add("RNase-H-mediated degradation", -3,
          "Rule5: uniform MOE designs have low support for RNase-H of SHANK3 mRNA")
    }
    # Rule 6: 5m-dC reported but unconfirmed
    if (isTRUE(chem$five_methyl_dC_reported_unconfirmed)) {
      add("unresolved", 1,
          "Rule6: 5m-dC reported_but_needs_confirmation — caveat only; not treated as confirmed chemistry")
    } else if (isTRUE(chem$five_methyl_dC_confirmed)) {
      add("unresolved", 1, "Rule6: confirmed 5m-dC may alter recognition vs unmodified models")
    }
    # Rule 7: structure
    if (has_struct) {
      add("RNA-structure remodeling", 2, "Rule7: site structure context available (still hypothetical)")
    } else {
      add("RNA-structure remodeling", 0, "Rule7: limited structure evidence")
    }
    # Rule 8
    add("unresolved", 2, "Rule8: sequence/DB evidence alone cannot confirm mechanism")
    # Rule 9: FAM label
    if (isTRUE(chem$fam_label_present) && !isTRUE(chem$fam_label_in_unlabeled_mechanism)) {
      add("unresolved", 0,
          "Rule9: FAM present in supplier construct but excluded from unlabeled mechanism interpretation")
    }

    if (has_repressor) {
      add("RBP occlusion", 3, "Top RBPs include repressive class")
      add("translation regulation", 1, "Repressor occlusion may affect translation/stability")
    }
    if (any(top_rbps$rna_process == "splicing", na.rm = TRUE)) {
      add("splicing or processing", 3, "Splicing-related RBP prioritized")
    }
    if (any(top_rbps$rna_process == "translation", na.rm = TRUE)) {
      add("translation regulation", 2, "Translation-related RBP prioritized")
    }
    if (regulatory_tf) {
      add("transcriptional/enhancer regulation", 2, "Regulatory genomic context + TF motifs")
    }
    if (!exact) {
      add("off-target mechanism", 2, "No/uncertain exact SHANK3 match — off-target possible")
      add("unresolved", 2, "Ambiguity statement applies")
    } else {
      add("off-target mechanism", 0, "Exact SHANK3 match present; off-target still possible")
    }

    # Stabilization as direct ASO effect (non-RNase-H)
    add("RNA stabilization", 1, "Occupancy mechanisms can stabilize RNA in some contexts (hypothesis)")

    rnase_h_status <- as.character(chem$rnase_h_direct_SHANK3_mRNA %||% "low_support")
    ord <- sort(scores, decreasing = TRUE)
    for (lab in names(ord)) {
      rows[[length(rows) + 1]] <- data.frame(
        aso_id = aso_id,
        mechanism_label = lab,
        hypothesis_score = as.numeric(ord[[lab]]),
        rank = which(names(ord) == lab),
        is_top = which(names(ord) == lab) == 1,
        rnase_h_status = rnase_h_status,
        rnase_h_supported = FALSE,
        chemistry_core = chem$core_chemistry,
        chemistry_architecture = chem$architecture,
        chemistry_gapmer = chem$gapmer,
        chemistry_moe = chem$moe,
        chemistry_all_bases_moe = chem$all_bases_in_moe_brackets,
        chemistry_ps = identical(chem$backbone, "phosphorothioate"),
        chemistry_ps_linkage_status = chem$ps_linkage_status,
        chemistry_5mdC = as.character(chem$five_methyl_dC),
        chemistry_fam_present = chem$fam_label_present,
        chemistry_fam_in_unlabeled_mechanism = chem$fam_label_in_unlabeled_mechanism,
        exact_shank3_match = exact,
        approximate_shank3_match = approx,
        top_rbps = paste(head(top_rbps$rbp, 5), collapse = ","),
        rationale = paste(rationale[[lab]], collapse = " || "),
        caveat = paste(
          "Hypothesis only. Do not treat as confirmed mechanism.",
          "FAM excluded from unlabeled mechanism interpretation.",
          if (isTRUE(chem$five_methyl_dC_reported_unconfirmed))
            "5m-dC reported but unconfirmed." else "",
          if (identical(chem$ps_linkage_status, "awaiting_supplier_notation"))
            "PS linkage map awaiting supplier notation (not invented)." else "",
          if (!exact) ambiguity_statement(cfg) else ""
        ),
        status = "scored",
        warning = "Chemistry-aware ranking; not experimental proof",
        stringsAsFactors = FALSE
      )
    }
  }

  # Provenance note for missing PS notation
  if (identical(chem$ps_linkage_status, "awaiting_supplier_notation")) {
    append_provenance(
      "PS_linkage_pattern",
      "parse_from_supplier_notation",
      "supplier_notation_string_not_provided",
      "GRCh38",
      notes = paste(
        "TODO: supplier PS linkage notation string required.",
        "No PS map invented. Status=awaiting_supplier_notation."
      )
    )
  }
  append_provenance(
    "ASO_chemistry_config",
    paste(chem$core_chemistry, chem$architecture, chem$gapmer, sep = "/"),
    "file://_config.yml#chemistry",
    "GRCh38",
    notes = chem$note
  )

  hyp <- bind_rows(rows) %>% arrange(aso_id, rank)
  out <- file.path(root, "results/07_mechanism_hypotheses.tsv")
  write_tsv(hyp, out)
  log_message("Step 07 mechanism ranking complete (revised chemistry)")
  list(hypotheses = hyp, tsv = out, chemistry = chem)
}
