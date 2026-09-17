# R/01_sequence_qc.R — ASO sequence QC (Biostrings)
# Canonical bases for alignment/motifs; retain chemistry annotation (PS/MOE/5m-dC).

run_sequence_qc <- function(cfg) {
  suppressPackageStartupMessages({
    library(Biostrings)
    library(dplyr)
  })
  root <- cfg$project_root %||% project_root()
  ensure_dirs(root)
  chem <- chemistry_annotation(cfg)

  rows <- lapply(cfg$asos, function(aso) {
    raw <- aso$sequence %||% ""
    cleaned <- gsub("[^ACGTUacgtu]", "", raw)
    cleaned <- toupper(gsub("U", "T", cleaned))
    valid <- grepl("^[ACGT]+$", cleaned) && nchar(cleaned) > 0
    dna <- if (valid) DNAString(cleaned) else DNAString("")
    rc <- if (valid) as.character(reverseComplement(dna)) else NA_character_
    rna_compat <- if (valid) chartr("T", "U", cleaned) else NA_character_
    gc <- if (valid) {
      alphabetFrequency(dna, baseOnly = TRUE)
      round(100 * letterFrequency(dna, letters = c("G", "C"), as.prob = TRUE)[[1]], 2)
    } else NA_real_
    len <- nchar(cleaned)
    # Homopolymers
    homo <- if (valid) {
      m <- gregexpr("(A{4,}|C{4,}|G{4,}|T{4,})", cleaned, perl = TRUE)[[1]]
      if (m[[1]] == -1) "" else paste(regmatches(cleaned, list(m))[[1]], collapse = ",")
    } else NA_character_
    extreme_gc <- !is.na(gc) && (gc < 25 || gc > 75)
    # Self-complementarity: longest RC k-mer match within sequence (k=4..6)
    self_comp_score <- if (valid && len >= 6) {
      rc_s <- rc
      scores <- vapply(4:min(8, floor(len / 2)), function(k) {
        kmers <- unique(substring(cleaned, 1:(len - k + 1), k:len))
        sum(vapply(kmers, function(km) grepl(km, rc_s, fixed = TRUE), logical(1)))
      }, numeric(1))
      as.integer(sum(scores))
    } else 0L
    # Simple hairpin: stem of >=4 with loop 3–8
    hairpin_flag <- if (valid && len >= 11) {
      found <- FALSE
      for (stem in 4:6) {
        for (loop in 3:8) {
          if (2 * stem + loop > len) next
          for (i in 1:(len - 2 * stem - loop + 1)) {
            left <- substr(cleaned, i, i + stem - 1)
            right <- substr(cleaned, i + stem + loop, i + 2 * stem + loop - 1)
            if (identical(left, as.character(reverseComplement(DNAString(right))))) {
              found <- TRUE
              break
            }
          }
          if (found) break
        }
        if (found) break
      }
      found
    } else FALSE
    cpg_count <- 0L
    if (valid) {
      m <- gregexpr("CG", cleaned, fixed = TRUE)[[1]]
      cpg_count <- if (m[[1]] == -1) 0L else length(m)
    }
    # 5m-dC: reported but needs confirmation — annotate tentatively; retain caveat
    chem_seq <- if (valid && (isTRUE(chem$five_methyl_dC_confirmed) ||
                              isTRUE(chem$five_methyl_dC_reported_unconfirmed))) {
      gsub("C", "5mC", cleaned)
    } else if (valid) {
      cleaned
    } else NA_character_
    warnings <- c()
    if (!valid) warnings <- c(warnings, "invalid_bases")
    if (nzchar(homo %||% "")) warnings <- c(warnings, "homopolymer")
    if (isTRUE(extreme_gc)) warnings <- c(warnings, "extreme_gc")
    if (isTRUE(hairpin_flag)) warnings <- c(warnings, "hairpin_potential")
    if (self_comp_score >= 5) warnings <- c(warnings, "self_complementarity")
    if (cpg_count >= 2) warnings <- c(warnings, "cpg_rich")
    if (isTRUE(chem$five_methyl_dC_reported_unconfirmed)) {
      warnings <- c(warnings, "5m-dC_reported_unconfirmed")
    }
    if (identical(chem$ps_linkage_status, "awaiting_supplier_notation")) {
      warnings <- c(warnings, "PS_linkage_awaiting_supplier_notation")
    }
    if (isTRUE(chem$fam_label_present) && !isTRUE(chem$fam_label_in_unlabeled_mechanism)) {
      warnings <- c(warnings, "FAM_present_excluded_from_unlabeled_mechanism")
    }

    data.frame(
      aso_id = aso$id,
      input_sequence = raw,
      cleaned_sequence = cleaned,
      valid_acgt = valid,
      length_nt = len,
      gc_percent = gc,
      reverse_complement = rc,
      rna_compatible_sequence = rna_compat,
      chemistry_core = chem$core_chemistry,
      chemistry_backbone = chem$backbone,
      chemistry_moe = chem$moe,
      chemistry_all_bases_moe = chem$all_bases_in_moe_brackets,
      chemistry_five_methyl_dC = as.character(chem$five_methyl_dC),
      chemistry_five_methyl_dC_confirmed = chem$five_methyl_dC_confirmed,
      chemistry_architecture = chem$architecture,
      chemistry_gapmer = chem$gapmer,
      chemistry_rnase_h_compatible = chem$rnase_h_compatible,
      chemistry_rnase_h_direct_SHANK3_mRNA = chem$rnase_h_direct_SHANK3_mRNA,
      chemistry_ps_linkage_pattern = chem$ps_linkage_pattern,
      chemistry_ps_linkage_status = chem$ps_linkage_status,
      chemistry_fam_present = chem$fam_label_present,
      chemistry_fam_in_unlabeled_mechanism = chem$fam_label_in_unlabeled_mechanism,
      all_C_annotated_as_5mdC = isTRUE(chem$five_methyl_dC_confirmed) ||
        isTRUE(chem$five_methyl_dC_reported_unconfirmed),
      chemistry_annotated_sequence = chem_seq,
      canonical_bases_for_alignment = cleaned,
      homopolymers = homo %||% "",
      extreme_gc = extreme_gc,
      self_complementarity_score = self_comp_score,
      hairpin_potential = hairpin_flag,
      cpg_count = as.integer(cpg_count),
      qc_warnings = paste(warnings, collapse = ";"),
      status = if (valid) "ok" else "invalid_sequence",
      stringsAsFactors = FALSE
    )
  })

  qc <- bind_rows(rows)
  out_tsv <- file.path(root, "results/01_sequence_qc.tsv")
  out_fa <- file.path(root, "results/01_aso_sequences.fasta")
  write_tsv(qc, out_tsv)

  # FASTA with canonical bases
  valid_qc <- qc %>% filter(valid_acgt)
  if (nrow(valid_qc)) {
    dna_set <- DNAStringSet(valid_qc$cleaned_sequence)
    names(dna_set) <- paste0(
      valid_qc$aso_id,
      "|chemistry=", chem$core_chemistry,
      ";PS;MOE_brackets=", chem$all_bases_in_moe_brackets,
      "|architecture=", chem$architecture,
      "|gapmer=", chem$gapmer,
      "|rnaseH_SHANK3=", chem$rnase_h_direct_SHANK3_mRNA,
      "|5m-dC=", chem$five_methyl_dC,
      "|FAM_in_mech=", chem$fam_label_in_unlabeled_mechanism,
      "|PS_linkage=", chem$ps_linkage_status
    )
    writeXStringSet(dna_set, out_fa)
  } else {
    writeLines("", out_fa)
  }

  log_message("Step 01 sequence QC complete: ", nrow(qc), " ASOs")
  list(qc = qc, tsv = out_tsv, fasta = out_fa)
}
