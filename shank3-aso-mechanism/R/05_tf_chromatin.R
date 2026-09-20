# R/05_tf_chromatin.R — TF / chromatin (conditional)
# Only if target is nuclear regulatory / enhancer / promoter / antisense / intronic regulatory RNA.
# TF motifs on genomic DNA only, never RNA. Separate from RBP evidence.

.is_regulatory_context <- function(mapping) {
  if (!nrow(mapping)) return(FALSE)
  ft <- tolower(paste(c(mapping$feature_type, mapping$match_class), collapse = " "))
  # Intronic / genomic gene-body / flanking / antisense / promoter-like contexts
  any(grepl("intronic|promoter|enhancer|antisense|genomic_flanking|genomic_gene_body|genomic_|regulatory", ft))
}

.run_tf_motif_scan_dna <- function(cfg, map_result) {
  suppressPackageStartupMessages({
    library(Biostrings)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(dplyr)
  })
  # Minimal JASPAR-like IUPAC motifs (DNA) for demonstration of genomic scan only
  motifs <- data.frame(
    tf = c("REST", "CTCF", "MEF2A", "NEUROD1", "SP1"),
    motif_id = c("MA0138", "MA0139", "MA0052", "MA1109", "MA0079"),
    dna_regex = c("TTCAGCAC[CT]","CCGCG[ACGT]GG[AC]GG","CTA[AT]AAAATAG","CAGCTG[ACGT]{0,2}CAGCTG","GGGCGG"),
    stringsAsFactors = FALSE
  )
  append_provenance("JASPAR_curated_subset", "JASPAR-style curated regex subset",
                    cfg$resources$jaspar$url, "GRCh38",
                    notes = "Small DNA motif subset for genomic scan; not full JASPAR PWM library")
  mapping <- map_result$mapping %>%
    dplyr::filter(!is.na(chrom), !is.na(start), !is.na(end), grepl("^chr", chrom))
  rows <- list()
  if (!nrow(mapping)) return(dplyr::bind_rows(rows))

  for (i in seq_len(nrow(mapping))) {
    chrom <- mapping$chrom[[i]]
    # Expand to ±200 for TF context on DNA
    up <- as.integer(cfg$window_upstream %||% 200)
    down <- as.integer(cfg$window_downstream %||% 200)
    gr_start <- max(1L, as.integer(mapping$start[[i]]) - up)
    gr_end <- as.integer(mapping$end[[i]]) + down
    dna <- tryCatch({
      getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, gr_start, gr_end)
    }, error = function(e) NULL)
    if (is.null(dna)) next
    seq <- as.character(dna)
    for (j in seq_len(nrow(motifs))) {
      hits <- gregexpr(motifs$dna_regex[[j]], seq, ignore.case = TRUE, perl = TRUE)[[1]]
      if (hits[[1]] == -1) next
      ml <- attr(hits, "match.length")
      for (k in seq_along(hits)) {
        rows[[length(rows) + 1]] <- data.frame(
          aso_id = mapping$aso_id[[i]],
          tf = motifs$tf[[j]],
          motif_id = motifs$motif_id[[j]],
          source = "JASPAR_style_DNA_motif_scan",
          chrom = chrom,
          start = gr_start + hits[[k]] - 1L,
          end = gr_start + hits[[k]] + ml[[k]] - 2L,
          strand = "+",
          score = NA_real_,
          assay = "motif_scan_only",
          cell_type = paste(cfg$cell_types, collapse = ","),
          evidence_note = "DNA motif presence ≠ ChIP occupancy; not RBP evidence",
          genome_build = "GRCh38",
          resource_url = cfg$resources$jaspar$url,
          download_date = as.character(Sys.Date()),
          status = "motif_scan_ok",
          warning = "No ChIP peak confirmation in this automated run",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows)) dplyr::bind_rows(rows) else data.frame()
}

run_tf_chromatin <- function(cfg, map_result) {
  suppressPackageStartupMessages(library(dplyr))
  root <- cfg$project_root
  ensure_dirs(root)
  mapping <- map_result$mapping
  regulatory <- .is_regulatory_context(mapping)

  append_provenance("ENCODE_TF", "ENCODE Portal", "https://www.encodeproject.org/", "GRCh38",
                    notes = "TF/chromatin step conditional")
  append_provenance("ChIP-Atlas", cfg$resources$chip_atlas$version, cfg$resources$chip_atlas$url, "GRCh38",
                    notes = "Referenced; peaks not fabricated")
  append_provenance("SCREEN_cCREs", cfg$resources$screen$version, cfg$resources$screen$url, "GRCh38",
                    notes = "Referenced for enhancer/promoter context")

  out_tsv <- file.path(root, "results/05_tf_chromatin.tsv")
  out_bed <- file.path(root, "results/05_regulatory_locus.bed")

  if (!regulatory) {
    df <- data.frame(
      aso_id = unique(mapping$aso_id %||% NA_character_),
      tf = NA_character_, motif_id = NA_character_, source = NA_character_,
      chrom = NA_character_, start = NA_integer_, end = NA_integer_, strand = NA_character_,
      score = NA_real_, assay = NA_character_, cell_type = NA_character_,
      evidence_note = "TF/chromatin analysis skipped: target context not classified as nuclear regulatory/enhancer/promoter/antisense/intronic regulatory RNA",
      genome_build = "GRCh38", resource_url = NA_character_,
      download_date = as.character(Sys.Date()),
      status = "skipped_not_regulatory_context",
      warning = "No TF inference on RNA; step gated by locus class",
      stringsAsFactors = FALSE
    )
    write_tsv(df, out_tsv)
    writeLines("", out_bed)
    log_message("Step 05 TF/chromatin skipped (non-regulatory context)")
    return(list(tf = df, bed = out_bed, tsv = out_tsv, ran = FALSE))
  }

  # Regulatory locus BED from genomic mapping rows
  gen <- mapping %>% filter(!is.na(chrom), grepl("^chr", chrom), !is.na(start))
  if (nrow(gen)) {
    up <- as.integer(cfg$window_upstream %||% 200)
    down <- as.integer(cfg$window_downstream %||% 200)
    bed <- gen %>% transmute(
      chrom, start = pmax(0L, as.integer(start) - up - 1L), end = as.integer(end) + down,
      name = paste(aso_id, match_class, sep = "_"), score = 500, strand = ifelse(is.na(strand), ".", strand)
    )
    write.table(bed, out_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    writeLines("", out_bed)
  }

  tf_df <- .run_tf_motif_scan_dna(cfg, map_result)
  if (!nrow(tf_df)) {
    tf_df <- data.frame(
      aso_id = unique(mapping$aso_id),
      tf = NA_character_, motif_id = NA_character_, source = "JASPAR/ENCODE/ChIP-Atlas/SCREEN",
      chrom = NA_character_, start = NA_integer_, end = NA_integer_, strand = NA_character_,
      score = NA_real_, assay = "none", cell_type = paste(cfg$cell_types, collapse = ","),
      evidence_note = "Regulatory context flagged but no motif hits / ChIP peaks available",
      genome_build = "GRCh38", resource_url = cfg$resources$jaspar$url,
      download_date = as.character(Sys.Date()),
      status = "no_tf_hits", warning = "ChIP peaks not auto-fetched; motifs empty",
      stringsAsFactors = FALSE
    )
  }
  write_tsv(tf_df, out_tsv)
  log_message("Step 05 TF/chromatin complete. rows=", nrow(tf_df))
  list(tf = tf_df, bed = out_bed, tsv = out_tsv, ran = TRUE)
}
