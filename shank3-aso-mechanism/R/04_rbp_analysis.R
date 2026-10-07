# R/04_rbp_analysis.R — RBP CLIP / motif / expression / function evidence
# Sources: ENCODE eCLIP, UCSC ENCODE, POSTAR3, oRNAment, ATtRACT, RBPDB, RBP2GO
# Never invent peaks. Unreachable DBs → empty tables with status/warning.

# Evidence strength rules:
#   strong   = direct CLIP peak overlap (0 nt) at site in relevant cell type + motif optional
#   moderate = overlap within 10–30 nt OR motif hit + neuronal/NPC expression
#   weak     = nearest distance only OR motif without CLIP OR expression alone
#   none     = no supporting evidence

.rbp_required_clip_cols <- c(
  "aso_id", "rbp", "clip_source", "clip_dataset_id", "cell_type",
  "chrom", "peak_start", "peak_end", "peak_strand",
  "site_start", "site_end", "overlap_0nt", "overlap_10nt", "overlap_30nt",
  "nearest_distance_nt", "evidence_strength", "genome_build",
  "resource_version", "resource_url", "download_date", "status", "warning"
)

.rbp_required_motif_cols <- c(
  "aso_id", "rbp", "motif_source", "motif_id", "motif_sequence_or_regex",
  "window_id", "match_start", "match_end", "match_strand", "score",
  "evidence_strength", "resource_version", "resource_url", "download_date",
  "status", "warning"
)

.rbp_required_expr_cols <- c(
  "aso_id", "rbp", "expression_source", "cell_type", "expressed",
  "expression_value", "expression_units", "resource_version", "resource_url",
  "download_date", "status", "warning"
)

.rbp_required_func_cols <- c(
  "aso_id", "rbp", "function_source", "molecular_function", "rna_process",
  "regulatory_direction_hint", "shank3_relevance_note", "resource_version",
  "resource_url", "download_date", "status", "warning"
)

# Curated literature-backed RBP priors for SHANK3 / neuronal RNA (NOT invented CLIP peaks).
# Used only for motif/expression/function annotation layers with explicit status.
.neuronal_rbp_priors <- function() {
  data.frame(
    rbp = c("FMR1", "ELAVL2", "ELAVL4", "RBFOX1", "RBFOX3", "PTBP2", "CELF4",
            "ATXN2", "TDP43", "HNRNPA2B1", "LIN28A", "PUM2", "MOV10", "UPF1",
            "AGO2", "IGF2BP1", "CSTF2", "SRSF1"),
    molecular_function = c(
      "mRNA binding / translation regulation",
      "3'UTR binding / stabilization",
      "3'UTR binding / stabilization",
      "splicing regulation",
      "splicing regulation (neuronal)",
      "splicing repression (neuronal)",
      "splicing / translation",
      "stress granule / translation",
      "splicing / RNA stability",
      "splicing / processing",
      "miRNA/biogenesis / translation",
      "3'UTR repression",
      "helicase / RNAi-related",
      "NMD / surveillance",
      "RNAi / miRNA repression",
      "m6A reader / stabilization",
      "3' end processing",
      "splicing activation"
    ),
    rna_process = c(
      "translation", "stabilization", "stabilization", "splicing", "splicing",
      "splicing", "splicing", "translation", "stability", "processing",
      "translation", "repression", "remodeling", "surveillance", "repression",
      "stabilization", "processing", "splicing"
    ),
    regulatory_direction_hint = c(
      "context_dependent", "activating_stabilizing", "activating_stabilizing",
      "context_dependent", "context_dependent", "repressive_splicing",
      "context_dependent", "context_dependent", "context_dependent",
      "context_dependent", "context_dependent", "repressive",
      "structure_related", "repressive", "repressive",
      "activating_stabilizing", "processing", "activating_splicing"
    ),
    neuron_expressed = TRUE,
    npc_expressed = TRUE,
    stringsAsFactors = FALSE
  )
}

# Simple consensus motifs (IUPAC) for a subset — from public ATtRACT/RBPDB-style consensus.
# Explicitly marked as consensus approximations, not measured affinities.
.rbp_consensus_motifs <- function() {
  data.frame(
    rbp = c("FMR1", "ELAVL2", "ELAVL4", "RBFOX1", "RBFOX3", "PTBP2", "TDP43", "PUM2", "SRSF1", "HNRNPA2B1"),
    motif_id = c("FMR1_ACYK", "ELAVL_AUUU", "ELAVL_AUUU", "RBFOX_GCAUG", "RBFOX_GCAUG",
                 "PTBP_UCUU", "TDP43_UG", "PUM_UGUANAUA", "SRSF1_GGAGGA", "HNRNPA2_UAGG"),
    motif_regex = c("AC[CT][GT]", "ATTT+", "ATTT+", "GCATG", "GCATG",
                    "TCTT", "(TG){3,}", "TGTAAATA", "GGAGGA", "TAGG"),
    motif_source = "curated_public_consensus_approx",
    stringsAsFactors = FALSE
  )
}

.query_encode_eclip_metadata <- function(cfg) {
  root <- cfg$project_root
  url <- "https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP&assembly=GRCh38&format=json&limit=50"
  dest <- file.path(root, "data/raw/encode_eclip_search.json")
  res <- tryCatch({
    resp <- httr2::request(url) |>
      httr2::req_headers(Accept = "application/json") |>
      httr2::req_timeout(60) |>
      httr2::req_perform()
    txt <- httr2::resp_body_string(resp)
    writeLines(txt, dest)
    append_provenance("ENCODE_eCLIP_search", "ENCODE Portal API", url, "GRCh38",
                      notes = "Metadata search only; peak BED download attempted separately")
    jsonlite::fromJSON(txt, simplifyVector = FALSE)
  }, error = function(e) {
    append_provenance("ENCODE_eCLIP_search", "ENCODE Portal API", url, "GRCh38",
                      notes = paste("FAILED", conditionMessage(e)))
    NULL
  })
  res
}

.download_postar_hint <- function(cfg) {
  # POSTAR3 does not provide a simple bulk SHANK3 API without scraping;
  # record attempt and return empty with clear status.
  url <- cfg$resources$postar3$url %||% "http://postar.ncrnalab.org/"
  append_provenance("POSTAR3", cfg$resources$postar3$version %||% "POSTAR3", url, "GRCh38",
                    notes = "No automated bulk SHANK3 peak dump used; overlaps not fabricated")
  list(status = "not_downloaded_bulk", warning = "POSTAR3 requires interactive/query UI; no peaks invented")
}

.scan_motifs_in_windows <- function(windows_fa, motifs, aso_ids) {
  suppressPackageStartupMessages(library(Biostrings))
  rows <- list()
  if (!file.exists(windows_fa)) return(dplyr::bind_rows(rows))
  seqs <- tryCatch(readDNAStringSet(windows_fa), error = function(e) DNAStringSet())
  if (!length(seqs) || grepl("^no_windows", names(seqs)[1])) return(dplyr::bind_rows(rows))
  for (nm in names(seqs)) {
    aso_id <- sub("\\|.*", "", nm)
    seq <- as.character(seqs[[nm]])
    for (j in seq_len(nrow(motifs))) {
      # Convert simple regex on DNA
      pat <- motifs$motif_regex[[j]]
      hits <- tryCatch(gregexpr(pat, seq, ignore.case = TRUE, perl = TRUE)[[1]], error = function(e) -1)
      if (hits[[1]] == -1) next
      ml <- attr(hits, "match.length")
      for (k in seq_along(hits)) {
        rows[[length(rows) + 1]] <- data.frame(
          aso_id = aso_id,
          rbp = motifs$rbp[[j]],
          motif_source = motifs$motif_source[[j]],
          motif_id = motifs$motif_id[[j]],
          motif_sequence_or_regex = pat,
          window_id = nm,
          match_start = as.integer(hits[[k]]),
          match_end = as.integer(hits[[k]] + ml[[k]] - 1L),
          match_strand = "+",
          score = NA_real_,
          evidence_strength = "weak",
          resource_version = "curated_consensus_v1",
          resource_url = "ATtRACT/RBPDB-style consensus (curated)",
          download_date = as.character(Sys.Date()),
          status = "motif_scan_ok",
          warning = "Consensus motif approx; not CLIP-validated at this locus",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows)) dplyr::bind_rows(rows) else {
    df <- as.data.frame(matrix(nrow = 0, ncol = length(.rbp_required_motif_cols)))
    names(df) <- .rbp_required_motif_cols
    df
  }
}

run_rbp_analysis <- function(cfg, map_result, qc_result) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
  })
  root <- cfg$project_root
  ensure_dirs(root)
  today <- as.character(Sys.Date())
  mapping <- map_result$mapping
  priors <- .neuronal_rbp_priors()
  motifs <- .rbp_consensus_motifs()

  # --- CLIP overlaps ---
  encode_meta <- .query_encode_eclip_metadata(cfg)
  postar <- .download_postar_hint(cfg)
  for (src in c("oRNAment", "ATtRACT", "RBPDB", "RBP2GO", "UCSC_ENCODE")) {
    u <- cfg$resources[[tolower(src)]]$url %||% cfg$resources[[gsub("-", "", tolower(src))]]$url
    # map keys
  }
  append_provenance("oRNAment", cfg$resources$ornament$version, cfg$resources$ornament$url, "GRCh38",
                    notes = "Motif resource referenced; bulk matrices may be unreachable")
  append_provenance("ATtRACT", cfg$resources$attract$version, cfg$resources$attract$url, "GRCh38",
                    notes = "Motif resource referenced")
  append_provenance("RBPDB", cfg$resources$rbpdb$version, cfg$resources$rbpdb$url, "GRCh38",
                    notes = "Motif resource referenced")
  append_provenance("RBP2GO", cfg$resources$rbp2go$version, cfg$resources$rbp2go$url, "GRCh38",
                    notes = "Function annotation resource referenced")

  clip_rows <- list()
  # Attempt: if mapping has genomic coords, we cannot invent eCLIP peaks.
  # Write explicit unavailable rows per ASO documenting sources tried.
  for (aso_id in unique(qc_result$qc$aso_id)) {
    sites <- mapping %>% filter(aso_id == !!aso_id, !is.na(chrom), !is.na(start))
    if (!nrow(sites)) {
      clip_rows[[length(clip_rows) + 1]] <- data.frame(
        aso_id = aso_id, rbp = NA_character_, clip_source = "ENCODE_eCLIP;POSTAR3;UCSC_ENCODE",
        clip_dataset_id = NA_character_, cell_type = paste(cfg$cell_types, collapse = ","),
        chrom = NA_character_, peak_start = NA_integer_, peak_end = NA_integer_, peak_strand = NA_character_,
        site_start = NA_integer_, site_end = NA_integer_,
        overlap_0nt = FALSE, overlap_10nt = FALSE, overlap_30nt = FALSE,
        nearest_distance_nt = NA_integer_, evidence_strength = "none",
        genome_build = "GRCh38", resource_version = "n/a",
        resource_url = cfg$resources$encode_eclip$url,
        download_date = today, status = "no_site_coords",
        warning = "No genomic site for overlap; CLIP peaks not invented",
        stringsAsFactors = FALSE
      )
      next
    }
    for (r in seq_len(min(3, nrow(sites)))) {
      clip_rows[[length(clip_rows) + 1]] <- data.frame(
        aso_id = aso_id, rbp = NA_character_,
        clip_source = "ENCODE_eCLIP;POSTAR3;UCSC_ENCODE",
        clip_dataset_id = if (!is.null(encode_meta)) "encode_search_cached" else NA_character_,
        cell_type = paste(cfg$cell_types, collapse = ","),
        chrom = sites$chrom[[r]], peak_start = NA_integer_, peak_end = NA_integer_,
        peak_strand = NA_character_,
        site_start = as.integer(sites$start[[r]]), site_end = as.integer(sites$end[[r]]),
        overlap_0nt = FALSE, overlap_10nt = FALSE, overlap_30nt = FALSE,
        nearest_distance_nt = NA_integer_, evidence_strength = "none",
        genome_build = "GRCh38",
        resource_version = if (!is.null(encode_meta)) "ENCODE_API_live" else "unavailable",
        resource_url = cfg$resources$encode_eclip$url,
        download_date = today,
        status = if (!is.null(encode_meta)) "metadata_only_no_peak_overlap" else "encode_unreachable",
        warning = paste(
          "eCLIP peak BED files not auto-downloaded for SHANK3 locus in this run;",
          postar$warning, "; overlaps left empty (not fabricated)"
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  clip_df <- bind_rows(clip_rows)
  # ensure columns
  for (cc in .rbp_required_clip_cols) if (!cc %in% names(clip_df)) clip_df[[cc]] <- NA
  clip_df <- clip_df[, .rbp_required_clip_cols]

  # --- Motifs ---
  motif_df <- .scan_motifs_in_windows(map_result$windows_fasta, motifs, qc_result$qc$aso_id)
  for (cc in .rbp_required_motif_cols) if (!cc %in% names(motif_df)) motif_df[[cc]] <- NA
  if (!nrow(motif_df)) {
    motif_df <- data.frame(matrix(ncol = length(.rbp_required_motif_cols), nrow = 0))
    names(motif_df) <- .rbp_required_motif_cols
  } else {
    motif_df <- motif_df[, .rbp_required_motif_cols]
  }

  # --- Expression ---
  expr_rows <- list()
  for (aso_id in unique(qc_result$qc$aso_id)) {
    for (ct in cfg$cell_types) {
      for (j in seq_len(nrow(priors))) {
        expr_rows[[length(expr_rows) + 1]] <- data.frame(
          aso_id = aso_id,
          rbp = priors$rbp[[j]],
          expression_source = "curated_neuronal_prior_not_quant_RNAseq",
          cell_type = ct,
          expressed = if (ct == "neuron") priors$neuron_expressed[[j]] else priors$npc_expressed[[j]],
          expression_value = NA_real_,
          expression_units = "binary_prior",
          resource_version = "curated_v1",
          resource_url = "literature_curated_neuronal_RBP_set",
          download_date = today,
          status = "prior_annotation",
          warning = "Not a quantitative ENCODE/GTEx measurement; flag as prior only",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  expr_df <- bind_rows(expr_rows)
  expr_df <- expr_df[, .rbp_required_expr_cols]

  # --- Function ---
  func_rows <- list()
  for (aso_id in unique(qc_result$qc$aso_id)) {
    for (j in seq_len(nrow(priors))) {
      func_rows[[length(func_rows) + 1]] <- data.frame(
        aso_id = aso_id,
        rbp = priors$rbp[[j]],
        function_source = "RBP2GO_style_curated",
        molecular_function = priors$molecular_function[[j]],
        rna_process = priors$rna_process[[j]],
        regulatory_direction_hint = priors$regulatory_direction_hint[[j]],
        shank3_relevance_note = "Candidate neuronal RBP; locus-specific CLIP not confirmed in this run",
        resource_version = cfg$resources$rbp2go$version %||% "RBP2GO",
        resource_url = cfg$resources$rbp2go$url,
        download_date = today,
        status = "curated_annotation",
        warning = "Functional labels are priors; not proof of SHANK3 engagement",
        stringsAsFactors = FALSE
      )
    }
  }
  func_df <- bind_rows(func_rows)
  func_df <- func_df[, .rbp_required_func_cols]

  out_clip <- file.path(root, "results/04_rbp_clip_overlap.tsv")
  out_motif <- file.path(root, "results/04_rbp_motif_overlap.tsv")
  out_expr <- file.path(root, "results/04_rbp_expression.tsv")
  out_func <- file.path(root, "results/04_rbp_function.tsv")
  write_tsv(clip_df, out_clip)
  write_tsv(motif_df, out_motif)
  write_tsv(expr_df, out_expr)
  write_tsv(func_df, out_func)

  log_message("Step 04 RBP analysis complete. motif_hits=", nrow(motif_df))
  list(clip = clip_df, motif = motif_df, expression = expr_df, function_annot = func_df,
       clip_tsv = out_clip, motif_tsv = out_motif, expr_tsv = out_expr, func_tsv = out_func)
}
