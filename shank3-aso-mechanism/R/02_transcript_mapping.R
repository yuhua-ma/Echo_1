# R/02_transcript_mapping.R — SHANK3 transcript / genomic mapping (GRCh38 only)
# Exact match of reverse complement against coding/noncoding SHANK3 transcripts,
# antisense/intronic, nearby GENCODE, genomic locus both strands.
# If no exact SHANK3 transcript match: do NOT invent a transcript.

.gencode_cache_path <- function(root) {
  file.path(root, "data/raw/gencode.v45.annotation.gtf.gz")
}

download_gencode <- function(cfg) {
  root <- cfg$project_root
  dest <- .gencode_cache_path(root)
  url <- cfg$resources$gencode$url
  version <- cfg$resources$gencode$version
  if (file.exists(dest) && file.info(dest)$size > 1e6) {
    append_provenance("GENCODE_GTF", version, url, "GRCh38", notes = "cached")
    return(list(ok = TRUE, path = dest, status = "cached", warning = "", url = url, version = version))
  }
  res <- safe_download(url, dest, "GENCODE_GTF", version, "GRCh38",
                       notes = "Full GENCODE annotation GTF", timeout = 600)
  res$version <- version
  res
}

.load_shank3_from_bsgenome <- function(cfg) {
  suppressPackageStartupMessages({
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(Biostrings)
  })
  # SHANK3 Entrez: 85358
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  eg <- "85358"
  tx_ids <- tryCatch({
    suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = eg, columns = c("SYMBOL", "ENSEMBL", "ENSEMBLTRANS"),
                         keytype = "ENTREZID"))
  }, error = function(e) NULL)

  # Gene locus from TxDb
  genes_gr <- genes(txdb)
  shank3_gene <- genes_gr[genes_gr$gene_id == eg]
  if (length(shank3_gene) == 0) {
    # Fallback known GRCh38 locus for SHANK3 (chr22)
    shank3_gene <- GRanges("chr22", IRanges(50674415, 50733212), strand = "+", gene_id = eg)
  }
  seqlevelsStyle(shank3_gene) <- "UCSC"

  # Transcripts overlapping gene
  txs <- transcriptsBy(txdb, by = "gene")[[eg]]
  if (is.null(txs) || length(txs) == 0) {
    txs <- transcripts(txdb, filter = list(gene_id = eg))
  }

  # Extract transcript sequences (exonic, spliced)
  tx_seqs <- list()
  tx_meta <- list()
  if (length(txs)) {
    exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
    # Map tx names
    tx_names <- names(txs)
    if (is.null(tx_names) || all(is.na(tx_names))) {
      # use tx_id / tx_name from mcols
      if ("tx_name" %in% names(mcols(txs))) tx_names <- mcols(txs)$tx_name
      else tx_names <- as.character(mcols(txs)$tx_id)
    }
    for (i in seq_along(txs)) {
      tn <- tx_names[[i]]
      if (is.null(tn) || is.na(tn) || !tn %in% names(exons_by_tx)) next
      ex <- exons_by_tx[[tn]]
      if (is.null(ex) || length(ex) == 0) next
      seq <- tryCatch(extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, exons_by_tx[tn])[[1]],
                      error = function(e) NULL)
      if (is.null(seq)) next
      tx_seqs[[tn]] <- seq
      tx_meta[[tn]] <- data.frame(
        transcript_id = tn,
        gene_symbol = "SHANK3",
        gene_id = eg,
        seqnames = as.character(seqnames(txs[i])),
        start = start(txs[i]),
        end = end(txs[i]),
        strand = as.character(strand(txs[i])),
        biotype = "protein_coding_or_knownGene",
        source = "TxDb.Hsapiens.UCSC.hg38.knownGene",
        stringsAsFactors = FALSE
      )
    }
  }

  # Also extract full gene locus genomic DNA (both strands handled at match time)
  pad <- 5000L
  locus <- flank(shank3_gene, width = pad, both = TRUE)
  locus <- trim(locus)
  locus_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, locus)[[1]]

  list(
    gene = shank3_gene,
    locus = locus,
    locus_seq = locus_seq,
    tx_seqs = DNAStringSet(tx_seqs),
    tx_meta = if (length(tx_meta)) dplyr::bind_rows(tx_meta) else data.frame(),
    annotation_source = "TxDb.Hsapiens.UCSC.hg38.knownGene + BSgenome.Hsapiens.UCSC.hg38",
    annotation_version = as.character(packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene")),
    genome_build = "GRCh38/hg38"
  )
}

.find_matches_in_subject <- function(query_dna, subject_dna, max_mismatch = 0L) {
  if (length(subject_dna) == 0 || nchar(as.character(query_dna)) == 0) {
    return(data.frame())
  }
  hits <- matchPattern(query_dna, subject_dna, max.mismatch = max_mismatch, with.indels = FALSE)
  if (length(hits) == 0) return(data.frame())
  data.frame(
    start = start(hits),
    end = end(hits),
    width = width(hits),
    stringsAsFactors = FALSE
  )
}

.run_offtarget_screen <- function(aso_rc, genome, max_hits = 50) {
  # Lightweight off-target: search SHANK3 chromosome only (chr22) for exact RC
  # Full transcriptome screen is expensive; document limitation.
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  chr <- BSgenome.Hsapiens.UCSC.hg38[["chr22"]]
  hits_plus <- matchPattern(DNAString(aso_rc), chr, max.mismatch = 0)
  hits_minus <- matchPattern(reverseComplement(DNAString(aso_rc)), chr, max.mismatch = 0)
  rows <- list()
  if (length(hits_plus)) {
    n <- min(length(hits_plus), max_hits)
    rows[[1]] <- data.frame(
      chrom = "chr22", start = start(hits_plus)[1:n], end = end(hits_plus)[1:n],
      strand = "+", match_type = "genomic_exact", stringsAsFactors = FALSE
    )
  }
  if (length(hits_minus)) {
    n <- min(length(hits_minus), max_hits)
    rows[[2]] <- data.frame(
      chrom = "chr22", start = start(hits_minus)[1:n], end = end(hits_minus)[1:n],
      strand = "-", match_type = "genomic_exact_opposite_query", stringsAsFactors = FALSE
    )
  }
  if (length(rows)) dplyr::bind_rows(rows) else data.frame()
}

run_transcript_mapping <- function(cfg, qc_result) {
  suppressPackageStartupMessages({
    library(Biostrings)
    library(GenomicRanges)
    library(dplyr)
    library(readr)
  })
  root <- cfg$project_root
  ensure_dirs(root)
  up <- as.integer(cfg$window_upstream %||% 200)
  down <- as.integer(cfg$window_downstream %||% 200)

  gencode_dl <- download_gencode(cfg)
  append_provenance(
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    as.character(utils::packageVersion("TxDb.Hsapiens.UCSC.hg38.knownGene")),
    "Bioconductor AnnotationPackage",
    "GRCh38",
    notes = "Primary transcript source for SHANK3 mapping"
  )
  append_provenance(
    "BSgenome.Hsapiens.UCSC.hg38",
    as.character(utils::packageVersion("BSgenome.Hsapiens.UCSC.hg38")),
    "Bioconductor BSgenome",
    "GRCh38",
    notes = "Genome sequence GRCh38/hg38 only"
  )

  annot <- tryCatch(.load_shank3_from_bsgenome(cfg), error = function(e) {
    log_message("Annotation load failed: ", conditionMessage(e))
    NULL
  })

  qc <- qc_result$qc
  map_rows <- list()
  window_seqs <- list()
  window_beds <- list()
  offtarget_rows <- list()

  for (i in seq_len(nrow(qc))) {
    aso_id <- qc$aso_id[[i]]
    if (!isTRUE(qc$valid_acgt[[i]])) {
      map_rows[[length(map_rows) + 1]] <- data.frame(
        aso_id = aso_id, match_class = "invalid_aso", exact_shank3_match = FALSE,
        approximate_shank3_match = FALSE, feature_type = NA_character_,
        transcript_id = NA_character_, gene_symbol = "SHANK3",
        chrom = NA_character_, start = NA_integer_, end = NA_integer_,
        strand = NA_character_, mismatches = NA_integer_,
        rank = NA_integer_, support_note = "ASO failed QC",
        neuronal_expression_support = "unknown",
        status = "invalid_aso", warning = "skip_mapping",
        stringsAsFactors = FALSE
      )
      next
    }
    sense <- qc$cleaned_sequence[[i]]
    rc <- qc$reverse_complement[[i]]
    found_exact <- FALSE
    found_approx <- FALSE
    candidates <- list()

    if (!is.null(annot) && length(annot$tx_seqs)) {
      for (tn in names(annot$tx_seqs)) {
        subj <- annot$tx_seqs[[tn]]
        # Exact RC (ASO binds RNA via Watson-Crick to transcript)
        m0 <- .find_matches_in_subject(DNAString(rc), subj, 0L)
        if (nrow(m0)) {
          found_exact <- TRUE
          for (j in seq_len(nrow(m0))) {
            candidates[[length(candidates) + 1]] <- data.frame(
              aso_id = aso_id, match_class = "exact_transcript_rc",
              exact_shank3_match = TRUE, approximate_shank3_match = FALSE,
              feature_type = "transcript_exon_spliced",
              transcript_id = tn, gene_symbol = "SHANK3",
              chrom = annot$tx_meta$seqnames[annot$tx_meta$transcript_id == tn][1],
              start = m0$start[[j]], end = m0$end[[j]],
              strand = annot$tx_meta$strand[annot$tx_meta$transcript_id == tn][1],
              mismatches = 0L, query_coord_space = "transcript",
              support_note = "Exact reverse-complement match in SHANK3 spliced transcript",
              neuronal_expression_support = "SHANK3_neuronal_gene",
              isoform_shared = NA, rank_score = 100,
              status = "ok", warning = "",
              stringsAsFactors = FALSE
            )
            # Extract window from transcript
            w_start <- max(1, m0$start[[j]] - up)
            w_end <- min(length(subj), m0$end[[j]] + down)
            win <- subseq(subj, w_start, w_end)
            wid <- paste0(aso_id, "|", tn, "|tx", w_start, "-", w_end)
            window_seqs[[wid]] <- win
            window_beds[[length(window_beds) + 1]] <- data.frame(
              chrom = paste0("transcript:", tn),
              start = w_start - 1L, end = w_end,
              name = wid, score = 1000, strand = "+",
              aso_id = aso_id, match_class = "exact_transcript_rc",
              stringsAsFactors = FALSE
            )
          }
        }
        # Approximate 1–2 mismatch
        for (mm in 1:2) {
          mA <- .find_matches_in_subject(DNAString(rc), subj, as.integer(mm))
          if (nrow(mA)) {
            found_approx <- TRUE
            # exclude exact already counted: keep only if not exact for this position
            for (j in seq_len(nrow(mA))) {
              if (found_exact && any(abs(mA$start[[j]] - vapply(candidates, function(x) x$start[[1]], 1)) < 1)) next
              candidates[[length(candidates) + 1]] <- data.frame(
                aso_id = aso_id, match_class = paste0("approx_transcript_rc_mm", mm),
                exact_shank3_match = FALSE, approximate_shank3_match = TRUE,
                feature_type = "transcript_exon_spliced",
                transcript_id = tn, gene_symbol = "SHANK3",
                chrom = annot$tx_meta$seqnames[annot$tx_meta$transcript_id == tn][1],
                start = mA$start[[j]], end = mA$end[[j]],
                strand = annot$tx_meta$strand[annot$tx_meta$transcript_id == tn][1],
                mismatches = as.integer(mm), query_coord_space = "transcript",
                support_note = paste0("Approximate RC match (", mm, " mismatch) in SHANK3 transcript"),
                neuronal_expression_support = "SHANK3_neuronal_gene",
                isoform_shared = NA, rank_score = 80 - 10 * mm,
                status = "ok", warning = "approximate_only",
                stringsAsFactors = FALSE
              )
            }
          }
        }
        # Sense exact (less likely binding mode for ASO→RNA; record)
        mS <- .find_matches_in_subject(DNAString(sense), subj, 0L)
        if (nrow(mS)) {
          for (j in seq_len(nrow(mS))) {
            candidates[[length(candidates) + 1]] <- data.frame(
              aso_id = aso_id, match_class = "exact_transcript_sense",
              exact_shank3_match = TRUE, approximate_shank3_match = FALSE,
              feature_type = "transcript_sense_identical",
              transcript_id = tn, gene_symbol = "SHANK3",
              chrom = annot$tx_meta$seqnames[annot$tx_meta$transcript_id == tn][1],
              start = mS$start[[j]], end = mS$end[[j]],
              strand = annot$tx_meta$strand[annot$tx_meta$transcript_id == tn][1],
              mismatches = 0L, query_coord_space = "transcript",
              support_note = "Exact sense match (same sequence as transcript; atypical ASO binding geometry)",
              neuronal_expression_support = "SHANK3_neuronal_gene",
              isoform_shared = NA, rank_score = 60,
              status = "ok", warning = "sense_match_not_rc",
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    # Genomic locus search (both strands) for RC and sense
    if (!is.null(annot)) {
      loc <- annot$locus_seq
      gene_gr <- annot$gene
      locus_gr <- annot$locus
      for (qname in c("rc", "sense")) {
        q <- if (qname == "rc") rc else sense
        mG <- .find_matches_in_subject(DNAString(q), loc, 0L)
        if (nrow(mG)) {
          if (qname == "rc") found_exact <- TRUE
          for (j in seq_len(nrow(mG))) {
            # Map local coords to genome
            g_start <- start(locus_gr) + mG$start[[j]] - 1L
            g_end <- start(locus_gr) + mG$end[[j]] - 1L
            # Feature type relative to gene
            in_gene <- g_start >= start(gene_gr) && g_end <= end(gene_gr)
            feat <- if (in_gene) "genomic_gene_body_or_intronic" else "genomic_flanking"
            candidates[[length(candidates) + 1]] <- data.frame(
              aso_id = aso_id,
              match_class = paste0("exact_genomic_", qname),
              exact_shank3_match = in_gene || qname == "rc",
              approximate_shank3_match = FALSE,
              feature_type = feat,
              transcript_id = NA_character_, gene_symbol = "SHANK3",
              chrom = as.character(seqnames(gene_gr)),
              start = as.integer(g_start), end = as.integer(g_end),
              strand = as.character(strand(gene_gr)),
              mismatches = 0L, query_coord_space = "genomic",
              support_note = "Exact match in SHANK3 genomic locus window (±5kb)",
              neuronal_expression_support = "SHANK3_neuronal_gene",
              isoform_shared = NA, rank_score = if (qname == "rc") 90 else 50,
              status = "ok", warning = "",
              stringsAsFactors = FALSE
            )
            w_start <- max(1L, mG$start[[j]] - up)
            w_end <- min(length(loc), mG$end[[j]] + down)
            win <- subseq(loc, w_start, w_end)
            wid <- paste0(aso_id, "|genomic|", as.character(seqnames(gene_gr)), ":",
                          start(locus_gr) + w_start - 1, "-", start(locus_gr) + w_end - 1)
            window_seqs[[wid]] <- win
            window_beds[[length(window_beds) + 1]] <- data.frame(
              chrom = as.character(seqnames(gene_gr)),
              start = as.integer(start(locus_gr) + w_start - 1L - 1L),
              end = as.integer(start(locus_gr) + w_end - 1L),
              name = wid, score = 900, strand = as.character(strand(gene_gr)),
              aso_id = aso_id, match_class = paste0("exact_genomic_", qname),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    # Off-target screen on chr22
    ot <- tryCatch(.run_offtarget_screen(rc, NULL), error = function(e) data.frame())
    if (nrow(ot)) {
      ot$aso_id <- aso_id
      offtarget_rows[[length(offtarget_rows) + 1]] <- ot
    }

    if (!length(candidates)) {
      map_rows[[length(map_rows) + 1]] <- data.frame(
        aso_id = aso_id, match_class = "no_exact_shank3_match",
        exact_shank3_match = FALSE, approximate_shank3_match = found_approx,
        feature_type = NA_character_, transcript_id = NA_character_,
        gene_symbol = "SHANK3", chrom = if (!is.null(annot)) as.character(seqnames(annot$gene)) else NA_character_,
        start = NA_integer_, end = NA_integer_, strand = NA_character_,
        mismatches = NA_integer_, query_coord_space = NA_character_,
        support_note = ambiguity_statement(cfg),
        neuronal_expression_support = "unknown",
        isoform_shared = NA, rank_score = 0,
        status = "no_exact_match", warning = "ambiguity_statement_required",
        stringsAsFactors = FALSE
      )
    } else {
      cdf <- bind_rows(candidates)
      # Rank: exactness, support, neuronal expression, feature type
      cdf <- cdf %>%
        mutate(
          feature_rank = case_when(
            grepl("transcript", feature_type) ~ 3,
            grepl("gene_body", feature_type) ~ 2,
            TRUE ~ 1
          ),
          rank_score = rank_score + feature_rank
        ) %>%
        arrange(desc(rank_score), mismatches, match_class) %>%
        mutate(rank = row_number())
      # isoform sharing: same start across transcripts
      if (sum(cdf$exact_shank3_match, na.rm = TRUE) > 1) {
        cdf$isoform_shared <- duplicated(paste(cdf$start, cdf$end)) | duplicated(paste(cdf$start, cdf$end), fromLast = TRUE)
      }
      map_rows[[length(map_rows) + 1]] <- cdf
    }
  }

  mapping <- bind_rows(map_rows)
  out_map <- file.path(root, "results/02_shank3_mapping.tsv")
  write_tsv(mapping, out_map)

  out_fa <- file.path(root, "results/02_target_windows.fasta")
  if (length(window_seqs)) {
    ss <- DNAStringSet(window_seqs)
    writeXStringSet(ss, out_fa)
  } else {
    writeLines(">no_windows\n", out_fa)
  }

  out_bed <- file.path(root, "results/02_target_windows.bed")
  if (length(window_beds)) {
    bed <- bind_rows(window_beds)
    # BED: chrom start end name score strand
    write.table(bed[, c("chrom", "start", "end", "name", "score", "strand")],
                out_bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    writeLines("", out_bed)
  }

  out_ot <- file.path(root, "results/02_offtargets_chr22.tsv")
  if (length(offtarget_rows)) write_tsv(bind_rows(offtarget_rows), out_ot) else {
    write_tsv(data.frame(
      aso_id = character(), chrom = character(), start = integer(), end = integer(),
      strand = character(), match_type = character(),
      status = character(), warning = character()
    ), out_ot)
  }

  # Annotation provenance summary
  write_tsv(data.frame(
    resource = c("GENCODE", "TxDb", "BSgenome"),
    version = c(
      gencode_dl$version %||% cfg$resources$gencode$version,
      if (!is.null(annot)) annot$annotation_version else NA,
      as.character(utils::packageVersion("BSgenome.Hsapiens.UCSC.hg38"))
    ),
    url = c(gencode_dl$url %||% cfg$resources$gencode$url, "Bioconductor", "Bioconductor"),
    genome_build = "GRCh38",
    download_date = as.character(Sys.Date()),
    gencode_status = gencode_dl$status,
    gencode_warning = gencode_dl$warning %||% "",
    stringsAsFactors = FALSE
  ), file.path(root, "data/processed/02_annotation_provenance.tsv"))

  n_exact <- sum(mapping$exact_shank3_match, na.rm = TRUE)
  if (n_exact == 0) {
    log_message("AMBIGUITY: ", ambiguity_statement(cfg))
  }
  log_message("Step 02 mapping complete. exact rows=", n_exact)

  list(
    mapping = mapping,
    windows_fasta = out_fa,
    windows_bed = out_bed,
    offtargets = out_ot,
    tsv = out_map,
    gencode = gencode_dl,
    annot_ok = !is.null(annot)
  )
}
