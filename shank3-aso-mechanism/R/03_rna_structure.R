# R/03_rna_structure.R — ViennaRNA via system2() only
# If RNAfold/RNAplfold missing: record command, warn, continue — never fabricate.

.vienna_available <- function() {
  fold <- nzchar(Sys.which("RNAfold"))
  plfold <- nzchar(Sys.which("RNAplfold"))
  list(
    RNAfold = fold,
    RNAplfold = plfold,
    RNAfold_path = Sys.which("RNAfold"),
    RNAplfold_path = Sys.which("RNAplfold")
  )
}

.run_rnafold <- function(seq_chr, workdir, prefix) {
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  fa <- file.path(workdir, paste0(prefix, ".fa"))
  # RNA alphabet
  rna <- chartr("T", "U", toupper(seq_chr))
  writeLines(c(paste0(">", prefix), rna), fa)
  cmd <- "RNAfold"
  args <- c("--noPS", "-i", fa)
  recorded <- paste(cmd, paste(args, collapse = " "))
  out <- tryCatch({
    # system2() has no wd=; use absolute -i path
    system2(cmd, args = args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  list(command = recorded, stdout = out, sequence = rna)
}

.run_rnaplfold <- function(seq_chr, workdir, prefix, winsize = 70L, span = 40L) {
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  rna <- chartr("T", "U", toupper(seq_chr))
  fa <- file.path(workdir, paste0(prefix, "_pl.fa"))
  writeLines(c(paste0(">", prefix), rna), fa)
  cmd <- "RNAplfold"
  args <- c("-W", as.character(winsize), "-L", as.character(span), "-u", "1")
  recorded <- paste(cmd, paste(args, collapse = " "), "<", fa)
  # RNAplfold writes lunp into the working directory
  out <- tryCatch({
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(workdir)
    system2(cmd, args = args, stdout = TRUE, stderr = TRUE, stdin = basename(fa))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  lunp <- list.files(workdir, pattern = paste0("^", prefix, ".*_lunp$"), full.names = TRUE)
  if (!length(lunp)) lunp <- list.files(workdir, pattern = "_lunp$", full.names = TRUE)
  unpaired <- NA_real_
  if (length(lunp)) {
    # Prefer newest lunp for this prefix
    lunp <- lunp[order(file.info(lunp)$mtime, decreasing = TRUE)]
    tab <- tryCatch(utils::read.table(lunp[[1]], header = FALSE, comment.char = "#"),
                    error = function(e) NULL)
    if (!is.null(tab) && ncol(tab) >= 2) {
      unpaired <- mean(as.numeric(tab[[2]]), na.rm = TRUE)
    }
  }
  list(command = recorded, stdout = out, mean_unpaired_prob = unpaired, lunp_files = lunp)
}

.classify_site_context <- function(dotbracket, site_start, site_end) {
  if (is.na(dotbracket) || !nzchar(dotbracket)) return("unknown")
  db <- strsplit(dotbracket, "")[[1]]
  s <- max(1, site_start)
  e <- min(length(db), site_end)
  region <- db[s:e]
  if (all(region == ".")) return("loop_or_unpaired")
  if (all(region %in% c("(", ")"))) return("stem")
  if (any(region == ".") && any(region %in% c("(", ")"))) return("bulge_or_junction")
  "mixed"
}

.parse_rnafold_output <- function(stdout_lines) {
  # Typical: sequence line then structure + energy, e.g. .((...)). ( -0.20)
  lines <- as.character(stdout_lines)
  lines <- lines[!grepl("^ERROR", lines)]
  # Find lines that look like dot-bracket + energy
  is_struct <- grepl("^[.()]+[[:space:]]+\\([[:space:]]*[-0-9.]+[[:space:]]*\\)[[:space:]]*$", lines)
  # RNAcofold may include & in structure token region — handle trailing energy only
  if (!any(is_struct)) {
    is_struct <- grepl("\\([[:space:]]*[-0-9.]+[[:space:]]*\\)[[:space:]]*$", lines) &
      grepl("[.()]", lines)
  }
  struct_line <- lines[is_struct]
  if (!length(struct_line)) {
    return(list(structure = NA_character_, mfe = NA_real_))
  }
  sl <- trimws(struct_line[[length(struct_line)]])
  energy_match <- regexpr("\\([[:space:]]*([-0-9.]+)[[:space:]]*\\)[[:space:]]*$", sl, perl = TRUE)
  if (energy_match[[1]] == -1) {
    return(list(structure = NA_character_, mfe = NA_real_))
  }
  cap <- attr(energy_match, "capture.start")
  clen <- attr(energy_match, "capture.length")
  mfe <- as.numeric(substr(sl, cap[[1]], cap[[1]] + clen[[1]] - 1L))
  structure <- trimws(substr(sl, 1, energy_match[[1]] - 1L))
  # Drop trailing hybrid ampersand artifacts if present
  structure <- sub("[[:space:]]*&[[:space:]]*$", "", structure)
  list(structure = structure, mfe = mfe)
}

.run_hybridization_approx <- function(aso_rna, target_rna, workdir, prefix) {
  # Approximate duplex via RNAcofold; NOT exact with MOE/PS/5m-dC
  if (nzchar(Sys.which("RNAcofold"))) {
    fa <- file.path(workdir, paste0(prefix, "_co.fa"))
    writeLines(c(paste0(">", prefix), paste0(aso_rna, "&", target_rna)), fa)
    out <- tryCatch(
      system2("RNAcofold", c("--noPS"), stdin = fa, stdout = TRUE, stderr = TRUE),
      error = function(e) paste("ERROR", conditionMessage(e))
    )
    parsed <- .parse_rnafold_output(out)
    return(list(
      command = paste("RNAcofold --noPS <", fa),
      energy = parsed$mfe,
      note = "approx_RNAcofold_canonical_RNA"
    ))
  }
  res <- .run_rnafold(paste0(aso_rna, "AAAACCCC", target_rna), workdir, paste0(prefix, "_hyb"))
  parsed <- .parse_rnafold_output(res$stdout)
  list(command = res$command, energy = parsed$mfe, note = "approx_concat_fold_not_true_duplex")
}

run_rna_structure <- function(cfg, map_result, qc_result) {
  suppressPackageStartupMessages({
    library(Biostrings)
    library(dplyr)
    library(ggplot2)
  })
  root <- cfg$project_root
  ensure_dirs(root)
  figdir <- file.path(root, "figures/03_structure")
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
  work <- file.path(root, "data/processed/vienna")
  dir.create(work, recursive = TRUE, showWarnings = FALSE)

  vien <- .vienna_available()
  warn_file <- file.path(root, "logs/vienna_status.txt")
  writeLines(c(
    paste("RNAfold available:", vien$RNAfold, vien$RNAfold_path),
    paste("RNAplfold available:", vien$RNAplfold, vien$RNAplfold_path),
    paste("checked_at:", format(Sys.time(), tz = "UTC", usetz = TRUE))
  ), warn_file)

  if (!vien$RNAfold) {
    log_message("WARNING: RNAfold missing — structure results will be empty; commands recorded.")
  }

  windows_fa <- map_result$windows_fasta
  struct_rows <- list()
  acc_rows <- list()
  commands <- list()

  has_windows <- file.exists(windows_fa) && length(readLines(windows_fa, warn = FALSE)) > 2
  seqs <- if (has_windows) {
    tryCatch(readDNAStringSet(windows_fa), error = function(e) DNAStringSet())
  } else DNAStringSet()
  # drop placeholder
  if (length(seqs) && grepl("^no_windows", names(seqs)[1])) seqs <- DNAStringSet()

  qc <- qc_result$qc

  # ASO self-structure
  for (i in seq_len(nrow(qc))) {
    if (!isTRUE(qc$valid_acgt[[i]])) next
    aso_id <- qc$aso_id[[i]]
    aso_rna <- chartr("T", "U", qc$cleaned_sequence[[i]])
    prefix <- paste0("aso_self_", aso_id)
    if (vien$RNAfold) {
      res <- .run_rnafold(aso_rna, work, prefix)
      parsed <- .parse_rnafold_output(res$stdout)
      commands[[length(commands) + 1]] <- data.frame(
        aso_id = aso_id, context = "aso_self", command = res$command,
        status = "ran", stringsAsFactors = FALSE
      )
      struct_rows[[length(struct_rows) + 1]] <- data.frame(
        aso_id = aso_id, window_id = prefix, context = "aso_self",
        sequence_rna = aso_rna,
        structure_dotbracket = parsed$structure %||% NA_character_,
        mfe_kcal_mol = parsed$mfe %||% NA_real_,
        site_structure_class = NA_character_,
        mean_unpaired_probability = NA_real_,
        approx_hybridization_energy = NA_real_,
        local_structural_change_note = "ASO self-fold only; MOE/PS/5m-dC not modeled",
        chemistry_caveat = "Values are canonical-RNA approximations; not exact biochemical measurements for MOE/PS/5m-dC ASOs",
        status = "ok", warning = "",
        stringsAsFactors = FALSE
      )
    } else {
      cmd <- paste("RNAfold --noPS -i", file.path(work, paste0(prefix, ".fa")))
      commands[[length(commands) + 1]] <- data.frame(
        aso_id = aso_id, context = "aso_self", command = cmd,
        status = "skipped_missing_RNAfold", stringsAsFactors = FALSE
      )
      struct_rows[[length(struct_rows) + 1]] <- data.frame(
        aso_id = aso_id, window_id = prefix, context = "aso_self",
        sequence_rna = aso_rna,
        structure_dotbracket = NA_character_,
        mfe_kcal_mol = NA_real_,
        site_structure_class = NA_character_,
        mean_unpaired_probability = NA_real_,
        approx_hybridization_energy = NA_real_,
        local_structural_change_note = NA_character_,
        chemistry_caveat = "ViennaRNA missing; no fabricated structure",
        status = "vienna_missing", warning = "RNAfold not found",
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(seqs)) {
    for (nm in names(seqs)) {
      aso_id <- sub("\\|.*", "", nm)
      seq_chr <- as.character(seqs[[nm]])
      prefix <- gsub("[^A-Za-z0-9_]+", "_", nm)
      site_len <- nchar(qc$cleaned_sequence[qc$aso_id == aso_id][1] %||% "18")
      # site roughly centered: window = up + site + down; site starts at up+1
      up <- as.integer(cfg$window_upstream %||% 200)
      site_start <- up + 1L
      site_end <- up + as.integer(site_len)

      if (vien$RNAfold) {
        res <- .run_rnafold(seq_chr, work, prefix)
        parsed <- .parse_rnafold_output(res$stdout)
        commands[[length(commands) + 1]] <- data.frame(
          aso_id = aso_id, context = "window", command = res$command,
          status = "ran", stringsAsFactors = FALSE
        )
        site_class <- .classify_site_context(parsed$structure, site_start, site_end)
        # hybridization approx
        aso_rna <- chartr("T", "U", qc$cleaned_sequence[qc$aso_id == aso_id][[1]])
        site_seq <- substr(chartr("T", "U", seq_chr), site_start, min(nchar(seq_chr), site_end))
        hyb <- .run_hybridization_approx(aso_rna, site_seq, work, paste0(prefix, "_hyb"))
        commands[[length(commands) + 1]] <- data.frame(
          aso_id = aso_id, context = "hybridization_approx", command = hyb$command,
          status = "ran", stringsAsFactors = FALSE
        )
        unpaired <- NA_real_
        if (vien$RNAplfold) {
          pl <- .run_rnaplfold(seq_chr, work, prefix)
          unpaired <- pl$mean_unpaired_prob
          commands[[length(commands) + 1]] <- data.frame(
            aso_id = aso_id, context = "RNAplfold", command = pl$command,
            status = "ran", stringsAsFactors = FALSE
          )
        } else {
          commands[[length(commands) + 1]] <- data.frame(
            aso_id = aso_id, context = "RNAplfold", command = "RNAplfold -W 70 -L 40 -u 1",
            status = "skipped_missing_RNAplfold", stringsAsFactors = FALSE
          )
        }
        struct_rows[[length(struct_rows) + 1]] <- data.frame(
          aso_id = aso_id, window_id = nm, context = "target_window",
          sequence_rna = chartr("T", "U", seq_chr),
          structure_dotbracket = parsed$structure %||% NA_character_,
          mfe_kcal_mol = parsed$mfe %||% NA_real_,
          site_structure_class = site_class,
          mean_unpaired_probability = unpaired,
          approx_hybridization_energy = hyb$energy %||% NA_real_,
          local_structural_change_note = paste(
            "Compare window MFE to ASO-bound approx (", hyb$note, ").",
            "Not a precise delta-delta-G for modified ASO chemistry."
          ),
          chemistry_caveat = "Canonical RNA ViennaRNA model; MOE/PS/5m-dC not represented",
          status = "ok", warning = "",
          stringsAsFactors = FALSE
        )
        acc_rows[[length(acc_rows) + 1]] <- data.frame(
          aso_id = aso_id, window_id = nm,
          site_start = site_start, site_end = site_end,
          site_structure_class = site_class,
          mean_unpaired_probability = unpaired,
          window_mfe = parsed$mfe %||% NA_real_,
          status = "ok", warning = if (vien$RNAplfold) "" else "RNAplfold_missing",
          stringsAsFactors = FALSE
        )
        # Simple plot of unpaired if available / structure length
        if (!is.null(parsed$structure) && !is.na(parsed$structure)) {
          db <- strsplit(parsed$structure, "")[[1]]
          dfp <- data.frame(pos = seq_along(db), paired = as.integer(db != "."))
          p <- ggplot(dfp, aes(pos, paired)) +
            geom_col(width = 1) +
            annotate("rect", xmin = site_start, xmax = site_end, ymin = -0.05, ymax = 1.05,
                     alpha = 0.15) +
            labs(title = paste("Structure pairing:", aso_id),
                 subtitle = "Shaded = approximate ASO site; canonical RNA only",
                 x = "Position", y = "Paired (1) / unpaired (0)") +
            theme_bw()
          ggsave(file.path(figdir, paste0(prefix, "_pairing.png")), p, width = 8, height = 3, dpi = 120)
        }
      } else {
        cmd <- paste("RNAfold --noPS -i <window.fa>")
        commands[[length(commands) + 1]] <- data.frame(
          aso_id = aso_id, context = "window", command = cmd,
          status = "skipped_missing_RNAfold", stringsAsFactors = FALSE
        )
        struct_rows[[length(struct_rows) + 1]] <- data.frame(
          aso_id = aso_id, window_id = nm, context = "target_window",
          sequence_rna = chartr("T", "U", seq_chr),
          structure_dotbracket = NA_character_,
          mfe_kcal_mol = NA_real_,
          site_structure_class = NA_character_,
          mean_unpaired_probability = NA_real_,
          approx_hybridization_energy = NA_real_,
          local_structural_change_note = NA_character_,
          chemistry_caveat = "ViennaRNA missing; results not fabricated",
          status = "vienna_missing", warning = "RNAfold not found",
          stringsAsFactors = FALSE
        )
        acc_rows[[length(acc_rows) + 1]] <- data.frame(
          aso_id = aso_id, window_id = nm,
          site_start = NA_integer_, site_end = NA_integer_,
          site_structure_class = NA_character_,
          mean_unpaired_probability = NA_real_,
          window_mfe = NA_real_,
          status = "vienna_missing", warning = "RNAfold not found",
          stringsAsFactors = FALSE
        )
      }
    }
  } else {
    log_message("No target windows for structure analysis")
    struct_rows[[length(struct_rows) + 1]] <- data.frame(
      aso_id = NA_character_, window_id = NA_character_, context = "none",
      sequence_rna = NA_character_, structure_dotbracket = NA_character_,
      mfe_kcal_mol = NA_real_, site_structure_class = NA_character_,
      mean_unpaired_probability = NA_real_, approx_hybridization_energy = NA_real_,
      local_structural_change_note = NA_character_,
      chemistry_caveat = "No windows",
      status = "no_windows", warning = "No target windows extracted",
      stringsAsFactors = FALSE
    )
  }

  struct_df <- bind_rows(struct_rows)
  acc_df <- if (length(acc_rows)) bind_rows(acc_rows) else data.frame(
    aso_id = character(), window_id = character(), site_start = integer(),
    site_end = integer(), site_structure_class = character(),
    mean_unpaired_probability = numeric(), window_mfe = numeric(),
    status = character(), warning = character()
  )
  cmd_df <- if (length(commands)) bind_rows(commands) else data.frame()

  out_s <- file.path(root, "results/03_rna_structure.tsv")
  out_a <- file.path(root, "results/03_accessibility.tsv")
  out_c <- file.path(root, "logs/03_vienna_commands.tsv")
  write_tsv(struct_df, out_s)
  write_tsv(acc_df, out_a)
  write_tsv(cmd_df, out_c)

  log_message("Step 03 structure complete. vienna_fold=", vien$RNAfold)
  list(structure = struct_df, accessibility = acc_df, commands = cmd_df,
       tsv = out_s, accessibility_tsv = out_a, vienna = vien)
}
