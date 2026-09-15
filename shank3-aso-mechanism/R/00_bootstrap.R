# R/00_bootstrap.R — package bootstrap + shared helpers
# Chemistry-aware SHANK3 ASO MoA workflow (R only)

suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Run bootstrap_packages().")
  }
})

bootstrap_packages <- function(lib = Sys.getenv("R_LIBS_USER", unset = .libPaths()[[1]])) {
  cran <- c(
    "targets", "tarchetypes", "renv", "yaml", "tidyverse", "data.table",
    "httr2", "jsonlite", "readr", "glue", "ggplot2", "patchwork", "gt",
    "knitr", "rmarkdown", "BiocManager", "here", "fs", "withr", "dplyr",
    "tidyr", "purrr", "tibble", "stringr"
  )
  bioc <- c(
    "Biostrings", "GenomicRanges", "IRanges", "rtracklayer", "GenomicFeatures",
    "Rsamtools", "BiocFileCache", "AnnotationHub", "biomaRt",
    "BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db", "GenomeInfoDb"
  )
  inst <- rownames(installed.packages(lib.loc = lib))
  need_cran <- setdiff(cran, inst)
  if (length(need_cran)) {
    install.packages(need_cran, repos = "https://cloud.r-project.org", lib = lib, Ncpus = 2)
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = lib)
  }
  inst <- rownames(installed.packages(lib.loc = lib))
  need_bioc <- setdiff(bioc, inst)
  if (length(need_bioc)) {
    BiocManager::install(need_bioc, ask = FALSE, update = FALSE, lib = lib, Ncpus = 2)
  }
  invisible(TRUE)
}

project_root <- function() {
  # Prefer directory containing _targets.R
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(wd, "_targets.R"))) return(wd)
  if (file.exists(file.path(wd, "shank3-aso-mechanism", "_targets.R"))) {
    return(file.path(wd, "shank3-aso-mechanism"))
  }
  # Walk up
  cur <- wd
  for (i in 1:6) {
    if (file.exists(file.path(cur, "_targets.R"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  wd
}

load_config <- function(path = NULL) {
  root <- project_root()
  path <- path %||% file.path(root, "_config.yml")
  cfg <- yaml::read_yaml(path)
  stopifnot(!is.null(cfg$asos), length(cfg$asos) >= 1)
  # Guard against accidental key typo
  if (!is.null(cfg[["as os"]])) {
    stop("Config key 'as os' is invalid; use 'asos'.")
  }
  cfg$project_root <- root
  cfg$download_date <- as.character(Sys.Date())
  cfg
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (is.character(a) && !nzchar(a[[1]]))) b else a

ensure_dirs <- function(root = project_root()) {
  dirs <- c(
    "data/raw", "data/processed", "results", "figures", "figures/03_structure",
    "logs", "report"
  )
  for (d in dirs) dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
  invisible(dirs)
}

append_provenance <- function(resource_name, version, url, genome_build = "GRCh38",
                              notes = "", root = project_root()) {
  ensure_dirs(root)
  path <- file.path(root, "data/raw/PROVENANCE.tsv")
  row <- data.frame(
    resource = resource_name,
    version = version %||% "unknown",
    url = url %||% "",
    genome_build = genome_build,
    download_date = as.character(Sys.Date()),
    download_datetime_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = notes,
    stringsAsFactors = FALSE
  )
  if (file.exists(path)) {
    old <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    # replace same resource+url if re-downloaded
    keep <- !(old$resource == row$resource & old$url == row$url)
    utils::write.table(rbind(old[keep, , drop = FALSE], row), path,
                       sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    utils::write.table(row, path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  path
}

safe_download <- function(url, destfile, resource_name, version, genome_build = "GRCh38",
                          notes = "", timeout = 120) {
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  status <- "ok"
  warning_msg <- ""
  ok <- FALSE
  tryCatch({
    old <- getOption("timeout")
    on.exit(options(timeout = old), add = TRUE)
    options(timeout = timeout)
    utils::download.file(url, destfile = destfile, mode = "wb", quiet = TRUE)
    ok <- file.exists(destfile) && file.info(destfile)$size > 0
    if (!ok) {
      status <- "empty_download"
      warning_msg <- "Downloaded file missing or empty"
    }
  }, error = function(e) {
    status <<- "download_failed"
    warning_msg <<- conditionMessage(e)
  })
  append_provenance(resource_name, version, url, genome_build,
                    notes = paste(notes, status, warning_msg))
  list(ok = ok, path = destfile, status = status, warning = warning_msg, url = url)
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(df, path)
  path
}

log_message <- function(..., root = project_root()) {
  ensure_dirs(root)
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(root, "logs/pipeline.log"), append = TRUE)
  invisible(msg)
}

chemistry_annotation <- function(cfg) {
  list(
    backbone = cfg$chemistry$backbone %||% "phosphorothioate",
    moe = isTRUE(cfg$chemistry$moe),
    five_methyl_dC = isTRUE(cfg$chemistry$five_methyl_dC),
    architecture = cfg$chemistry$architecture %||% "unknown",
    rnase_h_compatible = cfg$chemistry$rnase_h_compatible %||% "unknown",
    note = paste(
      "All cytosines annotated as 5-methyl-dC.",
      "Architecture unknown — do not assume gapmer design.",
      "RNase-H compatibility unknown — do not assume RNase-H activity."
    )
  )
}

ambiguity_statement <- function(cfg) {
  cfg$ambiguity_statement %||% paste(
    "No exact SHANK3 transcript match was identified.",
    "Do not invent a transcript; mechanism interpretation remains provisional."
  )
}

empty_result <- function(cols, status = "unavailable", warning = "") {
  df <- as.data.frame(matrix(nrow = 0, ncol = length(cols)))
  names(df) <- cols
  if ("status" %in% cols) {
    # leave empty; callers may rbind rows with status
  }
  attr(df, "status") <- status
  attr(df, "warning") <- warning
  df
}
