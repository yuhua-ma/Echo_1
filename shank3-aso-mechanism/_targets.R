# _targets.R — SHANK3 ASO MoA workflow orchestration (R / targets)
# All analysis in R. ViennaRNA only via system2().

library(targets)

tar_option_set(
  packages = c(
    "yaml", "dplyr", "tidyr", "readr", "ggplot2", "glue", "Biostrings",
    "GenomicRanges", "IRanges", "jsonlite"
  ),
  format = "rds",
  error = "continue"
)

# Source modules relative to project root (directory containing _targets.R)
lapply(list.files("R", pattern = "[.]R$", full.names = TRUE), source)

list(
  tar_target(config, {
    cfg <- load_config()
    ensure_dirs(cfg$project_root)
    cfg
  }),

  tar_target(step01_sequence_qc, run_sequence_qc(config)),

  tar_target(step02_mapping, run_transcript_mapping(config, step01_sequence_qc)),

  tar_target(step03_structure, run_rna_structure(config, step02_mapping, step01_sequence_qc)),

  tar_target(step04_rbp, run_rbp_analysis(config, step02_mapping, step01_sequence_qc)),

  tar_target(step05_tf, run_tf_chromatin(config, step02_mapping)),

  # Structure is optional input to ranking (motif/CLIP driven); still sequenced after step03 when available
  tar_target(step06_rbp_rank, {
    struct <- tryCatch(step03_structure, error = function(e) NULL)
    rank_rbps(config, step04_rbp, struct)
  }),

  tar_target(step07_mechanisms, rank_mechanisms(
    config, step01_sequence_qc, step02_mapping, step03_structure,
    step04_rbp, step06_rbp_rank, step05_tf
  )),

  tar_target(step07_validation, build_validation_table(config, step06_rbp_rank)),

  tar_target(step07_summary, print_pipeline_summary(
    config, step01_sequence_qc, step02_mapping, step06_rbp_rank,
    step07_mechanisms, step03_structure
  )),

  tar_target(report_index, write_report_index(config)),

  tar_target(render_reports, {
    # Depend on upstream results
    list(step07_mechanisms, step07_validation, step07_summary, report_index)
    root <- config$project_root
    rmds <- c(
      "report/target_mapping_report.Rmd",
      "report/structure_report.Rmd",
      "report/rbp_evidence_report.Rmd",
      "report/tf_chromatin_report.Rmd",
      "report/final_report.Rmd",
      "report/FINAL_SHANK3_ASO_MECHANISM_REPORT.Rmd"
    )
    outs <- list()
    for (rmd in rmds) {
      path <- file.path(root, rmd)
      if (!file.exists(path)) next
      outs[[rmd]] <- tryCatch({
        rmarkdown::render(
          path,
          output_format = "html_document",
          quiet = TRUE,
          envir = new.env(parent = globalenv())
        )
      }, error = function(e) {
        log_message("Render failed for ", rmd, ": ", conditionMessage(e))
        conditionMessage(e)
      })
    }
    # Final report also PDF if possible
    final <- file.path(root, "report/FINAL_SHANK3_ASO_MECHANISM_REPORT.Rmd")
    pdf_out <- tryCatch({
      rmarkdown::render(
        final,
        output_format = "pdf_document",
        quiet = TRUE,
        envir = new.env(parent = globalenv())
      )
    }, error = function(e) {
      log_message("PDF render failed: ", conditionMessage(e))
      conditionMessage(e)
    })
    list(html = outs, pdf = pdf_out)
  })
)
