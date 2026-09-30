#!/usr/bin/env Rscript
# qPCR bar plots: Cortex (CTX) + Striatum (STR)
# Bars = mean; error bars = SEM; points = individual replicates (jitter).
# Values transcribed from screenshot tables (European commas → decimal points).
# Units: relative / normalized expression as shown in source tables (assay scale unknown).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# Paths: resolve relative to this script when possible
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

csv_path <- file.path(script_dir, "qpcr_ctx_str_long.csv")
# Prefer Project-store media/qpcr when present; else write PNGs next to this script.
media_candidates <- c(
  file.path(script_dir, "..", "..", "..", "media", "qpcr"),
  file.path(dirname(script_dir), "..", "..", "media", "qpcr"),
  file.path(getwd(), "media", "qpcr")
)
out_dir <- NULL
for (cand in media_candidates) {
  cand_norm <- suppressWarnings(normalizePath(cand, mustWork = FALSE))
  # Only use a candidate if its parent media/ (or media itself) already exists
  if (dir.exists(dirname(cand_norm)) && basename(dirname(cand_norm)) == "media") {
    out_dir <- cand_norm
    break
  }
  if (dir.exists(cand_norm)) {
    out_dir <- cand_norm
    break
  }
}
if (is.null(out_dir)) out_dir <- script_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("HT-Veh", "HT-New SCR", "HT-2,6", "WT-Veh", "WT-New SCR", "WT-2,6")

dat <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(
    group = factor(group, levels = group_levels),
    region = factor(region, levels = c("CTX", "STR")),
    region_label = factor(
      ifelse(region == "CTX", "Cortex (CTX)", "Striatum (STR)"),
      levels = c("Cortex (CTX)", "Striatum (STR)")
    )
  )

stopifnot(!any(is.na(dat$expression)))
stopifnot(all(levels(dat$group) %in% group_levels))

# Mean ± SEM by region × group
summary_tbl <- dat %>%
  group_by(region, region_label, group) %>%
  summarise(
    n = n(),
    mean = mean(expression),
    sd = sd(expression),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

# Optional summary CSV beside the long table
summary_path <- file.path(script_dir, "qpcr_ctx_str_mean_sem.csv")
write_csv(summary_tbl %>% select(region, group, n, mean, sd, sem), summary_path)

# Fill colors: slate / teal / ochre (avoid purple-default AI look)
fill_cols <- c(
  "HT-Veh" = "#4A5568",
  "HT-New SCR" = "#718096",
  "HT-2,6" = "#2C7A7B",
  "WT-Veh" = "#C05621",
  "WT-New SCR" = "#DD6B20",
  "WT-2,6" = "#B7791F"
)

theme_qpcr <- function() {
  theme_classic(base_size = 13, base_family = "sans") %+replace%
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(size = 10, color = "#4A5568", hjust = 0, margin = margin(b = 10)),
      plot.caption = element_text(size = 8, color = "#718096", hjust = 0, margin = margin(t = 10)),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, margin = margin(r = 8)),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 10, color = "#1A202C"),
      axis.text.y = element_text(size = 10, color = "#1A202C"),
      axis.line = element_line(color = "#2D3748", linewidth = 0.4),
      axis.ticks = element_line(color = "#2D3748", linewidth = 0.35),
      panel.grid.major.y = element_line(color = "#EDF2F7", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(12, 14, 10, 12)
    )
}

make_region_plot <- function(region_code, title_text) {
  d <- dat %>% filter(region == region_code)
  s <- summary_tbl %>% filter(region == region_code)
  y_max <- max(d$expression, s$mean + s$sem, na.rm = TRUE) * 1.12

  ggplot() +
    geom_col(
      data = s,
      aes(x = group, y = mean, fill = group),
      width = 0.72,
      color = "#1A202C",
      linewidth = 0.25,
      alpha = 0.92
    ) +
    geom_errorbar(
      data = s,
      aes(x = group, ymin = mean - sem, ymax = mean + sem),
      width = 0.22,
      linewidth = 0.55,
      color = "#1A202C"
    ) +
    geom_jitter(
      data = d,
      aes(x = group, y = expression),
      width = 0.12,
      height = 0,
      size = 2.1,
      shape = 21,
      fill = "white",
      color = "#1A202C",
      stroke = 0.55,
      alpha = 0.95
    ) +
    scale_fill_manual(values = fill_cols, breaks = group_levels) +
    scale_x_discrete(limits = group_levels, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, y_max)) +
    labs(
      title = title_text,
      subtitle = "Bars = mean; error bars = SEM; points = individual animals (unequal n)",
      y = "Normalized expression",
      caption = paste0(
        "Assay units as in source spreadsheet (relative/normalized scale; absolute units unknown). ",
        "Transcribed from ", unique(d$source_file), "."
      )
    ) +
    theme_qpcr()
}

p_ctx <- make_region_plot("CTX", "qPCR — Cortex (CTX)")
p_str <- make_region_plot("STR", "qPCR — Striatum (STR)")

# Combined faceted panel
p_combined <- ggplot() +
  geom_col(
    data = summary_tbl,
    aes(x = group, y = mean, fill = group),
    width = 0.72,
    color = "#1A202C",
    linewidth = 0.25,
    alpha = 0.92
  ) +
  geom_errorbar(
    data = summary_tbl,
    aes(x = group, ymin = mean - sem, ymax = mean + sem),
    width = 0.22,
    linewidth = 0.55,
    color = "#1A202C"
  ) +
  geom_jitter(
    data = dat,
    aes(x = group, y = expression),
    width = 0.12,
    height = 0,
    size = 1.9,
    shape = 21,
    fill = "white",
    color = "#1A202C",
    stroke = 0.5,
    alpha = 0.95
  ) +
  facet_wrap(~region_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = fill_cols, breaks = group_levels) +
  scale_x_discrete(limits = group_levels, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "qPCR — Cortex (CTX) vs Striatum (STR)",
    subtitle = "Bars = mean; error bars = SEM; points = individual animals (unequal n)",
    y = "Normalized expression",
    caption = "Assay units as in source spreadsheet (relative/normalized scale; absolute units unknown)."
  ) +
  theme_qpcr() +
  theme(
    strip.background = element_rect(fill = "#F7FAFC", color = NA),
    strip.text = element_text(face = "bold", size = 12, margin = margin(4, 4, 4, 4)),
    panel.spacing = unit(1.2, "lines")
  )

# Also stack via patchwork as an alternate combined layout
p_stack <- p_ctx / p_str + plot_annotation(
  caption = "Assay units as in source spreadsheet (relative/normalized scale; absolute units unknown)."
)

out_ctx <- file.path(out_dir, "qpcr_CTX_barplot.png")
out_str <- file.path(out_dir, "qpcr_STR_barplot.png")
out_comb <- file.path(out_dir, "qpcr_CTX_STR_combined.png")

ggsave(out_ctx, p_ctx, width = 8.2, height = 5.4, dpi = 300, bg = "white")
ggsave(out_str, p_str, width = 8.2, height = 5.4, dpi = 300, bg = "white")
ggsave(out_comb, p_combined, width = 12.5, height = 5.6, dpi = 300, bg = "white")

message("Wrote: ", out_ctx)
message("Wrote: ", out_str)
message("Wrote: ", out_comb)
message("Summary: ", summary_path)
print(as.data.frame(summary_tbl %>% select(region, group, n, mean, sem)))
