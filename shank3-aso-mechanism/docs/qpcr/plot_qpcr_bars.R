#!/usr/bin/env Rscript
# qPCR bar plots: Cortex (CTX) + Striatum (STR)
# Bars = mean; error bars = SEM; points = individual replicates (jitter).
# Stats: within each region × genotype, one-way ANOVA + Tukey HSD across
# treatments (Veh, New SCR, 2.6). Significant pairs only annotated on plots.
# Y-axis: Absolute expression. No bottom caption / “vs” titles.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

csv_path <- file.path(script_dir, "qpcr_ctx_str_long.csv")
media_candidates <- c(
  file.path(script_dir, "..", "..", "..", "media", "qpcr"),
  file.path(dirname(script_dir), "..", "..", "media", "qpcr"),
  file.path(getwd(), "media", "qpcr")
)
out_dir <- NULL
for (cand in media_candidates) {
  cand_norm <- suppressWarnings(normalizePath(cand, mustWork = FALSE))
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
treat_levels <- c("Veh", "New SCR", "2.6")

p_star <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

treat_to_suffix <- function(t) ifelse(as.character(t) == "2.6", "2,6", as.character(t))

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
dat <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(
    group = factor(group, levels = group_levels),
    genotype = factor(genotype, levels = c("HT", "WT")),
    treatment = factor(treatment, levels = treat_levels),
    region = factor(region, levels = c("CTX", "STR")),
    region_label = factor(
      ifelse(region == "CTX", "Cortex (CTX)", "Striatum (STR)"),
      levels = c("Cortex (CTX)", "Striatum (STR)")
    )
  )

stopifnot(!any(is.na(dat$expression)))
stopifnot(all(as.character(unique(dat$group)) %in% group_levels))

# Mean ± SEM
summary_tbl <- dat %>%
  group_by(region, region_label, group, genotype, treatment) %>%
  summarise(
    n = n(),
    mean = mean(expression),
    sd = sd(expression),
    sem = sd / sqrt(n),
    .groups = "drop"
  )

summary_path <- file.path(script_dir, "qpcr_ctx_str_mean_sem.csv")
write_csv(summary_tbl %>% select(region, group, n, mean, sd, sem), summary_path)

# ---------------------------------------------------------------------------
# Stats: region × genotype → one-way ANOVA + Tukey HSD on treatment
# ---------------------------------------------------------------------------
stat_rows <- list()
for (reg in levels(dat$region)) {
  for (gen in levels(dat$genotype)) {
    sub <- dat %>% filter(region == reg, genotype == gen)
    ns <- sub %>% count(treatment, name = "n")
    min_n <- min(ns$n)
    caution <- if (min_n < 3) {
      sprintf(
        "Low power / caution: min n=%d in this region×genotype block (WT-Veh n=2); do not over-interpret ns or borderline p",
        min_n
      )
    } else {
      ""
    }

    fit <- aov(expression ~ treatment, data = sub)
    aov_tab <- summary(fit)[[1]]
    aov_p <- aov_tab[["Pr(>F)"]][1]
    aov_F <- aov_tab[["F value"]][1]
    df1 <- aov_tab[["Df"]][1]
    df2 <- aov_tab[["Df"]][2]

    tuk <- TukeyHSD(fit, "treatment")
    tk <- as.data.frame(tuk$treatment)
    tk$contrast <- rownames(tk)

    for (i in seq_len(nrow(tk))) {
      contrast <- tk$contrast[i]
      t1 <- NA_character_
      t2 <- NA_character_
      for (a in treat_levels) {
        for (b in treat_levels) {
          if (a != b && identical(contrast, paste0(a, "-", b))) {
            t1 <- as.character(a)
            t2 <- as.character(b)
          }
        }
      }
      g1 <- paste0(gen, "-", treat_to_suffix(t1))
      g2 <- paste0(gen, "-", treat_to_suffix(t2))
      # Order groups left-to-right on plot axis for brackets
      pos1 <- match(g1, group_levels)
      pos2 <- match(g2, group_levels)
      if (!is.na(pos1) && !is.na(pos2) && pos1 > pos2) {
        tmp <- g1; g1 <- g2; g2 <- tmp
        tmp <- t1; t1 <- t2; t2 <- tmp
        tmp_n1 <- ns$n[ns$treatment == t1]
        tmp_n2 <- ns$n[ns$treatment == t2]
      } else {
        tmp_n1 <- ns$n[ns$treatment == t1]
        tmp_n2 <- ns$n[ns$treatment == t2]
      }

      stat_rows[[length(stat_rows) + 1]] <- data.frame(
        region = reg,
        region_label = as.character(
          unique(dat$region_label[dat$region == reg])
        ),
        genotype = as.character(gen),
        group1 = g1,
        group2 = g2,
        treatment1 = t1,
        treatment2 = t2,
        contrast = contrast,
        n_group1 = as.integer(tmp_n1),
        n_group2 = as.integer(tmp_n2),
        mean_diff = tk$diff[i],
        anova_df1 = df1,
        anova_df2 = df2,
        anova_F = aov_F,
        anova_p = aov_p,
        test = "one-way ANOVA + Tukey HSD",
        p_adj = tk$`p adj`[i],
        p_adj_method = "Tukey HSD",
        significance = p_star(tk$`p adj`[i]),
        note = caution,
        stringsAsFactors = FALSE
      )
    }
  }
}

stats_tbl <- bind_rows(stat_rows)
stats_path <- file.path(script_dir, "qpcr_ctx_str_stats.csv")
write_csv(stats_tbl, stats_path)

sig_tbl <- stats_tbl %>% filter(significance != "ns", !is.na(significance))

# Bracket y-positions stacked within each region (and genotype block)
build_annotations <- function(sig_df, summary_df, tip_pad = 0.04, step = 0.055) {
  if (nrow(sig_df) == 0) {
    return(data.frame(
      region = character(),
      region_label = character(),
      group1 = character(),
      group2 = character(),
      y_position = numeric(),
      annotations = character(),
      stringsAsFactors = FALSE
    ))
  }
  # Max of mean+SEM and raw points per group for vertical clearance
  tip <- summary_df %>%
    mutate(tip = mean + sem) %>%
    select(region, group, tip)
  raw_max <- dat %>%
    group_by(region, group) %>%
    summarise(raw_max = max(expression), .groups = "drop")
  tips <- tip %>%
    left_join(raw_max, by = c("region", "group")) %>%
    mutate(ceiling = pmax(tip, raw_max, na.rm = TRUE))

  out <- list()
  for (reg in unique(sig_df$region)) {
    sreg <- sig_df %>% filter(region == reg)
    # Sort: shorter spans first, then left position
    sreg <- sreg %>%
      mutate(
        span = match(group2, group_levels) - match(group1, group_levels),
        left = match(group1, group_levels)
      ) %>%
      arrange(left, span)
    base_y <- max(tips$ceiling[tips$region == reg], na.rm = TRUE) + tip_pad
    for (i in seq_len(nrow(sreg))) {
      g1 <- sreg$group1[i]
      g2 <- sreg$group2[i]
      pair_ceil <- max(
        tips$ceiling[tips$region == reg & tips$group %in% c(g1, g2)],
        na.rm = TRUE
      )
      y_pos <- max(base_y, pair_ceil + tip_pad) + (i - 1) * step
      out[[length(out) + 1]] <- data.frame(
        region = reg,
        region_label = sreg$region_label[i],
        group1 = g1,
        group2 = g2,
        y_position = y_pos,
        annotations = sreg$significance[i],
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(out)
}

ann_all <- build_annotations(sig_tbl, summary_tbl)

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
      plot.title = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 8)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
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

# Build bracket geometry from annotation table (reliable with facets)
bracket_layers_data <- function(ann_df, tip_frac = 0.025) {
  if (nrow(ann_df) == 0) {
    return(list(
      segs = data.frame(),
      labs = data.frame()
    ))
  }
  # tip length relative to y span of annotations in each region
  ann_df <- ann_df %>%
    group_by(region) %>%
    mutate(tip = pmax(y_position * tip_frac, 0.012)) %>%
    ungroup() %>%
    mutate(
      x1 = match(group1, group_levels),
      x2 = match(group2, group_levels),
      region_label = factor(region_label, levels = levels(dat$region_label))
    )

  segs <- bind_rows(
    # horizontal bar
    ann_df %>% transmute(region, region_label, x = x1, xend = x2, y = y_position, yend = y_position),
    # left tip
    ann_df %>% transmute(region, region_label, x = x1, xend = x1, y = y_position - tip, yend = y_position),
    # right tip
    ann_df %>% transmute(region, region_label, x = x2, xend = x2, y = y_position - tip, yend = y_position)
  )
  labs <- ann_df %>%
    transmute(
      region, region_label,
      x = (x1 + x2) / 2,
      y = y_position,
      label = annotations
    )
  list(segs = segs, labs = labs)
}

add_bracket_geoms <- function(p, br, use_facet_aes = FALSE) {
  if (nrow(br$segs) == 0) return(p)
  p +
    geom_segment(
      data = br$segs,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.4,
      color = "#1A202C",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = br$labs,
      aes(x = x, y = y, label = label),
      vjust = -0.35,
      size = 4.5,
      color = "#1A202C",
      inherit.aes = FALSE
    )
}

make_region_plot <- function(region_code, title_text) {
  d <- dat %>% filter(region == region_code)
  s <- summary_tbl %>% filter(region == region_code)
  ann <- ann_all %>% filter(region == region_code)
  br <- bracket_layers_data(ann)

  y_data_max <- max(d$expression, s$mean + s$sem, na.rm = TRUE)
  y_ann_max <- if (nrow(ann)) max(ann$y_position) else y_data_max
  y_max <- max(y_data_max * 1.10, y_ann_max * 1.12)

  p <- ggplot() +
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
    labs(title = title_text, y = "Absolute expression") +
    theme_qpcr()

  add_bracket_geoms(p, br)
}

p_ctx <- make_region_plot("CTX", "qPCR — Cortex (CTX)")
p_str <- make_region_plot("STR", "qPCR — Striatum (STR)")

br_all <- bracket_layers_data(ann_all)

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
  {
    if (nrow(br_all$segs) > 0) {
      list(
        geom_segment(
          data = br_all$segs,
          aes(x = x, xend = xend, y = y, yend = yend),
          linewidth = 0.4,
          color = "#1A202C",
          inherit.aes = FALSE
        ),
        geom_text(
          data = br_all$labs,
          aes(x = x, y = y, label = label),
          vjust = -0.35,
          size = 4.2,
          color = "#1A202C",
          inherit.aes = FALSE
        )
      )
    } else {
      list()
    }
  } +
  facet_wrap(~region_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = fill_cols, breaks = group_levels) +
  scale_x_discrete(limits = group_levels, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = "qPCR — Cortex (CTX), Striatum (STR)",
    y = "Absolute expression"
  ) +
  theme_qpcr() +
  theme(
    strip.background = element_rect(fill = "#F7FAFC", color = NA),
    strip.text = element_text(face = "bold", size = 12, margin = margin(4, 4, 4, 4)),
    panel.spacing = unit(1.2, "lines")
  )

out_ctx <- file.path(out_dir, "qpcr_CTX_barplot.png")
out_str <- file.path(out_dir, "qpcr_STR_barplot.png")
out_comb <- file.path(out_dir, "qpcr_CTX_STR_combined.png")

ggsave(out_ctx, p_ctx, width = 8.2, height = 5.6, dpi = 300, bg = "white")
ggsave(out_str, p_str, width = 8.2, height = 5.6, dpi = 300, bg = "white")
ggsave(out_comb, p_combined, width = 12.5, height = 5.8, dpi = 300, bg = "white")

message("Wrote: ", out_ctx)
message("Wrote: ", out_str)
message("Wrote: ", out_comb)
message("Stats: ", stats_path)
message("Summary: ", summary_path)
message("Significant pairs:")
print(as.data.frame(sig_tbl %>% select(region, genotype, group1, group2, p_adj, significance, note)))
message("Underpowered blocks (note nonempty):")
print(as.data.frame(
  stats_tbl %>%
    filter(nzchar(note)) %>%
    distinct(region, genotype, note)
))
