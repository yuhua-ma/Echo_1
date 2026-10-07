#!/usr/bin/env Rscript
# qPCR bar plots: Cortex (CTX) + Striatum (STR), Day 7 + Day 21 (HT)
# Bars = mean; error bars = SEM; points = individual replicates (jitter).
# Day 21 is HT-only (no WT Day 21 invented).
# Stats: within each region × day × genotype, one-way ANOVA + Tukey HSD
# across treatments (Veh, New SCR, 2.6). Significant pairs only annotated.
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
day_levels <- c("Day 7", "Day 21")

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
    day_num = as.integer(day),
    day = factor(paste("Day", day_num), levels = day_levels),
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
# Day 21 must be HT-only
stopifnot(all(dat$genotype[dat$day_num == 21] == "HT"))

summary_tbl <- dat %>%
  group_by(day, day_num, region, region_label, group, genotype, treatment) %>%
  summarise(
    n = n(),
    mean = mean(expression),
    sd = ifelse(n() > 1, sd(expression), NA_real_),
    sem = ifelse(n() > 1, sd(expression) / sqrt(n()), NA_real_),
    .groups = "drop"
  )

summary_path <- file.path(script_dir, "qpcr_ctx_str_mean_sem.csv")
write_csv(
  summary_tbl %>% select(day_num, region, group, genotype, treatment, n, mean, sd, sem),
  summary_path
)

# ---------------------------------------------------------------------------
# Stats: region × day × genotype → ANOVA + Tukey on treatment
# ---------------------------------------------------------------------------
stat_rows <- list()
for (reg in levels(dat$region)) {
  for (dnum in sort(unique(dat$day_num))) {
    for (gen in levels(dat$genotype)) {
      sub <- dat %>% filter(region == reg, day_num == dnum, genotype == gen)
      if (nrow(sub) == 0) next
      n_treats <- n_distinct(sub$treatment)
      if (n_treats < 2) next

      ns <- sub %>% count(treatment, name = "n")
      min_n <- min(ns$n)
      caution <- if (min_n < 3) {
        sprintf(
          "Low power / caution: min n=%d in region×day×genotype (day %d, %s %s)",
          min_n, dnum, reg, as.character(gen)
        )
      } else {
        ""
      }

      # Need ≥2 levels with data and residual df possible
      if (n_treats < 3 && nrow(sub) <= n_treats) {
        # still try if possible
      }

      fit <- tryCatch(
        aov(expression ~ treatment, data = sub),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      aov_tab <- summary(fit)[[1]]
      aov_p <- aov_tab[["Pr(>F)"]][1]
      aov_F <- aov_tab[["F value"]][1]
      df1 <- aov_tab[["Df"]][1]
      df2 <- aov_tab[["Df"]][2]
      if (is.na(aov_F) || is.na(df2) || df2 < 1) {
        # cannot run Tukey; record ANOVA failure / insufficient df
        treats_present <- as.character(ns$treatment)
        pairs <- utils::combn(treats_present, 2, simplify = FALSE)
        for (pr in pairs) {
          g1 <- paste0(gen, "-", treat_to_suffix(pr[1]))
          g2 <- paste0(gen, "-", treat_to_suffix(pr[2]))
          pos1 <- match(g1, group_levels)
          pos2 <- match(g2, group_levels)
          if (!is.na(pos1) && !is.na(pos2) && pos1 > pos2) {
            tmp <- g1; g1 <- g2; g2 <- tmp
            tmp <- pr[1]; pr[1] <- pr[2]; pr[2] <- tmp
          }
          stat_rows[[length(stat_rows) + 1]] <- data.frame(
            day = dnum,
            region = reg,
            region_label = as.character(unique(dat$region_label[dat$region == reg])),
            genotype = as.character(gen),
            group1 = g1,
            group2 = g2,
            treatment1 = pr[1],
            treatment2 = pr[2],
            contrast = paste0(pr[2], "-", pr[1]),
            n_group1 = ns$n[ns$treatment == pr[1]],
            n_group2 = ns$n[ns$treatment == pr[2]],
            mean_diff = NA_real_,
            anova_df1 = df1,
            anova_df2 = df2,
            anova_F = aov_F,
            anova_p = aov_p,
            test = "one-way ANOVA (Tukey not run; insufficient residual df)",
            p_adj = NA_real_,
            p_adj_method = NA_character_,
            significance = NA_character_,
            note = paste(caution, "Insufficient residual df for pairwise Tukey."),
            stringsAsFactors = FALSE
          )
        }
        next
      }

      tuk <- TukeyHSD(fit, "treatment")
      tk <- as.data.frame(tuk$treatment)
      tk$contrast <- rownames(tk)
      treats <- treat_levels

      for (i in seq_len(nrow(tk))) {
        contrast <- tk$contrast[i]
        t1 <- NA_character_
        t2 <- NA_character_
        for (a in treats) {
          for (b in treats) {
            if (a != b && identical(contrast, paste0(a, "-", b))) {
              t1 <- as.character(a)
              t2 <- as.character(b)
            }
          }
        }
        if (is.na(t1) || is.na(t2)) next
        # Skip if either treatment absent in this block
        if (!(t1 %in% ns$treatment) || !(t2 %in% ns$treatment)) next

        g1 <- paste0(gen, "-", treat_to_suffix(t1))
        g2 <- paste0(gen, "-", treat_to_suffix(t2))
        pos1 <- match(g1, group_levels)
        pos2 <- match(g2, group_levels)
        if (!is.na(pos1) && !is.na(pos2) && pos1 > pos2) {
          tmp <- g1; g1 <- g2; g2 <- tmp
          tmp <- t1; t1 <- t2; t2 <- tmp
        }
        n1 <- ns$n[ns$treatment == t1]
        n2 <- ns$n[ns$treatment == t2]

        stat_rows[[length(stat_rows) + 1]] <- data.frame(
          day = dnum,
          region = reg,
          region_label = as.character(unique(dat$region_label[dat$region == reg])),
          genotype = as.character(gen),
          group1 = g1,
          group2 = g2,
          treatment1 = t1,
          treatment2 = t2,
          contrast = contrast,
          n_group1 = as.integer(n1),
          n_group2 = as.integer(n2),
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
}

stats_tbl <- bind_rows(stat_rows)
stats_path <- file.path(script_dir, "qpcr_ctx_str_stats.csv")
write_csv(stats_tbl, stats_path)

sig_tbl <- stats_tbl %>%
  filter(!is.na(significance), significance != "ns")

# ---------------------------------------------------------------------------
# Plot theme / colors
# ---------------------------------------------------------------------------
# Distinguish days by fill; keep group identity via x-axis order
day_fills <- c("Day 7" = "#4A5568", "Day 21" = "#2B6CB0")
day_shapes <- c("Day 7" = 21, "Day 21" = 24)

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
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      plot.margin = margin(12, 14, 10, 12)
    )
}

# Bracket annotations for dodged bars: x at group center ± dodge offset by day
# Only annotate when both compared groups exist on that day (same day block).
build_day_annotations <- function(sig_df, summary_df, tip_pad = 0.03, step = 0.05) {
  # sig_df$day is integer day number from stats table
  if (nrow(sig_df) == 0) {
    return(data.frame(
      day = factor(levels = day_levels),
      day_num = integer(),
      region = character(),
      region_label = character(),
      group1 = character(),
      group2 = character(),
      x1 = numeric(),
      x2 = numeric(),
      y_position = numeric(),
      annotations = character(),
      stringsAsFactors = FALSE
    ))
  }

  tips <- summary_df %>%
    mutate(ceiling = mean + ifelse(is.na(sem), 0, sem)) %>%
    select(day, day_num, region, group, ceiling)
  raw_max <- dat %>%
    group_by(day, day_num, region, group) %>%
    summarise(raw_max = max(expression), .groups = "drop")
  tips <- tips %>%
    left_join(raw_max, by = c("day", "day_num", "region", "group")) %>%
    mutate(ceiling = pmax(ceiling, raw_max, na.rm = TRUE))

  # Dodge width 0.75 with 2 day levels → offset ≈ ±0.75/4
  dodge_w <- 0.75
  day_offsets <- c("Day 7" = -dodge_w / 4, "Day 21" = dodge_w / 4)

  out <- list()
  keys <- sig_df %>% distinct(region, day)
  for (k in seq_len(nrow(keys))) {
    reg <- keys$region[k]
    dnum <- keys$day[k]
    dlab <- paste("Day", dnum)
    sreg <- sig_df %>% filter(region == reg, day == dnum) %>%
      mutate(
        span = match(group2, group_levels) - match(group1, group_levels),
        left = match(group1, group_levels)
      ) %>%
      arrange(left, span)

    tip_sub <- tips %>% filter(region == reg, day_num == dnum)
    base_y <- if (nrow(tip_sub)) max(tip_sub$ceiling, na.rm = TRUE) + tip_pad else tip_pad
    off <- unname(day_offsets[dlab])
    if (is.na(off)) off <- 0

    for (i in seq_len(nrow(sreg))) {
      g1 <- sreg$group1[i]
      g2 <- sreg$group2[i]
      pair_ceil <- max(
        tip_sub$ceiling[tip_sub$group %in% c(g1, g2)],
        na.rm = TRUE
      )
      if (!is.finite(pair_ceil)) pair_ceil <- base_y
      y_pos <- max(base_y, pair_ceil + tip_pad) + (i - 1) * step
      out[[length(out) + 1]] <- data.frame(
        day = factor(dlab, levels = day_levels),
        day_num = as.integer(dnum),
        region = reg,
        region_label = sreg$region_label[i],
        group1 = g1,
        group2 = g2,
        x1 = match(g1, group_levels) + off,
        x2 = match(g2, group_levels) + off,
        y_position = y_pos,
        annotations = sreg$significance[i],
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(out)
}

ann_all <- build_day_annotations(sig_tbl, summary_tbl)

bracket_layers_data <- function(ann_df, tip_frac = 0.03) {
  if (nrow(ann_df) == 0) {
    return(list(segs = data.frame(), labs = data.frame()))
  }
  ann_df <- ann_df %>%
    group_by(region, day_num) %>%
    mutate(tip = pmax(y_position * tip_frac, 0.01)) %>%
    ungroup() %>%
    mutate(region_label = factor(region_label, levels = levels(dat$region_label)))

  segs <- bind_rows(
    ann_df %>% transmute(region, region_label, day, day_num, x = x1, xend = x2, y = y_position, yend = y_position),
    ann_df %>% transmute(region, region_label, day, day_num, x = x1, xend = x1, y = y_position - tip, yend = y_position),
    ann_df %>% transmute(region, region_label, day, day_num, x = x2, xend = x2, y = y_position - tip, yend = y_position)
  )
  labs <- ann_df %>%
    transmute(
      region, region_label, day, day_num,
      x = (x1 + x2) / 2,
      y = y_position,
      label = annotations
    )
  list(segs = segs, labs = labs)
}

make_region_plot <- function(region_code, title_text) {
  d <- dat %>% filter(region == region_code)
  s <- summary_tbl %>% filter(region == region_code)
  ann <- ann_all %>% filter(region == region_code)
  br <- bracket_layers_data(ann)

  y_data_max <- max(d$expression, s$mean + ifelse(is.na(s$sem), 0, s$sem), na.rm = TRUE)
  y_ann_max <- if (nrow(ann)) max(ann$y_position) else y_data_max
  y_max <- max(y_data_max * 1.12, y_ann_max * 1.14)

  pd <- position_dodge(width = 0.75)

  p <- ggplot() +
    geom_col(
      data = s,
      aes(x = group, y = mean, fill = day),
      position = pd,
      width = 0.68,
      color = "#1A202C",
      linewidth = 0.25,
      alpha = 0.92
    ) +
    geom_errorbar(
      data = s,
      aes(x = group, ymin = mean - sem, ymax = mean + sem, group = day),
      position = pd,
      width = 0.2,
      linewidth = 0.5,
      color = "#1A202C",
      na.rm = TRUE
    ) +
    geom_point(
      data = d,
      aes(x = group, y = expression, fill = day, shape = day),
      position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
      size = 2.0,
      color = "#1A202C",
      stroke = 0.45,
      alpha = 0.95
    ) +
    scale_fill_manual(values = day_fills, breaks = day_levels, drop = FALSE) +
    scale_shape_manual(values = day_shapes, breaks = day_levels, drop = FALSE) +
    scale_x_discrete(limits = group_levels, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, y_max)) +
    labs(title = title_text, y = "Absolute expression") +
    theme_qpcr()

  if (nrow(br$segs) > 0) {
    p <- p +
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
        size = 4.2,
        color = "#1A202C",
        inherit.aes = FALSE
      )
  }
  p
}

p_ctx <- make_region_plot("CTX", "qPCR — Cortex (CTX), Day 7 and Day 21")
p_str <- make_region_plot("STR", "qPCR — Striatum (STR), Day 7 and Day 21")

# Combined facet
br_all <- bracket_layers_data(ann_all)
pd <- position_dodge(width = 0.75)
p_combined <- ggplot() +
  geom_col(
    data = summary_tbl,
    aes(x = group, y = mean, fill = day),
    position = pd,
    width = 0.68,
    color = "#1A202C",
    linewidth = 0.25,
    alpha = 0.92
  ) +
  geom_errorbar(
    data = summary_tbl,
    aes(x = group, ymin = mean - sem, ymax = mean + sem, group = day),
    position = pd,
    width = 0.2,
    linewidth = 0.5,
    color = "#1A202C",
    na.rm = TRUE
  ) +
  geom_point(
    data = dat,
    aes(x = group, y = expression, fill = day, shape = day),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    size = 1.8,
    color = "#1A202C",
    stroke = 0.4,
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
          size = 4.0,
          color = "#1A202C",
          inherit.aes = FALSE
        )
      )
    } else list()
  } +
  facet_wrap(~region_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = day_fills, breaks = day_levels, drop = FALSE) +
  scale_shape_manual(values = day_shapes, breaks = day_levels, drop = FALSE) +
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

# Day-21-only HT panels (helpful)
make_day21_ht_plot <- function(region_code, title_text) {
  d <- dat %>% filter(region == region_code, day_num == 21, genotype == "HT")
  s <- summary_tbl %>% filter(region == region_code, day_num == 21, genotype == "HT")
  ht_groups <- c("HT-Veh", "HT-New SCR", "HT-2,6")
  ann <- ann_all %>% filter(region == region_code, day_num == 21)
  # Remap x to 1..3 for HT-only axis
  if (nrow(ann) > 0) {
    ann2 <- ann %>%
      mutate(
        x1 = match(group1, ht_groups),
        x2 = match(group2, ht_groups)
      )
  } else {
    ann2 <- ann
  }
  br <- bracket_layers_data(ann2)

  d <- d %>% mutate(group = factor(as.character(group), levels = ht_groups))
  s <- s %>% mutate(group = factor(as.character(group), levels = ht_groups))

  y_data_max <- max(d$expression, s$mean + ifelse(is.na(s$sem), 0, s$sem), na.rm = TRUE)
  y_ann_max <- if (nrow(ann2)) max(ann2$y_position) else y_data_max
  y_max <- max(y_data_max * 1.15, y_ann_max * 1.18)

  p <- ggplot() +
    geom_col(
      data = s,
      aes(x = group, y = mean),
      fill = day_fills[["Day 21"]],
      width = 0.65,
      color = "#1A202C",
      linewidth = 0.25,
      alpha = 0.92
    ) +
    geom_errorbar(
      data = s,
      aes(x = group, ymin = mean - sem, ymax = mean + sem),
      width = 0.18,
      linewidth = 0.5,
      color = "#1A202C",
      na.rm = TRUE
    ) +
    geom_jitter(
      data = d,
      aes(x = group, y = expression),
      width = 0.1,
      height = 0,
      size = 2.2,
      shape = 24,
      fill = "white",
      color = "#1A202C",
      stroke = 0.5
    ) +
    scale_x_discrete(limits = ht_groups, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, y_max)) +
    labs(title = title_text, y = "Absolute expression") +
    theme_qpcr() +
    theme(legend.position = "none")

  if (nrow(br$segs) > 0) {
    p <- p +
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
        size = 4.2,
        color = "#1A202C",
        inherit.aes = FALSE
      )
  }
  p
}

p_ctx_d21 <- make_day21_ht_plot("CTX", "qPCR — Cortex (CTX), Day 21 HT")
p_str_d21 <- make_day21_ht_plot("STR", "qPCR — Striatum (STR), Day 21 HT")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_ctx <- file.path(out_dir, "qpcr_CTX_Day7_Day21_barplot.png")
out_str <- file.path(out_dir, "qpcr_STR_Day7_Day21_barplot.png")
out_comb <- file.path(out_dir, "qpcr_CTX_STR_Day7_Day21_combined.png")
out_ctx_d21 <- file.path(out_dir, "qpcr_CTX_Day21_HT_barplot.png")
out_str_d21 <- file.path(out_dir, "qpcr_STR_Day21_HT_barplot.png")

# Also refresh legacy filenames used previously (Day7+Day21 combined view)
out_ctx_legacy <- file.path(out_dir, "qpcr_CTX_barplot.png")
out_str_legacy <- file.path(out_dir, "qpcr_STR_barplot.png")
out_comb_legacy <- file.path(out_dir, "qpcr_CTX_STR_combined.png")

ggsave(out_ctx, p_ctx, width = 9.0, height = 5.8, dpi = 300, bg = "white")
ggsave(out_str, p_str, width = 9.0, height = 5.8, dpi = 300, bg = "white")
ggsave(out_comb, p_combined, width = 13.0, height = 6.0, dpi = 300, bg = "white")
ggsave(out_ctx_d21, p_ctx_d21, width = 6.5, height = 5.2, dpi = 300, bg = "white")
ggsave(out_str_d21, p_str_d21, width = 6.5, height = 5.2, dpi = 300, bg = "white")
ggsave(out_ctx_legacy, p_ctx, width = 9.0, height = 5.8, dpi = 300, bg = "white")
ggsave(out_str_legacy, p_str, width = 9.0, height = 5.8, dpi = 300, bg = "white")
ggsave(out_comb_legacy, p_combined, width = 13.0, height = 6.0, dpi = 300, bg = "white")

message("Wrote: ", out_ctx)
message("Wrote: ", out_str)
message("Wrote: ", out_comb)
message("Wrote: ", out_ctx_d21)
message("Wrote: ", out_str_d21)
message("Stats: ", stats_path)
message("Summary: ", summary_path)

message("Significant pairs:")
print(as.data.frame(sig_tbl %>%
  select(day, region, genotype, group1, group2, p_adj, significance, note)))

message("Day 21 significant:")
print(as.data.frame(sig_tbl %>%
  filter(day == 21) %>%
  select(region, genotype, group1, group2, p_adj, significance, note)))

message("Underpowered blocks:")
print(as.data.frame(
  stats_tbl %>%
    filter(nzchar(note)) %>%
    distinct(day, region, genotype, note)
))
