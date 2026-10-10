#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse); library(metafor)
})

args <- commandArgs(trailingOnly = TRUE)
level <- if (length(args) > 0) args[1] else "species"

eff_file <- sprintf('out/effects_long_%s.csv', level)
if (!file.exists(eff_file)) {
  message(sprintf('Input missing for level %s', level))
  quit(save = "no", status = 0)
}
eff <- read_csv(eff_file, show_col_types = FALSE)

pool <- function(df) {
  if (nrow(df) < 2) return(tibble())
  m <- tryCatch(rma(yi = df$yi, vi = df$vi, method = 'REML'), error = function(e) NULL)
  if (is.null(m)) return(tibble())
  est_val <- as.numeric(m$b)
  tibble(est = est_val, se = m$se, pval = m$pval, I2 = m$I2, k = m$k, n_agree = sum(sign(df$yi) == sign(est_val), na.rm = TRUE), n_both = sum(df$both_methods, na.rm = TRUE))
}

scorecard <- eff |> group_by(species) |> filter(n() >= 2) |> group_modify(~ pool(.x)) |> ungroup()
if (nrow(scorecard) > 0) {
  scorecard <- scorecard |> mutate(padj = p.adjust(pval, 'BH'), replicates = !is.na(padj) & padj < 0.05 & n_agree >= ceiling(k * 0.6)) |> arrange(padj)
  write_csv(scorecard, sprintf('out/replication_scorecard_%s.csv', level))
  message(sprintf('Meta-analysis complete for %s. Found %d replicated features.', level, sum(scorecard$replicates)))
} else {
  message(sprintf('No replicated features found for %s.', level))
}
