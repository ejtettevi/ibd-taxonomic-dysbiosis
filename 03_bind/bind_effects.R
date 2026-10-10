#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
level <- if (length(args) > 0) args[1] else "species"

files <- list.files(sprintf('out/%s', level), pattern = 'coefs\\.csv$', recursive = TRUE, full.names = TRUE)
if (length(files) == 0) {
  message(sprintf('No coefs.csv files found for level %s', level))
  quit(save = "no", status = 0)
}

all_eff <- map_df(files, function(f) {
  cohort <- basename(dirname(f))
  read_csv(f, show_col_types = FALSE) |> mutate(cohort = cohort)
})

long_eff <- all_eff |> 
  pivot_longer(cols = c(yi, vi, q_maaslin, lfc_anc, q_anc), names_to = "metric", values_to = "val") |> 
  pivot_wider(names_from = metric, values_from = val)

write_csv(long_eff, sprintf('out/effects_long_%s.csv', level))
message(sprintf('Effects bound successfully for level %s', level))
