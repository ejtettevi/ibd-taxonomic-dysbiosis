#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse); library(phyloseq); library(Maaslin2); library(ANCOMBC)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Usage: Rscript 02_da/da_cohort.R <cohort> <level>')
cohort <- args[1]
level <- args[2]

if (level == "species") {
  abund_file <- sprintf('data/taxon/%s/combined_counts.csv', cohort)
} else {
  abund_file <- sprintf('data/taxon/%s/genus_counts.csv', cohort)
}
meta_file  <- sprintf('data/metadata/%s.csv', cohort)

if (!file.exists(abund_file) || !file.exists(meta_file)) {
  message(sprintf('Files missing for %s at level %s. Skipping.', cohort, level))
  quit(save = "no", status = 0)
}

outdir <- sprintf('out/%s/%s', level, cohort)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

abund_df <- read_csv(abund_file, show_col_types = FALSE)
if (nrow(abund_df) == 0) {
  message(sprintf('Empty table for %s. Skipping.', cohort))
  quit(save = "no", status = 0)
}

id_col_name <- colnames(abund_df)[1]
abund <- abund_df |> 
  mutate(across(1, as.character)) |> 
  column_to_rownames(id_col_name) |> 
  as.matrix()

meta <- read_csv(meta_file, show_col_types = FALSE) |> 
  column_to_rownames('sample_name')

meta$status <- factor(meta$status, levels = c('control', 'PD'))
if (length(unique(meta$status)) < 2) {
  message(sprintf('Cohort %s only has one status group. Skipping.', cohort))
  quit(save = "no", status = 0)
}

common <- intersect(colnames(abund), rownames(meta))
if (length(common) < 3) {
  message(sprintf('Too few common samples (%d) for %s. Skipping.', length(common), cohort))
  quit(save = "no", status = 0)
}
abund <- abund[, common, drop = FALSE]
meta  <- meta[common, , drop = FALSE]
abund <- abund[!rownames(abund) %in% c('unclassified', 'UNFRAC'), , drop = FALSE]

covariates <- c('status')
if ('age' %in% colnames(meta)) covariates <- c(covariates, 'age')
if ('sex' %in% colnames(meta)) covariates <- c(covariates, 'sex')

meta_clean <- meta |> drop_na(all_of(covariates))
abund_clean <- abund[, rownames(meta_clean), drop = FALSE]

mas_out <- file.path(outdir, 'maaslin2')
try({
  Maaslin2(input_data = as.data.frame(t(abund_clean)), input_metadata = meta_clean, output = mas_out, 
           fixed_effects = covariates, reference = 'status,control', normalization = 'TSS', transform = 'LOG', 
           min_prevalence = 0.05, plot_heatmap = FALSE, plot_scatter = FALSE)
})

mas_res_file <- file.path(mas_out, 'all_results.tsv')
mas <- if (file.exists(mas_res_file)) {
  read_tsv(mas_res_file, show_col_types = FALSE) |> filter(metadata == 'status', value == 'PD') |> 
  transmute(species = feature, yi = coef, vi = stderr^2, q_maaslin = qval)
} else { tibble() }

ps <- phyloseq(otu_table(abund_clean, taxa_are_rows = TRUE), sample_data(meta_clean))
anc <- tryCatch({
  ancombc2(data = ps, fix_formula = paste(covariates, collapse = ' + '), group = 'status', p_adj_method = 'BH', prv_cut = 0.05)
}, error = function(e) NULL)

anc_res <- if (!is.null(anc)) {
  anc_df <- anc$res
  lfc_col <- grep('lfc_statusPD', colnames(anc_df), value = TRUE)
  q_col <- grep('q_statusPD', colnames(anc_df), value = TRUE)
  anc_df |> transmute(species = taxon, lfc_anc = .data[[lfc_col]], q_anc = .data[[q_col]])
} else { tibble() }

if (nrow(mas) > 0) {
  coefs <- mas |> left_join(anc_res, by = 'species') |>
    mutate(both_methods = !is.na(q_anc) & ((q_maaslin < 0.05 & q_anc <= 0.1) | (q_maaslin <= 0.1 & q_anc < 0.05)) & sign(yi) == sign(lfc_anc))
  write_csv(coefs, file.path(outdir, 'coefs.csv'))
}
