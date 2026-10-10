#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

taxon_in <- 'data/taxon'
meta_in  <- 'data/metadata'
taxon_out <- 'taxon'
meta_out  <- 'metadata'

dir.create(taxon_out, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_out, recursive = TRUE, showWarnings = FALSE)

taxon_files <- list.files(taxon_in, pattern = '.csv|.tsv', full.names = TRUE)
if (length(taxon_files) == 0) stop('No abundance tables found in data/taxon/.')

for (f in taxon_files) {
  cohort <- tools::file_path_sans_ext(basename(f))
  message(sprintf('[%s] Harmonizing taxon table...', cohort))
  
  # Use read_csv but explicitly handle the first column
  df <- read_csv(f, show_col_types = FALSE)
  if (nrow(df) == 0) {
    message(sprintf('Skipping empty file %s', f))
    next
  }
  
  col_names <- colnames(df)
  # Try to find a name column
  target_col <- col_names[grep('species|taxon|genus|feature|clade|motu', col_names, ignore.case = TRUE)][1]
  
  if (is.na(target_col)) {
    target_col <- col_names[1]
    message(sprintf('Using first column as ID for %s: %s', f, target_col))
  }
  
  # 1. Rename and select
  df <- df |> rename(species = !!sym(target_col)) |> select(species, everything())
  
  # 2. STRICKT CLEANING: Remove any row where species is NA, empty, or looks like a number
  # This is the critical fix for the duplicate row.names error
  df <- df |> 
    filter(!is.na(species)) |> 
    filter(species != "") |>
    filter(!grepl("^[0-9\\.]+$", species)) |> # Remove rows where the "species" is just a number
    mutate(species = as.character(species)) |> 
    mutate(species = str_trim(species)) |> 
    mutate(species = str_replace_all(species, ' ', '_'))
  
  # 3. Remove duplicates
  df <- df |> distinct(species, .keep_all = TRUE)
  
  if (nrow(df) == 0) {
    message(sprintf('No valid species rows left in %s after cleaning. Skipping.', f))
    next
  }

  dir.create(file.path(taxon_out, cohort), recursive = TRUE, showWarnings = FALSE)
  write_csv(df, file.path(taxon_out, cohort, 'combined_counts.csv'))
  
  # Genus level aggregation
  df_gen <- df |> 
    mutate(genus = str_extract(species, '.*\\.g__[^\\.]*')) |> 
    filter(!is.na(genus)) |>
    select(genus, everything(), -species) |>
    group_by(genus) |>
    summarise(across(where(is.numeric), sum), .groups = 'drop')
  
  write_csv(df_gen, file.path(taxon_out, cohort, 'genus_counts.csv'))
}

# METADATA
meta_files <- list.files(meta_in, pattern = '.csv|.tsv', full.names = TRUE)
for (f in meta_files) {
  raw_name <- tools::file_path_sans_ext(basename(f))
  cohort <- raw_name
  if (grepl('metadata.smg', raw_name)) cohort <- 'mOTUs'
  if (grepl('metadata.16S', raw_name)) cohort <- 'wallen2022'
  
  message(sprintf('Harmonizing metadata: %s -> %s', raw_name, cohort))
  df <- read_csv(f, show_col_types = FALSE)
  col_names <- colnames(df)
  target_id <- col_names[grep('sample|id|name', col_names, ignore.case = TRUE)][1]
  if (!is.na(target_id)) { df <- df |> rename(sample_name = !!sym(target_id)) |> select(sample_name, everything()) }
  
  target_status <- NA
  if ("PD" %in% col_names) { target_status <- "PD" } 
  else { target_status <- col_names[grep('status|diag|group|cond|disease|case', col_names, ignore.case = TRUE)][1] }
  
  if (!is.na(target_status)) {
    df <- df |> mutate(status_raw = as.character(!!sym(target_status)))
    df <- df |> mutate(status = case_when(
      grepl('PD|Case|Disease|1|Patient', status_raw, ignore.case = TRUE) ~ 'PD',
      grepl('Control|HC|Healthy|0|Normal', status_raw, ignore.case = TRUE) ~ 'control',
      TRUE ~ 'unknown'
    ))
    df <- df |> select(-status_raw)
  } else { df$status <- 'unknown' }
  write_csv(df, file.path(meta_out, sprintf('%s.csv', cohort)))
}
message('Harmonization complete.')
