#!/usr/bin/env Rscript
#
# Fetch REAL processed metagenomic tables for three IBD cohorts from
# curatedMetagenomicData (Bioconductor 3.18, DOI 10.18129/B9.bioc.curatedMetagenomicData).
#
# These are shotgun metagenomes processed upstream with MetaPhlAn3, distributed as
# TreeSummarizedExperiment objects via ExperimentHub. Only the processed abundance
# tables and curated metadata are downloaded; no raw reads are fetched or stored.
#
# Output layout matches what the existing pipeline stages consume:
#   data/taxon/<cohort>.csv                 flat species table (consumed by 01_harmonize)
#   data/taxon/<cohort>/combined_counts.csv species table (consumed at "species" level)
#   data/taxon/<cohort>/genus_counts.csv    genus rollup (consumed at "genus" level)
#   data/metadata/<cohort>.csv              sample_name / status / covariates
#
# NOTE on naming: curatedMetagenomicData rownames are MetaPhlAn-style and
# pipe-delimited, e.g. "k__Bacteria|p__Firmicutes|...|g__Blautia|s__Blautia_vulpis".
# The downstream harmonizer and any Greengenes-style tooling expects dot-delimited
# lineages, so lineages are converted to dots here.

suppressPackageStartupMessages({
  library(tidyverse)
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
})

# --- cohort definitions -------------------------------------------------------
# baseline_only: keep visit_number == 1 so longitudinal studies contribute a
# single cross-sectional timepoint instead of repeated measures per subject.
# NOTE: PMIDs below are the curated study's own citations, read from the
# colData(PMID) field of each TreeSummarizedExperiment at fetch time and
# cross-checked against pmid_expected. Do not hand-edit these into provenance:
# earlier hardcoded guesses pointed at completely unrelated papers.
studies <- list(
  list(
    cohort        = "Nielsen2014",
    study         = "NielsenHB_2014",
    pmid_expected = "24997787",    # Nielsen HB et al., Nat Biotechnol 2014
    description   = "Shotgun metagenomes, IBD vs healthy controls, cross-sectional",
    baseline_only = FALSE,
    age_field     = NULL,          # only a coarse age_category is available
    sex_field     = "gender"       # complete
  ),
  list(
    cohort        = "Ijaz2017",
    study         = "IjazUZ_2017",
    pmid_expected = "28222161",    # Ijaz UZ et al., PLoS One 2017
    description   = paste("Crohn's disease dysbiosis; controls are unaffected",
                          "genetically-linked kindred, cross-sectional"),
    baseline_only = FALSE,
    age_field     = "age",         # numeric years
    sex_field     = NULL           # gender missing for 83/94 samples; excluded
  ),
  list(
    cohort        = "HMP2",
    study         = "HMP_2019_ibdmdb",
    pmid_expected = c("29311644", "31142855"), # Schirmer 2018; Lloyd-Price 2019
    description   = "HMP2 IBD multi-omics, baseline visit only",
    baseline_only = TRUE,
    age_field     = "age",         # numeric years
    sex_field     = "gender"       # complete
  )
)

taxon_dir <- "data/taxon"
meta_dir  <- "data/metadata"
dir.create(taxon_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir,  recursive = TRUE, showWarnings = FALSE)

# --- helpers ------------------------------------------------------------------
# "k__A|p__B|g__C|s__D" -> "k__A.p__B.g__C.s__D"  (drop empty / taxon-rank entries)
# Vectorised: strsplit() returns a list, so each lineage is collapsed separately.
dots <- function(x) {
  vapply(strsplit(x, "|", fixed = TRUE), function(p) {
    p <- p[nzchar(p) & !grepl("^t__", p)]
    paste(p, collapse = ".")
  }, character(1), USE.NAMES = FALSE)
}

fetch_one <- function(spec) {
  message(sprintf("=== %s (%s) ===", spec$cohort, spec$study))

  se <- curatedMetagenomicData(
    paste0(spec$study, ".relative_abundance"),
    dryrun = FALSE, counts = FALSE, rownames = "long"
  )[[1]]

  cd <- as.data.frame(colData(se))
  message(sprintf("  raw: %d species x %d samples", nrow(se), ncol(se)))

  keep <- rep(TRUE, nrow(cd))
  if ("body_site" %in% colnames(cd)) {
    keep <- keep & grepl("stool|fece|feces|colonic", cd$body_site, ignore.case = TRUE)
    message(sprintf("  stool only: %d samples", sum(keep)))
  }
  if (spec$baseline_only && "visit_number" %in% colnames(cd)) {
    baseline <- !is.na(cd$visit_number) & cd$visit_number == 1
    message(sprintf("  baseline visit only: %d samples", sum(baseline)))
    keep <- keep & baseline
  }

  # Only true case/control contrasts are usable for the disease-vs-control design.
  disease <- as.character(cd$disease)
  keep <- keep & disease %in% c("healthy", "IBD")
  cd <- cd[keep, , drop = FALSE]
  disease <- disease[keep]

  # One sample per subject: keep the earliest visit among the survivors.
  if ("subject_id" %in% colnames(cd)) {
    ord <- order(cd$subject_id)
    cd <- cd[ord, , drop = FALSE]
    disease <- disease[ord]
    first <- !duplicated(cd$subject_id)
    cd <- cd[first, , drop = FALSE]
    disease <- disease[first]
    message(sprintf("  unique subjects: %d", nrow(cd)))
  }

  abund <- assay(se, "relative_abundance")[, rownames(cd), drop = FALSE]

  # --- species level ---
  # assay() is already features (rows) x samples (columns); keep that orientation.
  sp <- as.data.frame(abund, check.names = FALSE)
  sp_names <- rownames(se)

  # Row labels come from the FULL lineage, not the bare species annotation.
  # rowData(se)$species holds only e.g. "Roseburia_hominis", which is ambiguous
  # (the same epithet can occur in more than one genus) and collides across cohorts
  # once make.unique() starts appending disambiguating suffixes. The pipe-delimited
  # lineage is unique and consistent, so it is both the label and the grouping key.
  species_lin <- dots(sp_names)
  rownames(sp) <- make.unique(species_lin)

  # --- genus rollup, summed over species within each genus lineage ---
  # The genus key must come from the FULL lineage with the species suffix removed,
  # not from species_col (which holds only e.g. "Blautia_vulpis" and so carries no
  # g__ rank). Genus identity is the lineage truncated after its g__ rank, which
  # keeps same-named genera in different families as separate rows.
  genus_of <- function(lin) {
    has <- grepl("\\.g__", lin)
    out <- rep(NA_character_, length(lin))
    out[has] <- sub("^(([^.]*\\.)*g__[^.]*).*$", "\\1", lin[has])
    out
  }
  genus_key <- genus_of(dots(sub("\\|s__.*$", "", sp_names)))

  keep_sp <- rowSums(sp) > 0
  sp <- sp[keep_sp, , drop = FALSE]
  genus_key <- genus_key[keep_sp]

  keep_g <- !is.na(genus_key)
  genus_tab <- rowsum(sp[keep_g, , drop = FALSE],
                      group = genus_key[keep_g], reorder = FALSE)

  message(sprintf("  cases: %d | controls: %d", sum(disease == "IBD"), sum(disease == "healthy")))
  message(sprintf("  species retained: %d | genera retained: %d", nrow(sp), nrow(genus_tab)))

  # --- write taxon tables ---
  cdir <- file.path(taxon_dir, spec$cohort)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)

  species_out <- sp
  colnames(species_out) <- rownames(cd)

  # --- write taxon tables ---
  cdir <- file.path(taxon_dir, spec$cohort)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)

  # readr::write_csv() does not emit row names, so an explicit identifier column is
  # required: 01_harmonize reads the flat file's first column as the taxon ID, and
  # the analysis stages read the per-cohort files with index_col = 0.
  # Built with cbind on plain vectors: data.frame() does not support `:=`.
  write_table <- function(mat, id_col, path) {
    id <- data.frame(rownames(mat), stringsAsFactors = FALSE)
    names(id) <- id_col
    df <- cbind(id, as.data.frame(mat, check.names = FALSE))
    readr::write_csv(df, path)
  }
  write_table(species_out, "species", file.path(taxon_dir, paste0(spec$cohort, ".csv")))
  write_table(species_out, "species", file.path(cdir, "combined_counts.csv"))
  write_table(genus_tab, "genus", file.path(cdir, "genus_counts.csv"))

  # --- write metadata ---
  meta <- tibble(
    sample_name = rownames(cd),
    status      = ifelse(disease == "IBD", "PD", "control")
  )
  if (!is.null(spec$age_field) && spec$age_field %in% colnames(cd)) {
    meta$age <- suppressWarnings(as.numeric(cd[[spec$age_field]]))
  }
  if (!is.null(spec$sex_field) && spec$sex_field %in% colnames(cd)) {
    meta$sex <- as.character(cd[[spec$sex_field]])
  }
  if ("disease_subtype" %in% colnames(cd)) {
    meta$disease_subtype <- as.character(cd$disease_subtype)
  }
  readr::write_csv(meta, file.path(meta_dir, paste0(spec$cohort, ".csv")))

  # Read the citation from the object's own curated metadata rather than trusting
  # a hardcoded string, and warn loudly if it disagrees with the expected PMID.
  pmids <- character(0)
  if ("PMID" %in% colnames(cd)) {
    pmids <- unique(as.character(cd$PMID))
    pmids <- pmids[!is.na(pmids) & nzchar(pmids)]
  }
  expected <- spec$pmid_expected
  if (length(pmids) && !any(expected %in% pmids)) {
    warning(sprintf("%s: curated PMID(s) %s do not include expected %s",
                    spec$cohort, paste(pmids, collapse = ","),
                    paste(expected, collapse = ",")))
  }
  message(sprintf("  PMID(s): %s", paste(pmids, collapse = ", ")))

  invisible(list(cohort = spec$cohort, n_samples = nrow(cd),
                 n_species = nrow(species_out), n_genera = nrow(genus_tab),
                 pmids = pmids, study = spec$study,
                 description = spec$description))
}

results <- lapply(studies, fetch_one)

# --- record provenance in cohorts.tsv ----------------------------------------
# This is a REPLACE, not an append: earlier versions used bind_rows(rows, tsv),
# so every re-run duplicated all three rows and stale placeholder cohorts
# (PRJNA123456 etc.) accumulated alongside the real ones.
managed <- vapply(studies, `[[`, character(1), "cohort")
source_tag <- "curatedMetagenomicData (MetaPhlAn3)"

rows <- tibble(
  cohort    = vapply(studies, `[[`, character(1), "cohort"),
  source_id = vapply(studies, `[[`, character(1), "study"),
  source    = source_tag,
  note      = vapply(results, function(r) {
    paste0("PMID ", paste(r$pmids, collapse = ", "), "; ", r$description)
  }, character(1))
)

tsv <- if (file.exists("cohorts.tsv")) {
  readr::read_tsv("cohorts.tsv", show_col_types = FALSE)
} else {
  tibble(cohort = character(0), source_id = character(0),
         source = character(0), note = character(0))
}

# Drop any previous row this script owns, then append the freshly written ones.
# The old file used a 'bioproject' column for curated study names; those are not
# BioProject accessions, so the column is renamed to source_id and stale rows
# from that schema are dropped too.
if ("bioproject" %in% colnames(tsv)) {
  tsv <- tsv %>%
    transmute(cohort = as.character(cohort),
              source_id = as.character(bioproject),
              source = as.character(source),
              note = as.character(note))
}
keep <- !(tsv$cohort %in% managed) & !(tsv$source == source_tag)
combined <- bind_rows(rows, tsv[keep, , drop = FALSE])
combined <- combined[order(factor(combined$cohort, levels = managed),
                           combined$cohort), , drop = FALSE]
readr::write_tsv(combined, "cohorts.tsv")
message(sprintf("\ncohorts.tsv: %d rows (%d managed by this script)",
                nrow(combined), length(managed)))

message("\n=== Summary ===")
for (r in results) {
  message(sprintf("  %-12s samples=%4d species=%4d genera=%4d",
                  r$cohort, r$n_samples, r$n_species, r$n_genera))
}
message("Done.")