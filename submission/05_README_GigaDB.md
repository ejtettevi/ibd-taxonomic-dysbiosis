# GigaDB README — ibd-taxonomic-dysbiosis

Supporting data for: **Cross-cohort reproducibility and transportability of
the taxonomic dysbiosis signature of inflammatory bowel disease**.

Linked manuscript: GigaByte (Technical Release).
Authors: Edward Jenner Tettevi (ORCID: 0000-0002-2448-2679),
Mike Y. Osei-Atweneboana — Biomedical and Public Health Research Unit,
CSIR – Water Research Institute (CSIR-WRI), Ghana.

---

## 1. Purpose

This dataset contains every output file that underlies the statistics,
tables, and figures of the linked manuscript. It is a re-analysis of three
publicly available 16S rRNA metagenomic cohorts of inflammatory bowel
disease (IBD): Nielsen et al. 2014 (n=318), Ijaz et al. 2017 (n=50), and
HMP2/IBDMDB (n=130) — 498 samples in total (197 IBD, 301 non-IBD controls).

The analysis pipeline (1) estimates differential abundance per cohort with
ANCOM-BC2 (and MaAsLin2 as a secondary method), (2) pools cohort effects
with random-effects meta-analysis (metafor, REML), (3) defines a replicated
"consensus" dysbiosis signature, and (4) benchmarks transportability with
nested leave-one-cohort-out random-forest classifiers (fixed seed).

## 2. Pipeline provenance

All files in this dataset were generated end-to-end by the script
`run_all.sh` (v1.0, commit 915e6a7 of
https://github.com/ejtettevi/ibd-taxonomic-dysbiosis) from the processed
tables in `curatedMetagenomicData` (MetaPhlAn 3 profiles). The exact
dependency set is reproducible with the included `environment.yml`
(conda env `pd-microbiome`; Python 3.11, R 4.3). No raw sequencing reads
were downloaded or stored.

## 3. File dictionary

| File path | Format | Contents | Manuscript mapping |
|-----------|--------|----------|--------------------|
| `cohorts.tsv` | TSV | Cohort/manifest metadata for the 3 cohorts (498 samples) | Table 1 |
| `environment.yml` | YAML | Conda environment specification | Availability |
| `LICENSE` | text | BSD-3-Clause | Availability |
| `out/species/<cohort>/coefs.csv` | CSV | Species differential-abundance coefficients (ANCOM-BC2 + MaAsLin2) | §2.3, §3.1 |
| `out/genus/<cohort>/coefs.csv` | CSV | Genus differential-abundance coefficients | §2.3, §3.1 |
| `out/effects_long_species.csv` | CSV | Long-form per-study species effect sizes | §3.2 |
| `out/effects_long_genus.csv` | CSV | Long-form per-study genus effect sizes | §3.2 |
| `out/replication_scorecard_species.csv` | CSV | Species meta-analysis + replication summary (27 replicated taxa) | Table 2 |
| `out/replication_scorecard_genus.csv` | CSV | Genus meta-analysis + replication summary (15 replicated taxa) | Table 3 |
| `out/ml_results/roc_auc_summary.csv` | CSV | Six pairwise train→test AUC values (mean 0.746) | Table 4, Fig 9 |
| `out/ml_results/predictions_<trn>__<tst>__genus.csv` | CSV | Predicted scores for all held-out samples per fold | Fig 9 |
| `figure_1_composition.png` | PNG | Per-cohort mean relative abundance by status | Figure 1 |
| `figure_2_pcoa.png` | PNG | PCA projection of consensus-taxon profiles | Figure 2 |
| `figure_3_heatmap.png` | PNG | Spearman rank correlation matrix | Figure 3 |
| `figure_4_diversity.png` | PNG | Within-cohort Shannon diversity | Figure 4 |
| `figure_5_importance.png` | PNG | Genus-level random-forest feature importance | Figure 5 |
| `figure_6_network.png` | PNG | Case/control co-occurrence networks | Figure 6 |
| `figure_7_flowchart.png` | PNG | Analysis workflow and data provenance | Figure 7 |
| `figure_8_forest_plot.png` | PNG | Meta-analytic forest plot (species) | Figure 8 |
| `figure_9_roc_curves.png` | PNG | Leave-one-cohort-out ROC curves | Figure 9 |

## 4. Reproduction instructions

```bash
# 1. Create the conda environment
conda env create -f environment.yml
conda activate pd-microbiome

# 2. Regenerate ALL outputs from curatedMetagenomicData
bash run_all.sh fetch      # download processed tables (network required)
bash run_all.sh all        # harmonize -> analyze -> ml -> figures -> manuscript
```

`run_all.sh all` regenerates every CSV in this dataset, the figures, and the
manuscript draft markdown. All random-forest models use a fixed seed, so a
clean run reproduces identical outputs.

## 5. Compute note

Estimating the complete pipeline for this 498-sample, three-cohort analysis
takes under 20 minutes on a single 8-core laptop with < 8 GB RAM.

## 6. License terms

The analysis code and dataset outputs are released under the BSD-3-Clause
license (see `LICENSE`). Source data remain subject to the terms of
`curatedMetagenomicData` and the original cohort data providers.

## 7. Known issues / limitations

- The three cohorts use different targeted 16S rRNA gene regions and library
  preparation protocols, which contributes to between-study heterogeneity;
  this is quantified (I²) in the scorecards.
- Shannon diversity is computed on the eight consensus taxa only, so the
  resulting values are not comparable in magnitude to whole-community
  diversity.
- Co-occurrence networks in Figure 6 are descriptive (no multiple-testing
  correction), as stated in the manuscript.

---

*Reviewed for submission. Please verify file hashes and byte sizes at export
time if the curation team requires them.*