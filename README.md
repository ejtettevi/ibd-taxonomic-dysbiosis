<div align="center">

# IBD Taxonomic Dysbiosis — Cross-Cohort Reproducibility

[![Language: R](https://img.shields.io/badge/R-4.3+-276DC3.svg)](https://www.r-project.org/)
[![Language: Python](https://img.shields.io/badge/Python-3.11+-3776AB.svg)](https://www.python.org/)
[![Status: Open Source](https://img.shields.io/badge/status-open%20source-brightgreen.svg)](LICENSE)

**Do taxonomic dysbiosis signatures of inflammatory bowel disease (IBD) reproduce across independent cohorts?**

A transparent, reproducible pipeline that asks whether case–control taxonomic signatures discovered in one shotgun-metagenomic IBD cohort are consistent (directionally reproducible) and transportable (predictive) in other cohorts — at both **species** and **genus** resolution.

</div>

---

## Highlights

- **3 independent IBD cohorts** — Nielsen 2014, Ijaz 2017, HMP2 — all from [`curatedMetagenomicData`](https://bioconductor.org/packages/release/data/experiment/html/curatedMetagenomicData.html) (MetaPhlAn3 profiles; no raw reads fetched).
- **Two differential-abundance (DA) analyses per cohort** — [MaAsLin2](https://github.com/biobakery/Maaslin2) and [ANCOM-BC2](https://github.com/FrederickHuangLin/ANCOMBC) run independently, with a concordance rule requiring both methods to agree on direction.
- **Random-effects meta-analysis** across cohorts (`metafor::rma`, REML) with a replication scorecard: **27 species** and **15 genera** replicate after Benjamini–Hochberg FDR control.
- **Cross-cohort ML validation** — Random Forest trained on one cohort, tested on another (leave-one-cohort-out), **mean AUC 0.746** on genus-level taxonomy.
- **Reproducibility built-in** — one `run_all.sh` entry point, curated inputs, deterministic seeds, and a fixed-format replication scorecard.

---

## Repository layout

```
.
├── 00_setup/                 Conda environment bootstrap
├── 01_fetch_processed/       curatedMetagenomicData fetch (MetaPhlAn3 tables + metadata)
├── 01_fetch_upload/          (legacy) raw-read SRA download + upload scaffolding
├── 01_harmonize/             standardize lineages + build species/genus count tables
├── 02_da/                    per-cohort differential abundance (MaAsLin2 + ANCOM-BC2)
├── 03_bind/                  bind per-cohort effects into one long-format table
├── 04_meta/                  random-effects meta-analysis + replication scorecard
├── 05_ml/                    cross-cohort Random Forest transportability (LOCO)
├── data/
│   ├── metadata/             sample_name / status / covariates per cohort
│   └── taxon/                species & genus count tables (with per-cohort dirs)
├── figures/                  intermediate figure renders
├── figures_final/            publication-quality figures 1–9
├── generate_*figures*.py     figure pipelines (final = publication versions)
├── environment.yml           full conda environment (R + Python)
├── cohorts.tsv               cohort provenance manifest (source, PMIDs)
└── run_all.sh                one-command pipeline driver
```

> `06_manuscript/`, `out/`, and `docx/` are build artifacts and excluded from version control (see `.gitignore`).

---

## Pipeline

```mermaid
graph LR
    A[fetch_curated_metagenomic.R<br/>curatedMetagenomicData] --> B[harmonize_tables.R<br/>species / genus counts]
    B --> C[da_cohort.R × 2 levels<br/>MaAsLin2 + ANCOM-BC2]
    C --> D[bind_effects.R<br/>long-format effects]
    D --> E[meta.R (metafor REML)<br/>replication scorecard]
    B --> F[ml_crosscohort.py<br/>Random Forest, LOCO]
    E --> G[generate_final_publication_figures.py]
    F --> G
```

**Stage 2 — DA concordance rule.** A species counts as a cohort-level signal only when MaAsLin2 and ANCOM-BC2 agree on effect *direction* and both pass a permissive significance threshold (`q < 0.05`/`≤ 0.1` in either method under the rule — see `02_da/da_cohort.R`). This guards against single-method artifacts before any meta-analysis.

**Stage 4 — replication definition.** A taxon is *replicated* when the pooled effect is FDR-significant (`padj < 0.05`) **and** at least 60% of contributing cohort estimates agree in sign (`n_agree ≥ ceil(0.6·k)`).

**Stage 5 — transportability.** A Random Forest (100 trees, `class_weight='balanced'`, features ≥ 20% prevalence in the *training* cohort only) is fit per train/test cohort pair; AUC is computed on the held-out cohort. The pooled-frame model's AUC is the headline transportability statistic.

---

## Quick start

### 1. Environment

```bash
# miniforge/miniconda with mamba recommended
./00_setup/setup_env.sh          # creates 'pd-microbiome' from environment.yml
conda activate pd-microbiome
```

### 2. Fetch real data

```bash
./run_all.sh fetch
# Downloads MetaPhlAn3 abundance tables + curated metadata for all 3 cohorts
# (processed only — no raw reads are stored)
```

### 3. Run the analysis

```bash
./run_all.sh harmonize   # Stage 1: standardize + build count tables
./run_all.sh analyze     # Stage 2+3: DA per cohort, bind, meta-analysis
./run_all.sh ml          # Stage 4: cross-cohort ML (LOCO)
./run_all.sh figures     # Stage 5: publication figures 1–9
# or, end-to-end from a clean tree:
./run_all.sh all
```

Individual stages can be re-run independently; each stage clears its own outputs so a stale result can never silently mix with a fresh run.

---

## Outputs

| Artifact | Produced by | Notes |
|---|---|---|
| `out/replication_scorecard_species.csv` / `_genus.csv` | `04_meta/meta.R` | **Primary result** — pooled estimate, SE, p, FDR, I², k, sign-agreement, replicated flag |
| `out/effects_long_<level>.csv` | `03_bind/bind_effects.R` | Long-format per-cohort effects for pooling |
| `out/<level>/<cohort>/coefs.csv` | `02_da/da_cohort.R` | Per-cohort DA concordance results |
| `out/ml_results/auc_results.txt` | `05_ml/ml_crosscohort.py` | All train→test AUCs |
| `out/ml_results/roc_auc_summary.csv` | figure pipeline | Aggregated LOCO AUC summary |
| `figures_final/figure_1..9_*.png` | figure pipeline | Flowchart, forest, ROC, composition, PCoA, heatmap, diversity, importance, network |

---

## Cohorts

| Cohort | Source study | PMID(s) | Design | Controls |
|---|---|---|---|---|
| `Nielsen2014` | Nielsen HB, *Nat Biotechnol* 2014 | [24997787](https://pubmed.ncbi.nlm.nih.gov/24997787/) | IBD vs healthy, cross-sectional | unrelated healthy |
| `Ijaz2017` | Ijaz UZ, *PLoS One* 2017 | [28222161](https://pubmed.ncbi.nlm.nih.gov/28222161/) | Crohn's disease, cross-sectional | unaffected genetically-linked kindred |
| `HMP2` | Schirmer 2018 / Lloyd-Price 2019 | [29311644](https://pubmed.ncbi.nlm.nih.gov/29311644/), [31142855](https://pubmed.ncbi.nlm.nih.gov/31142855/) | IBD multi-omics, baseline visit only | non-IBD |

Provenance is pinned in `cohorts.tsv` and read from `colData(PMID)` at fetch time — PMIDs are not hand-edited, so a mis-attributed paper cannot silently enter the record.

---

## Requirements

**R (4.3+):** `tidyverse`, `phyloseq`, `metafor`, `Maaslin2`, `ANCOMBC`, `curatedMetagenomicData` (Bioconductor 3.18), `SummarizedExperiment`

**Python (3.11+):** `pandas`, `numpy`, `scikit-learn`, `matplotlib`, `networkx`

All pinned in `environment.yml`.

---

## License

This repository is provided for transparency and reproducibility of the analysis. Input data remain property of the respective cohort studies and are redistributed by `curatedMetagenomicData` under its own terms. See the linked study licenses before any broader redistribution.