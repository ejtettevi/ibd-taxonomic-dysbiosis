#!/usr/bin/env bash
set -euo pipefail

# Resolve Rscript and Python portably: prefer executables already on PATH
# (e.g. an activated conda env), otherwise fall back to the environment named
# in environment.yml. No absolute machine-specific paths are hardcoded.
find_cmd() {
  local name="$1"
  if command -v "$name" >/dev/null 2>&1; then
    command -v "$name"
  elif [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/$name" ]; then
    echo "$CONDA_PREFIX/bin/$name"
  else
    echo "Error: $name not found on PATH and not under CONDA_PREFIX" >&2
    echo "Hint: run 'conda activate pd-microbiome' first." >&2
    exit 1
  fi
}
RSCRIPT="$(find_cmd Rscript)"
PYTHON="$(find_cmd python)"

stage="${1:-all}"
COHORTS_FILE="cohorts.tsv"

if [ ! -f "$COHORTS_FILE" ]; then
  echo "Error: $COHORTS_FILE not found!" >&2
  echo "Hint: run './run_all.sh fetch' once to download the processed tables." >&2
  exit 1
fi

run_fetch() {
  echo "==> [Stage 0: Fetch] Downloading processed tables from curatedMetagenomicData..."
  # Requires network access. Downloads curated MetaPhlAn3 profiles + metadata only
  # (no raw reads) and rewrites cohorts.tsv with verified provenance.
  $RSCRIPT ./01_fetch_processed/fetch_curated_metagenomic.R
}

run_harmonize() {
  echo "==> [Stage 1: Harmonize] Standardizing processed tables..."
  $RSCRIPT ./01_harmonize/harmonize_tables.R
}

run_analyze() {
  echo "==> [Stage 2: Analyze] Running DA & Meta-Analysis..."
  # Clear previous per-level results so a re-analysis cannot silently mix outputs
  # from a previous taxon-labelling scheme with the current one.
  rm -f out/replication_scorecard_species.csv out/replication_scorecard_genus.csv \
        out/effects_long_species.csv out/effects_long_genus.csv
  
  for level in "species" "genus"; do
    echo "Analyzing at $level level..."
    # Get list of cohorts by finding directories in data/taxon (without .csv extension)
    cohorts=($(find data/taxon -maxdepth 1 -mindepth 1 -type d -exec basename {} \;))
    
    for cohort in "${cohorts[@]}"; do
      $RSCRIPT ./02_da/da_cohort.R "$cohort" "$level"
    done
    
    $RSCRIPT ./03_bind/bind_effects.R "$level"
    $RSCRIPT ./04_meta/meta.R "$level"
  done
}

run_ml() {
  echo "==> [Stage 3: ML] Running Cross-Cohort Validation (Genus Level)..."
  rm -f out/ml_results/auc_results.txt
  mkdir -p out/ml_results
  
  # Use directories for cohorts to avoid .csv extension issues
  cohorts=($(find data/taxon -maxdepth 1 -mindepth 1 -type d -exec basename {} \;))
  n=${#cohorts[@]}
  for ((i=0; i<n; i++)); do
    for ((j=0; j<n; j++)); do
      if [ $i -ne $j ]; then
        $PYTHON ./05_ml/ml_crosscohort.py "${cohorts[$i]}" "${cohorts[$j]}" "genus"
      fi
    done
  done
}

run_figures() {
  echo "==> [Stage 4: Figures] Rendering publication figures..."
  $PYTHON ./generate_final_publication_figures.py
}

run_ms() {
  echo "==> [Stage 5: Manuscript] Generating draft..."
  $PYTHON ./06_manuscript/generate_ms.py
}

run_docx() {
  echo "==> [Stage 6: DOCX] Converting manuscript to Word..."
  if ! command -v pandoc >/dev/null 2>&1; then
    echo "Error: pandoc not found; cannot produce .docx" >&2
    exit 1
  fi
  # Regenerate the markdown if it is missing or older than the generator, otherwise
  # the .docx would silently reflect a stale draft.
  if [ ! -f out/manuscript_submission.md ] || \
     [ 06_manuscript/generate_ms.py -nt out/manuscript_submission.md ]; then
    $PYTHON ./06_manuscript/generate_ms.py
  fi
  # Reference template supplies Candara 11 pt and justified body text.
  [ -f docx/reference.docx ] || $PYTHON ./06_manuscript/make_reference_docx.py
  pandoc out/manuscript_submission.md -o out/manuscript_draft.docx \
    --from markdown+implicit_figures \
    --reference-doc docx/reference.docx
  echo "Wrote out/manuscript_draft.docx"
}

case "$stage" in
  fetch) run_fetch ;;
  harmonize) run_harmonize ;;
  analyze) run_analyze ;;
  ml) run_ml ;;
  figures) run_figures ;;
  ms) run_ms ;;
  docx) run_docx ;;
  all)
    run_harmonize
    run_analyze
    run_ml
    run_figures
    run_ms
    run_docx
    ;;
  *) echo "Usage: $0 [fetch|harmonize|analyze|ml|figures|ms|docx|all]"; exit 1 ;;
esac
