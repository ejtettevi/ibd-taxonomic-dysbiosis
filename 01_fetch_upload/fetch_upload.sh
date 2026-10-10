#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <cohort> <bioproject> <project_id>"
  exit 1
fi

cohort="$1"
bioproject="$2"
project_id="$3"
outdir="raw/${cohort}"
mkdir -p "$outdir"

for cmd in pysradb prefetch fasterq-dump gzip seqtoid; do
  if ! command -v "$cmd" &>/dev/null; then
    echo "Error: $cmd not found." >&2
    exit 1
  fi
done

pysradb metadata "$bioproject" --detailed --saveto "${outdir}/sra_meta.tsv"
col=$(head -n 1 "${outdir}/sra_meta.tsv" | tr '\t' '\n' | grep -nxi "run_accession" | cut -d: -f1 || true)
awk -v c="$col" -F'\t' 'NR>1 {print $c}' "${outdir}/sra_meta.tsv" | tr -d '\r' | sort -u | grep -v '^$' > "${outdir}/runs.txt"

while read -r acc; do
  [ -z "$acc" ] && continue
  echo "Processing $acc..."
  
  # Cleanup lock
  rm -f "${outdir}/${acc}/${acc}.sra.lock" 2>/dev/null || true
  
  prefetch "$acc" -O "$outdir"
  if [ -d "${outdir}/${acc}" ]; then
    fasterq-dump "${outdir}/${acc}" -O "$outdir" --split-files -e 4
  else
    fasterq-dump "$acc" -O "$outdir" --split-files -e 4
  fi
  find "$outdir" -name "${acc}*.fastq" -exec gzip -f {} +
  
  files=$(find "$outdir" -name "${acc}*.fastq.gz")
  for f in $files; do
    seqtoid metagenomics upload-sample --project "$project_id" --sample-name "$acc" --metadata-csv "metadata/${cohort}.csv" --sequencing-platform 'Illumina' "$f"
  done
  
  rm -rf "${outdir}/${acc}" "${outdir}/${acc}.sra"
  find "$outdir" -name "${acc}*.fastq.gz" -delete
done < "${outdir}/runs.txt"
