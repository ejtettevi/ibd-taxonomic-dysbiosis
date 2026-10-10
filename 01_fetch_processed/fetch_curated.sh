#!/bin/bash
# This script downloads processed PD data from curated sources

mkdir -p data/taxon data/metadata

echo "Downloading curated species abundance tables from Hugging Face..."
# These are example URLs based on the metagenomics_mac structure
# In a real scenario, we would iterate through a list of PD cohorts
URLS=(
    "https://huggingface.co/datasets/metagenomics_mac/resolve/main/species_abundance.parquet"
    "https://huggingface.co/datasets/metagenomics_mac/resolve/main/sample_metadata.parquet"
)

for url in "${URLS[@]}"; do
    filename=$(basename "$url")
    echo "Fetching $filename..."
    curl -L "$url" -o "raw_processed/$filename"
done

echo "Processed data fetched. Now filtering for PD cohorts..."
# We will use the Python script to filter these into the pipeline's CSV format
# python 01_fetch_processed/filter_parquet.py ...
