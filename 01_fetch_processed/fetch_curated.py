import pandas as pd
import requests
import os
import zipfile
import io
import shutil

ZENODO_URLS = {
    "profiles": "https://zenodo.org/record/14261087/files/SMG_profiles.zip?download=1",
    "metadata": "https://zenodo.org/record/14261087/files/Metadata_All.zip?download=1"
}

def download_and_extract(url, folder):
    print(f"Downloading {url}...")
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            z.extractall(folder)
        return True
    except Exception as e:
        print(f"Download failed: {e}")
        return False

def main():
    os.makedirs("data/taxon", exist_ok=True)
    os.makedirs("data/metadata", exist_ok=True)
    
    temp_taxon = "temp_taxon"
    temp_meta = "temp_meta"
    
    if os.path.exists(temp_taxon): shutil.rmtree(temp_taxon)
    if os.path.exists(temp_meta): shutil.rmtree(temp_meta)
    
    if not download_and_extract(ZENODO_URLS["profiles"], temp_taxon): return
    if not download_and_extract(ZENODO_URLS["metadata"], temp_meta): return
    
    # --- Taxon Processing ---
    print("\nProcessing Taxon files...")
    found_taxon = 0
    for root, dirs, files in os.walk(temp_taxon):
        for file in files:
            if file.endswith(('.csv', '.tsv', '.txt')):
                path = os.path.join(root, file)
                try:
                    # Read with python engine to handle mixed delimiters
                    df = pd.read_csv(path, sep=None, engine='python')
                    
                    # DETECT TRANSPOSE: If columns look like sample IDs (e.g., 'PD.01', 'SRR...')
                    # then the species are the index/first column.
                    first_col_name = df.columns[0]
                    sample_keywords = ['PD', 'CO', 'SRR', 'ERR']
                    is_transposed = any(any(kw in str(col) for kw in sample_keywords) for col in df.columns)
                    
                    if is_transposed:
                        # The first column contains the species names
                        # Set it as index and keep the samples as columns
                        df = df.set_index(df.columns[0])
                        # Rename index to 'species'
                        df.index.name = 'species'
                        # Reset index to make 'species' a column again
                        df = df.reset_index()
                    else:
                        # Normal format: check if first column is species
                        if not any(kw in first_col_name.lower() for kw in ['species', 'taxon', 'motu']):
                            df = df.rename(columns={first_col_name: 'species'})
                    
                    cohort = os.path.splitext(file)[0]
                    df.to_csv(f"data/taxon/{cohort}.csv", index=False)
                    found_taxon += 1
                    print(f"✅ Saved {file} -> data/taxon/{cohort}.csv")
                except Exception as e:
                    print(f"Error reading {file}: {e}")

    # --- Metadata Processing ---
    print("\nProcessing Metadata files...")
    found_meta = 0
    for root, dirs, files in os.walk(temp_meta):
        for file in files:
            if file.endswith(('.csv', '.tsv', '.txt')):
                path = os.path.join(root, file)
                try:
                    df = pd.read_csv(path, sep=None, engine='python')
                    
                    # Sample ID
                    sample_col = None
                    for c in df.columns:
                        if any(kw in c.lower() for kw in ['sample', 'id', 'name']):
                            sample_col = c
                            break
                    if sample_col:
                        df = df.rename(columns={sample_col: 'sample_name'})
                    
                    # Status
                    status_col = None
                    for c in df.columns:
                        if any(kw in c.lower() for kw in ['status', 'diag', 'group', 'cond', 'disease', 'case']):
                            status_col = c
                            break
                    if status_col:
                        df = df.rename(columns={status_col: 'status'})
                    else:
                        df['status'] = 'unknown'
                    
                    cohort = os.path.splitext(file)[0]
                    df.to_csv(f"data/metadata/{cohort}.csv", index=False)
                    found_meta += 1
                    print(f"✅ Saved {file} -> data/metadata/{cohort}.csv")
                except Exception as e:
                    print(f"Error reading {file}: {e}")

    shutil.rmtree(temp_taxon, ignore_errors=True)
    shutil.rmtree(temp_meta, ignore_errors=True)
    print(f"\nDone. Found {found_taxon} taxon tables and {found_meta} metadata tables.")
    if found_taxon > 0:
        print("You can now run: ./run_all.sh all")

if __name__ == "__main__":
    main()
