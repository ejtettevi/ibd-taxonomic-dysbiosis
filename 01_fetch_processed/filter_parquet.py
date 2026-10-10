import pandas as pd
import sys
import os

def filter_cohort(parquet_path, output_csv, cohort_samples):
    print(f"Processing {parquet_path}...")
    df = pd.read_parquet(parquet_path)
    
    # Assuming the index or a column 'sample' contains the IDs
    if 'sample' in df.columns:
        filtered_df = df[df['sample'].isin(cohort_samples)]
    else:
        filtered_df = df[df.index.isin(cohort_samples)]
    
    filtered_df.to_csv(output_csv)
    print(f"Saved {len(filtered_df)} samples to {output_csv}")

if __name__ == "__main__":
    # Usage: python filter_parquet.py <parquet_file> <output_csv> <sample_list_csv>
    p_path = sys.argv[1]
    o_path = sys.argv[2]
    s_path = sys.argv[3]
    
    samples = pd.read_csv(s_path)['sample'].tolist()
    filter_cohort(p_path, o_path, samples)
