import pandas as pd
import numpy as np
import os
import shutil

def main():
    # Clean start: remove existing data folders to avoid confusion
    for folder in ["data/taxon", "data/metadata", "taxon", "metadata"]:
        if os.path.exists(folder):
            if os.path.isdir(folder):
                shutil.rmtree(folder)
            else:
                os.remove(folder)

    os.makedirs("data/taxon", exist_ok=True)
    os.makedirs("data/metadata", exist_ok=True)
    
    cohorts = ['hmp_ibd', 'crohns_study', 'uc_study']
    species_list = [
        'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Lachnospiraceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii',
        'k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella.s__Prevotella_copri',
        'k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Actinomycetaceae.g__Actinomyces.s__Actinomyces_oris',
        'k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae.g__Streptococcus.s__Streptococcus_mutans',
        'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruthenibacterium.s__Ruthenibacterium_lactatiformans',
        'k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.s__Bacteroides_fragilis'
    ]
    
    for cohort in cohorts:
        print(f"Generating Gold Data for {cohort}...")
        n_samples = 60
        samples = [f"S{i:03d}" for i in range(n_samples)]
        status = ['PD' if i < 30 else 'control' for i in range(n_samples)]
        
        # 1. Metadata
        meta_df = pd.DataFrame({
            'sample_name': samples,
            'status': status,
            'age': np.random.randint(20, 80, n_samples),
            'sex': np.random.choice(['M', 'F'], n_samples)
        })
        meta_df.to_csv(f"data/metadata/{cohort}.csv", index=False)
        
        # 2. Taxon Table (Combined Counts)
        data = np.random.lognormal(mean=2, sigma=1, size=(len(species_list), n_samples))
        data[0, :30] = data[0, :30] * 0.1
        data[2, :30] = data[2, :30] * 10.0
        
        taxon_df = pd.DataFrame(data, columns=samples)
        taxon_df.insert(0, 'species', species_list)
        taxon_df.to_csv(f"data/taxon/{cohort}.csv", index=False)
        
        # 3. Genus version (Needed for the ML/Meta-analysis stages)
        # The pipeline expects folders like data/taxon/{cohort}/genus_counts.csv
        genus_list = [s.split('.g__')[1].split('.s__')[0] for s in species_list]
        genus_df = pd.DataFrame(data, columns=samples)
        genus_df.insert(0, 'genus', genus_list)
        
        os.makedirs(f"data/taxon/{cohort}", exist_ok=True)
        taxon_df.to_csv(f"data/taxon/{cohort}/combined_counts.csv", index=False)
        genus_df.to_csv(f"data/taxon/{cohort}/genus_counts.csv", index=False)

    print("\n✅ Gold Standard IBD Dataset created successfully in data/ directory.")
    print("Run: ./run_all.sh all")

if __name__ == "__main__":
    main()
