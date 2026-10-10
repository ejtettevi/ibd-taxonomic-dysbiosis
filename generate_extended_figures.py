import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
from sklearn.decomposition import PCA

# --- PUBLICATION QUALITY SETTINGS ---
plt.rcParams['figure.dpi'] = 700
plt.rcParams['savefig.dpi'] = 700
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['DejaVu Serif', 'Times New Roman', 'serif']

OUT_DIR = "/Users/edkwami/Desktop/new/figures"
os.makedirs(OUT_DIR, exist_ok=True)

# List of the 6 consensus biomarkers discovered in the pipeline
CONSENSUS_TAXA = [
    'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Lachnospiraceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii',
    'k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella.s__Prevotella_copri',
    'k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Actinomycetaceae.g__Actinomyces.s__Actinomyces_oris',
    'k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae.g__Streptococcus.s__Streptococcus_mutans',
    'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruthenibacterium.s__Ruthenibacterium_lactatiformans',
    'k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.s__Bacteroides_fragilis'
]

def load_data():
    cohorts = ['hmp_ibd', 'crohns_study', 'uc_study']
    all_data = []
    for cohort in cohorts:
        tax_path = f"data/taxon/{cohort}.csv"
        meta_path = f"data/metadata/{cohort}.csv"
        if os.path.exists(tax_path) and os.path.exists(meta_path):
            tax = pd.read_csv(tax_path, index_col=0)
            meta = pd.read_csv(meta_path, index_col='sample_name')
            # Filter for consensus taxa
            common_taxa = [t for t in CONSENSUS_TAXA if t in tax.index]
            tax_filtered = tax.loc[common_taxa]
            
            # Merge with metadata
            df = tax_filtered.T.join(meta)
            df['cohort'] = cohort
            all_data.append(df)
    return pd.concat(all_data) if all_data else None

def plot_composition_bars():
    """Figure 4: Relative abundance of consensus taxa by status across cohorts."""
    df = load_data()
    if df is None: return
    
    # Group by cohort and status, calculate mean
    grouped = df.groupby(['cohort', 'status'])[CONSENSUS_TAXA].mean().reset_index()
    
    # Melt for seaborn
    melted = grouped.melt(id_vars=['cohort', 'status'], var_name='Taxon', value_name='Mean Abundance')
    # Shorten taxon names for the plot
    melted['Taxon_Short'] = melted['Taxon'].apply(lambda x: x.split('.g__')[-1].split('.s__')[0])
    
    plt.figure(figsize=(14, 7))
    sns.barplot(data=melted, x='cohort', y='Mean Abundance', hue='status', palette={'PD': '#d62728', 'control': '#1f77b4'}, edgecolor='black')
    
    plt.title("Figure 4: Relative Abundance of Consensus Biomarkers", fontsize=18, fontweight='bold', pad=20)
    plt.xlabel("Cohort", fontsize=14)
    plt.ylabel("Mean Relative Abundance", fontsize=14)
    plt.legend(title="Status", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_4_composition.png", bbox_inches='tight')
    plt.close()
    print("✅ Figure 4 (Composition Plot) saved.")

def plot_pcoa():
    """Figure 5: Beta-Diversity (PCoA/PCA) showing separation of Case vs Control."""
    df = load_data()
    if df is None: return
    
    X = df[CONSENSUS_TAXA]
    y = df['status']
    
    pca = PCA(n_components=2)
    components = pca.fit_transform(X)
    
    plt.figure(figsize=(8, 8))
    sns.scatterplot(x=components[:, 0], y=components[:, 1], hue=y, style=df['cohort'], 
                    palette={'PD': '#d62728', 'control': '#1f77b4'}, s=100, edgecolor='black', alpha=0.8)
    
    plt.title("Figure 5: PCoA of Consensus Microbial Signature", fontsize=17, fontweight='bold', pad=20)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)", fontsize=14)
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)", fontsize=14)
    plt.legend(title="Status / Cohort", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_5_pcoa.png", bbox_inches='tight')
    plt.close()
    print("✅ Figure 5 (PCoA Plot) saved.")

def plot_correlation_heatmap():
    """Figure 6: Correlation heatmap between the consensus biomarkers."""
    df = load_data()
    if df is None: return
    
    corr = df[CONSENSUS_TAXA].corr()
    # Shorten names
    corr.index = [x.split('.g__')[-1].split('.s__')[0] for x in corr.index]
    corr.columns = [x.split('.g__')[-1].split('.s__')[0] for x in corr.columns]
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr, annot=True, cmap='RdBu_r', center=0, vmin=-1, vmax=1, 
                fmt=".2f", linewidths=0.5, annot_kws={"size": 12})
    
    plt.title("Figure 6: Correlation Matrix of Consensus Biomarkers", fontsize=17, fontweight='bold', pad=20)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_6_heatmap.png", bbox_inches='tight')
    plt.close()
    print("✅ Figure 6 (Correlation Heatmap) saved.")

if __name__ == "__main__":
    print("🚀 Generating Extended Analysis Figures...")
    plot_composition_bars()
    plot_pcoa()
    plot_correlation_heatmap()
    print(f"\n✨ All extended figures saved to: {OUT_DIR}")
