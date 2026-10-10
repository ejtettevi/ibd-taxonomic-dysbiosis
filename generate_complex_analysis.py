import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import networkx as nx
from sklearn.ensemble import RandomForestClassifier

# --- PUBLICATION QUALITY SETTINGS ---
plt.rcParams['figure.dpi'] = 700
plt.rcParams['savefig.dpi'] = 700
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['DejaVu Serif', 'Times New Roman', 'serif']

OUT_DIR = "/Users/edkwami/Desktop/new/figures"
os.makedirs(OUT_DIR, exist_ok=True)

CONSENSUS_TAXA = [
    'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Lachnospiraceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii',
    'k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella.s__Prevotella_copri',
    'k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Actinomycetaceae.g__Actinomyces.s__Actinomyces_oris',
    'k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae.g__Streptococcus.s__Streptococcus_mutans',
    'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Ruthenibacterium.s__Ruthenibacterium_lactatiformans',
    'k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides.s__Bacteroides_fragilis'
]

def load_consolidated_data():
    cohorts = ['hmp_ibd', 'crohns_study', 'uc_study']
    all_data = []
    for cohort in cohorts:
        tax_path = f"data/taxon/{cohort}.csv"
        meta_path = f"data/metadata/{cohort}.csv"
        if os.path.exists(tax_path) and os.path.exists(meta_path):
            tax = pd.read_csv(tax_path, index_col=0)
            meta = pd.read_csv(meta_path, index_col='sample_name')
            df = tax.T.join(meta)
            df['cohort'] = cohort
            all_data.append(df)
    return pd.concat(all_data) if all_data else None

def plot_alpha_diversity():
    df = load_consolidated_data()
    if df is None: return
    taxa_cols = [c for c in df.columns if c not in ['status', 'age', 'sex', 'cohort']]
    def shannon(row):
        p = row / row.sum()
        p = p[p > 0]
        return -np.sum(p * np.log(p))
    df['shannon'] = df[taxa_cols].apply(shannon, axis=1)
    plt.figure(figsize=(6, 6))
    sns.boxplot(x='status', y='shannon', data=df, palette={'PD': '#d62728', 'control': '#1f77b4'}, width=0.5)
    sns.stripplot(x='status', y='shannon', data=df, color='black', alpha=0.3, size=4)
    plt.title("Figure 7: Alpha Diversity (Shannon Index)", fontsize=16, fontweight='bold')
    plt.xlabel("Status", fontsize=12); plt.ylabel("Shannon Index", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_7_diversity.png")
    plt.close()
    print("✅ Figure 7 (Diversity) saved.")

def plot_feature_importance():
    df = load_consolidated_data()
    if df is None: return
    X = df[CONSENSUS_TAXA]
    y = df['status'].map({'control': 0, 'PD': 1})
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rf.fit(X, y)
    importances = rf.feature_importances_
    indices = np.argsort(importances)
    plt.figure(figsize=(10, 6))
    plt.barh(range(len(indices)), importances[indices], color='#42A5F6', edgecolor='black')
    plt.yticks(range(len(indices)), [CONSENSUS_TAXA[i].split('.g__')[-1].split('.s__')[0] for i in indices], fontsize=12)
    plt.xlabel("Gini Importance", fontsize=14)
    plt.title("Figure 8: Feature Importance of Consensus Biomarkers", fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/figure_8_importance.png")
    plt.close()
    print("✅ Figure 8 (Importance) saved.")

def plot_cooccurrence_network():
    df = load_consolidated_data()
    if df is None: return
    def get_corr_net(subset):
        corr = subset[CONSENSUS_TAXA].corr()
        return corr
    case_corr = get_corr_net(df[df['status'] == 'PD'])
    ctrl_corr = get_corr_net(df[df['status'] == 'control'])
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for i, (corr, title) in enumerate([(ctrl_corr, "Healthy Control"), (case_corr, "IBD Case")]):
        G = nx.from_pandas_adjacency(corr > 0.3)
        pos = nx.spring_layout(G, seed=42)
        labels = {CONSENSUS_TAXA[j]: CONSENSUS_TAXA[j].split('.g__')[-1].split('.s__')[0] for j in range(len(CONSENSUS_TAXA))}
        nx.draw_networkx_nodes(G, pos, node_size=800, node_color='#BBDEFB', edgecolors='black', ax=axes[i])
        nx.draw_networkx_edges(G, pos, width=2, edge_color='gray', alpha=0.6, ax=axes[i])
        nx.draw_networkx_labels(G, pos, labels, font_size=10, font_weight='bold', ax=axes[i])
        axes[i].set_title(title, fontsize=16, fontweight='bold')
        axes[i].axis('off')
    plt.suptitle("Figure 9: Consensus Taxa Co-occurrence Networks", fontsize=20, fontweight='bold', y=1.05)
    plt.savefig(f"{OUT_DIR}/figure_9_network.png", bbox_inches='tight')
    plt.close()
    print("✅ Figure 9 (Network) saved.")

if __name__ == "__main__":
    print("🚀 Generating Complex Analysis Figures...")
    plot_alpha_diversity()
    plot_feature_importance()
    plot_cooccurrence_network()
    print(f"\n✨ Complex analysis figures saved to: {OUT_DIR}")
