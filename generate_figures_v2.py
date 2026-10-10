import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import os

# --- PUBLICATION QUALITY SETTINGS ---
plt.rcParams['figure.dpi'] = 700
plt.rcParams['savefig.dpi'] = 700
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['DejaVu Serif', 'Times New Roman', 'serif']

OUT_DIR = "/Users/edkwami/Desktop/new/figures"
os.makedirs(OUT_DIR, exist_ok=True)

def generate_figure_1_process_map():
    """Generates a high-end Process Map with inputs and outputs."""
    fig, ax = plt.subplots(figsize=(22, 10))
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 8)
    ax.axis('off')

    # Colors: Gradient from soft grey to deep blue
    colors = ['#F5F5F5', '#E3F2FD', '#BBDEFB', '#64B5F6', '#1565C0']
    text_colors = ['#333333', '#333333', '#222222', '#111111', '#FFFFFF']
    
    # Structure: (Title, Detail, Input, Output, x, y, color)
    stages = [
        ("Data Acquisition", "Public Repositories\n(HMP, Zenodo)", "Raw Tables", "Harmonized CSVs", 1.0, 4, colors[0]),
        ("Harmonization", "Taxonomic Mapping\n& Schema Alignment", "Harmonized CSVs", "Standardized Tables", 5.5, 4, colors[1]),
        ("Differential\nAbundance", "Maaslin2 & ANCOM-BC2\n(Dual-Algorithm)", "Standardized Tables", "Coef Tables", 10.0, 4, colors[2]),
        ("Meta-Analysis", "REML Aggregation\n(Weighted Effect)", "Coef Tables", "Replication Scorecard", 14.5, 4, colors[3]),
        ("Validation", "Random Forest\nCross-Cohort ML", "Scorecard", "Generalizable Signature", 19.0, 4, colors[4])
    ]

    for i, (title, detail, inp, outp, x, y, color) in enumerate(stages):
        # 1. Main Process Box
        box = patches.FancyBboxPatch(
            (x, y-0.8), 3.0, 1.6, 
            boxstyle="round,pad=0.1", 
            facecolor=color, edgecolor='#444444', linewidth=2, zorder=2
        )
        ax.add_patch(box)
        
        # 2. Title and Detail
        ax.text(x + 1.5, y + 0.3, title, ha='center', va='center', fontsize=15, fontweight='bold', color=text_colors[i], zorder=3)
        ax.text(x + 1.5, y - 0.1, detail, ha='center', va='center', fontsize=11, style='italic', color=text_colors[i], zorder=3)
        
        # 3. Input Label (below box)
        ax.text(x + 1.5, y - 1.1, f"Input:\n{inp}", ha='center', va='center', fontsize=10, color='#666666', fontweight='bold', zorder=3)
        
        # 4. Output Label (below box)
        ax.text(x + 1.5, y - 1. same_as_input := 1. same_as_input := 1.5, # adjusted manually
               # I'll just use simple spacing
               )
        # Re-adjusting output labels
        ax.text(x + 1.5, y - 1.8, f"Output:\n{outp}", ha='center', va='center', fontsize=10, color='#444444', fontweight='bold', zorder=3)

    # 5. Professional Arrows
    for i in range(len(stages)-1):
        start_x = stages[i][4] + 3.0 + 0.1
        end_x = stages[i+1][4] - 0.1
        ax.annotate('', xy=(end_x, y), xytext=(start_x, y),
                    arrowprops=dict(arrowstyle='-|>', color='#555555', lw=3, mutation_scale=25), zorder=1)

    plt.title("Figure 1: Integrated Meta-Analytical Pipeline for Microbiome Biomarker Discovery", 
              fontsize=24, fontweight='bold', pad=60, color='#222222')
    
    plt.savefig(f"{OUT_DIR}/figure_1_flowchart.png", bbox_inches='tight', transparent=False, facecolor='white')
    plt.close()
    print("✅ Figure 1 (High-End Process Map) saved.")

if __name__ == "__main__":
    generate_figure_1_process_map()
