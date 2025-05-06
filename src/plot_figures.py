import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from generate_plots.plotting_utils import (
    plot_figure1,
    plot_figure2
)

def main():
    # Create output directory for figures
    figures_dir = "outputs/figures"
    os.makedirs(figures_dir, exist_ok=True)
    
    # Read data - main dataset
    processed_data_matrix = pd.read_csv("outputs/processed_data_matrix.csv", index_col=0)
    design_matrix = pd.read_csv("outputs/design_matrix.csv", index_col=0)

    # Read data - rerun dataset
    rerun_data_matrix = pd.read_csv("outputs/processed_rerun_data_matrix.csv", index_col=0)
    rerun_design_matrix = pd.read_csv("outputs/rerun_design_matrix.csv", index_col=0)
    
    print(f"Main data shape: {processed_data_matrix.shape}")
    print(f"Rerun data shape: {rerun_data_matrix.shape}")
    
    # Create Figure 1 - Main data analysis
    print("\nGenerating Figure 1...")
    print("This figure will include:")
    print("- Venn diagram of S. aureus, P. aeruginosa vs Ctrl peptide overlap")
    print("- Swarmplot of peptide counts by group and day")
    print("- KDE plot of peptide length distribution by group")
    print("- UMAP visualization with group colors and day shapes")
    
    fig1 = plot_figure1(processed_data_matrix, design_matrix, figures_dir)
    
    # Create Figure 2 - Rerun data comparison
    print("\nGenerating Figure 2...")
    print("This figure will include:")
    print("- Venn diagrams for reruns vs corresponding main groups")
    print("- Venn diagram showing overlap between rerun groups")
    print("- Multiple small plots showing per-sample intensity correlations")
    
    fig2 = plot_figure2(processed_data_matrix, design_matrix, 
                      rerun_data_matrix, rerun_design_matrix, figures_dir)
    
    print(f"\nAll figures saved to {figures_dir}")
    print("Figure 1: figure1.png")
    print("Figure 2: figure2.png")

if __name__ == "__main__":
    main()


