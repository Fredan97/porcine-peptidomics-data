import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional
from matplotlib_venn import venn3, venn3_circles, venn2
import numpy as np
import numpy as np
from sklearn.preprocessing import StandardScaler
from umap import UMAP
import math
import matplotlib.gridspec as gridspec

COLOR_PALETTE = {
    "S. aureus": "#ff9e6d",  # Light orange
    "P. aeruginosa": "#5b8ff9",  # Medium blue
    "Ctrl": "#61c0bf",  # Teal
    "Double infection": "#b264c5",  # Purple
    "Accidental double infection": "#ffcc5c",  # Yellow
}

LIGHT_BLUE = "#c6dbef"
MEDIUM_BLUE = "#6baed6"
DARK_BLUE = "#2171b5"

LIGHT_ORANGE = "#fdd0a2"
MEDIUM_ORANGE = "#fd8d3c"
DARK_ORANGE = "#d94801"


sns.set_context("paper")


def plot_peptide_counts_swarm(
    data_matrix: pd.DataFrame,
    design_matrix: pd.DataFrame,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Create swarmplot of peptide counts by group and day.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))

    common_samples = list(
        set(data_matrix.columns).intersection(set(design_matrix.index))
    )
    data_subset = data_matrix[common_samples]
    design_subset = design_matrix.loc[common_samples]

    peptide_stats = []
    for sample in common_samples:
        peptide_count = data_subset[sample].dropna().shape[0]
        group = design_subset.loc[sample, "group"]
        day = design_subset.loc[sample, "day"]

        peptide_stats.append(
            {
                "Sample": sample,
                "Peptide Count": peptide_count,
                "Group": group,
                "Day": day,
            }
        )

    stats_df = pd.DataFrame(peptide_stats)
    sns.swarmplot(
        data=stats_df,
        x="Group",
        y="Peptide Count",
        hue="Day",
        palette="Set2",
        size=8,
        ax=ax,
    )

    ax.set_xlabel("Group")
    ax.set_ylabel("Number of Peptides")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.legend(frameon=False, loc="upper right")
    sns.despine(ax=ax)

    return ax


def plot_length_distribution_kdeplot(
    data_matrix: pd.DataFrame,
    design_matrix: pd.DataFrame,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Create KDE plots of peptide length distribution by group.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))

    common_samples = list(
        set(data_matrix.columns).intersection(set(design_matrix.index))
    )
    data_subset = data_matrix[common_samples]
    design_subset = design_matrix.loc[common_samples]

    group_peptides = {}
    for group in design_subset["group"].unique():
        group_samples = design_subset[design_subset["group"] == group].index
        peptides = set()
        for sample in group_samples:
            if sample in data_subset.columns:
                sample_peptides = data_subset[sample].dropna().index.tolist()
                peptides.update(sample_peptides)

        group_peptides[group] = list(peptides)
    group_lengths = {}
    for group, peptides in group_peptides.items():
        lengths = [len(peptide) for peptide in peptides]
        group_lengths[group] = lengths

    for group, lengths in group_lengths.items():
        sns.kdeplot(
            lengths,
            label=group,
            ax=ax,
            color=COLOR_PALETTE.get(group, "gray"),
            fill=True,
            alpha=0.3,
            linewidth=2,
        )

    ax.set_xlabel("Peptide Length (amino acids)")
    ax.set_ylabel("Density")
    ax.legend(frameon=False)
    sns.despine(ax=ax)

    return ax


def plot_rerun_group_overlap_venn(
    rerun_data: pd.DataFrame,
    rerun_design: pd.DataFrame,
    main_design: pd.DataFrame,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Create Venn diagram showing overlap between rerun groups.
    """

    rerun_groups = {}
    for sample in rerun_design.index:
        if sample in main_design.index:
            group = main_design.loc[sample, "group"]
            if group not in rerun_groups:
                rerun_groups[group] = []
            rerun_groups[group].append(sample)

    group_peptides = {}
    for group, samples in rerun_groups.items():
        peptides = set()
        for sample in samples:
            matching_rerun_cols = [col for col in rerun_data.columns if sample in col]
            for col in matching_rerun_cols:
                sample_peptides = rerun_data[col].dropna().index
                peptides.update(sample_peptides)

        group_peptides[group] = peptides

    if len(group_peptides) == 3:
        group_names = list(group_peptides.keys())
        venn = venn3(
            [
                group_peptides[group_names[0]],
                group_peptides[group_names[1]],
                group_peptides[group_names[2]],
            ],
            set_labels=group_names,
            ax=ax,
        )

        if venn:
            for i, group in enumerate(group_names):
                if group in COLOR_PALETTE:
                    idx = ["100", "010", "001"][i]
                    if venn.get_patch_by_id(idx):
                        venn.get_patch_by_id(idx).set_color(COLOR_PALETTE[group])

    elif len(group_peptides) == 2:
        group_names = list(group_peptides.keys())
        venn = venn2(
            [group_peptides[group_names[0]], group_peptides[group_names[1]]],
            set_labels=group_names,
            ax=ax,
        )

        if venn:
            for i, group in enumerate(group_names):
                if group in COLOR_PALETTE:
                    idx = ["10", "01"][i]
                    if venn.get_patch_by_id(idx):
                        venn.get_patch_by_id(idx).set_color(COLOR_PALETTE[group])

    sns.despine(ax=ax, left=True, bottom=True)

    return ax


def plot_rerun_intensity_per_sample(
    main_data: pd.DataFrame,
    rerun_data: pd.DataFrame,
    rerun_design: pd.DataFrame,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Create small scatter plots for each rerun sample vs its main counterpart.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 10))

    rerun_to_main = {}
    for rerun_sample in rerun_design.index:
        if rerun_sample in main_data.columns:
            rerun_to_main[rerun_sample] = rerun_sample

    num_samples = len(rerun_to_main)
    if num_samples == 0:
        ax.text(
            0.5,
            0.5,
            "No matching samples found between main and rerun data",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        return ax

    grid_cols = min(3, num_samples)
    grid_rows = math.ceil(num_samples / grid_cols)

    ax.clear()
    ax.axis("off")

    gs = gridspec.GridSpecFromSubplotSpec(
        grid_rows, grid_cols, subplot_spec=ax.get_subplotspec(), wspace=0.4, hspace=0.4
    )

    for i, (rerun_sample, main_sample) in enumerate(rerun_to_main.items()):
        rerun_cols = [col for col in rerun_data.columns if rerun_sample in col]

        if not rerun_cols or main_sample not in main_data.columns:
            continue

        rerun_col = rerun_cols[0]

        rerun_peptides = set(rerun_data[rerun_col].dropna().index)
        main_peptides = set(main_data[main_sample].dropna().index)
        common_peptides = rerun_peptides.intersection(main_peptides)

        # Create subplot
        row = i // grid_cols
        col = i % grid_cols

        # Create subplot with GridSpec
        subax = plt.Subplot(ax.figure, gs[row, col])
        ax.figure.add_subplot(subax)

        # Collect intensity pairs
        main_intensities = []
        rerun_intensities = []

        for peptide in common_peptides:
            rerun_intensity = rerun_data.loc[peptide, rerun_col]
            main_intensity = main_data.loc[peptide, main_sample]

            main_intensities.append(main_intensity)
            rerun_intensities.append(rerun_intensity)

        subax.scatter(
            main_intensities,
            rerun_intensities,
            alpha=0.6,
            s=10,
            c=COLOR_PALETTE.get("P. aeruginosa", "#5b8ff9"),
            edgecolor="none",
        )

        min_int = min(min(main_intensities), min(rerun_intensities))
        max_int = max(max(main_intensities), max(rerun_intensities))
        subax.plot(
            [min_int, max_int],
            [min_int, max_int],
            color="#cccccc",
            linestyle="--",
            linewidth=1,
            alpha=0.7,
        )

        sample_name = rerun_sample
        if len(sample_name) > 10:
            sample_name = sample_name.split("_")[-1]

        subax.text(
            0.5,
            0.95,
            sample_name,
            ha="center",
            va="top",
            transform=subax.transAxes,
            fontsize=8,
        )

        subax.set_xlim(min_int * 0.8, max_int * 1.2)
        subax.set_ylim(min_int * 0.8, max_int * 1.2)

        if row == grid_rows - 1:  # Bottom row
            subax.set_xlabel("log(I)", fontsize=8)
        else:
            subax.set_xticklabels([])

        if col == 0:  # First column
            subax.set_ylabel("log(I)", fontsize=8)
        else:
            subax.set_yticklabels([])

        subax.tick_params(axis="both", which="major", labelsize=6)

        # Apply despine to each subplot
        sns.despine(ax=subax)
    ax.set_xlabel("Main")
    ax.set_ylabel("Rerun")

    return ax


def plot_group_peptide_venn(
    data_matrix: pd.DataFrame,
    design_matrix: pd.DataFrame,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Plot Venn diagram of peptide overlap between different treatment groups (S. aureus, P. aeruginosa, Ctrl).
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))

    common_samples = list(
        set(data_matrix.columns).intersection(set(design_matrix.index))
    )
    data_subset = data_matrix[common_samples]
    design_subset = design_matrix.loc[common_samples]

    group_peptides = {}
    for group in ["S. aureus", "P. aeruginosa", "Ctrl"]:
        group_samples = design_subset[design_subset["group"] == group].index
        peptides = set()
        for sample in group_samples:
            if sample in data_subset.columns:

                sample_peptides = data_subset[sample].dropna().index.tolist()
                peptides.update(sample_peptides)

        group_peptides[group] = peptides

    if len(group_peptides) >= 3:
        venn = venn3(
            [
                group_peptides["S. aureus"],
                group_peptides["P. aeruginosa"],
                group_peptides["Ctrl"],
            ],
            ("S. aureus", "P. aeruginosa", "Ctrl"),
            ax=ax,
        )

        venn.get_patch_by_id("100").set_color(COLOR_PALETTE["S. aureus"])
        venn.get_patch_by_id("010").set_color(COLOR_PALETTE["P. aeruginosa"])
        venn.get_patch_by_id("001").set_color(COLOR_PALETTE["Ctrl"])

        venn.get_label_by_id("100").set_color(DARK_ORANGE)
        venn.get_label_by_id("010").set_color(DARK_BLUE)
        venn.get_label_by_id("001").set_color("#2c7c7c")

        venn3_circles(
            [
                group_peptides["S. aureus"],
                group_peptides["P. aeruginosa"],
                group_peptides["Ctrl"],
            ],
            ax=ax,
            linewidth=0.5,
        )

    ax.set_title("")
    sns.despine(ax=ax, left=True, bottom=True)

    return ax


def plot_umap_with_day_shape(
    data_matrix: pd.DataFrame,
    design_matrix: pd.DataFrame,
    highlight_samples: Optional[list] = None,
    ax: Optional[plt.Axes] = None,
    legend_above: bool = False,
) -> plt.Axes:
    """
    Create UMAP plot with samples colored by group and shaped by day.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))

    matrix = data_matrix.copy().T

    matrix = matrix.dropna(axis=1, how="all")

    matrix = matrix.fillna(0)

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(matrix)

    reducer = UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
    embedding = reducer.fit_transform(scaled_data)

    common_samples = list(set(matrix.index).intersection(set(design_matrix.index)))
    sample_indices = [matrix.index.get_loc(sample) for sample in common_samples]
    embedding_subset = embedding[sample_indices]
    design_subset = design_matrix.loc[common_samples]

    day_markers = {"Day 1": "o", "Day 2": "^", "Day 3": "s"}

    group_colors = COLOR_PALETTE

    for group in design_subset["group"].unique():
        for day in design_subset["day"].unique():
            mask = (design_subset["group"] == group) & (design_subset["day"] == day)

            points_indices = np.where(mask)[0]

            ax.scatter(
                embedding_subset[points_indices, 0],
                embedding_subset[points_indices, 1],
                c=group_colors.get(group, "gray"),
                marker=day_markers.get(day, "x"),
                s=70,
                alpha=0.7,
                label=f"{group} - {day}",
            )

    if highlight_samples:
        highlight_indices = [
            i for i, sample in enumerate(common_samples) if sample in highlight_samples
        ]

        if highlight_indices:
            ax.scatter(
                embedding_subset[highlight_indices, 0],
                embedding_subset[highlight_indices, 1],
                s=120,
                facecolors="none",
                edgecolors="red",
                linewidth=1.5,
                label="Rerun samples",
            )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    legend_elements = []

    for group, color in group_colors.items():
        if group in design_subset["group"].unique():
            legend_elements.append(
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=color,
                    markersize=8,
                    label=group,
                )
            )

    for day, marker in day_markers.items():
        if day in design_subset["day"].unique():
            legend_elements.append(
                plt.Line2D(
                    [0], [0], marker=marker, color="black", markersize=8, label=day
                )
            )

    if highlight_samples:
        legend_elements.append(
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="w",
                markeredgecolor="red",
                markeredgewidth=1.5,
                markersize=8,
                label="Rerun Samples",
            )
        )

    if legend_above:
        ax.legend(
            handles=legend_elements,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.15),
            ncol=3,
        )
    else:
        ax.legend(handles=legend_elements, frameon=False, loc="best", ncol=3)

    sns.despine(ax=ax)

    return ax


def plot_rerun_venn_diagrams(
    main_data: pd.DataFrame,
    main_design: pd.DataFrame,
    rerun_data: pd.DataFrame,
    rerun_design: pd.DataFrame,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Create Venn diagrams comparing rerun samples to their corresponding main samples by group.
    """

    rerun_samples = rerun_design.index

    rerun_groups = {}
    for sample in rerun_samples:
        if sample in main_design.index:
            group = main_design.loc[sample, "group"]
            if group not in rerun_groups:
                rerun_groups[group] = []
            rerun_groups[group].append(sample)

    num_groups = len(rerun_groups)
    for i, (group, samples) in enumerate(rerun_groups.items()):
        spacing = 0.05
        width = (1.0 - spacing * (num_groups - 1)) / num_groups
        x_pos = i * (width + spacing)

        subax = ax.inset_axes([x_pos, 0, width, 1])

        main_group_samples = main_design[main_design["group"] == group].index
        main_peptides = set()
        for sample in main_group_samples:
            if sample in main_data.columns:
                sample_peptides = main_data[sample].dropna().index
                main_peptides.update(sample_peptides)

        rerun_peptides = set()
        for sample in samples:
            matching_rerun_cols = [col for col in rerun_data.columns if sample in col]
            for col in matching_rerun_cols:
                sample_peptides = rerun_data[col].dropna().index
                rerun_peptides.update(sample_peptides)

        venn = venn2([main_peptides, rerun_peptides], (f"Main", f"Rerun"), ax=subax)

        if venn:
            main_color = COLOR_PALETTE.get(group, "#aaaaaa")
            rerun_color = "#9ecae1"  # Light blue for all reruns

            venn.get_patch_by_id("10").set_color(main_color)
            venn.get_patch_by_id("01").set_color(rerun_color)
            venn.get_patch_by_id("11").set_color("#b3b3b3")

        sns.despine(ax=subax, left=True, bottom=True)

    ax.axis("off")

    return ax


def plot_figure1(
    data_matrix: pd.DataFrame, design_matrix: pd.DataFrame, output_dir: str
) -> plt.Figure:
    """
    Create the first composite figure as requested.
    """

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(
        2, 2, figure=fig, height_ratios=[1, 1.2], width_ratios=[1, 1.1]
    )

    # Top left: Venn diagram of peptide overlap between groups
    print("Creating group peptide Venn diagram...")
    ax1 = fig.add_subplot(gs[0, 0])
    plot_group_peptide_venn(data_matrix, design_matrix, ax=ax1)

    # Top right: Swarmplot of peptide counts
    print("Creating peptide counts swarmplot...")
    ax2 = fig.add_subplot(gs[0, 1])
    plot_peptide_counts_swarm(data_matrix, design_matrix, ax=ax2)

    # Bottom left: KDE plot of peptide length distribution
    print("Creating peptide length KDE plot...")
    ax3 = fig.add_subplot(gs[1, 0])
    plot_length_distribution_kdeplot(data_matrix, design_matrix, ax=ax3)

    # Bottom right: UMAP with group colors and day shapes
    print("Creating UMAP plot with group colors and day shapes...")
    ax4 = fig.add_subplot(gs[1, 1])
    plot_umap_with_day_shape(data_matrix, design_matrix, ax=ax4, legend_above=True)

    plt.tight_layout()

    output_file = os.path.join(output_dir, "figure1.svg")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")

    return fig


def plot_figure2(
    main_data: pd.DataFrame,
    main_design: pd.DataFrame,
    rerun_data: pd.DataFrame,
    rerun_design: pd.DataFrame,
    output_dir: str,
) -> plt.Figure:
    """
    Create the second composite figure as requested.
    """
    fig = plt.figure(figsize=(12, 10))

    gs = gridspec.GridSpec(
        2, 2, figure=fig, height_ratios=[1, 1.2], width_ratios=[1, 1.1]
    )

    # Top left: Venn diagrams for reruns vs corresponding main groups
    print("Creating rerun vs main Venn diagrams...")
    ax1 = fig.add_subplot(gs[0, 0])
    plot_rerun_venn_diagrams(main_data, main_design, rerun_data, rerun_design, ax=ax1)

    # Top right: Venn diagram showing overlap between rerun groups
    print("Creating rerun group overlap Venn diagram...")
    ax2 = fig.add_subplot(gs[0, 1])
    plot_rerun_group_overlap_venn(rerun_data, rerun_design, main_design, ax=ax2)

    bottom_gs = gridspec.GridSpecFromSubplotSpec(
        1,
        2,
        subplot_spec=gs[1, :],
        width_ratios=[1, 1.1],  # Give slightly more width to the UMAP
        wspace=0.3,  # Add more space between the plots
    )

    # Bottom left: Per-sample intensity correlation plots (smaller)
    print("Creating per-sample intensity correlation plots...")
    ax3 = fig.add_subplot(bottom_gs[0, 0])
    plot_rerun_intensity_per_sample(main_data, rerun_data, rerun_design, ax=ax3)

    # Bottom right: UMAP with rerun highlights and day shapes
    print("Creating UMAP with rerun highlights and day shapes...")
    ax4 = fig.add_subplot(bottom_gs[0, 1])

    rerun_samples = list(rerun_design.index)

    combined_data = pd.concat([main_data, rerun_data], axis=1)
    combined_data = combined_data.loc[:, ~combined_data.columns.duplicated()]

    plot_umap_with_day_shape(
        combined_data,
        main_design,
        highlight_samples=rerun_samples,
        ax=ax4,
        legend_above=True,
    )

    ax4.yaxis.set_label_coords(-0.2, 0.5)
    plt.tight_layout(pad=2.0)

    output_file = os.path.join(output_dir, "figure2.svg")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")

    return fig
