"""
Visualization Module
====================

Publication-quality figures for plant gene expression analysis:
- PCA scatter plots colored by condition
- Volcano plots for differential expression
- Heatmaps of differentially expressed genes
- Expression boxplots for marker genes
- GO enrichment bar plots
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for script usage
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from pathlib import Path

# Global style settings
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
    "figure.figsize": (8, 6),
})

# Color palette for conditions
CONDITION_COLORS = {
    "Control": "#2E86AB",
    "Drought": "#A23B72",
    "Heat": "#F18F01",
    "Cold": "#1B998B",
}

# Custom diverging colormap
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "diverging_custom",
    ["#2166AC", "#F7F7F7", "#B2182B"],
    N=256,
)


def plot_pca(
    pca_scores: pd.DataFrame,
    metadata_df: pd.DataFrame,
    variance: pd.DataFrame,
    output_path: str,
) -> str:
    """
    Create a PCA scatter plot with PC1 vs PC2, colored by condition.

    Parameters
    ----------
    pca_scores : pd.DataFrame
        PCA scores (samples x components).
    metadata_df : pd.DataFrame
        Sample metadata with "Condition" column.
    variance : pd.DataFrame
        Explained variance ratios.
    output_path : str
        Path to save the figure.

    Returns
    -------
    Path to the saved figure.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Merge scores with metadata
    plot_df = pca_scores.copy()
    plot_df["Condition"] = metadata_df.loc[
        plot_df.index.intersection(metadata_df.index), "Condition"
    ]

    var_pc1 = variance.loc[variance["PC"] == "PC1", "Explained_Variance"].values[0]
    var_pc2 = variance.loc[variance["PC"] == "PC2", "Explained_Variance"].values[0]

    # Scatter plot
    for condition in CONDITION_COLORS:
        subset = plot_df[plot_df["Condition"] == condition]
        ax.scatter(
            subset["PC1"], subset["PC2"],
            c=CONDITION_COLORS[condition],
            label=condition,
            s=120,
            edgecolors="white",
            linewidth=0.8,
            alpha=0.9,
            zorder=3,
        )

        # Label each point
        for idx, row in subset.iterrows():
            ax.annotate(
                idx,
                (row["PC1"], row["PC2"]),
                fontsize=7,
                ha="center",
                va="bottom",
                xytext=(0, 7),
                textcoords="offset points",
            )

    ax.set_xlabel(f"PC1 ({var_pc1:.1%} variance)", fontweight="bold")
    ax.set_ylabel(f"PC2 ({var_pc2:.1%} variance)", fontweight="bold")
    ax.set_title("Principal Component Analysis of Gene Expression", fontweight="bold")
    ax.legend(frameon=True, fancybox=True, shadow=True)
    ax.axhline(y=0, color="grey", linestyle="--", alpha=0.3)
    ax.axvline(x=0, color="grey", linestyle="--", alpha=0.3)
    ax.grid(True, alpha=0.15)

    sns.despine()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

    print(f"[Viz] PCA plot saved to: {output_path}")
    return output_path


def plot_volcano(
    de_results: pd.DataFrame,
    comparison_name: str,
    output_path: str,
    log2fc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
    n_label: int = 15,
) -> str:
    """
    Create a volcano plot for differential expression results.

    Parameters
    ----------
    de_results : pd.DataFrame
        Results from differential_expression().
    comparison_name : str
        Name of the comparison (e.g., "Drought_vs_Control").
    output_path : str
        Path to save the figure.
    log2fc_threshold : float
        Log2 fold-change threshold line.
    pval_threshold : float
        Adjusted p-value threshold line.
    n_label : int
        Number of top DE genes to label.

    Returns
    -------
    Path to the saved figure.
    """
    df = de_results.copy()
    df["-log10(p_adj)"] = -np.log10(df["Adj_P_Value"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(9, 7))

    # Color by regulation status
    colors = {"Up": "#B2182B", "Down": "#2166AC", "NS": "#888888"}
    for reg_status in ["Up", "Down", "NS"]:
        subset = df[df["Regulation"] == reg_status]
        if len(subset) > 0:
            ax.scatter(
                subset["Log2FC"],
                subset["-log10(p_adj)"],
                c=colors[reg_status],
                label=f"{reg_status}regulated ({len(subset)} genes)"
                      if reg_status != "NS" else f"NS ({len(subset)} genes)",
                s=8,
                alpha=0.6 if reg_status == "NS" else 0.9,
                rasterized=True,
            )

    # Threshold lines
    ax.axhline(
        -np.log10(pval_threshold), color="grey", linestyle="--", alpha=0.5,
        linewidth=0.8,
    )
    ax.axvline(log2fc_threshold, color="grey", linestyle="--", alpha=0.5, linewidth=0.8)
    ax.axvline(-log2fc_threshold, color="grey", linestyle="--", alpha=0.5, linewidth=0.8)

    # Label top genes
    top_genes = df[df["Significant"]].head(n_label)
    for _, row in top_genes.iterrows():
        ax.annotate(
            row["Gene"],
            (row["Log2FC"], row["-log10(p_adj)"]),
            fontsize=6,
            ha="center",
            va="bottom",
            xytext=(0, 4),
            textcoords="offset points",
            arrowprops=dict(arrowstyle="-", color="grey", alpha=0.3),
        )

    # Clean up comparison name for title
    display_name = comparison_name.replace("_vs_", " vs ")
    ax.set_xlabel("Log₂ Fold Change", fontweight="bold")
    ax.set_ylabel("-Log₁₀ Adjusted P-value", fontweight="bold")
    ax.set_title(f"Volcano Plot: {display_name}", fontweight="bold")
    ax.legend(frameon=True, fancybox=True, shadow=True, markerscale=3)
    ax.grid(True, alpha=0.15)

    sns.despine()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

    print(f"[Viz] Volcano plot saved to: {output_path}")
    return output_path


def plot_heatmap(
    expression_df: pd.DataFrame,
    de_results: pd.DataFrame,
    metadata_df: pd.DataFrame,
    comparison_name: str,
    output_path: str,
    max_genes: int = 50,
) -> str:
    """
    Create a clustered heatmap of the top differentially expressed genes.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Normalized expression matrix.
    de_results : pd.DataFrame
        Differential expression results.
    metadata_df : pd.DataFrame
        Sample metadata with "Condition" column.
    comparison_name : str
        Name of the comparison.
    output_path : str
        Path to save the figure.
    max_genes : int
        Maximum number of genes to show.

    Returns
    -------
    Path to the saved figure.
    """
    # Get top DE genes (by adjusted p-value)
    top_genes = de_results[
        de_results["Significant"]
    ].head(max_genes)

    if len(top_genes) == 0:
        print(f"[Viz] No significant genes to plot for {comparison_name}")
        return output_path

    # Get expression for these genes
    plot_data = expression_df.loc[
        [g for g in top_genes["Gene"] if g in expression_df.index]
    ]

    if plot_data.empty:
        return output_path

    # Sort samples by condition
    sample_order = metadata_df.sort_values("Condition").index
    sample_order = [s for s in sample_order if s in plot_data.columns]
    plot_data = plot_data[sample_order]

    # Create color annotation for conditions
    condition_colors = [CONDITION_COLORS[metadata_df.loc[s, "Condition"]]
                        for s in plot_data.columns]

    # Z-score normalize rows
    plot_data_z = plot_data.subtract(plot_data.mean(axis=1), axis=0)
    plot_data_z = plot_data_z.div(plot_data.std(axis=1), axis=0)

    # Determine figure height based on number of genes
    n_genes = plot_data_z.shape[0]
    fig_height = max(4, min(16, n_genes * 0.3 + 2))

    fig_width = max(6, len(sample_order) * 0.4 + 3)
    g = sns.clustermap(
        plot_data_z,
        col_colors=[condition_colors],
        cmap=DIVERGING_CMAP,
        vmin=-2.5,
        vmax=2.5,
        figsize=(fig_width, fig_height),
        row_cluster=True,
        col_cluster=True,
        dendrogram_ratio=(0.08, 0.08),
        cbar_pos=(0.02, 0.85, 0.03, 0.15),
        xticklabels=True,
        yticklabels=True if n_genes <= 40 else False,
        linewidths=0,
        method="ward",
        metric="euclidean",
    )

    # Add condition legend
    display_name = comparison_name.replace("_vs_", " vs ")
    g.ax_heatmap.set_ylabel("Genes")
    g.ax_heatmap.set_xlabel("Samples")
    g.fig.suptitle(
        f"Expression Heatmap: Differentially Expressed Genes\n{display_name}",
        fontweight="bold", y=1.02,
    )

    # Create legend for conditions
    handles = [mpatches.Patch(color=color, label=cond)
               for cond, color in CONDITION_COLORS.items()
               if cond in metadata_df["Condition"].values]

    g.ax_col_dendrogram.legend(
        handles=handles, loc="upper right", frameon=True,
        fancybox=True, shadow=True, ncol=2, fontsize=8,
    )

    plt.savefig(output_path, dpi=200)
    plt.close()

    print(f"[Viz] Heatmap saved to: {output_path}")
    return output_path


def plot_marker_expression(
    expression_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    marker_genes: list,
    output_path: str,
) -> str:
    """
    Create boxplots showing expression of known marker genes across conditions.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Normalized expression matrix.
    metadata_df : pd.DataFrame
        Sample metadata with "Condition" column.
    marker_genes : list
        List of marker gene names to plot.
    output_path : str
        Path to save the figure.

    Returns
    -------
    Path to the saved figure.
    """
    available = [g for g in marker_genes if g in expression_df.index]
    if not available:
        print("[Viz] No marker genes found in expression data")
        return output_path

    n_genes = len(available)
    n_cols = min(4, n_genes)
    n_rows = (n_genes + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
    axes = axes.flatten() if n_genes > 1 else [axes]

    for i, gene in enumerate(available):
        ax = axes[i]
        gene_expr = expression_df.loc[gene]
        plot_data = gene_expr.to_frame("Expression")
        plot_data["Condition"] = metadata_df.loc[plot_data.index, "Condition"]

        sns.boxplot(
            data=plot_data, x="Condition", y="Expression",
            palette=CONDITION_COLORS, ax=ax, width=0.5,
            order=["Control", "Drought", "Heat", "Cold"],
        )
        sns.stripplot(
            data=plot_data, x="Condition", y="Expression",
            palette=CONDITION_COLORS, ax=ax,
            size=8, edgecolor="black", linewidth=0.5, jitter=True,
            order=["Control", "Drought", "Heat", "Cold"],
        )

        ax.set_title(gene, fontweight="bold", fontsize=11)
        ax.set_xlabel("")
        ax.set_ylabel("Expression (log₂ TPM+1)")
        ax.grid(True, alpha=0.15, axis="y")

    # Hide unused axes
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        "Expression of Known Stress Marker Genes Across Conditions",
        fontweight="bold", fontsize=13, y=1.02,
    )
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

    print(f"[Viz] Marker gene expression plot saved to: {output_path}")
    return output_path


def plot_go_enrichment(
    go_results: pd.DataFrame,
    comparison_name: str,
    output_path: str,
    max_terms: int = 15,
) -> str:
    """
    Create a bar plot of the top enriched GO terms.

    Parameters
    ----------
    go_results : pd.DataFrame
        Results from go_enrichment().
    comparison_name : str
        Name of the comparison.
    output_path : str
        Path to save the figure.
    max_terms : int
        Maximum number of terms to display.

    Returns
    -------
    Path to the saved figure.
    """
    if go_results.empty:
        print(f"[Viz] No GO results to plot for {comparison_name}")
        return output_path

    top_terms = go_results[go_results["Significant"]].head(max_terms)

    if top_terms.empty:
        # Plot non-significant top terms instead
        top_terms = go_results.head(max_terms)
        title_suffix = " (top terms)"
    else:
        title_suffix = ""

    fig, ax = plt.subplots(figsize=(8, max(4, len(top_terms) * 0.35)))

    y_pos = range(len(top_terms))
    colors = plt.cm.RdYlBu_r(
        np.linspace(0.3, 0.8, len(top_terms))
    )

    bars = ax.barh(y_pos, -np.log10(top_terms["Adj_P_Value"].values.clip(1e-15)),
                   color=colors, edgecolor="white", linewidth=0.5)

    # Add number of DE genes per term
    for i, (_, row) in enumerate(top_terms.iterrows()):
        adj_p_val = max(row["Adj_P_Value"], 1e-15)
        ax.text(
            -np.log10(adj_p_val) + 0.1,
            i,
            f"({int(row['DE_Genes_in_Term'])})",
            va="center",
            fontsize=8,
            color="grey",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_terms["GO_Term"].str.replace("_", " "), fontsize=9)
    ax.set_xlabel("-Log₁₀ Adjusted P-value", fontweight="bold")
    display_name = comparison_name.replace("_vs_", " vs ").replace("_", " ")
    ax.set_title(
        f"GO Enrichment Analysis{title_suffix}\n{display_name}",
        fontweight="bold",
    )
    ax.axvline(-np.log10(0.05), color="red", linestyle="--", alpha=0.5,
               linewidth=0.8, label="p = 0.05")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.15, axis="x")

    sns.despine(left=False, bottom=False)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

    print(f"[Viz] GO enrichment plot saved to: {output_path}")
    return output_path


def generate_all_figures(
    expression_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    analysis_results: dict,
    output_dir: str,
) -> dict:
    """
    Generate all figures for the analysis report.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Normalized expression matrix.
    metadata_df : pd.DataFrame
        Sample metadata.
    analysis_results : dict
        Results from analysis.run_full_analysis().
    output_dir : str
        Directory to save figures.

    Returns
    -------
    dict mapping figure names to file paths.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    figures = {}

    print("\n" + "=" * 60)
    print("GENERATING FIGURES")
    print("=" * 60)

    # 1. PCA plot
    print("\n[1/5] PCA plot...")
    figures["pca"] = plot_pca(
        analysis_results["pca"]["scores"],
        metadata_df,
        analysis_results["pca"]["variance"],
        str(output_path / "01_pca_plot.png"),
    )

    # 2. Volcano plots for each comparison
    print("\n[2/5] Volcano plots...")
    for comp_name, de_df in analysis_results["de_results"].items():
        figures[f"volcano_{comp_name}"] = plot_volcano(
            de_df, comp_name,
            str(output_path / f"02_volcano_{comp_name}.png"),
        )

    # 3. Heatmaps
    print("\n[3/5] Expression heatmaps...")
    for comp_name, de_df in analysis_results["de_results"].items():
        figures[f"heatmap_{comp_name}"] = plot_heatmap(
            expression_df, de_df, metadata_df, comp_name,
            str(output_path / f"03_heatmap_{comp_name}.png"),
        )

    # 4. Marker gene expression
    print("\n[4/5] Marker gene expression...")
    marker_genes = [
        "RD29A", "HSP70", "CBF3", "NCED3",
        "COR15A", "HSP101", "WRKY33", "RBCS1A",
    ]
    figures["markers"] = plot_marker_expression(
        expression_df, metadata_df, marker_genes,
        str(output_path / "04_marker_expression.png"),
    )

    # 5. GO enrichment plots
    print("\n[5/5] GO enrichment plots...")
    for comp_name, go_df in analysis_results["go_results"].items():
        figures[f"go_{comp_name}"] = plot_go_enrichment(
            go_df, comp_name,
            str(output_path / f"05_go_enrichment_{comp_name}.png"),
        )

    print("\n" + "=" * 60)
    print(f"All figures saved to: {output_path.resolve()}")
    print("=" * 60)

    return figures
