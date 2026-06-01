"""
Analysis Module
===============

Statistical analysis of plant gene expression data including:
- Principal Component Analysis (PCA) for dimensionality reduction
- Differential Expression Analysis using t-tests with FDR correction
- Gene Set Enrichment Analysis (simulated GO terms)
- Summary statistics and result tables
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from itertools import combinations


def run_pca(
    expression_df: pd.DataFrame,
    n_components: int = 5,
) -> dict:
    """
    Perform PCA on the expression matrix.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Normalized expression matrix (genes x samples).
    n_components : int
        Number of principal components to compute.

    Returns
    -------
    dict with keys:
        'loadings' : PC loadings (genes x components)
        'scores'   : PC scores (samples x components)
        'variance' : Explained variance ratio
    """
    # Scale data (genes as features, samples as observations → transpose)
    X = expression_df.T.values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=min(n_components, X.shape[0], X.shape[1]))
    scores = pca.fit_transform(X_scaled)

    # Loadings: contribution of each gene to each PC
    loadings = pd.DataFrame(
        pca.components_.T,
        index=expression_df.index,
        columns=[f"PC{i+1}" for i in range(pca.n_components_)],
    )

    scores_df = pd.DataFrame(
        scores,
        index=expression_df.columns,
        columns=[f"PC{i+1}" for i in range(pca.n_components_)],
    )

    variance_df = pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(pca.n_components_)],
        "Explained_Variance": pca.explained_variance_ratio_,
        "Cumulative_Variance": np.cumsum(pca.explained_variance_ratio_),
    })

    print(f"[Analysis] PCA complete: {pca.n_components_} components "
          f"explain {variance_df['Cumulative_Variance'].iloc[-1]:.1%} of variance")

    return {
        "loadings": loadings,
        "scores": scores_df,
        "variance": variance_df,
    }


def differential_expression(
    expression_df: pd.DataFrame,
    condition_a: str,
    condition_b: str,
    metadata_df: pd.DataFrame,
    log2fc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
) -> pd.DataFrame:
    """
    Perform differential expression analysis between two conditions.

    Uses Welch's t-test (unequal variance) with Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Normalized expression matrix (genes x samples).
    condition_a : str
        Name of the first condition (e.g., "Control").
    condition_b : str
        Name of the second condition (e.g., "Drought").
    metadata_df : pd.DataFrame
        Sample metadata with "Condition" column.
    log2fc_threshold : float
        Minimum absolute log2 fold-change for calling a gene significant.
    pval_threshold : float
        Adjusted p-value threshold for significance.

    Returns
    -------
    pd.DataFrame with columns:
        Gene, Mean_A, Mean_B, Log2FC, P_Value, Adj_P_Value, Significant, Regulation
    """
    samples_a = metadata_df[metadata_df["Condition"] == condition_a].index
    samples_b = metadata_df[metadata_df["Condition"] == condition_b].index

    # Filter to matching samples
    common_samples_a = [s for s in samples_a if s in expression_df.columns]
    common_samples_b = [s for s in samples_b if s in expression_df.columns]

    if len(common_samples_a) < 2 or len(common_samples_b) < 2:
        raise ValueError(
            f"Need at least 2 samples per condition. Got: "
            f"{condition_a}={len(common_samples_a)}, "
            f"{condition_b}={len(common_samples_b)}"
        )

    results = []
    for gene in expression_df.index:
        expr_a = expression_df.loc[gene, common_samples_a].values
        expr_b = expression_df.loc[gene, common_samples_b].values

        mean_a = expr_a.mean()
        mean_b = expr_b.mean()

        # Log2 fold-change (condition_b vs condition_a)
        log2fc = mean_b - mean_a

        # Welch's t-test
        t_stat, p_val = stats.ttest_ind(expr_b, expr_a, equal_var=False)

        results.append({
            "Gene": gene,
            f"Mean_{condition_a}": round(mean_a, 3),
            f"Mean_{condition_b}": round(mean_b, 3),
            "Log2FC": round(log2fc, 4),
            "P_Value": p_val,
        })

    result_df = pd.DataFrame(results)

    # Multiple testing correction (Benjamini-Hochberg)
    result_df["Adj_P_Value"] = _benjamini_hochberg(result_df["P_Value"].values)
    result_df["Significant"] = (
        (np.abs(result_df["Log2FC"]) >= log2fc_threshold)
        & (result_df["Adj_P_Value"] <= pval_threshold)
    )
    result_df["Regulation"] = "NS"
    result_df.loc[
        (result_df["Significant"]) & (result_df["Log2FC"] > 0), "Regulation"
    ] = "Up"
    result_df.loc[
        (result_df["Significant"]) & (result_df["Log2FC"] < 0), "Regulation"
    ] = "Down"

    # Sort by adjusted p-value
    result_df = result_df.sort_values("Adj_P_Value")

    n_up = (result_df["Regulation"] == "Up").sum()
    n_down = (result_df["Regulation"] == "Down").sum()
    print(f"[Analysis] Differential expression: {condition_b} vs {condition_a}")
    print(f"    {n_up} genes upregulated, {n_down} genes downregulated "
          f"(|Log2FC| >= {log2fc_threshold}, adj.P < {pval_threshold})")

    return result_df


def _benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(p_values)
    sorted_idx = np.argsort(p_values)
    sorted_p = p_values[sorted_idx]
    adjusted = np.minimum(1, sorted_p * n / (np.arange(1, n + 1)))
    # Ensure monotonicity
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])
    # Restore original order
    result = np.zeros(n)
    result[sorted_idx] = adjusted
    return result


# Simulated GO terms for enrichment analysis
SIMULATED_GO_TERMS = {
    "GO:0006950": "response to stress",
    "GO:0009409": "response to cold",
    "GO:0009414": "response to water deprivation",
    "GO:0009408": "response to heat",
    "GO:0006970": "response to osmotic stress",
    "GO:0009266": "response to temperature stimulus",
    "GO:0009651": "response to salt stress",
    "GO:0009737": "response to abscisic acid",
    "GO:0009415": "response to light stimulus",
    "GO:0015979": "photosynthesis",
    "GO:0006952": "defense response",
    "GO:0009607": "response to biotic stimulus",
    "GO:0009719": "response to endogenous stimulus",
    "GO:0010200": "response to chitin",
    "GO:0009628": "response to abiotic stimulus",
    "GO:0006355": "regulation of transcription, DNA-templated",
    "GO:0006412": "translation",
    "GO:0005975": "carbohydrate metabolic process",
    "GO:0006629": "lipid metabolic process",
    "GO:0006811": "ion transport",
}


def _assign_go_terms(gene_list: list) -> dict:
    """
    Simulate GO term assignments for enrichment testing.
    In a real project, this would use actual annotations from TAIR/Ensembl.
    """
    np.random.seed(42)
    term_map = {}
    go_ids = list(SIMULATED_GO_TERMS.keys())

    # Assign GO terms with realistic biases
    for gene in gene_list:
        # Each gene gets 1-5 random GO terms
        n_terms = np.random.randint(1, 6)
        selected = np.random.choice(go_ids, min(n_terms, len(go_ids)), replace=False)
        term_map[gene] = list(selected)

    # Enrich stress-related terms for known marker genes
    stress_go_ids = ["GO:0006950", "GO:0009409", "GO:0009414",
                     "GO:0009408", "GO:0009266", "GO:0009651"]
    for gene in gene_list:
        if gene in [
            "RD29A", "RD29B", "NCED3", "CBF1", "CBF2", "CBF3",
            "HSP70", "HSP90.1", "HSP101", "COR15A", "COR47",
        ]:
            term_map[gene] = list(set(term_map.get(gene, []) + stress_go_ids))

    return term_map


def go_enrichment(
    de_results: pd.DataFrame,
    gene_list: list,
    pval_threshold: float = 0.05,
) -> pd.DataFrame:
    """
    Perform GO term enrichment analysis on differentially expressed genes.

    Uses Fisher's exact test to identify overrepresented GO terms.

    Parameters
    ----------
    de_results : pd.DataFrame
        Results from differential_expression().
    gene_list : list
        Full list of genes in the background (for the universe).
    pval_threshold : float
        Adjusted p-value threshold for significance.

    Returns
    -------
    pd.DataFrame with enriched GO terms.
    """
    significant_genes = de_results[
        de_results["Significant"]
    ]["Gene"].tolist()

    if len(significant_genes) == 0:
        print("[Analysis] No significant genes found for GO enrichment")
        return pd.DataFrame()

    # Simulate GO term assignments
    go_assignments = _assign_go_terms(gene_list)

    # Build GO term -> gene mapping
    term_to_genes = {}
    for gene, terms in go_assignments.items():
        for term in terms:
            if term not in term_to_genes:
                term_to_genes[term] = []
            term_to_genes[term].append(gene)

    sig_set = set(significant_genes)
    bg_set = set(gene_list)

    results = []
    for go_id, go_genes in term_to_genes.items():
        go_set = set(go_genes)
        # Remove genes not in our dataset
        go_set = go_set & bg_set

        if len(go_set) < 3:
            continue

        # 2x2 contingency table
        a = len(sig_set & go_set)  # DE + in GO term
        b = len(go_set - sig_set)  # Not DE + in GO term
        c = len(sig_set - go_set)  # DE + not in GO term
        d = len(bg_set - go_set - sig_set)  # Neither

        # Fisher's exact test (right-tailed: overrepresentation)
        _, p_val = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

        results.append({
            "GO_ID": go_id,
            "GO_Term": SIMULATED_GO_TERMS.get(go_id, "unknown"),
            "Genes_in_Term": len(go_set),
            "DE_Genes_in_Term": a,
            "P_Value": p_val,
            "Representation": "Enriched" if a > 0 else "",
        })

    if not results:
        return pd.DataFrame()

    result_df = pd.DataFrame(results)
    result_df["Adj_P_Value"] = _benjamini_hochberg(result_df["P_Value"].values)
    result_df["Significant"] = result_df["Adj_P_Value"] <= pval_threshold
    result_df = result_df[result_df["DE_Genes_in_Term"] >= 2].sort_values("Adj_P_Value")

    n_sig = result_df["Significant"].sum()
    print(f"[Analysis] GO enrichment: {n_sig} significant terms found "
          f"(adj.P < {pval_threshold})")

    return result_df


def run_full_analysis(
    expression_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    control_condition: str = "Control",
    log2fc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
) -> dict:
    """
    Run the complete analysis pipeline: PCA + all pairwise DE + GO enrichment.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Normalized expression matrix.
    metadata_df : pd.DataFrame
        Sample metadata with "Condition" column.
    control_condition : str
        Name of the control/baseline condition.
    log2fc_threshold : float
        |Log2FC| threshold for calling DE genes.
    pval_threshold : float
        Adjusted p-value threshold.

    Returns
    -------
    dict with keys: 'pca', 'de_results', 'go_results'
    """
    print("\n" + "=" * 60)
    print("ANALYSIS PIPELINE")
    print("=" * 60)

    # PCA
    print("\n[1/3] Running PCA...")
    pca_results = run_pca(expression_df)

    # Differential expression (each stress condition vs control)
    print("\n[2/3] Running differential expression analyses...")
    conditions = [c for c in metadata_df["Condition"].unique()
                  if c != control_condition]
    de_results = {}
    for stress_cond in conditions:
        print(f"\n  >> {stress_cond} vs {control_condition}")
        de_results[f"{stress_cond}_vs_{control_condition}"] = differential_expression(
            expression_df, control_condition, stress_cond,
            metadata_df, log2fc_threshold, pval_threshold,
        )

    # GO enrichment
    print("\n[3/3] Running GO enrichment analysis...")
    go_results = {}
    for comp_name, de_df in de_results.items():
        print(f"\n  >> {comp_name}")
        go_results[comp_name] = go_enrichment(
            de_df, expression_df.index.tolist(), pval_threshold,
        )

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)

    return {
        "pca": pca_results,
        "de_results": de_results,
        "go_results": go_results,
    }
