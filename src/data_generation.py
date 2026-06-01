"""
Data Generation Module
======================

Generates simulated Arabidopsis thaliana gene expression data under four
abiotic stress conditions: control, drought, heat, and cold.

The simulation is designed to produce realistic patterns:
- Most genes show modest expression changes (log2 fold-change ~0)
- A subset of stress-responsive genes show significant up/down regulation
- Known stress-marker genes are included with realistic response patterns
- Technical variability (noise) follows a realistic distribution
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple, Optional

# Seed for reproducibility
RANDOM_SEED = 42

# Number of genes to simulate
N_GENES = 5000

# Number of samples per condition
SAMPLES_PER_CONDITION = 4

# Conditions
CONDITIONS = ["Control", "Drought", "Heat", "Cold"]

# Known stress-marker genes with their expected log2 fold-changes
# (gene_name: {condition: log2FC})
STRESS_MARKERS = {
    # Drought-responsive markers
    "RD29A":    {"Drought": 5.2, "Heat": 1.1, "Cold": 3.8},
    "RD29B":    {"Drought": 4.8, "Heat": 0.8, "Cold": 3.5},
    "NCED3":    {"Drought": 4.1, "Heat": 0.5, "Cold": 2.0},
    "DREB2A":   {"Drought": 3.5, "Heat": 2.8, "Cold": 1.5},
    "LEA14":    {"Drought": 3.9, "Heat": 0.3, "Cold": 2.1},
    "P5CS1":    {"Drought": 3.2, "Heat": 0.6, "Cold": 0.9},

    # Heat-responsive markers
    "HSP70":    {"Drought": 1.2, "Heat": 6.0, "Cold": -0.5},
    "HSP90.1":  {"Drought": 0.8, "Heat": 5.5, "Cold": -0.3},
    "HSP17.6A": {"Drought": 0.5, "Heat": 7.2, "Cold": -0.8},
    "HSFA1D":   {"Drought": 1.5, "Heat": 3.5, "Cold": 0.2},
    "HSP101":   {"Drought": 0.9, "Heat": 4.8, "Cold": -0.1},
    "MBF1C":    {"Drought": 0.6, "Heat": 3.2, "Cold": 0.3},

    # Cold-responsive markers
    "CBF1":     {"Drought": 0.3, "Heat": -1.0, "Cold": 6.5},
    "CBF2":     {"Drought": 0.4, "Heat": -0.8, "Cold": 6.8},
    "CBF3":     {"Drought": 0.2, "Heat": -1.2, "Cold": 7.0},
    "COR15A":   {"Drought": 0.9, "Heat": -0.5, "Cold": 5.2},
    "COR47":    {"Drought": 1.1, "Heat": -0.3, "Cold": 4.5},
    "ICE1":     {"Drought": 0.1, "Heat": -0.6, "Cold": 2.8},

    # Multi-stress & regulatory
    "ABF3":     {"Drought": 3.0, "Heat": 1.5, "Cold": 1.2},
    "ABF4":     {"Drought": 2.8, "Heat": 1.2, "Cold": 1.5},
    "MYB96":    {"Drought": 3.3, "Heat": 0.4, "Cold": 0.7},
    "WRKY33":   {"Drought": 2.0, "Heat": 2.5, "Cold": 1.8},
    "ZAT10":    {"Drought": 2.5, "Heat": 3.0, "Cold": 2.2},

    # Photosynthesis-related (downregulated under stress)
    "RBCS1A":   {"Drought": -2.5, "Heat": -3.0, "Cold": -1.8},
    "RBCS2B":   {"Drought": -2.2, "Heat": -2.8, "Cold": -1.5},
    "LHCB1.1":  {"Drought": -2.0, "Heat": -2.5, "Cold": -1.2},
    "LHCB2.1":  {"Drought": -1.8, "Heat": -2.2, "Cold": -1.0},
}

# Additional random stress-responsive genes (non-marker)
N_RESPONSIVE_GENES = 200


def _generate_gene_names(random_state: np.random.RandomState) -> list:
    """Generate realistic Arabidopsis-style gene locus identifiers."""
    genes = [f"AT{random_state.randint(1, 6)}G{random_state.randint(10000, 99999)}"
             for _ in range(N_GENES - len(STRESS_MARKERS) - N_RESPONSIVE_GENES)]
    # Add stress-responsive gene names
    responsive = [f"AT{random_state.randint(1, 6)}G{random_state.randint(10000, 99999)}"
                  for _ in range(N_RESPONSIVE_GENES)]
    return list(STRESS_MARKERS.keys()) + responsive + genes


def _generate_expression_profiles(
    n_genes: int,
    n_samples: int,
    conditions: list,
    gene_names: list,
    random_state: np.random.RandomState,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate expression data with condition-specific effects.

    Parameters
    ----------
    n_genes : int
    n_samples : int
    conditions : list
    gene_names : list
        Pre-generated gene names for matching stress markers.
    random_state : np.random.RandomState

    Returns
    -------
    expression_matrix : np.ndarray, shape (n_genes, n_samples)
    condition_labels : np.ndarray, shape (n_samples,)
    """
    samples_per_cond = n_samples // len(conditions)
    true_means = np.zeros((n_genes, len(conditions)))

    # Basal expression: log2(TPM+1) ~ Normal(4, 2)
    basal = random_state.normal(4, 2, n_genes)
    true_means[:, 0] = basal  # Control

    # Generate condition effects
    for j, cond in enumerate(conditions[1:], start=1):
        effects = random_state.normal(0, 0.3, n_genes)
        # Assign some genes to be stress-responsive in this condition
        n_resp = max(10, int(n_genes * 0.02))
        resp_idx = random_state.choice(n_genes, n_resp, replace=False)
        # Most responsive genes are mildly affected
        effects[resp_idx] = random_state.normal(0, 1.2, n_resp)
        # A few are strongly affected
        strong_idx = random_state.choice(resp_idx, max(3, n_resp // 10),
                                        replace=False)
        effects[strong_idx] = random_state.choice([-4, -3, -2, 2, 3, 4],
                                                  len(strong_idx))
        true_means[:, j] = basal + effects

    # Override with known stress markers using provided gene names
    for i, gene in enumerate(gene_names):
        if gene in STRESS_MARKERS:
            for j, cond in enumerate(conditions):
                if cond in STRESS_MARKERS[gene]:
                    true_means[i, j] = basal[i] + STRESS_MARKERS[gene][cond]

    # Generate samples with technical noise
    expression = np.zeros((n_genes, n_samples))
    condition_labels = []
    sample_names = []

    for j, cond in enumerate(conditions):
        for s in range(samples_per_cond):
            col = j * samples_per_cond + s
            # Technical noise: Negative Binomial-like (overdispersed)
            noise = random_state.normal(0, 0.4, n_genes)
            # Add some outlier noise to a few samples
            if s == 0:  # First sample of each condition has extra noise
                noise += random_state.normal(0, 0.2, n_genes)
            expression[:, col] = true_means[:, j] + noise
            condition_labels.append(cond)
            sample_names.append(f"{cond[:3]}_{s+1:02d}")

    return expression, np.array(condition_labels), sample_names


def generate_dataset(
    output_dir: Optional[str] = None,
    n_genes: int = N_GENES,
    samples_per_condition: int = SAMPLES_PER_CONDITION,
    seed: int = RANDOM_SEED,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate a complete simulated plant gene expression dataset.

    Parameters
    ----------
    output_dir : str or Path, optional
        Directory to save CSV files. If None, only returns DataFrames.
    n_genes : int
        Number of genes to simulate.
    samples_per_condition : int
        Number of biological replicates per condition.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    expression_df : pd.DataFrame
        Expression matrix (genes x samples) with log2(TPM+1) values.
    metadata_df : pd.DataFrame
        Sample metadata with condition labels.
    """
    np.random.seed(seed)
    random_state = np.random.RandomState(seed)

    total_samples = len(CONDITIONS) * samples_per_condition
    gene_names = _generate_gene_names(random_state)

    expression_matrix, conditions_arr, sample_names = _generate_expression_profiles(
        n_genes, total_samples, CONDITIONS, gene_names, random_state
    )

    # Ensure unique gene names
    seen = set()
    unique_names = []
    for g in gene_names:
        if g in seen:
            i = 1
            while f"{g}_{i}" in seen:
                i += 1
            unique_names.append(f"{g}_{i}")
        else:
            unique_names.append(g)
        seen.add(unique_names[-1])
    gene_names = unique_names

    # Build DataFrames
    expression_df = pd.DataFrame(
        expression_matrix,
        index=gene_names,
        columns=sample_names,
    )
    expression_df.index.name = "Gene"

    metadata_df = pd.DataFrame({
        "Sample": sample_names,
        "Condition": conditions_arr,
        "Replicate": [f"R{(i % samples_per_condition) + 1}"
                      for i in range(total_samples)],
    })
    metadata_df = metadata_df.set_index("Sample")

    # Save if output_dir provided
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        expression_df.to_csv(output_path / "expression_data.csv")
        metadata_df.to_csv(output_path / "metadata.csv")
        print(f"[✓] Data saved to {output_path.resolve()}")
        print(f"    - {len(gene_names)} genes x {total_samples} samples")
        print(f"    - Conditions: {', '.join(CONDITIONS)}")

    return expression_df, metadata_df


if __name__ == "__main__":
    # When run directly, generate sample data
    data_dir = Path(__file__).resolve().parent.parent / "data"
    expr, meta = generate_dataset(output_dir=str(data_dir))
    print("\nExpression data shape:", expr.shape)
    print("\nSample expression values (first 5 genes):")
    print(expr.head())
    print("\nSample metadata:")
    print(meta)
