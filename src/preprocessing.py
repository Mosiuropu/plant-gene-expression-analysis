"""
Preprocessing Module
====================

Quality control, filtering, and normalization of gene expression data.

Steps:
1. Filter lowly-expressed genes
2. Normalize using quantile normalization
3. Detect potential outlier samples
"""

import numpy as np
import pandas as pd


def filter_low_expression(
    expression_df: pd.DataFrame,
    min_mean: float = 1.0,
    min_samples_pct: float = 0.5,
) -> pd.DataFrame:
    """
    Remove genes with low expression across samples.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Expression matrix (genes x samples).
    min_mean : float
        Minimum mean expression value to retain a gene.
    min_samples_pct : float
        Minimum fraction of samples that must have expression > min_mean.

    Returns
    -------
    Filtered expression DataFrame.
    """
    n_min = max(1, int(expression_df.shape[1] * min_samples_pct))
    expressed_in_samples = (expression_df > min_mean).sum(axis=1)
    mask = expressed_in_samples >= n_min

    n_removed = (~mask).sum()
    if n_removed > 0:
        print(f"[Preprocessing] Removed {n_removed} lowly-expressed genes "
              f"(mean < {min_mean} in >= {1-min_samples_pct:.0%} of samples)")

    return expression_df.loc[mask].copy()


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Quantile normalize expression values across samples.

    This ensures the distribution of expression values is the same
    across all samples, which is critical for comparing expression
    levels between different conditions.

    Parameters
    ----------
    df : pd.DataFrame
        Expression matrix (genes x samples).

    Returns
    -------
    Quantile-normalized DataFrame.
    """
    rank_mean = df.stack().groupby(
        df.rank(method="first").stack().astype(int)
    ).mean()

    normalized = df.rank(method="min").stack().astype(int).map(rank_mean).unstack()
    normalized.index = df.index
    normalized.columns = df.columns

    print(f"[Preprocessing] Applied quantile normalization "
          f"({df.shape[0]} genes x {df.shape[1]} samples)")
    return normalized


def detect_outlier_samples(
    expression_df: pd.DataFrame,
    z_score_threshold: float = 3.0,
) -> tuple:
    """
    Detect potential outlier samples using median absolute deviation (MAD).

    Parameters
    ----------
    expression_df : pd.DataFrame
        Expression matrix (genes x samples).
    z_score_threshold : float
        Z-score threshold for flagging an outlier.

    Returns
    -------
    outlier_mask : pd.Series
        Boolean mask indicating outlier samples (True = outlier).
    outlier_stats : pd.DataFrame
        Summary statistics for each sample.
    """
    sample_means = expression_df.mean()
    sample_mads = expression_df.sub(expression_df.median(axis=0), axis=1).abs().median(axis=0)

    # Compute robust Z-scores
    median_mean = sample_means.median()
    mad_mean = np.median(np.abs(sample_means - median_mean))
    if mad_mean == 0:
        mad_mean = sample_means.std()

    z_scores = (sample_means - median_mean) / (mad_mean * 1.4826)
    outlier_mask = np.abs(z_scores) > z_score_threshold

    outlier_stats = pd.DataFrame({
        "Mean_Expression": sample_means,
        "MAD": sample_mads,
        "Z_Score": z_scores,
        "Outlier": outlier_mask,
    })

    n_outliers = outlier_mask.sum()
    if n_outliers > 0:
        print(f"[Preprocessing] Detected {n_outliers} potential outlier sample(s):")
        for sample in outlier_mask[outlier_mask].index:
            print(f"    - {sample} (Z-score: {z_scores[sample]:.2f})")
    else:
        print("[Preprocessing] No outlier samples detected")

    return outlier_mask, outlier_stats


def run_preprocessing(
    expression_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    min_mean: float = 1.0,
    min_samples_pct: float = 0.5,
    normalize: bool = True,
) -> dict:
    """
    Run the full preprocessing pipeline.

    Parameters
    ----------
    expression_df : pd.DataFrame
        Raw expression matrix.
    metadata_df : pd.DataFrame
        Sample metadata.
    min_mean : float
        Minimum mean expression filter.
    min_samples_pct : float
        Minimum fraction of samples with expression > min_mean.
    normalize : bool
        Whether to apply quantile normalization.

    Returns
    -------
    dict with keys: 'expression', 'metadata', 'outlier_mask', 'outlier_stats'
    """
    print("=" * 60)
    print("PREPROCESSING PIPELINE")
    print("=" * 60)

    # Step 1: Filter lowly expressed genes
    print("\n[1/3] Filtering lowly-expressed genes...")
    filtered = filter_low_expression(expression_df, min_mean, min_samples_pct)

    # Step 2: Normalize
    print(f"\n[2/3] Normalizing expression values...")
    if normalize:
        normalized = quantile_normalize(filtered)
    else:
        normalized = filtered
        print("[Preprocessing] Skipping quantile normalization")

    # Step 3: Outlier detection
    print(f"\n[3/3] Detecting outlier samples...")
    outlier_mask, outlier_stats = detect_outlier_samples(normalized)

    print("\n" + "=" * 60)
    print(f"Preprocessing complete!")
    print(f"  Genes retained: {expression_df.shape[0]} -> {normalized.shape[0]}")
    print(f"  Samples: {normalized.shape[1]}")
    print("=" * 60)

    return {
        "expression": normalized,
        "metadata": metadata_df.copy(),
        "outlier_mask": outlier_mask,
        "outlier_stats": outlier_stats,
    }
