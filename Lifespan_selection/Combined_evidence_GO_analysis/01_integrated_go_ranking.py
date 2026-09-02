#!/usr/bin/env python
# Author: Athanasios Toros
# Email: athanasiostoros@gmail.com

import os
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import yaml


def _to_namespace(obj):
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: _to_namespace(v) for k, v in obj.items()})
    return obj


def load_config(path=Path(__file__).parent / "config.yaml"):
    """Load config.yaml into a nested SimpleNamespace."""
    with open(path) as fh:
        raw = yaml.safe_load(fh)
    return _to_namespace(raw)


cfg = load_config()

# Must be set before importing polars.
os.environ["POLARS_MAX_THREADS"] = str(cfg.params.threads)
import polars as pl
from polars import col as c
import polars.selectors as cs

from goatools.obo_parser import GODag
from statsmodels.stats.multitest import multipletests

np.random.seed(cfg.params.seed)

# ---------------------------------------------------------------------------
# 1. Load inputs
# ---------------------------------------------------------------------------
incidence_df = pl.read_parquet(cfg.inputs.incidence_matrix)
go_annotations_df = pl.read_csv(cfg.inputs.go_annotations, separator="\t")
# I don't know why, but there are duplicate rows in this go file. Remove them.
go_annotations_df = go_annotations_df.unique()
godag = GODag(cfg.inputs.go_dag)

if not set(go_annotations_df.select("GID").unique()["GID"].to_list()).issubset(
    set(incidence_df.select("feature").unique()["feature"].to_list())
):
    raise ValueError("Some GO terms are missing in the incidence matrix.")


# ---------------------------------------------------------------------------
# 2. Merge CNV gains + deletions per regime
# ---------------------------------------------------------------------------
# Evidence layers with a small amount of "hits" lead to high variance in their
# scoring, which translates to a single evidence category dominating the final signal. CNV
# gains and deletions are merged into a single per-regime CNV category.
incidence_df = (
    incidence_df.with_columns(
        cnv_early=(c("cnv_early_deletions") + c("cnv_early_gains")).clip(upper_bound=1),
        cnv_late=(c("cnv_late_deletions") + c("cnv_late_gains")).clip(upper_bound=1),
    )
    .drop(cs.contains("gains"))
    .drop(cs.contains("deletions"))
)

# ---------------------------------------------------------------------------
# 3. Hypergeometric Z-scores per GO term x evidence layer
# ---------------------------------------------------------------------------
# X ~ Hypergeometric(N, K, n); N = filtered gene universe size, K = GO gene
# set size, n = hits per layer, x = observed hits per GO term.
# mu = n * K/N ; sigma^2 = n * K/N * (1 - K/N) * (N-n)/(N-1) ; Z = (x-mu)/sigma


def create_GO_hypergeometric_Zscores_matrix(
    incidence_df: pl.DataFrame,
    go_term_annotations_df: pl.DataFrame,
    column_universe_sizes: dict[str, int],  # N_c
    column_hits_nr: dict[str, int],  # n_c
) -> pl.DataFrame:
    # Join and drop original feature column
    merged_df = incidence_df.join(
        go_term_annotations_df,
        left_on="feature",
        right_on="GID",
        how="inner",
    ).drop("feature")

    # Evidence columns
    value_cols = [col for col in incidence_df.columns if col != "feature"]

    # Aggregate per GO term:
    # - x_tc = count of hits
    # - K_tc = size of term-specific gene set (respect structural zeros)
    initial_aggregation = merged_df.group_by("GO").agg(
        *( 
            [c(col).sum().alias(col) for col in value_cols]
            + [
                c(col).is_not_null().sum().alias(f"{col}_geneset_size")
                for col in value_cols
            ]
        )
    )

    # Compute:
    # - Hypergeometric expected value μ = n * (K/N)
    # - Variance σ² = μ * (1 - K/N) * (N - n)/(N - 1)
    intermediate_df = initial_aggregation.select(
        [
            c("GO"),
            *[c(col) for col in value_cols],
            *[
                (
                    column_hits_nr[col]
                    * (c(f"{col}_geneset_size") / column_universe_sizes[col])
                ).alias(f"{col}_hypergeom_expected")
                for col in value_cols
            ],
            *[
                (
                    column_hits_nr[col]
                    * (
                        c(f"{col}_geneset_size") / column_universe_sizes[col]
                    )  # mean term
                    * (
                        1 - (c(f"{col}_geneset_size") / column_universe_sizes[col])
                    )  # (1 - K/N)
                    * (
                        (column_universe_sizes[col] - column_hits_nr[col])
                        / (column_universe_sizes[col] - 1)
                    )  # finite pop correction
                ).alias(f"{col}_hypergeom_variance")
                for col in value_cols
            ],
        ]
    )

    # Convert counts to Z-scores
    # Z = (x - μ) / sqrt(σ²)
    # Return 0 when variance = 0
    zscore_df = intermediate_df.select(
        [
            c("GO"),
            *[
                pl.when(c(f"{col}_hypergeom_variance") > 0)
                .then(
                    (c(col) - c(f"{col}_hypergeom_expected"))
                    / c(f"{col}_hypergeom_variance").sqrt()
                )
                .otherwise(
                    pl.lit(0.0)
                )  # changed None to 0, to avoid null integrated scores. This is debatable.
                .alias(col)
                for col in value_cols
            ],
        ]
    )

    return zscore_df


column_universe_sizes = incidence_df.select(
    [
        (c(col).is_not_null().sum()).alias(col)
        for col in incidence_df.columns
        if col != "feature"
    ]
).to_dicts()[0]

column_hits_nr = incidence_df.select(
    [(c(col).sum()).alias(col) for col in incidence_df.columns if col != "feature"]
).to_dicts()[0]

observed_zscore_matrix = create_GO_hypergeometric_Zscores_matrix(
    incidence_df,
    go_annotations_df,
    column_universe_sizes,
    column_hits_nr,
)

# ---------------------------------------------------------------------------
# 4. p-norm integration within Selection / DEG categories, per regime
# ---------------------------------------------------------------------------
# Within each evidence category, Z-scores are aggregated with a normalized
# p-norm. Regime score = Selection +
# shared Expression (DEG).


def power_mean_expr(category_columns, power_factor: float):
    cols = list(category_columns)  # avoid passing a set directly

    def _helper_expr(col_name: str) -> pl.Expr:
        # Cap the negative values to zero, we only care about enrichment. 
        # Underrepresentation might be due to lack of power, bias
        return (
            pl.when(c(col_name) >= 0)
            .then(c(col_name) ** power_factor)
            .otherwise(pl.lit(0.0))
        )

    powered = [_helper_expr(col) for col in cols]
    # (mean of powered values) ** (1/p)
    return pl.mean_horizontal(powered) ** (1.0 / power_factor)


def compute_integrated_go_scores(
    go_count_matrix: pl.DataFrame,
    column_categories: dict,
    power_factor: float = 2.0,
    include_go_term: bool = False,
) -> pl.DataFrame:

    # adjust selection columns to what's actually present
    selection_columns = set(column_categories["selection"]).intersection(
        set(go_count_matrix.columns)
    )

    output_cols = [
        c("selection_score"),
        c("DEG_score"),
    ]

    if include_go_term:
        output_cols = [c("GO")] + output_cols

    integrated_scores_df = (
        go_count_matrix.select(
            pl.all(),
            power_mean_expr(selection_columns, power_factor).alias("selection_score"),
            power_mean_expr(column_categories["DEG"], power_factor).alias("DEG_score"),
        )
        .select(
            *output_cols,
        )
        .with_columns((c("selection_score") + c("DEG_score")).alias("integrated_score"))
    )

    return integrated_scores_df


column_categories_early = {
    "selection": cfg.column_categories.early.selection,
    "DEG": cfg.column_categories.early.DEG,
}
column_categories_late = {
    "selection": cfg.column_categories.late.selection,
    "DEG": cfg.column_categories.late.DEG,
}

# ---------------------------------------------------------------------------
# 5. Permutation null
# ---------------------------------------------------------------------------
# Gene-set labels (GID<->GO) are shuffled to break the gene-GO association,
# giving an empirical null for the integrated score per regime.


def create_permuted_zscore_matrix(
    incidence_df: pl.DataFrame,
    go_term_annotations_df: pl.DataFrame,
    column_universe_sizes: dict[str, int],  # N_c
    column_hits_nr: dict[str, int],  # n_c
):

    # Permute gene set labels
    permuted_go_terms = go_term_annotations_df.select(
        pl.col("GID").shuffle().alias("GID"),
        pl.col("GO"),
    )

    # Create GO count matrix based on permuted labels
    permuted_go_zscore_matrix = create_GO_hypergeometric_Zscores_matrix(
        incidence_df, permuted_go_terms, column_universe_sizes, column_hits_nr
    )
    return permuted_go_zscore_matrix


def compute_integrated_go_scores_for_permutations(
    go_zscore_matrix: pl.DataFrame,
    column_categories: dict,
    permutation_index: int,
    power_factor: float = 2.0,
) -> pl.DataFrame:

    integrated_scores_df = (
        go_zscore_matrix.select(
            pl.all(),
            power_mean_expr(column_categories["selection"], power_factor).alias(
                f"selection_score_{permutation_index}"
            ),
            power_mean_expr(column_categories["DEG"], power_factor).alias(
                f"DEG_score_{permutation_index}"
            ),
        )
        .select(
            c(f"selection_score_{permutation_index}"),
            c(f"DEG_score_{permutation_index}"),
        )
        .with_columns(
            (
                c(f"selection_score_{permutation_index}")
                + c(f"DEG_score_{permutation_index}")
            ).alias(f"integrated_score_{permutation_index}")
        )
    )
    return integrated_scores_df


def run_permutation_analysis(
    incidence_df: pl.DataFrame,
    go_term_annotations_df: pl.DataFrame,
    column_categories_early: dict,
    column_categories_late: dict,
    column_universe_sizes: dict[str, int],
    column_hits_nr: dict[str, int],
    n_permutations: int,
    power_factor: float = 2.0,
):
    early_scores = []
    late_scores = []

    base_go_terms = (
        create_GO_hypergeometric_Zscores_matrix(
            incidence_df, go_term_annotations_df, column_universe_sizes, column_hits_nr
        )
        .select("GO")
        .sort("GO")
        .to_series()
    )

    if n_permutations <= 0:
        raise ValueError("n_permutations must be a positive integer")

    t_start = time.time()
    # Report every 10%; ensure at least step 1
    step = max(n_permutations // 10, 1)
    next_report = step

    for i in range(n_permutations):
        randomized_go_count_matrix = create_permuted_zscore_matrix(
            incidence_df=incidence_df,
            go_term_annotations_df=go_term_annotations_df,
            column_universe_sizes=column_universe_sizes,
            column_hits_nr=column_hits_nr,
        ).sort(
            "GO"
        )  # ensure consistent GO term order

        integrated_score_early = compute_integrated_go_scores_for_permutations(
            randomized_go_count_matrix,
            column_categories_early,
            power_factor=power_factor,
            permutation_index=i,
        )
        early_scores.append(integrated_score_early)

        integrated_score_late = compute_integrated_go_scores_for_permutations(
            randomized_go_count_matrix,
            column_categories_late,
            power_factor=power_factor,
            permutation_index=i,
        )
        late_scores.append(integrated_score_late)

        done = i + 1
        if done >= next_report or done == n_permutations:
            elapsed = time.time() - t_start
            frac_done = done / n_permutations
            # avoid division by zero: done > 0 always here
            est_total = elapsed / frac_done
            eta = est_total - elapsed

            percent = int(round(frac_done * 100))
            print(
                f"{percent}% of permutations finished "
                f"({done}/{n_permutations}), "
                f"took {elapsed:.1f} s so far, "
                f"ETA: {eta:.1f} s"
            )

            next_report += step

    early_scores_df = pl.concat(early_scores, how="horizontal").with_columns(
        base_go_terms.alias("GO")
    )
    late_scores_df = pl.concat(late_scores, how="horizontal").with_columns(
        base_go_terms.alias("GO")
    )

    return early_scores_df, late_scores_df


# Small GO terms introduce noise, since the Z-score estimator becomes
# unstable. Large terms tend to be uninformative, leading to reduced power.
go_annotations_df_filtered = go_annotations_df.join(
    go_annotations_df.group_by("GO")
    .agg(pl.len().alias("gene_set_size"))
    .filter(
        (pl.col("gene_set_size") >= cfg.params.minimum_genes_in_go_term)
        & (pl.col("gene_set_size") <= cfg.params.maximum_genes_in_go_term)
    )
    .select("GO"),
    on="GO",
    how="inner",
).unique()

# Also filter for GO terms with no representation in our dataset (structural
# zeros in all), since they do not provide any valuable information.
go_annotations_to_keep = (
    incidence_df.join(
        go_annotations_df_filtered, left_on="feature", right_on="GID", how="inner"
    )
    .drop("feature")
    .group_by("GO")
    .agg(
        # Nulls are structural zeros here, so count non-nulls as K
        cs.all()
        .exclude("GO")
        .is_not_null()
        .sum()
        .name.prefix("K_")
    )
    # If there is any representation in any evidence type, keep the GO term
    .filter(pl.max_horizontal(cs.starts_with("K_")) > 0)
).select("GO")

go_annotations_df_final = go_annotations_df_filtered.join(
    go_annotations_to_keep, on="GO", how="inner"
)

print(f"Running {cfg.params.n_permutations} permutations...")
early_scores_permuted, late_scores_permuted = run_permutation_analysis(
    incidence_df=incidence_df,
    go_term_annotations_df=go_annotations_df_final,
    column_categories_early=column_categories_early,
    column_categories_late=column_categories_late,
    column_universe_sizes=column_universe_sizes,
    column_hits_nr=column_hits_nr,
    n_permutations=cfg.params.n_permutations,
    power_factor=cfg.params.power_factor,
)

permuted_scores_early_path = Path(
    str(cfg.outputs.permuted_scores_early).format(
        n_permutations=cfg.params.n_permutations, power_factor=cfg.params.power_factor
    )
)
permuted_scores_late_path = Path(
    str(cfg.outputs.permuted_scores_late).format(
        n_permutations=cfg.params.n_permutations, power_factor=cfg.params.power_factor
    )
)
permuted_scores_early_path.parent.mkdir(parents=True, exist_ok=True)
early_scores_permuted.write_parquet(permuted_scores_early_path)
late_scores_permuted.write_parquet(permuted_scores_late_path)

# ---------------------------------------------------------------------------
# 6. Empirical p-values
# ---------------------------------------------------------------------------


def calculate_integrated_scores_and_pvalues(
    permuted_scores_df: pl.DataFrame,
    observed_zscore_matrix: pl.DataFrame,
    go_term_annotations_df: pl.DataFrame,
    column_categories: dict,
    n_permutations: int,
    power_factor: float = 2.0,
    score_components: list[str] = ["selection", "DEG"],
) -> pl.DataFrame:
    """
    Calculates integrated scores for both permuted and observed data,
    and computes unadjusted empirical p-values.
    """

    # 1. Infer number of permutations from the columns
    # Assumes columns are named like "selection_score_0", "DEG_score_0", etc.
    # We look for the max index in columns starting with the first component.
    perm_cols = [
        c
        for c in permuted_scores_df.columns
        if c.startswith(f"{score_components[0]}_score_")
    ]

    if len(perm_cols) != n_permutations:
        raise ValueError(
            f"Number of permutation columns ({len(perm_cols)}) does not match n_permutations ({n_permutations})."
        )

    # 2. Compute Integrated Scores for Permutations
    integrated_exprs = [
        pl.sum_horizontal(
            [pl.col(f"{comp}_score_{i}") for comp in score_components]
        ).alias(f"integration_score_{i}")
        for i in range(n_permutations)
    ]

    permuted_integrated = permuted_scores_df.select([pl.col("GO")] + integrated_exprs)

    # 3. Compute Observed Integrated Scores
    # We assume 'compute_integrated_go_scores' is available in your scope
    observed_scores = compute_integrated_go_scores(
        observed_zscore_matrix,
        column_categories,
        power_factor=power_factor,
        include_go_term=True,
    )

    # Prepare observed columns for joining
    obs_select_cols = [pl.col("GO"), pl.col("integrated_score").alias("observed_score")]
    for comp in score_components:
        obs_select_cols.append(pl.col(f"{comp}_score").alias(f"observed_{comp}_score"))

    observed_subset = observed_scores.select(obs_select_cols)

    # 4. Join Permuted and Observed
    # We use inner join to ensure we only test terms present in both
    # Might lead to silent loss of some GO terms if not present in permutations
    scores_combined = permuted_integrated.join(observed_subset, on="GO", how="inner")

    # 5. Calculate Empirical P-values
    # (Count(Permuted >= Observed) + 1) / (N_perm + 1)
    # Adding 1 to numerator and denominator for finite sample correction
    perm_score_cols = [
        c for c in scores_combined.columns if c.startswith("integration_score_")
    ]

    p_value_expr = (
        (
            pl.sum_horizontal(
                [pl.col("observed_score") <= pl.col(c) for c in perm_score_cols]
            )
            + 1
        )
        / (n_permutations + 1)
    ).alias("p_value")

    result_df = scores_combined.select(
        pl.col("GO"),
        pl.col("observed_score"),
        *[pl.col(f"observed_{comp}_score") for comp in score_components],
        p_value_expr,
    )

    # 6. Add Neg Log10 P-value and Gene Set Size
    # Calculate gene set sizes from annotations
    gene_set_sizes = go_term_annotations_df.group_by("GO").agg(
        pl.len().alias("gene_set_size")
    )

    final_df = (
        result_df.with_columns((-pl.col("p_value").log(10)).alias("neg_log10_p_value"))
        .join(gene_set_sizes, on="GO", how="left")
        .join(observed_zscore_matrix, on="GO", how="left")
    )

    return final_df


pvals_early = calculate_integrated_scores_and_pvalues(
    permuted_scores_df=early_scores_permuted,
    observed_zscore_matrix=observed_zscore_matrix,
    go_term_annotations_df=go_annotations_df_final,
    column_categories=column_categories_early,
    n_permutations=cfg.params.n_permutations,
    power_factor=cfg.params.power_factor,
)
pvals_late = calculate_integrated_scores_and_pvalues(
    permuted_scores_df=late_scores_permuted,
    observed_zscore_matrix=observed_zscore_matrix,
    go_term_annotations_df=go_annotations_df_final,
    column_categories=column_categories_late,
    n_permutations=cfg.params.n_permutations,
    power_factor=cfg.params.power_factor,
)

# ---------------------------------------------------------------------------
# 7. Per-category BH FDR
# ---------------------------------------------------------------------------
# Results are split by ontological category and corrected for multiple
# testing separately, since BP terms are the primary interest.


def annotate_and_adjust_pvalues(
    pvals_df: pl.DataFrame,
    go_lookup_df: pl.DataFrame,
    alpha: float = 0.05,
    method: str = "fdr_bh"
) -> dict[str, pl.DataFrame]:
    """
    Adds GO metadata, splits data by GO category (BP, MF, CC), 
    and applies FDR correction per category.

    Returns:
        dict[str, pl.DataFrame]: Dictionary with keys 'BP', 'MF', 'CC' 

        containing the processed DataFrames.
    """

    # 1. Join with GO Lookup metadata (Descriptions and Categories)
    annotated_df = pvals_df.join(
        go_lookup_df, on="GO", how="left"
    ).filter(
        ~pl.col("GO_category").is_null()
    )

    # 2. Define categories to process
    categories = {
        "BP": "biological_process", # or just "BP" depending on your lookup values
        "MF": "molecular_function", # or "MF"
        "CC": "cellular_component"  # or "CC"
    }

    # Check if lookup uses short codes (BP) or full names
    # Based on your snippet, your lookup maps to "BP", "MF", "CC"
    unique_cats = annotated_df.select("GO_category").unique().to_series().to_list()
    results = {}

    for cat_code in ["BP", "MF", "CC"]:
        if cat_code not in unique_cats:
            continue

        # Filter for specific category
        cat_df = annotated_df.filter(pl.col("GO_category") == cat_code)

        if cat_df.height == 0:
            continue

        # 3. Apply Multiple Testing Correction
        pvals = cat_df.select("p_value").to_series().to_list()

        _, adjusted_pvals, _, _ = multipletests(
            pvals,
            alpha=alpha,
            method=method,
        )

        # Add adjusted p-value column
        cat_df_adjusted = cat_df.with_columns(
            pl.Series("adjusted_p_value", adjusted_pvals)
        )

        results[cat_code] = cat_df_adjusted

    return results


namespace_map = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC",
}

go_lookup = pl.DataFrame(
    {
        "GO": list(godag.keys()),
        "GO_description": [term.name for term in godag.values()],
        "GO_category": [
            namespace_map.get(term.namespace, None) for term in godag.values()
        ],
    }
)

annotated_pvals_early = annotate_and_adjust_pvalues(
    pvals_early,
    go_lookup,
    alpha=cfg.params.alpha,
    method="fdr_bh",
)
annotated_pvals_late = annotate_and_adjust_pvalues(
    pvals_late,
    go_lookup,
    alpha=cfg.params.alpha,
    method="fdr_bh",
)

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
results_dir = Path(cfg.outputs.results_dir)
results_dir.mkdir(parents=True, exist_ok=True)

for regime, annotated_pvals in [("early", annotated_pvals_early), ("late", annotated_pvals_late)]:
    for category, df in annotated_pvals.items():
        out_path = results_dir / cfg.outputs.results_filename_template.format(
            regime=regime,
            category=category,
            power_factor=cfg.params.power_factor,
            n_permutations=cfg.params.n_permutations,
        )
        df.write_csv(out_path, separator="\t")
