#!/usr/bin/env python
"""
Author: Athanasios Toros
Email: athanasiostoros@gmail.com
"""

import multiprocessing as mp
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import polars as pl
import pywt
import yaml
from pybedtools import BedTool as Bed
from skmisc.loess import loess
from tqdm import tqdm


def load_config(path: Path = Path(__file__).parent / "config.yaml") -> SimpleNamespace:
    with open(path) as fh:
        raw = yaml.safe_load(fh)

    cfg = SimpleNamespace()
    cfg.allowed_chromosomes = raw["genome"]["allowed_chromosomes"]
    cfg.chromosome_x = raw["genome"]["chromosome_x"]

    paths = raw["pipeline"]["paths"]
    cfg.coverage_dir = Path(paths["coverage_dir"])
    cfg.working_dir = Path(paths["working_dir"])
    cfg.gc_content_file = Path(paths["gc_content_file"])
    cfg.mappability_file = Path(paths["mappability_file"])
    cfg.exclusion_low_coverage_regions = Path(paths["exclusion_low_coverage_regions"])
    cfg.exclusion_high_coverage_regions = Path(paths["exclusion_high_coverage_regions"])
    cfg.segmentation_indices_file = Path(paths["segmentation_indices_file"])

    params = raw["pipeline"]["params"]
    cfg.gc_lower_quantile = params["gc_lower_quantile"]
    cfg.gc_upper_quantile = params["gc_upper_quantile"]
    cfg.mappability_lower_bound = params["mappability_lower_bound"]
    cfg.median_coverage_threshold = params["median_coverage_threshold"]
    cfg.alpha_threshold = params["alpha_threshold"]
    cfg.fold_change_threshold = params["fold_change_threshold"]
    cfg.modwt_level = params["modwt_level"]
    cfg.modwt_wavelet = params["modwt_wavelet"]
    cfg.save_results = params["save_results"]
    return cfg


Config = load_config()

Working_100bp = SimpleNamespace()
Metadata = SimpleNamespace()
Filtering = SimpleNamespace()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def get_sample_name(file_path: Path) -> str:
    return file_path.stem.split("_")[0]


def load_coverage_sample_matrix(
    coverage_files: list[Path],
    chromosomes: list[str],
) -> pl.DataFrame:
    """
    Load the coverage files and return a DataFrame with the coverage values.
    """
    # Create an empty DataFrame to store the coverage values
    combined_df = None

    for coverage_file in coverage_files:
        # Extract the sample name from the file name
        sample_name = get_sample_name(coverage_file)

        # Read the coverage file into a DataFrame
        tmp_df = pl.scan_csv(
            coverage_file,
            separator="\t",
            has_header=False,
            new_columns=[
                "chrom",
                "start",
                "end",
                "coverage_" + sample_name,
            ],
        ).filter(pl.col("chrom").is_in(chromosomes))

        if combined_df is None:
            combined_df = tmp_df
        else:
            combined_df = combined_df.join(
                tmp_df,
                on=["chrom", "start", "end"],
                how="inner",
                suffix=sample_name,
            )

    return combined_df.collect()


Working_100bp.coverage_files = Config.coverage_dir.glob("*regions.bed.gz")
Working_100bp.all_samples_coverages = load_coverage_sample_matrix(
    Working_100bp.coverage_files,
    Config.allowed_chromosomes,
)


# ---------------------------------------------------------------------------
# Filtering: aberrant-coverage blacklist
# ---------------------------------------------------------------------------
def extract_bed_file_from_df(
    df: pl.DataFrame,
    output_file=None,
) -> Bed:
    """
    Extract a bed file from the DataFrame.
    """

    extracted_bed_data = df.select(
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
    )

    if output_file is None:
        output_file = tempfile.NamedTemporaryFile(
            suffix=".bed", delete=False
        ).name
        # Write the DataFrame to a bed file
        if Config.save_results == True:
            extracted_bed_data.write_csv(
                output_file,
                separator="\t",
                include_header=False,
            )
    else:
        # Write the DataFrame to a bed file
        if Config.save_results == True:
            extracted_bed_data.write_csv(
                output_file,
                separator="\t",
                include_header=False,
            )

    return Bed(output_file)


def find_overlapping_regions(
    original_regions: Bed,
    regions_to_compare: Bed,
) -> Bed:
    """
    Find overlapping regions between the bed file and the regions.
    """
    overlapping_regions = original_regions.intersect(
        regions_to_compare,
        u=True,
        wb=True,
    )

    return overlapping_regions


def black_list_filter_coverage_df(
    coverage_df: pl.DataFrame,
    exclusion_bed_file: Path,
) -> pl.DataFrame:
    """
    Filter the coverage DataFrame to remove outliers based on quantiles.
    """

    # Extract the bed file from the DataFrame
    extracted_bed_file = extract_bed_file_from_df(
        coverage_df,
        output_file=None,
    )

    # Find overlapping regions
    overlapping_regions_bed = find_overlapping_regions(
        extracted_bed_file,
        Bed(exclusion_bed_file),
    )

    # Turn into polars dataframe
    overlapping_regions_df = pl.read_csv(
        overlapping_regions_bed.fn,
        separator="\t",
        has_header=False,
        new_columns=[
            "chrom",
            "start",
            "end",
        ],
    )

    # Filter the initial dataframe based on antimatching these regions
    filtered_df = coverage_df.join(
        overlapping_regions_df,
        on=["chrom", "start", "end"],
        how="anti",
    )

    return filtered_df


# Consecutively filter the coverage DataFrame
Working_100bp.filtered_aberrant_coverage = black_list_filter_coverage_df(
    Working_100bp.all_samples_coverages,
    Config.exclusion_low_coverage_regions,
)
Working_100bp.filtered_aberrant_coverage = black_list_filter_coverage_df(
    Working_100bp.filtered_aberrant_coverage,
    Config.exclusion_high_coverage_regions,
)

# No separate low-quality filter is needed: coverage was computed with MAPQ >= 20.


# ---------------------------------------------------------------------------
# Filtering: GC content
# ---------------------------------------------------------------------------
Metadata.gc_content_df = pl.read_csv(
    Config.gc_content_file,
    separator="\t",
    has_header=False,
    new_columns=[
        "chrom",
        "start",
        "end",
        "gc_content",
    ],
).filter(pl.col("chrom").is_in(Config.allowed_chromosomes))


def filter_gc_content(
    coverage_df: pl.DataFrame,
    gc_content_df: pl.DataFrame,
    lower_quantile_bound: float = Config.gc_lower_quantile,
    upper_quantile_bound: float = Config.gc_upper_quantile,
) -> pl.DataFrame:
    """
    Filter the coverage DataFrame based on GC content.
    """

    # Join the coverage DataFrame with the GC content DataFrame
    merged_df = coverage_df.join(
        gc_content_df,
        on=["chrom", "start", "end"],
        how="inner",
    )
    # Calculate the quantiles for GC content
    lower_quantile_bound = merged_df.select(
        pl.col("gc_content").quantile(lower_quantile_bound)
    ).item()

    upper_quantile_bound = merged_df.select(
        pl.col("gc_content").quantile(upper_quantile_bound)
    ).item()

    # Filter based on the GC content
    filtered_df = merged_df.filter(
        (pl.col("gc_content") > lower_quantile_bound)
        & (pl.col("gc_content") < upper_quantile_bound)
    ).drop(
        "gc_content",
    )

    return filtered_df


Working_100bp.filtered_gc_content = filter_gc_content(
    Working_100bp.filtered_aberrant_coverage,
    Metadata.gc_content_df,
    Config.gc_lower_quantile,
    Config.gc_upper_quantile,
)


# ---------------------------------------------------------------------------
# Filtering: mappability
# ---------------------------------------------------------------------------
Metadata.mappability_df = pl.read_csv(
    Config.mappability_file,
    separator="\t",
    has_header=False,
    new_columns=[
        "chrom",
        "start",
        "end",
        "coverage",
    ],
).filter(pl.col("chrom").is_in(Config.allowed_chromosomes))

Metadata.mappability_df = Metadata.mappability_df.with_columns(
    (pl.col("coverage") / pl.col("coverage").max()).alias("mappability")
).drop("coverage")


def filter_mappability(
    coverage_df: pl.DataFrame,
    mappability_df: pl.DataFrame,
    lower_bound: float,
) -> pl.DataFrame:
    """
    Filter the coverage DataFrame based on mappability.
    """

    # Join the coverage DataFrame with the mappability DataFrame
    merged_df = coverage_df.join(
        mappability_df,
        on=["chrom", "start", "end"],
        how="inner",
    )

    # Filter based on the mappability
    filtered_df = merged_df.filter(
        (pl.col("mappability") > lower_bound)
    ).drop(
        "mappability",
    )

    return filtered_df


Working_100bp.filtered_mappability = filter_mappability(
    Working_100bp.filtered_gc_content,
    Metadata.mappability_df,
    Config.mappability_lower_bound,
)


# ---------------------------------------------------------------------------
# GC-content correction (LOESS)
# ---------------------------------------------------------------------------
Working_100bp.filtered_mappability = Working_100bp.filtered_mappability.join(
    Metadata.gc_content_df,
    on=["chrom", "start", "end"],
    how="inner",
)


def summarise_by_gc(
    df: pl.DataFrame,
    gc_col: str,
    depth_col: str,
    min_coverage: float = 1,
    resolution: int = 100,
) -> pl.DataFrame:
    """
    Collapse window-level depths to one row per GC%
    retaining:
      gc_pct       - integer GC content (0..100)
      mu_depth     - mean depth at that GC
      n_windows    - #windows at that GC  (=> Loess weight)
    """
    df = df.select(gc_col, depth_col).filter(
        (pl.col(depth_col) > min_coverage)
    )

    return (
        df.with_columns(
            (pl.col(gc_col) * resolution)
            .round(0)
            .cast(pl.Int8)
            .alias("gc_pct")
        )
        .group_by("gc_pct")
        .agg(
            [
                pl.mean(depth_col).alias("mu_depth"),
                pl.len().alias("n_windows"),
            ]
        )
        .sort("gc_pct")
    )


def gc_correct_array(
    input_array,
    gc_array,
    weights=None,
    span: float = 0.35,
    degree: int = 1,
    family: str = "symmetric",
):
    assert isinstance(
        input_array, np.ndarray
    ), "input_array must be a numpy array"
    assert isinstance(
        gc_array, np.ndarray
    ), "gc_array must be a numpy array"
    assert (
        input_array.shape == gc_array.shape
    ), "input_array and gc_array must have the same shape"

    fit = loess(
        x=gc_array,
        y=input_array,
        weights=weights,
        span=span,
        degree=degree,
        family=family,
        surface="direct",
    )
    fit.fit()

    return fit


def gc_correct_dataframe_col_binned(
    df: pl.DataFrame,
    columns_to_correct: list[str],
    gc_col: str,
    min_thr: float = 0.0,
    max_thr: float = 1.0,
    min_coverage=1,
    resolution: int = 100,
    span: float = 0.35,
):
    working = df.select(gc_col, *columns_to_correct)
    masked = working.filter(
        (pl.col(gc_col) > min_thr) & (pl.col(gc_col) < max_thr)
    )

    # --- scale original GC to match the training scale (!!!) -------------
    orig_gc_scaled = (df[gc_col].to_numpy() * resolution).round(0)

    corrected = {}
    for c in columns_to_correct:
        bdf = summarise_by_gc(
            masked,
            gc_col,
            c,
            min_coverage=min_coverage,
            resolution=resolution,
        )

        x = bdf["gc_pct"].to_numpy().astype(float)
        y = bdf["mu_depth"].to_numpy().astype(float)
        w = bdf["n_windows"].to_numpy().astype(float)

        fit = gc_correct_array(y, x, weights=w, span=span)

        # Scale based on the predicted gc. This is a multiplicative scaling approach
        #  I thought it would be best since GC effects are multiplicative
        y_hat = fit.predict(orig_gc_scaled).values

        # Avoid division by zero
        eps = 1e-6
        y_hat = np.maximum(y_hat, eps)

        m_global = df[c].mean()

        scale = m_global / y_hat

        corrected[c] = df[c].to_numpy() * scale

    return df.with_columns(
        [
            pl.Series(
                f"gc_corrected_{c.replace('coverage_', '')}",
                corrected[c],
            )
            for c in columns_to_correct
        ]
    )


coverage_cols = [
    col
    for col in Working_100bp.filtered_mappability.columns
    if col.startswith("coverage_")
]

Working_100bp.gc_corrected_df = gc_correct_dataframe_col_binned(
    Working_100bp.filtered_mappability,
    coverage_cols,
    "gc_content",
)

Working_100bp.gc_corrected_df = Working_100bp.gc_corrected_df.drop(
    pl.col(r"^coverage_.+$")
)


# ---------------------------------------------------------------------------
# Replicate merge and median normalization
# ---------------------------------------------------------------------------
# We are working on the gc corrected vals
Working_100bp.normalized_df = Working_100bp.gc_corrected_df.select(
    pl.col("chrom"),
    pl.col("start"),
    pl.col("end"),
    *[
        # Use the median to avoid outliers
        (pl.col(col) / pl.col(col).median()).alias(
            col.replace("gc_corrected_", "median_normalized_")
        )
        for col in Working_100bp.gc_corrected_df.columns
        if col.startswith("gc_corrected")
    ],
)
# Sort for visualization
sorted_columns = sorted(Working_100bp.normalized_df.columns)
Working_100bp.normalized_df = Working_100bp.normalized_df.select(sorted_columns)


# ---------------------------------------------------------------------------
# Mappability-bias correction (LOESS, reuses the GC-correction helpers)
# ---------------------------------------------------------------------------
Working_100bp.mappability_correction_df = Working_100bp.normalized_df.join(
    Metadata.mappability_df,
    on=["chrom", "start", "end"],
    how="inner",
)

coverage_cols = [
    col
    for col in Working_100bp.mappability_correction_df.columns
    if col.startswith("median_normalized_")
]

Working_100bp.mappability_correction_df = gc_correct_dataframe_col_binned(
    Working_100bp.mappability_correction_df,
    coverage_cols,
    "mappability",
    min_coverage=0.01,
    max_thr=0.95,
    span=0.3,
    resolution=100,
)

Working_100bp.mappability_correction_df = (
    Working_100bp.mappability_correction_df.select(
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        *[
            col
            for col in Working_100bp.mappability_correction_df.columns
            if col.startswith("gc_corrected")
        ],
        pl.col("mappability"),
    )
)
rename_map = {
    name: name.replace("gc_corrected_median_normalized_", "corrected_")
    for name in Working_100bp.mappability_correction_df.columns
    if name.startswith("gc_corrected_median_normalized_")
}
Working_100bp.mappability_correction_df = (
    Working_100bp.mappability_correction_df.rename(rename_map)
)


# ---------------------------------------------------------------------------
# Replicate collapse and log-ratio of medians
# ---------------------------------------------------------------------------
coverage_columns = [
    col
    for col in Working_100bp.mappability_correction_df.columns
    if col.startswith("corrected_")
]

column_sample_groups = {}
for col in coverage_columns:
    sample_name = col.replace("corrected_", "").split("-")[0]
    column_sample_groups.setdefault(sample_name, []).append(col)

Working_100bp.dereplicated_df = Working_100bp.mappability_correction_df.select(
    # keep your coordinate columns
    pl.col("chrom"),
    pl.col("start"),
    pl.col("end"),
    # for each sample, do a single horizontal sum
    *[
        pl.sum_horizontal(pl.col(c) for c in cols).alias(sample)
        for sample, cols in column_sample_groups.items()
    ],
)

# Polars does not provide a median horizontal function, so we need to do a list trick
Working_100bp.final_df = Working_100bp.dereplicated_df.clone()
Working_100bp.final_df = Working_100bp.final_df.select(
    pl.col("chrom"),
    pl.col("start"),
    pl.col("end"),
    pl.concat_list(pl.col(r"^E.+$")).list.median().alias("median_early"),
    pl.concat_list(pl.col(r"^L.+$")).list.median().alias("median_late"),
)

# Remember to also clip the values
eps = 1e-6
Working_100bp.final_df = Working_100bp.final_df.with_columns(
    pl.col("chrom"),
    pl.col("start"),
    pl.col("end"),
    log2_ratio_late_early=(
        ((pl.col("median_late") + eps) / (pl.col("median_early") + eps))
    )
    .log(2)
    .clip(-8, 8),
)


# ---------------------------------------------------------------------------
# Filter low-coverage bins (similar to CNVkit), per regime median
# ---------------------------------------------------------------------------
Working_100bp.final_df = Working_100bp.final_df.filter(
    (pl.col("median_late") > Config.median_coverage_threshold)
    & (pl.col("median_early") > Config.median_coverage_threshold)
)


# ---------------------------------------------------------------------------
# Denoising: stationary wavelet transform (MODWT) with SureShrink soft threshold
# ---------------------------------------------------------------------------
def sure_threshold(detail_coeffs: np.ndarray, sigma: float) -> float:
    """
    Donoho-Johnstone SureShrink (soft) threshold for one wavelet scale, with the SureShrink method.
    """
    n = detail_coeffs.size
    if n == 0 or sigma == 0:
        return 0.0

    d = detail_coeffs / sigma
    y = np.sort(d**2)  # ascending
    cumsum = np.add.accumulate(y)
    k = np.arange(1, n + 1)

    risk = (n - 2 * k + cumsum + (n - k) * y) / n
    idx = np.argmin(risk)
    t_sure = np.sqrt(y[idx]) * sigma

    # In SureShrink, the minimum threshold is capped by the universal threshold
    t_univ = sigma * np.sqrt(2.0 * np.log(n))
    return min(t_sure, t_univ)


def pad_signal_for_modwt(signal, level):

    # See pywt.swt for details on padding requirements
    n = signal.shape[0]
    pad_length = (2**level - (n % 2**level)) % 2**level
    if pad_length == 0:
        return signal
    return np.pad(signal, (0, pad_length), mode="wrap")


def modwt_sure_denoise(signal, wavelet="haar", level=None):

    # Set to  max level if not specified
    if level is None:
        level = pywt.swt_max_level(signal.size)
    padded_signal = pad_signal_for_modwt(signal, level)

    # This is equivalent to modwt see https://github.com/PyWavelets/pywt/pull/476
    coeffs = pywt.swt(
        padded_signal,
        wavelet=wavelet,
        level=level,
        trim_approx=True,
        norm=True,
    )

    # Estimate noise sigma via MAD on the finest detail band
    finest = coeffs[-1]
    sigma = np.median(np.abs(finest)) / 0.6745
    if not np.isfinite(sigma):

        print("Warning: sigma is not finite, setting to 0.0")
        sigma = 0.0

    # Apply SURE-shrink to each detail level
    new_coeffs = [coeffs[0]]

    for d in coeffs[1:]:
        thr = max((0.0, sure_threshold(d, sigma)))
        new_coeffs.append(pywt.threshold(d, thr, mode="soft"))

    # Norm MUST be set to True to ensure the inverse transform is correct
    denoised = pywt.iswt(new_coeffs, wavelet=wavelet, norm=True)

    # Trim padding
    return denoised[: len(signal)]


array_to_denoise = (
    Working_100bp.final_df["log2_ratio_late_early"].to_numpy().flatten()
)

denoised_values = modwt_sure_denoise(
    array_to_denoise,
    wavelet=Config.modwt_wavelet,
    level=Config.modwt_level,
)
Working_100bp.denoised_df = Working_100bp.final_df.with_columns(
    pl.Series(
        "log2_ratio_late_early_denoised",
        denoised_values,
    )
)

if Config.save_results:
    Working_100bp.denoised_df.write_csv(
        Config.working_dir / "for_plotting" / "denoised_coverage.tsv",
        separator="\t",
        include_header=True,
    )

# Save as a bed file
if Config.save_results == True:
    Working_100bp.denoised_df.select(
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col("log2_ratio_late_early_denoised"),
    ).write_csv(
        Config.working_dir
        / "A_obtectus_denoised_modwt_log2_ratio_late_early_100bp.bedgraph",
        separator="\t",
        include_header=False,
    )

# Save other tracks for manual inspection
if Config.save_results == True:
    Working_100bp.denoised_df.select(
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col("median_early"),
    ).write_csv(
        Config.working_dir / "A_obtectus_median_early_coverage_100bp.bedgraph",
        separator="\t",
        include_header=False,
    )
    Working_100bp.denoised_df.select(
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col("median_late"),
    ).write_csv(
        Config.working_dir / "A_obtectus_median_late_coverage_100bp.bedgraph",
        separator="\t",
        include_header=False,
    )

Working_100bp.segments_tmp_modwt = pl.read_csv(
    Config.segmentation_indices_file,
    separator="\t",
    has_header=False,
    new_columns=[
        "chrom",
        "start",
        "end",
        "segment_index",
    ],
).filter(pl.col("chrom").is_in(Config.allowed_chromosomes))

Working_100bp.segmented_df = Working_100bp.denoised_df.join(
    Working_100bp.segments_tmp_modwt,
    on=["chrom", "start", "end"],
    how="full",
    coalesce=True,
).sort(
    [
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
    ]
)


# ---------------------------------------------------------------------------
# Segment filtering: null fraction and length
# ---------------------------------------------------------------------------
Working_100bp.segmented_df = Working_100bp.segmented_df.with_columns(
    pl.col("segment_index")
    .interpolate(method="nearest")
    .forward_fill()
    .backward_fill()
    .alias("segment")
)

Filtering.null_proportion = Working_100bp.segmented_df.group_by("segment").agg(
    [
        (
            pl.col("log2_ratio_late_early_denoised").is_null().sum()
            / pl.len()
        ).alias("null_fraction")
    ]
)

# Rejoin the segment proportion info, and filter based on the null fraction
Working_100bp.segmented_df = (
    Working_100bp.segmented_df.join(
        Filtering.null_proportion,
        on="segment",
        how="inner",
    )
    .filter(pl.col("null_fraction") < 0.6)
    .drop(
        "null_fraction",
    )
)

# Length filtering: drop segments shorter than 1 kb (10 bins of 100 bp)
Filtering.segment_lengths = Working_100bp.segmented_df.group_by("segment").agg(
    [(pl.len()).alias("segment_length")]
)

# Rejoin the segment length info, and filter based on the length
Working_100bp.segmented_df = (
    Working_100bp.segmented_df.join(
        Filtering.segment_lengths,
        on="segment",
        how="inner",
    )
    .filter(pl.col("segment_length") > 10)
    .drop(
        "segment_length",
    )
)

# ---------------------------------------------------------------------------
# Permutation test for segment differentiation between regimes
# ---------------------------------------------------------------------------
def permute_segment_pval_optimized(X, n_perm=None, rng=None):
    """
    Non-parametric line-label permutation test for one CBS segment.

    Parameters
    ----------
    X : (n_lines, n_bins) float
        Coverage matrix for the segment.  Each row = one technical replicate / line.
        Rows 0-3 are the observed 'L', rows 4-8 the observed 'E'.
    n_perm : int or None
        If None or >= C(n_lines, n_E), enumerate all permutations.
        Otherwise draw `n_perm` random permutations.
    rng : np.random.Generator or None
        PRNG to use.

    Returns
    -------
    p_val : float
        Two-sided permutation P-value.
    stat_obs : float
        Observed E-minus-L segment statistic.
    """

    rng = np.random.default_rng(rng)

    n_lines, n_bins = X.shape
    if n_lines < 2 or n_lines % 2 != 0:
        raise ValueError(
            f"Need an even number of lines >= 2; got n_lines={n_lines}"
        )

    middle = n_lines // 2

    
    # Statistic is E-minus-L.
    idx_L_obs_true = np.arange(
        middle
    )  # Indices for observed L (rows 0 to middle-1)
    idx_E_obs_true = np.arange(
        middle, n_lines
    )  # Indices for observed E (rows middle to n_lines-1)

    def segment_stat_single(idx_E, idx_L):
        seg_E = X[idx_E]
        seg_L = X[idx_L]

        seg_E_median = np.median(seg_E, axis=1)
        seg_L_median = np.median(seg_L, axis=1)

        # For each coverage matrix, take the median of each bin (vertically)
        # Then the statistic is the mean of the median of 'E' minus the mean of the median of 'L'
        stat = seg_L_median.mean() - seg_E_median.mean()
        return stat

    stat_obs = segment_stat_single(idx_E_obs_true, idx_L_obs_true)

    # Generate or select permutations for the 'E' group
    possible_combinations = list(combinations(range(n_lines), middle))

    num_possible_combinations = len(possible_combinations)

    if n_perm is None or n_perm >= num_possible_combinations:
        selected_idx_E_tuples = possible_combinations
    else:
        chosen_indices = rng.choice(
            num_possible_combinations, size=n_perm, replace=False
        )
        selected_idx_E_tuples = [
            possible_combinations[i] for i in chosen_indices
        ]

    num_selected_perms = len(selected_idx_E_tuples)
    # Prepare index arrays for vectorized operations
    # idx_E_perm will store indices for the 'E' group for each permutation
    idx_E_perm = np.array(
        [list(idx_tuple) for idx_tuple in selected_idx_E_tuples],
        dtype=np.int_,
    )

    # idx_L_perm will store indices for the 'L' group
    # Size of L group is n_lines - middle
    idx_L_perm = np.empty(
        (num_selected_perms, n_lines - middle), dtype=np.int_
    )

    all_indices_arange = np.arange(n_lines)
    for i in range(num_selected_perms):
        idx_L_perm[i] = np.setdiff1d(
            all_indices_arange, idx_E_perm[i], assume_unique=True
        )

    # Vectorized calculation of statistics for all permutations
    # 1. Select data for all 'E' groups across permutations
    # X_E_permuted will have shape (num_selected_perms, middle, n_bins)
    X_E_permuted = X[idx_E_perm]

    # 2. Calculate median for each bin for 'E' groups
    # med_E_all will have shape (num_selected_perms, n_bins)
    med_E_all = np.median(
        X_E_permuted, axis=1
    )  # Median along the 'middle' dimension
    # 3. Calculate mean of medians for 'E' groups
    # mean_med_E_all will have shape (num_selected_perms,)
    mean_med_E_all = med_E_all.mean(
        axis=1
    )  # Mean along the 'n_bins' dimension
    # Repeat for 'L' groups
    # X_L_permuted will have shape (num_selected_perms, n_lines - middle, n_bins)
    X_L_permuted = X[idx_L_perm]

    # med_L_all will have shape (num_selected_perms, n_bins)
    med_L_all = np.median(
        X_L_permuted, axis=1
    )  # Median along the 'n_lines - middle' dimension
    # mean_med_L_all will have shape (num_selected_perms,)
    mean_med_L_all = med_L_all.mean(
        axis=1
    )  # Mean along the 'n_bins' dimension
    # Final statistics for all permutations
    stats_permuted = mean_med_E_all - mean_med_L_all
    extreme_count = np.sum(np.abs(stats_permuted) >= abs(stat_obs))
    p_val = (extreme_count + 1) / (num_selected_perms + 1)

    return p_val, stat_obs


def run_permutation_all_segments_per_chromosome(
    segments_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    n_perm: int = 70,
    rng: np.random.Generator = np.random.default_rng(42),
) -> pd.DataFrame:
    """
    Run permutation test on all segments in the segments_df for a single chromosome.

    Parameters
    ----------
    segments_df : pd.DataFrame
        DataFrame with at least ['chrom', 'start', 'end'] columns.
    coverage_df : pd.DataFrame
        DataFrame with at least ['chrom', 'start', 'end', ...] columns.
    n_perm : int
        Number of permutations to run.
    rng : np.random.Generator
        PRNG to use.

    Returns
    -------
    pd.DataFrame
        Copy of segments_df filtered to `chrom`, with a new 'p_value' column.
    """
    working_df = segments_df.copy()
    pvals = []

    # Iterate over segment rows
    for row in working_df.itertuples(index=False):

        seg_chrom = row.chrom
        start = row.start
        end = row.end

        # Subset coverage for this segment
        sub = coverage_df[
            (coverage_df["chrom"] == seg_chrom)
            & (coverage_df["start"] >= start)
            & (coverage_df["end"] <= end)
        ]

        # Drop unwanted cols, convert to numpy and transpose
        coverage_matrix = (
            sub.drop(
                columns=["chrom", "start", "end"],
                errors="ignore",
            )
            .to_numpy()
            .T
        )

        # Compute p-value for this segment
        p_val, stat_obs = permute_segment_pval_optimized(
            coverage_matrix, n_perm=n_perm, rng=rng
        )
        pvals.append(p_val)

    # Attach p-values
    working_df["p_value"] = pvals

    return working_df


def _process_chromosome(
    chrom_df_tuple: tuple,
    n_perm: int,
    seed: int,
):
    """
    Worker helper: runs permutations for one chromosome.
    """
    # create a private RNG for this worker
    rng = np.random.default_rng(seed)
    return run_permutation_all_segments_per_chromosome(
        segments_df=chrom_df_tuple[0],
        coverage_df=chrom_df_tuple[1],
        n_perm=n_perm,
        rng=rng,
    )


def run_permutation_all_segments_parallel(
    segments_df: pl.DataFrame,
    coverage_df: pl.DataFrame,
    n_perm: int = 70,
    rng: np.random.Generator = np.random.default_rng(42),
    n_workers: int = None,
) -> pl.DataFrame:
    """
    Parallelize run_permutation_all_segments_per_chromosome over all chromosomes.

    Parameters
    ----------
    segments_df : pl.DataFrame
        must have a "chrom" column
    coverage_df : pl.DataFrame
    n_perm : int
        Number of permutations to run per segment
    rng : np.random.Generator
        Base RNG; will be used to derive one seed per chromosome
    n_workers : int or None
        Number of parallel workers (defaults to os.cpu_count())

    Returns
    -------
    pl.DataFrame
        All segments with an added "p_value" column.
    """
    # get unique chromosomes in the order they appear
    chroms = segments_df["chrom"].unique().to_list()
    chroms.sort()

    # Get pandas chrom dataframes
    chrom_df_tuples_list = [
        (
            segments_df.filter(pl.col("chrom") == chrom).to_pandas(),
            coverage_df.filter(pl.col("chrom") == chrom).to_pandas(),
        )
        for chrom in chroms
    ]

    # derive one independent seed per chromosome
    seeds = rng.integers(low=0, high=2**32 - 1, size=len(chroms))

    results = []

    # The fork default context is not compatible with Polars
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as exe:
        # submit one job per chromosome
        futures = {
            exe.submit(
                _process_chromosome,
                chrom_df_tuple,
                n_perm,
                int(seed),
            ): chrom
            for chrom, chrom_df_tuple, seed in zip(
                chroms, chrom_df_tuples_list, seeds
            )
        }

        # collect results as they finish, with a top-level progress bar
        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Chromosomes",
        ):
            chrom = futures[future]
            df_chrom = future.result()
            results.append(df_chrom)

    results_df = pl.concat(
        [
            pl.from_pandas(chrom_df, include_index=True)
            for chrom_df in results
        ],
        how="vertical",
        rechunk=True,
    )

    # Rejoin with the original segments_df to keep all columns
    final_df = segments_df.join(
        results_df,
        on=["chrom", "start", "end", "segment"],
        how="left",
    )
    return final_df, results_df


Working_100bp.segmented_df = Working_100bp.segmented_df.drop("segment_index")

# Get the filtered segments
Working_100bp.filtered_segments_step_1_df = (
    Working_100bp.segmented_df.group_by("segment", maintain_order=True).agg(
        [
            pl.col("chrom").first().alias("chrom"),
            pl.col("start").first().alias("start"),
            pl.col("end").last().alias("end"),
            pl.col("log2_ratio_late_early_denoised")
            .mean()
            .alias("segment_mean"),
        ]
    )
)

# Actual run, warning this will take a while
Filtering.pval_segments_df, Filtering.tmp_df = (
    run_permutation_all_segments_parallel(
        Working_100bp.filtered_segments_step_1_df,
        Working_100bp.dereplicated_df,
        n_perm=70,
        rng=np.random.default_rng(42),
        n_workers=12,
    )
)

# Save the results cause I am not running this again
if Config.save_results == True:
    Filtering.pval_segments_df.write_csv(
        Config.working_dir / "A_obtectus_pval_segments_100bp_undo4.tsv",
        separator="\t",
        include_header=True,
    )

Filtering.pval_segments_df = pl.read_csv(
    Config.working_dir / "A_obtectus_pval_segments_100bp_undo4.tsv",
    separator="\t",
    has_header=True,
)


# ---------------------------------------------------------------------------
# Filtering on p-value and fold change
# ---------------------------------------------------------------------------
# rejoin the p-values with the original bin-level data
Filtering.segmented_df_pvals = Working_100bp.segmented_df.join(
    Filtering.pval_segments_df.select(pl.col("segment"), pl.col("p_value")),
    on="segment",
    how="inner",
)

Filtering.segments_df_only_significant = Filtering.segmented_df_pvals.filter(
    pl.col("p_value") < Config.alpha_threshold
)

# Aggregate data: calculate mean, and N (number)
Filtering.aggregated_df = Filtering.segments_df_only_significant.group_by(
    "segment"
).agg(
    [
        # don't use denoised here, get more ground truth data
        pl.mean("log2_ratio_late_early").alias("log2_ratio_late_early_mean"),
        pl.len().alias("n"),  # Calculate group size N
    ]
)

# A 1.5 fold change is the worst case for a real heterozygous CNV (2 -> 3 gain);
# deletions should exceed it, so it is a reasonable lower bound.

# reattach mean calculated above
Filtering.segments_df_only_significant = (
    Filtering.segments_df_only_significant.join(
        Filtering.aggregated_df.select(
            pl.col("segment"),
            pl.col("log2_ratio_late_early_mean"),
            pl.col("n"),
        ),
        on="segment",
        how="inner",
    )
)

lower_filter = pl.col("log2_ratio_late_early_mean") < -np.log2(
    Config.fold_change_threshold
)

upper_filter = pl.col("log2_ratio_late_early_mean") > np.log2(
    Config.fold_change_threshold
)

Filtering.segmented_df_log_fold_change = (
    Filtering.segments_df_only_significant.filter(
        (lower_filter | upper_filter)
    ).sort("segment")
)

# Create a bedgraph file with the segments that have a significant change
# Group segments, keeping the first start and last end
Working_100bp.filtered_segments_to_save = (
    Filtering.segmented_df_log_fold_change.group_by("segment")
    .agg(
        [
            # Last so that any errors show up in the end
            pl.last("chrom"),
            pl.first("start").alias("start"),
            pl.last("end").alias("end"),
            pl.mean("log2_ratio_late_early").alias(
                "log2_ratio_late_early_mean"
            ),
        ]
    )
    .sort(
        [
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end"),
        ]
    )
    .drop("segment")
)

Working_100bp.filtered_segments_to_save.write_csv(
    Config.working_dir
    / "A_obtectus_filtered_segments_log2_ratio_late_early_100bp.bedgraph",
    separator="\t",
    include_header=False,
)


# ---------------------------------------------------------------------------
# Calling deletions and gains from the per-regime medians
#
# A baseline is taken from the segments NOT called as differentiated. The X
# chromosome is scaled separately from the autosomes (it should be ~1.5x
# diploid), otherwise dividing by an all-segment mean would be wrong.
# ---------------------------------------------------------------------------
# Create a dataframe that does not have the significant segments
Working_100bp.non_significant_segments = Working_100bp.final_df.join(
    Filtering.segmented_df_log_fold_change.select("chrom", "start", "end"),
    on=["chrom", "start", "end"],
    how="anti",
)

Working_100bp.non_significant_segments_chrom_X = (
    Working_100bp.non_significant_segments.filter(
        pl.col("chrom") == Config.chromosome_x
    )
)
Working_100bp.non_significant_segments_excluding_X = (
    Working_100bp.non_significant_segments.filter(
        pl.col("chrom") != Config.chromosome_x
    )
)

# Calculate mean expectation for the non-significant segments
mean_coverage_non_significant_early_excluding_X = (
    Working_100bp.non_significant_segments_excluding_X["median_early"].mean()
)
mean_coverage_non_significant_late_excluding_X = (
    Working_100bp.non_significant_segments_excluding_X["median_late"].mean()
)

mean_coverage_non_significant_early_X = (
    Working_100bp.non_significant_segments_chrom_X["median_early"].mean()
)
mean_coverage_non_significant_late_X = (
    Working_100bp.non_significant_segments_chrom_X["median_late"].mean()
)

# 1) aggregate raw moments
agg = Filtering.segmented_df_log_fold_change.group_by("segment").agg(
    [
        pl.first("chrom").alias("chrom"),
        pl.first("start").alias("start"),
        pl.last("end").alias("end"),
        # raw means
        pl.mean("median_early").alias("mean_med_early"),
        pl.mean("median_late").alias("mean_med_late"),
        # sample size
        pl.len().alias("n"),
    ]
)

# 2) scale to CNV space, check for X
with_scaled = agg.with_columns(
    [
        # *2 so that it showcases the dizygotic CNV space
        # (i.e. 0.5 is a deletion, 1.5 is a duplication)
        pl.when(pl.col("chrom") == Config.chromosome_x)
        .then(
            ((pl.col("mean_med_early")) / mean_coverage_non_significant_early_X)
            * 2
        )
        .otherwise(
            (
                (pl.col("mean_med_early"))
                / mean_coverage_non_significant_early_excluding_X
            )
            * 2
        )
        .alias("cnv_early"),
        pl.when(pl.col("chrom") == Config.chromosome_x)
        .then(
            ((pl.col("mean_med_late")) / mean_coverage_non_significant_late_X)
            * 2
        )
        .otherwise(
            (
                (pl.col("mean_med_late"))
                / mean_coverage_non_significant_late_excluding_X
            )
            * 2
        )
        .alias("cnv_late"),
    ]
)

Working_100bp.cnv_calls_df = with_scaled

Working_100bp.cnv_calls_df = Working_100bp.cnv_calls_df.with_columns(
    difference_early_late=(pl.col("cnv_late") - pl.col("cnv_early")),
)

# Filter out cnv calls differences of at least 1
Working_100bp.cnv_calls_df_filtered = Working_100bp.cnv_calls_df.filter(
    pl.col("difference_early_late").abs() >= 1
)
Working_100bp.final_cnv_calls_to_write = (
    Working_100bp.cnv_calls_df_filtered.select(
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col("cnv_early").alias("cnv_early_scaled"),
        pl.col("cnv_late").alias("cnv_late_scaled"),
        pl.col("difference_early_late").alias("copy_number_early_minus_late"),
    )
)
Working_100bp.final_cnv_calls_to_write.write_csv(
    Config.working_dir
    / "A_obtectus_final_cnv_calls_100bp_undo4_corrected_for_X.tsv",
    separator="\t",
    include_header=True,
)

Working_100bp.final_cnv_bedgraph = Working_100bp.cnv_calls_df_filtered.select(
    pl.col("chrom"),
    pl.col("start"),
    pl.col("end"),
    pl.col("difference_early_late")
    .alias("copy_number_early_minus_late")
    .round(),
)
Working_100bp.final_cnv_bedgraph.write_csv(
    Config.working_dir
    / "A_obtectus_final_cnv_calls_100bp_undo4_corrected_for_X.bedgraph",
    separator="\t",
    include_header=False,
)
