#!/usr/bin/env python
"""
Author: Athanasios Toros
Email: athanasiostoros@gmail.com
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import polars as pl
import plotly.graph_objects as go
import yaml


def load_config(path: Path = Path(__file__).parent / "config.yaml") -> SimpleNamespace:
    with open(path) as fh:
        raw = yaml.safe_load(fh)

    cfg = SimpleNamespace()
    paths = raw["figures"]["paths"]
    cfg.working_dir = Path(paths["working_dir"])
    cfg.filtered_coverage = Path(paths["filtered_coverage"])
    cfg.mappability_corrected_coverage = Path(paths["mappability_corrected_coverage"])
    cfg.denoised_coverage = Path(paths["denoised_coverage"])
    return cfg


Config = load_config()
Plotting = SimpleNamespace()

# Representative region chosen for the figures (clear boundaries, some noise).
Plotting.region_to_analyze = {
    "chrom": "CAVLJG010000004.1",
    "start": 31969700,
    "end": 31977900,
}


def isolate_region(
    df: pl.DataFrame, region: dict, expand_ratio=1
) -> pl.DataFrame:
    """
    Isolate a specific region from the DataFrame.

    Args:
        df (pl.DataFrame): The DataFrame containing coverage data.
        region (dict): A dictionary with keys 'chrom', 'start', and 'end'.

    Returns:
        pl.DataFrame: Filtered DataFrame for the specified region.
    """

    if expand_ratio > 1:
        region = {
            "chrom": region["chrom"],
            "start": int(
                region["start"]
                - (region["end"] - region["start"])
                * (expand_ratio - 1)
                / 2
            ),
            "end": int(
                region["end"]
                + (region["end"] - region["start"])
                * (expand_ratio - 1)
                / 2
            ),
        }

    return df.filter(
        (pl.col("chrom") == region["chrom"])
        & (pl.col("start") >= region["start"])
        & (pl.col("end") <= region["end"])
    )


# ---------------------------------------------------------------------------
# Load pipeline-stage intermediates 
# ---------------------------------------------------------------------------
Plotting.original_filtered_coverage = pl.scan_parquet(Config.filtered_coverage)

coverage_columns = [
    col
    for col in Plotting.original_filtered_coverage.collect_schema().names()
    if col.startswith("coverage_")
]

Plotting.original_filtered_coverage = Plotting.original_filtered_coverage.select(
    [
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        *[
            (pl.col(col) / pl.col(col).mean()).alias(col)
            for col in coverage_columns
        ],
    ]
)

Plotting.corrected_df = pl.scan_parquet(Config.mappability_corrected_coverage)

Plotting.corrected_df = Plotting.corrected_df.select(
    pl.col("chrom"),
    pl.col("start"),
    pl.col("end"),
    *[
        (pl.col(col) / pl.col(col).mean()).alias(col)
        for col in Plotting.corrected_df.collect_schema().names()
        if col.startswith("corrected_")
    ],
)

Plotting.denoised_ratio_df = pl.scan_parquet(Config.denoised_coverage)

Plotting.denoised_ratio_df_region = isolate_region(
    Plotting.denoised_ratio_df,
    Plotting.region_to_analyze,
    expand_ratio=2,
).collect()

Plotting.denoised_ratio_df_region = Plotting.denoised_ratio_df_region.select(
    pl.col("chrom"),
    pl.col("start"),
    pl.col("end"),
    ((pl.col("start") + pl.col("end")) / 2).alias("midpoint"),
    pl.col("log2_ratio_late_early"),
    pl.col("log2_ratio_late_early_denoised"),
)


# ---------------------------------------------------------------------------
# Region-matrix plotter
# ---------------------------------------------------------------------------
def plot_region_matrix(
    df: pl.DataFrame,
    region: dict,
    columns_to_plot: list = None,
    expand_ratio: float = 1,
    y_axis_range: tuple = (0, None),
    color=None,
    title_text=None,
):
    # Isolate region and compute midpoints
    region_df = (
        isolate_region(df, region, expand_ratio=expand_ratio)
        .select(
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end"),
            *columns_to_plot,
        )
        .collect()
    )

    region_df = region_df.select(
        pl.col("chrom"),
        midpoint=(pl.col("start") + pl.col("end")) / 2,
        *(columns_to_plot or df.columns),
    )

    pdf = region_df.to_pandas()

    # Calculate median
    med = np.median(pdf[columns_to_plot], axis=1)
    pdf["median"] = med

    fig = go.Figure()

    # Plot all sample lines
    for col in columns_to_plot:
        fig.add_trace(
            go.Scatter(
                x=pdf["midpoint"],
                y=pdf[col],
                mode="lines",
                line=dict(color="lightgrey", width=1),
                hoverinfo="skip",
                connectgaps=False,
            )
        )

    # Plot median prominently
    fig.add_trace(
        go.Scatter(
            x=pdf["midpoint"],
            y=pdf["median"],
            mode="lines",
            name="Median",
            line=dict(color=color, width=2.5),
            connectgaps=False,
            opacity=0.7,
        )
    )

    fig.update_layout(
        template="plotly_white",
        width=1000,
        height=250,
        margin=dict(l=40, r=20, t=10, b=20),
        font=dict(family="Serif Bold", size=13),
        showlegend=False,
        yaxis=dict(
            range=y_axis_range,
            title=dict(
                text=title_text,
                font=dict(family="Serif Bold", size=13, color="black"),
            ),
            gridcolor="whitesmoke",
        ),
        xaxis=dict(
            showline=False,
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            title_text="",
        ),
    )

    return fig


# ---------------------------------------------------------------------------
# Normalized (filtered) coverage matrices
# ---------------------------------------------------------------------------
columns_early = [
    col
    for col in Plotting.original_filtered_coverage.collect_schema().names()
    if col.startswith("coverage_E")
]
columns_late = [
    col
    for col in Plotting.original_filtered_coverage.collect_schema().names()
    if col.startswith("coverage_L")
]

filtered_region_matrix_early_fig = plot_region_matrix(
    Plotting.original_filtered_coverage,
    Plotting.region_to_analyze,
    columns_to_plot=columns_early,
    expand_ratio=2,
    y_axis_range=(0, 2.3),
    color="#d16c06",
    title_text="Normalized Coverage",
)
filtered_region_matrix_early_fig.write_image(
    Config.working_dir / "filtered_region_matrix_early.svg"
)

filtered_region_matrix_late_fig = plot_region_matrix(
    Plotting.original_filtered_coverage,
    Plotting.region_to_analyze,
    columns_to_plot=columns_late,
    expand_ratio=2,
    y_axis_range=(0, 2.3),
    color="#067dd1",
    title_text="Normalized Coverage",
)
filtered_region_matrix_late_fig.write_image(
    Config.working_dir / "filtered_region_matrix_late.svg"
)


# ---------------------------------------------------------------------------
# GC + mappability corrected coverage matrices
# ---------------------------------------------------------------------------
columns_early = [
    col
    for col in Plotting.corrected_df.collect_schema().names()
    if col.startswith("corrected_E")
]
columns_late = [
    col
    for col in Plotting.corrected_df.collect_schema().names()
    if col.startswith("corrected_L")
]

filtered_region_matrix_corrected_early_fig = plot_region_matrix(
    Plotting.corrected_df,
    Plotting.region_to_analyze,
    columns_to_plot=columns_early,
    expand_ratio=2,
    y_axis_range=(0, 2.3),
    color="#d16c06",
    title_text="Corrected Coverage",
)
filtered_region_matrix_corrected_early_fig.write_image(
    Config.working_dir / "filtered_region_matrix_corrected_early.svg"
)

filtered_region_matrix_corrected_late_fig = plot_region_matrix(
    Plotting.corrected_df,
    Plotting.region_to_analyze,
    columns_to_plot=columns_late,
    expand_ratio=2,
    y_axis_range=(0, 2.3),
    color="#067dd1",
    title_text="Corrected Coverage",
)
filtered_region_matrix_corrected_late_fig.write_image(
    Config.working_dir / "filtered_region_matrix_corrected_late.svg"
)


# ---------------------------------------------------------------------------
# Raw vs. denoised log2(L/E) ratio
# ---------------------------------------------------------------------------
fig = go.Figure()

# Plot the raw log2 ratio
fig.add_trace(
    go.Scatter(
        x=Plotting.denoised_ratio_df_region["midpoint"],
        y=Plotting.denoised_ratio_df_region["log2_ratio_late_early"],
        mode="markers",  # important: switch to 'markers' if you want dots
        name="Raw Ratio",
        connectgaps=False,
        hoverinfo="skip",
        marker=dict(
            color="lightgrey",
            size=4,
            opacity=1,
            line=dict(
                color="darkgrey",  # or use "black" for higher contrast
                width=1,
            ),
        ),
    )
)

# Plot the denoised signal
fig.add_trace(
    go.Scatter(
        x=Plotting.denoised_ratio_df_region["midpoint"],
        y=Plotting.denoised_ratio_df_region["log2_ratio_late_early_denoised"],
        mode="lines",
        name="Denoised Ratio",
        line=dict(color="#a30522", width=2.5),  # A soft, clear red
        connectgaps=False,
    )
)

fig.update_layout(
    template="plotly_white",
    width=1000,
    height=250,
    margin=dict(l=40, r=20, t=10, b=20),
    font=dict(family="Serif Bold", size=12),
    showlegend=False,
    yaxis=dict(
        gridcolor="whitesmoke",
        range=[-0.5, 2.0],  # Adjusted for clarity
        title=dict(
            text="Log₂ Ratio<br>(Late/Early)",
            font=dict(family="Serif Bold", size=13, color="black"),
        ),
    ),
    xaxis=dict(
        showline=False,
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        title_text="",
    ),
)
fig.write_image(Config.working_dir / "log2_ratio_late_early_denoised.svg")


# ---------------------------------------------------------------------------
# Mean segment coverage ratio
# ---------------------------------------------------------------------------
# keep a copy of the per-bin table BEFORE the groupby
bins = Plotting.denoised_ratio_df_region.select(
    "chrom",
    "start",
    "end",
    pl.col("log2_ratio_late_early_denoised")
    .pow(2)
    .alias("original_scale_ratio"),
).with_columns(
    segment_id=(
        pl.when(pl.col("start") < Plotting.region_to_analyze["start"])
        .then(0)
        .when(pl.col("start") > Plotting.region_to_analyze["end"])
        .then(2)
        .otherwise(1)
    )
)

# aggregate once to get the collapsed segments
segments = bins.group_by("segment_id").agg(
    pl.first("chrom").alias("chrom"),
    pl.first("start").alias("seg_start"),
    pl.last("end").alias("seg_end"),
    pl.mean("original_scale_ratio").alias("mean_ratio"),
)

# re-attach the mean to every bin
plot_df = bins.join(
    segments.select(["segment_id", "mean_ratio"]),
    on="segment_id",
    how="left",  # keeps all bins
).sort(["segment_id", "start"])
plot_df = plot_df.with_columns(
    midpoint=(pl.col("start") + pl.col("end")) / 2,
)

fig = go.Figure()

# Plot the segment means
for group in plot_df.partition_by("segment_id"):
    fig.add_trace(
        go.Scatter(
            x=group["midpoint"].to_list(),
            y=group["mean_ratio"].to_list(),
            mode="lines",
            line=dict(color="black", width=2.5),
            connectgaps=False,
            showlegend=False,
        )
    )

fig.update_layout(
    template="plotly_white",
    width=1000,
    height=250,
    margin=dict(l=40, r=20, t=10, b=20),
    font=dict(family="Serif", size=13),
    showlegend=False,
    yaxis=dict(
        gridcolor="whitesmoke",
        title=dict(
            text="Mean Segment Coverage Ratio<br>(Late/Early)",
            font=dict(family="Serif Bold", size=13, color="black"),
        ),
    ),
    xaxis=dict(
        title=dict(
            text=f"Genomic Position on {Plotting.region_to_analyze['chrom']} (Mb)",
            font=dict(family="Serif Bold", size=13, color="black"),
        ),
        showline=True,
        showticklabels=True,
        showgrid=False,
        # Format ticks to be in Megabases (Mb)
        tickvals=np.linspace(
            plot_df["midpoint"].min(), plot_df["midpoint"].max(), 5
        ),
        ticktext=[
            f"{x/1_000_000:.4f} Mb"
            for x in np.linspace(
                plot_df["midpoint"].min(),
                plot_df["midpoint"].max(),
                5,
            )
        ],
    ),
)

# Add annotations for region boundaries
for pos, text in zip(
    [
        Plotting.region_to_analyze["start"],
        Plotting.region_to_analyze["end"],
    ],
    ["Start", "End"],
):
    fig.add_vline(
        x=pos,
        line_width=1.5,
        line_dash="dash",
        line_color="grey",
    )


# Add threshold for significant amplification
fig.add_shape(
    type="rect",
    xref="paper",
    yref="y",
    x0=0,
    x1=1,
    y0=1,
    y1=2.0,  # Using the y-axis max
    fillcolor="rgba(231, 76, 60, 0.2)",  # Light red tint
    line=dict(width=0),
    layer="below",
)
fig.write_image(Config.working_dir / "mean_segment_coverage_ratio.svg")
