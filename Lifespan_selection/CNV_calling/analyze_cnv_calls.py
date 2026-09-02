#!/usr/bin/env python
"""
Author: Athanasios Toros
Email: athanasiostoros@gmail.com
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import polars as pl
from plotly import express as px
from pybedtools import BedTool as bt
import yaml


def load_config(path: Path = Path(__file__).parent / "config.yaml") -> SimpleNamespace:
    with open(path) as fh:
        raw = yaml.safe_load(fh)

    cfg = SimpleNamespace()
    paths = raw["analysis"]["paths"]
    cfg.working_dir = Path(paths["working_dir"])
    cfg.data_dir = Path(paths["data_dir"])
    cfg.gtf_cds_annotation = Path(paths["gtf_cds_annotation"])
    cfg.length_copy_number_scatter = paths["length_copy_number_scatter"]

    params = raw["analysis"]["params"]
    cfg.call_deletion_threshold = params["call_deletion_threshold"]
    cfg.call_gain_threshold = params["call_gain_threshold"]
    cfg.gene_with_cnv_overlap_threshold = params["gene_with_cnv_overlap_threshold"]
    cfg.deletion_filter_threshold = params["deletion_filter_threshold"]
    cfg.gain_filter_threshold = params["gain_filter_threshold"]
    cfg.save_results_flag = params["save_results"]
    return cfg


Config = load_config()
Working = SimpleNamespace()


# ---------------------------------------------------------------------------
# Fix annotations: exon structure leaves duplicated gene_ids in the CDS BED,
# so collapse to one interval (min start, max end) per gene.
# ---------------------------------------------------------------------------
Working.duplicate_bed_annotations = pl.read_csv(
    Config.gtf_cds_annotation,
    has_header=False,
    new_columns=["chrom", "start", "end", "gene_id"],
    separator="\t",
)

Working.dereplicated_bed_annotations = (
    Working.duplicate_bed_annotations.group_by("gene_id")
    .agg(
        [
            pl.col("chrom").first(),
            pl.col("start").min().alias("start"),
            pl.col("end").max().alias("end"),
        ]
    )
    .sort(["chrom", "start", "end"])
    .select(
        [
            pl.col("chrom"),
            pl.col("start"),
            pl.col("end"),
            pl.col("gene_id"),
        ]
    )
)

if Config.save_results_flag:
    Working.dereplicated_bed_annotations.write_csv(
        Config.working_dir / "a_obtectus_gene_ids_deduplicated.bed",
        include_header=False,
        separator="\t",
    )


# ---------------------------------------------------------------------------
# CNV data
# ---------------------------------------------------------------------------
Working.cnv_calls = pl.read_csv(
    Config.data_dir
    / "A_obtectus_final_cnv_calls_100bp_undo4_corrected_for_X.tsv",
    has_header=True,
    separator="\t",
)


# ---------------------------------------------------------------------------
# Figure: copy number vs length of CNV calls
# ---------------------------------------------------------------------------
Plot = SimpleNamespace()

Plot.cnv_calls = Working.cnv_calls.with_columns(
    [
        (pl.col("end") - pl.col("start")).alias("length"),
        (pl.col("end") - pl.col("start")).log10().alias("log_length"),
        (pl.col("copy_number_early_minus_late").abs()).alias("abs_copy_number"),
    ]
)

# For the x-axis, we will use log10 of the length of CNV calls
tickvals = np.arange(3, 8)  # log10 values for 1kb to 10Mb

median_length = Plot.cnv_calls["length"].median()
mean_length = Plot.cnv_calls["length"].mean()


def format_bp(val):
    bp = 10**val
    if bp >= 1e6:
        return f"{bp/1e6:.0f} Mb"
    elif bp >= 1e3:
        return f"{bp/1e3:.0f} kb"
    else:
        return f"{int(bp)} bp"


ticktext = [format_bp(v) for v in tickvals]

fig = px.scatter(
    Plot.cnv_calls,
    x="log_length",
    y="copy_number_early_minus_late",
    marginal_x="histogram",
    marginal_y="histogram",
    opacity=0.5,
    template="plotly_white",
)
fig.update_layout(
    title="Copy Number vs Length of Segment calls",
    xaxis_title="Log-Length of Segment calls (bp)",
    yaxis_title="Population Copy Number-Early minus Late",
    xaxis=dict(
        tickmode="array",
        tickvals=tickvals,
        ticktext=ticktext,
    ),
    height=800,
)

fig.add_vline(
    x=np.log10(median_length),
    line_width=2,
    line_dash="dash",
    line_color="red",
    yref="paper",
    xanchor="center",
    yanchor="top",
    annotation_text=f"Median: {median_length:.0f} bp",
    annotation_position="top left",
    row=2,  # only draw in the main scatter's row
    col=1,  # only draw in the main scatter's column
)
fig.add_vline(
    x=np.log10(mean_length),
    line_width=2,
    line_dash="dash",
    line_color="purple",
    annotation_text=f"Mean: {mean_length:.0f} bp",
    annotation_position="top right",
    row=2,
    col=1,
)

fig.update_layout(
    annotations=[
        {**a, **{"y": 1.12}}
        for a in fig.to_dict()["layout"]["annotations"]
    ],
    width=800,
    height=800,
)
fig.write_image(Config.length_copy_number_scatter)


# ---------------------------------------------------------------------------
# Call deletions and gains per selection regime
#
# We take the scaled CNV calls (normalized median coverage in diploid space,
# 2N) and decide whether each segment is a deletion or a gain. These are
# population-level calls, so zygosity refers to the population-level average copy number.
# ---------------------------------------------------------------------------
def expr_call_deletion_gain_per_regime(
    cnv_scaled_col, selection_regime, Config=None
):
    """Polars expressions for per-regime deletion and gain magnitude.

    Deletion = (2 - scaled CN) below the deletion threshold; gain =
    (scaled CN - 2) above the gain threshold; otherwise 0. `selection_regime`
    must be "early" or "late" and only names the output columns.
    """
    if selection_regime not in ["early", "late"]:
        raise ValueError(
            "Selection regime must be either 'early' or 'late'."
        )

    return (
        pl.when(pl.col(cnv_scaled_col) < Config.call_deletion_threshold)
        .then((2 - pl.col(cnv_scaled_col)))
        .otherwise(0)
        .alias(f"{selection_regime}_deletion"),
        pl.when(pl.col(cnv_scaled_col) > Config.call_gain_threshold)
        .then((pl.col(cnv_scaled_col) - 2))
        .otherwise(0)
        .alias(f"{selection_regime}_gain"),
    )


Working.cnv_calls = Working.cnv_calls.with_columns(
    expr_call_deletion_gain_per_regime(
        "cnv_early_scaled", "early", Config=Config
    )
    + expr_call_deletion_gain_per_regime(
        "cnv_late_scaled", "late", Config=Config
    )
)

if Config.save_results_flag:
    Working.cnv_calls.write_csv(
        Config.working_dir
        / "A_obtectus_final_cnv_calls_100bp_undo4_with_deletion_gain.tsv",
        include_header=True,
        separator="\t",
    )


def expr_bankers_rounding_to_0_5(column):
    """
    Round a column to the nearest 0.5 increment using bankers' rounding.
    """
    return (
        ((pl.col(column) * 2).round() / 2)
        .cast(pl.Float64)
        .alias(f"{column}_rounded")
    )


# Round out the CNV calls to 0.5 increments
Working.cnv_calls = Working.cnv_calls.with_columns(
    [
        expr_bankers_rounding_to_0_5("early_deletion"),
        expr_bankers_rounding_to_0_5("early_gain"),
        expr_bankers_rounding_to_0_5("late_deletion"),
        expr_bankers_rounding_to_0_5("late_gain"),
    ]
)


# ---------------------------------------------------------------------------
# Save the rounded calls
# ---------------------------------------------------------------------------
Working.cnv_calls_rounded_to_save = Working.cnv_calls.select(
    [
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col("cnv_early_scaled"),
        pl.col("cnv_late_scaled"),
        pl.col("early_deletion_rounded"),
        pl.col("early_gain_rounded"),
        pl.col("late_deletion_rounded"),
        pl.col("late_gain_rounded"),
    ]
)
if Config.save_results_flag:
    Working.cnv_calls_rounded_to_save.write_csv(
        Config.working_dir
        / "A_obtectus_final_cnv_calls_100bp_undo4_with_deletion_gain_rounded.tsv",
        include_header=True,
        separator="\t",
    )


# ---------------------------------------------------------------------------
# Isolate deletion and gain calls per regime, for GO analysis
# ---------------------------------------------------------------------------
Working.early_deletions = Working.cnv_calls.select(
    [
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col("early_deletion"),
    ]
).filter(pl.col("early_deletion") != 0)

Working.late_deletions = Working.cnv_calls.select(
    [
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col("late_deletion"),
    ]
).filter(pl.col("late_deletion") != 0)
Working.early_gain = Working.cnv_calls.select(
    [
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col("early_gain"),
    ]
).filter(pl.col("early_gain") != 0)
Working.late_gain = Working.cnv_calls.select(
    [
        pl.col("chrom"),
        pl.col("start"),
        pl.col("end"),
        pl.col("late_gain"),
    ]
).filter(pl.col("late_gain") != 0)

if Config.save_results_flag:
    Working.early_deletions.write_csv(
        Config.working_dir / "A_obtectus_early_deletions_cnv_calls.bedgraph",
        include_header=False,
        separator="\t",
    )
    Working.late_deletions.write_csv(
        Config.working_dir / "A_obtectus_late_deletions_cnv_calls.bedgraph",
        include_header=False,
        separator="\t",
    )
    Working.early_gain.write_csv(
        Config.working_dir / "A_obtectus_early_gain_cnv_calls.bedgraph",
        include_header=False,
        separator="\t",
    )
    Working.late_gain.write_csv(
        Config.working_dir / "A_obtectus_late_gain_cnv_calls.bedgraph",
        include_header=False,
        separator="\t",
    )


# ---------------------------------------------------------------------------
# Intersect CNV calls with CDS regions
# ---------------------------------------------------------------------------
def intersect_cnv_with_gene_annotations(
    cnv_calls_bed_file,
    gene_annotations_bed_file,
    fraction_overlap,
    cnv_metadata_df,
    column_with_selection_regime_call,
    filter_by_cnv_type=None,
):

    cnv_calls_bed = bt(cnv_calls_bed_file)
    gene_annotations_bed = bt(gene_annotations_bed_file)

    # Intersect the CNV calls with the gene annotations
    intersection_bed = gene_annotations_bed.intersect(
        cnv_calls_bed,
        wa=True,  # Write the original gene annotations
        wb=True,  # Write the original CNV calls
        F=fraction_overlap,  # Fraction of overlap required to consider a hit
    )

    # Convert the intersection to a DataFrame
    intersections_df = pl.read_csv(
        intersection_bed.fn,
        has_header=False,
        separator="\t",
        new_columns=[
            "chrom",
            "start",
            "end",
            "gene_id",
            "cnv_chrom",
            "cnv_start",
            "cnv_end",
            "cnv_copy_number",
        ],
    )

    # Add metadata to the intersections DataFrame
    intersections_df = intersections_df.join(
        cnv_metadata_df,
        left_on=["cnv_chrom", "cnv_start", "cnv_end"],
        right_on=["chrom", "start", "end"],
        how="left",
    )

    # isolate actual selection regime call to the DataFrame
    intersections_df = intersections_df.select(
        "chrom",
        "start",
        "end",
        "gene_id",
        column_with_selection_regime_call,
    )

    # Some genes intersect with multiple CNV calls, so we need to deduplicate
    intersections_df = (
        intersections_df.group_by("gene_id")
        .agg(
            [
                pl.col("chrom").first(),  # Use the first chromosome
                pl.col("start").min().alias("start"),  # Use the minimum start
                pl.col("end").max().alias("end"),  # Use the maximum end
                pl.col(column_with_selection_regime_call)
                .mean()
                .alias(column_with_selection_regime_call),
            ]
        )
        .sort(["chrom", "start", "end"])
    )

    # Reorder the columns to match the original structure
    intersections_df = intersections_df.select(
        "chrom",
        "start",
        "end",
        "gene_id",
        column_with_selection_regime_call,
    )

    # We use this logic to filter out CNVs that are not relevant in a population level
    # Deletions and gains are handled differently based on the call type
    if filter_by_cnv_type is not None:
        intersections_df = intersections_df.filter(
            pl.col(column_with_selection_regime_call) >= filter_by_cnv_type
        )

    return intersections_df


Working.early_deletions_intersection = intersect_cnv_with_gene_annotations(
    Config.working_dir / "A_obtectus_early_deletions_cnv_calls.bedgraph",
    Config.working_dir / "a_obtectus_gene_ids_deduplicated.bed",
    fraction_overlap=Config.gene_with_cnv_overlap_threshold,
    cnv_metadata_df=Working.early_deletions,
    column_with_selection_regime_call="early_deletion",
    filter_by_cnv_type=Config.deletion_filter_threshold,
)
Working.late_deletions_intersection = intersect_cnv_with_gene_annotations(
    Config.working_dir / "A_obtectus_late_deletions_cnv_calls.bedgraph",
    Config.working_dir / "a_obtectus_gene_ids_deduplicated.bed",
    fraction_overlap=Config.gene_with_cnv_overlap_threshold,
    cnv_metadata_df=Working.late_deletions,
    column_with_selection_regime_call="late_deletion",
    filter_by_cnv_type=Config.deletion_filter_threshold,
)
Working.early_gain_intersection = intersect_cnv_with_gene_annotations(
    Config.working_dir / "A_obtectus_early_gain_cnv_calls.bedgraph",
    Config.working_dir / "a_obtectus_gene_ids_deduplicated.bed",
    fraction_overlap=Config.gene_with_cnv_overlap_threshold,
    cnv_metadata_df=Working.early_gain,
    column_with_selection_regime_call="early_gain",
    filter_by_cnv_type=Config.gain_filter_threshold,
)
Working.late_gain_intersection = intersect_cnv_with_gene_annotations(
    Config.working_dir / "A_obtectus_late_gain_cnv_calls.bedgraph",
    Config.working_dir / "a_obtectus_gene_ids_deduplicated.bed",
    fraction_overlap=Config.gene_with_cnv_overlap_threshold,
    cnv_metadata_df=Working.late_gain,
    column_with_selection_regime_call="late_gain",
    filter_by_cnv_type=Config.gain_filter_threshold,
)


# ---------------------------------------------------------------------------
# Save gene lists and annotated bedgraphs for GO analysis
# ---------------------------------------------------------------------------
def save_gene_list_from_intersection(
    intersection_df,
    save_path_bedgraph,
    save_path_gene_list,
    cnv_type=None,
    save_flag=False,
):

    if save_flag != True:
        return

    if (cnv_type is None) or (cnv_type not in ["deletion", "gain"]):
        raise ValueError(
            "Deletion type must be specified. Use 'Deletion' or 'gain'."
        )

    # First save an informative bedgraph file with gene ids and CNV
    intersection_df.write_csv(
        save_path_bedgraph,
        include_header=False,
        separator="\t",
    )

    # The columns are: chrom, start, end, gene_id, and the CNV type
    comment_row = "\t".join(
        [
            "#chrom",
            "start",
            "end",
            "gene_id",
            f"{"DEL" if cnv_type == "Deletion" else "GAIN"}_type",
        ]
    )
    content = []
    # Read the file and insert the comment row at the beginning
    with open(save_path_bedgraph, "r") as f:
        content = f.readlines()

    with open(save_path_bedgraph, "w") as f:
        # Write the comment row first
        f.write(comment_row + "\n")
        # Write the rest of the content
        f.writelines(content)

    # Now save the unique gene ids to a separate file
    intersection_df.select("gene_id").unique().write_csv(
        save_path_gene_list,
        include_header=False,
        separator="\t",
    )
    return


# Save the intersections to files
gene_list_save_dir = Config.working_dir / "gene_intersections"
gene_list_save_dir.mkdir(exist_ok=True)


save_gene_list_from_intersection(
    Working.early_deletions_intersection,
    gene_list_save_dir / "early_deletions_intersection_gene_ids.bedgraph",
    gene_list_save_dir / "early_deletions_intersection_gene_ids.list",
    cnv_type="deletion",
    save_flag=Config.save_results_flag,
)

save_gene_list_from_intersection(
    Working.late_deletions_intersection,
    gene_list_save_dir / "late_deletions_intersection_gene_ids.bedgraph",
    gene_list_save_dir / "late_deletions_intersection_gene_ids.list",
    cnv_type="deletion",
    save_flag=Config.save_results_flag,
)
save_gene_list_from_intersection(
    Working.early_gain_intersection,
    gene_list_save_dir / "early_gain_intersection_gene_ids.bedgraph",
    gene_list_save_dir / "early_gain_intersection_gene_ids.list",
    cnv_type="gain",
    save_flag=Config.save_results_flag,
)
save_gene_list_from_intersection(
    Working.late_gain_intersection,
    gene_list_save_dir / "late_gain_intersection_gene_ids.bedgraph",
    gene_list_save_dir / "late_gain_intersection_gene_ids.list",
    cnv_type="gain",
    save_flag=Config.save_results_flag,
)
