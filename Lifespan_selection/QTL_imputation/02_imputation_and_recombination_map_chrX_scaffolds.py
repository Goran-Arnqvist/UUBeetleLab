#!/usr/bin/env python
"""
Author: Athanasios Toros
Email: athanasiostoros@gmail.com
"""

import allel
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform

from qtl_common import (
    create_chromosome_dataframe_list,
    create_recombination_map_dataframe_for_chromosome,
    get_chromosome_label,
    impute_chromosome_df,
    isolate_pedigree_chromosome,
    kosambi_method,
    load_config,
)

config = load_config()

MODEL = config.segmentation.model
PENALTY = config.segmentation.penalty
NA_THRESHOLD = config.imputation.na_threshold
CHROMOSOME_X = config.genome.chromosome_x
SCAFFOLDS = list(
    range(config.genome.scaffold_first, config.genome.scaffold_last + 1)
)
SEGMENTATION_SCAFFOLD = config.genome.scaffold_segmentation_chromosome
ID_COLUMN = config.chrx_scaffolds.id_column
INDIVIDUAL_LABEL = config.chrx_scaffolds.individual_label

phased_regions_100kb = pd.read_csv(config.chrx_scaffolds.input_calls)


def create_combined_pedigree_recombination_map(
    position_count_AB_df,
    position_count_CD_df,
    total_individuals_AB,
    total_individuals_CD,
    chromosome_X=False,
):
    """
    Creates a combined recombination map from two pedigrees (AB and CD).

    This function merges recombination data from two different pedigrees, calculates combined
    recombination fractions, and generates interpolated genetic map distances and recombination rates.

    Parameters
    ----------
    position_count_AB_df : pandas.DataFrame
        DataFrame containing chromosome positions and recombination counts for pedigree AB.
        Must contain columns 'Chromosome_position' and 'count'.
    position_count_CD_df : pandas.DataFrame
        DataFrame containing chromosome positions and recombination counts for pedigree CD.
        Must contain columns 'Chromosome_position' and 'count'.
    total_individuals_AB : int
        Total number of individuals in pedigree AB.
    total_individuals_CD : int
        Total number of individuals in pedigree CD.

    Returns
    -------
    pandas.DataFrame
        Combined recombination map containing the following columns:
        - Chromosome_position: Position on chromosome
        - count_AB: Recombination count from pedigree AB
        - count_CD: Recombination count from pedigree CD
        - count_combined: Combined recombination counts
        - recombination_fraction: Combined recombination fraction
        - genetic_map_distance_cM: Genetic distance in centiMorgans (using Kosambi mapping function)
        - cumulative_genetic_map_distance_cM: Cumulative genetic distance
        - cumulative_genetic_map_distance_cM_interpolated: Interpolated cumulative genetic distance
        - delta_cM: Change in genetic distance between adjacent positions
        - delta_position: Change in physical position between adjacent positions
        - recombination_rate_cM_Mb: Recombination rate in cM/Mb

    Notes
    -----
    - Input DataFrames must be sorted by chromosome position
    - Uses linear interpolation for missing regions
    - Applies Kosambi mapping function to convert recombination fractions to genetic distances
    """

    assert position_count_AB_df[
        "Chromosome_position"
    ].is_monotonic_increasing, "AB dataframe is not sorted"
    assert position_count_CD_df[
        "Chromosome_position"
    ].is_monotonic_increasing, "CD dataframe is not sorted"

    # Step 1: Create combined count dataframe
    combined_df = pd.merge(
        position_count_AB_df[["Chromosome_position", "count"]],
        position_count_CD_df[["Chromosome_position", "count"]],
        on="Chromosome_position",
        how="outer",
        suffixes=("_AB", "_CD"),
    ).sort_values("Chromosome_position")

    # Step 2: Combine the counts, and calculate the recombination fraction for the number of gametes

    combined_df["count_combined"] = combined_df["count_AB"].fillna(
        0
    ) + combined_df["count_CD"].fillna(0)

    # set the 0 values to NaN, necessary for the interpolation
    combined_df["count_combined"] = combined_df[
        "count_combined"
    ].replace(0, np.nan)

    # Gamete pool. For an autosome each F2 contributes two informative gametes
    # (maternal + paternal). For the X only the maternal X recombines: the F1
    # father is hemizygous (XY) and transmits a single, non-recombinant X to his
    # daughters, so each F2 contributes ONE recombinant X gamete -- half the
    # autosomal pool. The published X map calls this with the default
    # chromosome_X=False, for which total_gametes reduces to
    # (total_individuals_AB + total_individuals_CD), i.e. one maternal gamete per
    # F2 -- the X-appropriate denominator. (The nr_individuals intermediate is
    # half the individual count and the branch labels do not read literally; see
    # NOTES.md, incl. the scaffold-12 reuse caveat.)
    nr_individuals = (total_individuals_AB + total_individuals_CD) / 2
    total_gametes = (
        nr_individuals if chromosome_X else 2 * nr_individuals
    )
    combined_df["recombination_fraction"] = (
        combined_df["count_combined"] / total_gametes
    )

    # Step 3: Calculate the genetic map distance, and create a cumulicative column

    combined_df["genetic_map_distance_cM"] = (
        combined_df["recombination_fraction"].apply(kosambi_method) * 100
    )

    combined_df["cumulative_genetic_map_distance_cM"] = combined_df[
        "genetic_map_distance_cM"
    ].cumsum()

    # Step 4: Interpolate the cumulicative genetic distance for missing regions

    combined_df["cumulative_genetic_map_distance_cM_interpolated"] = (
        combined_df["cumulative_genetic_map_distance_cM"].interpolate(
            method="linear"
        )
    )

    # Step 5: Calculate the delta cM and delta position
    combined_df["delta_cM"] = combined_df[
        "cumulative_genetic_map_distance_cM_interpolated"
    ].diff()

    combined_df["delta_position"] = combined_df[
        "Chromosome_position"
    ].diff()

    # Step 6: Calculate the interpolated recombination rate, cM/Mb
    combined_df["recombination_rate_cM_Mb"] = combined_df["delta_cM"] / (
        combined_df["delta_position"] / 1e6
    )

    return combined_df


def impute_chromosome_df_based_on_heuristics(
    chromosome_df, zygote_agreement_threshold, mode_strategy="skip"
):
    """
    Imputes rows in a DataFrame based on the mode value if a certain threshold of agreement is met.

    Parameters:
        chromosome_df (pd.DataFrame): The DataFrame containing chromosome segment data.
        zygote_agreement_threshold (float): A value between 0 and 1 representing the required fraction
            of entries in the row that must agree with the mode for imputation to occur.
        mode_strategy (str): How to handle rows with multiple equally frequent modes. Options:
            - 'first': Use the first mode.
            - 'random': Choose one of the modes randomly.
            - 'skip': Do not impute if there are multiple modes.

    Returns:
        pd.DataFrame: The DataFrame with imputed rows.
    """

    # Iterate over each row in the DataFrame
    for index, row in chromosome_df.iterrows():

        # Skip rows where all values are NaN
        if row.isna().all():
            continue

        # Calculate the mode(s) of the row, skipping NaNs
        modes = row.mode(dropna=True)

        # If there are no modes (empty Series), skip this row
        if modes.empty:
            continue

        # Check if there are multiple modes and handle according to mode_strategy
        if len(modes) > 1:
            if mode_strategy == "first":
                segment_mode = modes.iloc[0]
            elif mode_strategy == "random":
                segment_mode = np.random.choice(modes)
            elif mode_strategy == "skip":
                # Skip imputation for this row if there is ambiguity
                print(f"Skipping row {index} due to multiple modes.")
                continue
            else:
                raise ValueError(
                    "Invalid mode_strategy. Use 'first', 'random', or 'skip'."
                )
        else:
            # Only one mode exists
            segment_mode = modes.iloc[0]

        # Compute the fraction of values in the row that equal the chosen mode
        agreement_fraction = (
            row.dropna() == segment_mode
        ).sum() / row.shape[0]

        # If the fraction exceeds the threshold, impute the entire row with the chosen mode
        if agreement_fraction > zygote_agreement_threshold:
            chromosome_df.loc[index] = segment_mode

    return chromosome_df


# ===========================================================================
# Part 1: impute the X chromosome (chromosome 10)
# ===========================================================================
individuals_AB, chromX_pedigree_AB = isolate_pedigree_chromosome(
    phased_regions_100kb, "AB", [CHROMOSOME_X], id_column=ID_COLUMN
)
individuals_CD, chromX_pedigree_CD = isolate_pedigree_chromosome(
    phased_regions_100kb, "CD", [CHROMOSOME_X], id_column=ID_COLUMN
)

imputed_X_AB = impute_chromosome_df(
    chromX_pedigree_AB, PENALTY, NA_THRESHOLD, chromosome_X=True
)
imputed_X_CD = impute_chromosome_df(
    chromX_pedigree_CD, PENALTY, NA_THRESHOLD, chromosome_X=True
)

imputed_X_AB[ID_COLUMN] = individuals_AB
imputed_X_CD[ID_COLUMN] = individuals_CD

final_df_to_export = pd.concat([imputed_X_AB, imputed_X_CD], axis=0)
final_df_to_export.to_csv(
    config.chrx_scaffolds.x_imputed_output, index=False
)


# ===========================================================================
# Part 2: X recombination map (paternal-X recombination disallowed)
# ===========================================================================
genetic_map_chrom_X_AB = create_recombination_map_dataframe_for_chromosome(
    chromX_pedigree_AB,
    MODEL,
    PENALTY,
    filter_individuals_threshold=config.recombination_map.x_filter_individuals_threshold,
    chromosome_X=True,
)
genetic_map_chrom_X_CD = create_recombination_map_dataframe_for_chromosome(
    chromX_pedigree_CD,
    MODEL,
    PENALTY,
    filter_individuals_threshold=config.recombination_map.x_filter_individuals_threshold,
    chromosome_X=True,
)

position_count_AB = genetic_map_chrom_X_AB[
    ["Chromosome_position", "count"]
].copy()
position_count_CD = genetic_map_chrom_X_CD[
    ["Chromosome_position", "count"]
].copy()
total_individuals_AB = isolate_pedigree_chromosome(
    phased_regions_100kb, "AB", [CHROMOSOME_X], id_column=ID_COLUMN
)[1].shape[0]
total_individuals_CD = isolate_pedigree_chromosome(
    phased_regions_100kb, "CD", [CHROMOSOME_X], id_column=ID_COLUMN
)[1].shape[0]

combined_pedigree_recombination_map = (
    create_combined_pedigree_recombination_map(
        position_count_AB,
        position_count_CD,
        total_individuals_AB,
        total_individuals_CD,
    )
)
combined_pedigree_recombination_map.to_csv(
    config.chrx_scaffolds.x_recombination_map_output, index=False
)

# LD check for the X (data behind the published LD heatmap): combine both
# pedigrees and compute the Rogers-Huff r matrix (scikit-allel).
ld_chromosome_df = pd.concat(
    [chromX_pedigree_AB, chromX_pedigree_CD], axis=0
)
np_ld_chromosome = ld_chromosome_df.to_numpy().T
np_ld_chromosome = np.where(
    np.isnan(np_ld_chromosome), -1, np_ld_chromosome
)
ld_matrix = allel.rogers_huff_r(np_ld_chromosome)
square_matrix = squareform(ld_matrix)


# ===========================================================================
# Part 3: scaffolds (chromosomes 11-54)
# ===========================================================================

# --- Scaffold 12: segmentation map (kept as an exercise) + imputation ------
individuals_AB, scaffold_12_AB_df = isolate_pedigree_chromosome(
    phased_regions_100kb, "AB", [SEGMENTATION_SCAFFOLD], id_column=ID_COLUMN
)
individuals_CD, scaffold_12_CD_df = isolate_pedigree_chromosome(
    phased_regions_100kb, "CD", [SEGMENTATION_SCAFFOLD], id_column=ID_COLUMN
)

genetic_map_scaffold_12_AB = (
    create_recombination_map_dataframe_for_chromosome(
        scaffold_12_AB_df,
        MODEL,
        PENALTY,
        filter_individuals_threshold=config.recombination_map.scaffold12_filter_individuals_threshold_AB,
        chromosome_X=False,
    )
)
genetic_map_scaffold_12_CD = (
    create_recombination_map_dataframe_for_chromosome(
        scaffold_12_CD_df,
        MODEL,
        PENALTY,
        filter_individuals_threshold=config.recombination_map.scaffold12_filter_individuals_threshold_CD,
        chromosome_X=False,
    )
)

position_count_AB = genetic_map_scaffold_12_AB[
    ["Chromosome_position", "count"]
].copy()
position_count_CD = genetic_map_scaffold_12_CD[
    ["Chromosome_position", "count"]
].copy()
total_individuals_AB = isolate_pedigree_chromosome(
    phased_regions_100kb, "AB", [SEGMENTATION_SCAFFOLD], id_column=ID_COLUMN
)[1].shape[0]
total_individuals_CD = isolate_pedigree_chromosome(
    phased_regions_100kb, "CD", [SEGMENTATION_SCAFFOLD], id_column=ID_COLUMN
)[1].shape[0]

combined_pedigree_recombination_map_scaffold_12 = (
    create_combined_pedigree_recombination_map(
        position_count_AB,
        position_count_CD,
        total_individuals_AB,
        total_individuals_CD,
    )
)

imputed_scaffold_12_AB = impute_chromosome_df(
    scaffold_12_AB_df,
    PENALTY,
    config.imputation.scaffold12_na_threshold,
    chromosome_X=False,
)
imputed_scaffold_12_CD = impute_chromosome_df(
    scaffold_12_CD_df,
    PENALTY,
    config.imputation.scaffold12_na_threshold,
    chromosome_X=False,
)

# --- Remaining small scaffolds: modal heuristic imputation -----------------
heuristic_threshold = config.imputation.heuristic_agreement_threshold

imputed_small_scaffolds_AB = []
for scaffold_id, _, scaffold_df in create_chromosome_dataframe_list(
    phased_regions_100kb, pedigree="AB", chromosomes=SCAFFOLDS, id_column=ID_COLUMN
):
    if scaffold_id == f"Chromosome_{SEGMENTATION_SCAFFOLD}":
        continue
    imputed_scaffold = impute_chromosome_df_based_on_heuristics(
        scaffold_df, heuristic_threshold
    )
    imputed_small_scaffolds_AB.append((scaffold_id, imputed_scaffold))

imputed_small_scaffolds_CD = []
for scaffold_id, _, scaffold_df in create_chromosome_dataframe_list(
    phased_regions_100kb, pedigree="CD", chromosomes=SCAFFOLDS, id_column=ID_COLUMN
):
    if scaffold_id == f"Chromosome_{SEGMENTATION_SCAFFOLD}":
        continue
    imputed_scaffold = impute_chromosome_df_based_on_heuristics(
        scaffold_df, heuristic_threshold
    )
    imputed_small_scaffolds_CD.append((scaffold_id, imputed_scaffold))

dfs_to_concat_AB = (
    [df for _, df in imputed_small_scaffolds_AB]
    + [imputed_scaffold_12_AB]
    + [individuals_AB.to_frame()]
)
imputed_small_scaffolds_AB_df = pd.concat(dfs_to_concat_AB, axis=1)

dfs_to_concat_CD = (
    [df for _, df in imputed_small_scaffolds_CD]
    + [imputed_scaffold_12_CD]
    + [individuals_CD.to_frame()]
)
imputed_small_scaffolds_CD_df = pd.concat(dfs_to_concat_CD, axis=1)

imputed_small_scaffolds_total = pd.concat(
    [imputed_small_scaffolds_AB_df, imputed_small_scaffolds_CD_df], axis=0
)


# ===========================================================================
# Part 4: combine imputed X + scaffolds, restore the header row, write
# ===========================================================================
df_to_save = pd.concat([imputed_X_AB, imputed_X_CD], axis=0)
df_to_save = pd.concat([df_to_save, imputed_small_scaffolds_total], axis=1)

chromosome_row = df_to_save.columns.map(
    lambda col: get_chromosome_label(col, individual_label=INDIVIDUAL_LABEL)
)
chromosome_row = (
    pd.Series(chromosome_row, index=df_to_save.columns).to_frame().T
)

df_to_save = pd.concat([chromosome_row, df_to_save], ignore_index=True)
df_to_save = df_to_save.fillna("NA")

df_to_save.to_csv(
    config.chrx_scaffolds.x_scaffolds_imputed_output, index=False
)
