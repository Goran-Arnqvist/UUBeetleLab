#!/usr/bin/env python
"""
Author: Athanasios Toros
Email: athanasiostoros@gmail.com
"""

from concurrent.futures import ProcessPoolExecutor
from functools import partial

import allel
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform

from qtl_common import (
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
AUTOSOMES = config.genome.autosomes
ID_COLUMN = config.autosomes.id_column
INDIVIDUAL_LABEL = config.autosomes.individual_label
MAX_WORKERS = config.recombination_map.max_workers

phased_regions_100kb = pd.read_csv(config.autosomes.input_calls)


# ---------------------------------------------------------------------------
# Recombination-map runners
# ---------------------------------------------------------------------------
def helper_function(chromosome_df_tuple, filter_individuals_threshold=0.4):
    chrom, chromosome_df = chromosome_df_tuple
    try:
        tmp = create_recombination_map_dataframe_for_chromosome(
            chromosome_df,
            MODEL,
            PENALTY,
            filter_individuals_threshold=filter_individuals_threshold,
        )

        print(f"Finished processing for chromosome {chrom}")
        return chrom, tmp

    except Exception as e:
        print(f"Failed for chromosome {chrom}")
        print(e)
        return chrom, None


def process_all_chromosomes_pedigree(
    chromosome_dfs, filter_individuals_threshold=0.4
):

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        results = executor.map(
            helper_function,
            chromosome_dfs,
            [filter_individuals_threshold] * len(chromosome_dfs),
        )

    return list(results)


def create_combined_pedigree_recombination_map(
    position_count_AB_df,
    position_count_CD_df,
    total_individuals_AB,
    total_individuals_CD,
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

    total_gametes = (total_individuals_AB + total_individuals_CD) * 2
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


def calculate_recombination_fraction(linkage_disequillibrium):

    return (
        (1 - np.sqrt(linkage_disequillibrium))
        / 2
        * (1 + np.sqrt(linkage_disequillibrium))
    )


# ---------------------------------------------------------------------------
# Imputation runner
# ---------------------------------------------------------------------------
def helper_impute_all_chromosomes(
    chrom_id, pedigree, penalty, na_threshhold
):

    individuals, working_df = isolate_pedigree_chromosome(
        phased_regions_100kb, pedigree, [chrom_id], id_column=ID_COLUMN
    )
    working_df = working_df.copy()

    imputed_df = impute_chromosome_df(working_df, penalty, na_threshhold)
    print(f"Imputed chromosome {chrom_id}")
    return (chrom_id, individuals, imputed_df)


# ===========================================================================
# Part 1: recombination maps
# ===========================================================================

# --- Pedigree AB -----------------------------------------------------------
chromosome_dfs = []
for chrom in AUTOSOMES:
    _, tmp = isolate_pedigree_chromosome(
        phased_regions_100kb, "AB", [chrom], id_column=ID_COLUMN
    )
    chromosome_dfs.append((chrom, tmp))

pedigree_AB_genetic_maps = process_all_chromosomes_pedigree(
    chromosome_dfs,
    filter_individuals_threshold=config.recombination_map.autosome_filter_individuals_threshold_AB,
)

list_dfs_to_concat = []
for chrom, chrom_df in pedigree_AB_genetic_maps:
    tmp_df = chrom_df.copy()
    tmp_df["Chromosome"] = chrom
    list_dfs_to_concat.append(tmp_df)

all_genetic_maps_AB = pd.concat(list_dfs_to_concat, ignore_index=True)

# --- Pedigree CD -----------------------------------------------------------
chromosome_dfs = []
for chrom in AUTOSOMES:
    _, tmp = isolate_pedigree_chromosome(
        phased_regions_100kb, "CD", [chrom], id_column=ID_COLUMN
    )
    chromosome_dfs.append((chrom, tmp))

pedigree_CD_genetic_maps = process_all_chromosomes_pedigree(
    chromosome_dfs,
    filter_individuals_threshold=config.recombination_map.autosome_filter_individuals_threshold_CD,
)

list_dfs_to_concat_CD = []
for chrom, chrom_df in pedigree_CD_genetic_maps:
    if chrom_df is None:
        continue
    tmp_df = chrom_df.copy()
    tmp_df["Chromosome"] = chrom
    list_dfs_to_concat_CD.append(tmp_df)

all_genetic_maps_CD = pd.concat(list_dfs_to_concat_CD, ignore_index=True)

# --- Combine pedigrees into one map per chromosome -------------------------
# AB and CD belong to distinct populations, but a genetic map is a species-level
# property, so the evidence is combined from the raw per-region counts with a
# correction for the doubled gamete pool.
total_individuals_AB = isolate_pedigree_chromosome(
    phased_regions_100kb, "AB", [1], id_column=ID_COLUMN
)[1].shape[0]
total_individuals_CD = isolate_pedigree_chromosome(
    phased_regions_100kb, "CD", [1], id_column=ID_COLUMN
)[1].shape[0]

combined_genetic_maps_all_chromosomes = pd.DataFrame()
for chrom in AUTOSOMES:
    tmp_AB = all_genetic_maps_AB[all_genetic_maps_AB["Chromosome"] == chrom]
    tmp_CD = all_genetic_maps_CD[all_genetic_maps_CD["Chromosome"] == chrom]

    position_count_AB = tmp_AB[["Chromosome_position", "count"]].copy()
    position_count_CD = tmp_CD[["Chromosome_position", "count"]].copy()

    combined_genetic_map = create_combined_pedigree_recombination_map(
        position_count_AB,
        position_count_CD,
        total_individuals_AB,
        total_individuals_CD,
    )

    combined_genetic_map["Chromosome"] = chrom

    combined_genetic_maps_all_chromosomes = pd.concat(
        [combined_genetic_maps_all_chromosomes, combined_genetic_map],
        ignore_index=True,
    )

combined_genetic_maps_all_chromosomes.to_csv(
    config.autosomes.recombination_map_output, index=False
)

# --- LD validation (representative chromosome) -----------------------------
# Compare the recombination fraction implied by direct haplotype phasing with
# the one inferred from linkage disequilibrium (scikit-allel), as a check that
# the mapping assumptions fit the observed data.
ld_chromosome_df = isolate_pedigree_chromosome(
    phased_regions_100kb, "AB", [1], id_column=ID_COLUMN
)[1]

np_ld_chromosome = ld_chromosome_df.to_numpy().T
np_ld_chromosome = np.where(
    np.isnan(np_ld_chromosome), -1, np_ld_chromosome
)

ld_matrix = allel.rogers_huff_r(np_ld_chromosome)
square_matrix = squareform(ld_matrix)

vec_calculate_recombination_fraction = np.vectorize(
    calculate_recombination_fraction
)
recombination_fraction_matrix = vec_calculate_recombination_fraction(
    square_matrix
)


# ===========================================================================
# Part 2: imputation
# ===========================================================================

# --- Pedigree AB -----------------------------------------------------------
helper_partial = partial(
    helper_impute_all_chromosomes,
    pedigree="AB",
    penalty=PENALTY,
    na_threshhold=NA_THRESHOLD,
)

with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
    results = executor.map(helper_partial, AUTOSOMES)

results_final_AB = list(results)

dfs_to_concat = []
individuals_AB = None
for chrom_id, individuals, imputed_df in results_final_AB:
    if individuals_AB is None:
        individuals_AB = individuals
    dfs_to_concat.append(imputed_df)

imputed_pedigree_AB = pd.concat(dfs_to_concat, axis=1)
imputed_pedigree_AB[INDIVIDUAL_LABEL] = individuals_AB

# --- Pedigree CD -----------------------------------------------------------
helper_partial = partial(
    helper_impute_all_chromosomes,
    pedigree="CD",
    penalty=PENALTY,
    na_threshhold=NA_THRESHOLD,
)

with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
    results = executor.map(helper_partial, AUTOSOMES)

results_final_CD = list(results)

indificuals_CD = None
dfs_to_concat = []
for chrom_id, individuals, imputed_df in results_final_CD:
    if indificuals_CD is None:
        indificuals_CD = individuals
    dfs_to_concat.append(imputed_df)

imputed_pedigree_CD = pd.concat(dfs_to_concat, axis=1)
imputed_pedigree_CD[INDIVIDUAL_LABEL] = indificuals_CD

# --- Combine pedigrees, restore the chromosome header row, write -----------
imputed_pedigree = pd.concat(
    [imputed_pedigree_AB, imputed_pedigree_CD], axis=0
)

chromosome_row = imputed_pedigree.columns.map(
    lambda col: get_chromosome_label(col, individual_label=INDIVIDUAL_LABEL)
)
chromosome_row = (
    pd.Series(chromosome_row, index=imputed_pedigree.columns).to_frame().T
)

imputed_pedigree = pd.concat(
    [chromosome_row, imputed_pedigree], ignore_index=True
)

imputed_pedigree.to_excel(config.autosomes.imputed_output, index=False)
