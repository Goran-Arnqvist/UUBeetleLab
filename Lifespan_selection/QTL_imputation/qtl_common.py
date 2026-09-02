#!/usr/bin/env python
"""
Author: Athanasios Toros
Email: athanasiostoros@gmail.com
"""
from collections import namedtuple
from enum import Enum
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import ruptures as rpt
from scipy.stats import mode


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
def _to_namespace(obj):
    if isinstance(obj, dict):
        return SimpleNamespace(
            **{k: _to_namespace(v) for k, v in obj.items()}
        )
    return obj


def load_config(path=Path(__file__).parent / "config.yaml"):
    """Load config.yaml into a nested SimpleNamespace (e.g. cfg.segmentation.penalty)."""
    import yaml

    with open(path) as fh:
        raw = yaml.safe_load(fh)
    return _to_namespace(raw)


# ---------------------------------------------------------------------------
# Genotype model
# ---------------------------------------------------------------------------
class Marker(Enum):
    EE = 0
    EL = 1
    LL = 2


class Segment_transition(Enum):
    EE_to_EL = 0
    EL_to_LL = 1
    EL_to_EE = 2
    LL_to_EL = 3
    EE_to_LL = 4
    LL_to_EE = 5
    EE_to_EE = 6
    EL_to_EL = 7
    LL_to_LL = 8


# ---------------------------------------------------------------------------
# Import / formatting of the binned-genotype matrix
# ---------------------------------------------------------------------------
chromosome_tuple = namedtuple(
    "Chromosome_df", ["chromosome_id", "individuals", "dataframe"]
)


def isolate_pedigree_chromosome(
    total_dataframe,
    pedigree=None,
    chromosomes=None,
    id_column="Individual_name",
):

    working_dataframe = total_dataframe.copy()

    # Filter chromosomes row
    chrom_mask = working_dataframe[id_column] == "chrom"
    working_dataframe = working_dataframe.loc[~chrom_mask, :]

    # isolate the individual names, and change later in case of pedigree filtering
    individuals_series = working_dataframe[id_column]

    if pedigree == "AB":
        # filter rows where the id column contains AB or BA
        pedigree_mask = working_dataframe[id_column].str.contains(
            "AB"
        ) | working_dataframe[id_column].str.contains("BA")
        working_dataframe = working_dataframe.loc[pedigree_mask, :]

        # Reassign the individuals series
        individuals_series = individuals_series[pedigree_mask]

    elif pedigree == "CD":
        # filter rows where the id column contains CD or DC
        pedigree_mask = working_dataframe[id_column].str.contains(
            "CD"
        ) | working_dataframe[id_column].str.contains("DC")
        working_dataframe = working_dataframe.loc[pedigree_mask, :]

        # Reassign the individuals series
        individuals_series = individuals_series[pedigree_mask]

    if chromosomes is not None:

        # recreate chrom column
        chrom_series = working_dataframe.columns.str.extract(r"(\d+):")[
            0
        ].astype(float)

        # create chrom mask
        chrom_mask = chrom_series.isin(chromosomes)

        if not chrom_mask.any():
            raise ValueError(
                f"No columns match the chromosomes: {chromosomes}"
            )

        # Filter the DataFrame's columns using the mask.
        # I don't know why I have to convert the mask to a numpy array.
        working_dataframe = working_dataframe.iloc[
            :, chrom_mask.to_numpy()
        ]

    return individuals_series, working_dataframe


def create_chromosome_dataframe_list(
    total_dataframe,
    pedigree=None,
    chromosomes=None,
    id_column="Individual_name",
):

    for chrom in chromosomes:
        try:
            individuals_series, chrom_df = isolate_pedigree_chromosome(
                total_dataframe,
                pedigree=pedigree,
                chromosomes=[chrom],
                id_column=id_column,
            )
            yield chromosome_tuple(
                f"Chromosome_{chrom}", individuals_series, chrom_df
            )
        except ValueError as e:
            print(e)
            continue


def bin_to_midpoint_position(bin_name):
    """
    Translate a bin name into chromosome position

    Parameters:
    bin_name (str): The name of the bin, form "chr:start-end.1"

    Returns:
    int: The chromosome position
    """

    # remove the chromosome number
    bin_name = bin_name.split(":")[1]

    # remove .1 from the end
    bin_name = bin_name.split(".")[0]

    start = int(bin_name.split("-")[0])
    end = int(bin_name.split("-")[1])

    midpoint = (start + end) // 2

    return midpoint


def get_chromosome_label(col_name, individual_label="Individual"):

    if col_name == individual_label:
        return np.nan
    return int(col_name.split(":")[0])


# ---------------------------------------------------------------------------
# Segment model (segment -> modal genotype -> transition)
#
# Implements the segmentation model described in methods.md. In the retained
# bare path the maps/imputation call filter_recombination_events / impute_row
# directly; the Segment class was exercised by the (dropped) penalty-tuning
# sweep. Kept here as documentation of the method (see NOTES.md).
# ---------------------------------------------------------------------------
class Segment:
    def __init__(
        self,
        start_index,
        end_index,
        series,
        previous_segment,
        biological_length_threshold=10,
    ):

        # assign basic properties
        # those work based on the index, not the actual position, due to filtering
        self.start = start_index
        self.end = end_index
        self.segment_array = series[start_index:end_index].to_numpy()
        self.previous_segment = previous_segment

        # calculate basic statistics
        self.length = len(self.segment_array)
        # length noise - essentially a normalized absolute integral of the discrete derivative
        self.noise = np.sum(np.abs(np.diff(self.segment_array))) / (
            self.length - 1
        )

        # check if the segment is biologically improbable
        # due to interference
        self.biologically_improbable = False
        if self.length < biological_length_threshold:
            self.biologically_improbable = True

        # set zygote
        self.set_segment_zygote()
        self.set_segment_transition()

    def set_segment_zygote(self):

        # mode is assumed as the zygote
        result_list = mode(self.segment_array, nan_policy="omit")

        self.marker = Marker(int(result_list.mode))

        # set percentage of values that agree with the mode
        self.zygote_fraction = result_list.count / self.length

    def set_segment_transition(self):
        if self.previous_segment is None:
            self.transition = None
            return

        # Warning: the transition map is backwards
        transition_map = {
            (Marker.EE, Marker.EE): Segment_transition.EE_to_EE,
            (Marker.EE, Marker.EL): Segment_transition.EL_to_EE,
            (Marker.EE, Marker.LL): Segment_transition.LL_to_EE,
            (Marker.EL, Marker.EE): Segment_transition.EE_to_EL,
            (Marker.EL, Marker.EL): Segment_transition.EL_to_EL,
            (Marker.EL, Marker.LL): Segment_transition.LL_to_EL,
            (Marker.LL, Marker.EE): Segment_transition.EE_to_LL,
            (Marker.LL, Marker.EL): Segment_transition.EL_to_LL,
            (Marker.LL, Marker.LL): Segment_transition.LL_to_LL,
        }

        # Use the map to get the transition type
        self.transition = transition_map.get(
            (self.marker, self.previous_segment.marker), None
        )

    def __str__(self):
        return f"Segment from {self.start} to {self.end} with zygote {self.zygote}"


def get_segments_chromosome(row, recombination_events):

    # extend the segments to include the start and end of the chromosome
    recombination_events = [0] + recombination_events + [len(row) - 1]

    previous_segment = None
    for start, end in zip(
        recombination_events[:-1], recombination_events[1:]
    ):
        yield (
            previous_segment := Segment(
                start, end, row, previous_segment
            )
        )


# ---------------------------------------------------------------------------
# PELT changepoint detection + recombination-event validation heuristics
# ---------------------------------------------------------------------------
def get_recombination_events_chromosome(
    row, model, penalty, jump=1, min_size=1
):
    """
    Detects recombination events in genetic data using the PELT (Pruned Exact Linear Time) changepoint detection algorithm.

    Parameters:
    ----------
    row : pandas.Series
        A series containing genetic markers with possible NaN values.
    model : str
        The cost model to use for changepoint detection (e.g., 'l2', 'rbf', etc.).
    penalty : float
        The penalty value for the PELT algorithm. Controls the trade-off between
        the number of changepoints and the fit to the data.

    Returns:
    -------
    list
        Integer positions of detected recombination events in the original array
        (before NaN removal).

    Notes:
    -----
    - NaN values are removed before processing but their positions are tracked
      to return results in the original coordinate system
    - Returned positions correspond to the last position before a change occurs
    """
    # Drop NaN values but record their integer positions from the original index.
    non_na = row.dropna()
    positions = row.index.get_indexer(non_na.index)

    chromosome_array = non_na.to_numpy()

    # Get changepoints using ruptures
    # min_size == 1 by default, filter out later
    algo = rpt.Pelt(model=model, jump=jump, min_size=min_size).fit(
        chromosome_array
    )

    result = algo.predict(pen=penalty)

    # Note: p from result is 1-indexed
    # remove the last result to extend the segment to the end of the chromosome
    # Warning: this will return an empty list if there are no changepoints, should not be a problem
    events = [positions[p - 1] for p in result[:-1]]

    return events


def filter_recombination_events(
    recombination_indices,
    chromosome_row,
    min_distance_bp=1_000_000,
    chromosome_X=False,
):
    """
    Filter and validate recombination events in a chromosome row based on distance and segment modes.

    This function performs two-stage filtering of recombination events:
    1. Distance-based filtering: Removes events that are too close to each other, due to interference considerations
    2. Mode-based filtering: Validates events based on the mode of chromosome segments, duplicating double transitions

    Parameters
    ----------
    recombination_indices : list or array-like
        Indices where potential recombination events occur
    chromosome_row : pandas.Series
        A row from a chromosome matrix containing marker values
    min_distance : int, optional
        Minimum required distance between consecutive recombination events, by default 4

    Returns
    -------
    list
        Filtered indices where valid recombination events occur

    Notes
    -----
    - The function assumes the chromosome_row can be converted to a numpy array
    - For mode-based filtering, transitions involving 'EL' marker are treated as single events,
      while transitions between other markers are treated as double events
    - The first segment's mode is used as a reference for subsequent comparisons
    - The final indices are 0 indexed, corresponding to the original chromosome_row
    """

    # Assuming Pandas Series
    np_chromosome_row = chromosome_row.to_numpy()

    # Initialize filtering
    # turn the recombination indices into column names
    recombination_columns = chromosome_row.index[recombination_indices]
    # set as actual positions

    recombination_chromosome_positions = list(
        map(
            bin_to_midpoint_position,
            recombination_columns.to_list(),
        )
    )

    # WARNING: The last_selected approach will essentially make it so that
    # if we have multiple bad events in a row, only the middle ones will be removed
    filtered_indices = []
    # this selects the first index
    last_selected_position = -np.inf

    # Iterate and filter
    for index, actual_position in zip(
        recombination_indices,
        recombination_chromosome_positions,
    ):
        if actual_position - last_selected_position > min_distance_bp:
            filtered_indices.append(index)
            last_selected_position = actual_position

    # 2. Selection based on the mode of the segment, including transitions

    # create two lists for iterating the segments
    recombination_indices_plus_start = [0] + filtered_indices
    recombination_indices_plus_end = filtered_indices + [
        len(np_chromosome_row) - 1
    ]

    final_indices = []
    curr_mode = None
    last_mode = None

    for start, end in zip(
        recombination_indices_plus_start,
        recombination_indices_plus_end,
    ):
        segment = np_chromosome_row[
            start : end + 1
        ]  # +1 to include the end, check if problematic
        curr_mode = mode(segment, nan_policy="omit").mode

        # deal with first segment
        if last_mode is None:
            last_mode = curr_mode
            continue

        # Now let's filter concecutive segments with the same mode
        if curr_mode == last_mode:
            last_mode = curr_mode
            continue

        # Only for X!!!
        # If we assume that for F2 there is no recombination in the F1 male,
        # then only transitions with heterozygocity are valid
        # i.e. the male haplotype is E or L, in its entirety
        if chromosome_X and (
            Marker.EL.value
            not in {
                curr_mode,
                last_mode,
            }
        ):
            print("Warning: invalid transition for X")
            continue

        # if we reach this point, we have a valid recombinant segment

        # Now, since we know that the two modes are different,
        # we can work out the transition
        if Marker.EL.value in {
            curr_mode,
            last_mode,
        }:
            # here this might be wrong, but we work with the start of the segment
            # otherwise we would get a recombination at the end of the chroosome
            final_indices.append(start)
        else:
            # essentially, if none of the different modes is EL, we have a double transition
            final_indices.extend([start, start])

        last_mode = curr_mode

    return final_indices


def get_valid_recombination_events_chromosome(
    row,
    model,
    penalty,
    jump=1,
    min_size=1,
    min_distance_bp=1_000_000,
    chromosome_X=False,
):

    # Get recombination events
    recombination_events = get_recombination_events_chromosome(
        row,
        model,
        penalty=penalty,
        jump=jump,
        min_size=min_size,
    )

    # Filter the recombination events
    filtered_events = filter_recombination_events(
        recombination_events,
        row,
        min_distance_bp=min_distance_bp,
        chromosome_X=chromosome_X,
    )

    return filtered_events


def find_all_valid_recombination_events(
    chromosome_df,
    model,
    penalty,
    jump=1,
    min_size=1,
    min_distance_bp=1_000_000,
    chromosome_X=False,
):
    """
    Find all valid recombination events across individuals in a chromosome DataFrame.

    This function processes each individual in the chromosome data and identifies valid
    recombination events based on the specified model and parameters.

    Parameters
    ----------
    chromosome_df : pandas.DataFrame
        DataFrame containing chromosome data for multiple individuals.
    model : str
        Model used for identifying recombination events.
    penalty : float
        Penalty parameter for the recombination event detection.
    jump : int, optional
        Step size for scanning through positions. Default is 1.
    min_size : int, optional
        Minimum size requirement for a valid recombination event. Default is 1.
    min_distance_bp : int, optional
        Minimum distance in base pairs between recombination events. Default is 1,000,000.

    Returns
    -------
    list
        A list of indices where valid recombination events were detected across all individuals.

    Notes
    -----
    The function processes each individual separately and combines all detected
    recombination events into a single list.
    """

    # List to store total recombination indices
    recombination_indices = []

    # Iterate over each row in the DataFrame
    for i, row in chromosome_df.iterrows():
        # Get recombination events for the current individual
        recombination_events = (
            get_valid_recombination_events_chromosome(
                row,
                model,
                penalty=penalty,
                jump=jump,
                min_size=min_size,
                min_distance_bp=min_distance_bp,
                chromosome_X=chromosome_X,
            )
        )
        recombination_indices.extend(recombination_events)

    return recombination_indices


def recombination_counts_per_region(
    chromosome_df, validated_recombination_events
):
    """
    Calculate the recombination counts for each genomic region based on validated recombination events.
    This function takes a DataFrame containing chromosome data and a list of validated
    recombination event indices, and returns counts of recombination events per genomic region
    along with their corresponding chromosome positions.
    Parameters
    ----------
    chromosome_df : pandas.DataFrame
        DataFrame containing chromosome data where columns represent genomic regions
    validated_recombination_events : list
        List of indices corresponding to validated recombination events in the chromosome_df
    Returns
    -------
    pandas.DataFrame
        DataFrame containing:
            - count: Number of recombination events per region
            - Original_genomic_regions: Original genomic region identifiers
            - Chromosome_position: Midpoint position of each genomic region on the chromosome
    Notes
    -----
    The function assumes the existence of a `bin_to_midpoint_position` helper function
    to convert genomic regions to chromosome positions.
    """

    # Turn the indices into column names
    list_of_indices_columns = []
    for recombination_event in validated_recombination_events:
        list_of_indices_columns.append(
            chromosome_df.columns[recombination_event]
        )

    # Count the occurrences of each region
    series_of_indices_columns = pd.Series(list_of_indices_columns)
    recombination_counts = (
        series_of_indices_columns.value_counts().to_frame()
    )

    # Add the original genomic regions and chromosome positions
    recombination_counts["Original_genomic_regions"] = (
        recombination_counts.index
    )

    # Transform into real chromosome positions, and sort
    # This assumes the midpoint as the true position
    recombination_counts["Chromosome_position"] = recombination_counts[
        "Original_genomic_regions"
    ].apply(bin_to_midpoint_position)

    recombination_counts.sort_values(
        "Chromosome_position", inplace=True
    )
    recombination_counts = recombination_counts.reset_index(drop=True)

    return recombination_counts


# ---------------------------------------------------------------------------
# Genetic-map construction (Kosambi -> cumulative cM -> interpolation)
# ---------------------------------------------------------------------------
def kosambi_method(recombination_fraction):
    """
    Calculate genetic map distance using Kosambi's mapping function.

    This function converts recombination fraction to genetic map distance in Morgans
    using Kosambi's mapping function, which accounts for interference between crossovers.

    Parameters
    ----------
    recombination_fraction : float
        The recombination fraction between two loci, must be between 0 and 0.5.

    Returns
    -------
    float
        The genetic map distance in Morgans.

    Raises
    ------
    AssertionError
        If recombination_fraction is not between 0 and 0.5.

    Notes
    -----
    Kosambi's mapping function assumes some interference between crossovers and is
    defined as:
        d = 1/4 * ln((1 + 2r)/(1 - 2r))
    where d is the distance in Morgans and r is the recombination fraction.
    """

    # handle NaN values, usefull in combine_genetic_maps
    if np.isnan(recombination_fraction):
        return np.nan

    assert (
        recombination_fraction >= 0 and recombination_fraction < 0.5
    ), f"Invalid recombination fraction {recombination_fraction}, must be between 0 and 0.5."

    return 0.25 * (
        np.log(
            (1 + 2 * recombination_fraction)
            / (1 - 2 * recombination_fraction)
        )
    )


def calculate_genetic_map(
    recombination_counts, nr_individuals, chromosome_X=False
):
    """
    Calculate genetic map distances and recombination rates from recombination counts.

    This function takes recombination counts between markers and calculates genetic distances
    using Kosambi's mapping function, as well as various recombination metrics.

    Parameters
    ----------
    recombination_counts : pandas.DataFrame
        DataFrame containing columns:
        - count: Number of recombination events between adjacent markers
        - Original_genomic_regions: Marker identifiers
        - Chromosome_position: Physical position of markers on chromosome

    nr_individuals : int
        Number of individuals in the mapping population

    Returns
    -------
    pandas.DataFrame
        Input DataFrame with additional columns:
        - recombination_fraction: Proportion of recombinant gametes
        - genetic_map_distance: Map distance in Morgans
        - genetic_map_distance_cM: Map distance in centiMorgans
        - cumulative_genetic_map_distance_cM: Cumulative map distance in cM
        - recombination_rate_cM_Mb: Local recombination rate in cM/Mb

    Notes
    -----
    The recombination rate calculation for the first marker will result in NaN
    due to the position difference calculation.
    """

    working_df = recombination_counts.copy()
    # recombination counts is a dataframe with columns count, Original_genomic_regions, Chromosome_position

    # for each marker, the recombination fraction is the number of recombination events / total number of gametes
    # this is the same as the number of recombination events / (2 * number of individuals)
    # for chromosome X, the pool of possible recombinants is halved
    # (only the maternal X recombines: the hemizygous F1 father transmits a
    # single, non-recombinant X, so each F2 contributes one informative gamete)
    total_gametes = (
        nr_individuals if chromosome_X else 2 * nr_individuals
    )

    working_df["recombination_fraction"] = (
        working_df["count"] / total_gametes
    )

    # calculate the genetic map distance between each region

    working_df["genetic_map_distance"] = working_df[
        "recombination_fraction"
    ].apply(kosambi_method)

    # turn into cM
    working_df["genetic_map_distance_cM"] = (
        working_df["genetic_map_distance"] * 100
    )

    # calculate the cumulative genetic map distance
    working_df["cumulative_genetic_map_distance_cM"] = working_df[
        "genetic_map_distance_cM"
    ].cumsum()

    # calculate the rate of recombination, i.e. the  cM per Mbp

    distance_between_markers = working_df["Chromosome_position"].diff()

    # WARNING: The first value will be NaN, may lead to issues

    working_df["recombination_rate_cM_Mb"] = working_df[
        "genetic_map_distance_cM"
    ] / (
        distance_between_markers / 1e6
    )  # distance in Mbp

    return working_df


def find_closest_position(row, complete_positions_df):
    """
    Finds the index of the closest position in a reference dataframe to a given target position.

    Parameters
    ----------
    row : pandas.Series
        A row containing a 'Chromosome_position' value to find the closest match for
    complete_positions_df : pandas.DataFrame
        Reference dataframe containing 'Chromosome_position' column to search in

    Returns
    -------
    int
        Index of the position in complete_positions_df that is closest to the target position

    Notes
    -----
    Uses numpy's argmin to find the index with minimum absolute difference between positions
    """

    target_position = row["Chromosome_position"]

    # Calculate the absolute differences and find the index of the minimum difference
    idx = np.abs(
        complete_positions_df["Chromosome_position"] - target_position
    ).argmin()

    return idx


def interpolate_genetic_distances(
    recombination_counts_per_region,
    original_chromosome_df,
    step_size=100_000,
):

    # Step 1: get the original starting and ending positions

    original_columns = list(original_chromosome_df.columns)
    original_positions = list(
        map(bin_to_midpoint_position, original_columns)
    )

    chromosome_min_position = min(original_positions)
    chromosome_max_position = max(original_positions)

    complete_positions = np.arange(
        chromosome_min_position,
        chromosome_max_position + step_size,
        step_size,
    )
    # Step 2: Create a new dataframe for all of the regions, with step_size
    complete_positions_df = pd.DataFrame(
        {"Chromosome_position": complete_positions}
    )

    # Step 3: Find the closest position in the original dataframe for each position in the new dataframe
    for index, row in recombination_counts_per_region.iterrows():

        idx = find_closest_position(row, complete_positions_df)

        complete_positions_df.loc[
            idx, "cumulative_genetic_map_distance_cM"
        ] = row["cumulative_genetic_map_distance_cM"]

        # Extra, add the count column
        complete_positions_df.loc[idx, "count"] = row["count"]

    # Step 4: Remove any leading and trailing NaN values

    complete_positions_df = complete_positions_df.loc[
        complete_positions_df[
            "cumulative_genetic_map_distance_cM"
        ].first_valid_index() :
    ]
    complete_positions_df = complete_positions_df.loc[
        : complete_positions_df[
            "cumulative_genetic_map_distance_cM"
        ].last_valid_index()
    ]

    # Step 5: Interpolate the missing values
    complete_positions_df[
        "cumulative_genetic_map_distance_cM_interpolated"
    ] = complete_positions_df[
        "cumulative_genetic_map_distance_cM"
    ].interpolate(
        method="linear"
    )

    # Step 6: Calculate the delta cM and delta position
    complete_positions_df["delta_cM"] = complete_positions_df[
        "cumulative_genetic_map_distance_cM_interpolated"
    ].diff()
    complete_positions_df["delta_position"] = complete_positions_df[
        "Chromosome_position"
    ].diff()

    # Step 7: Calculate the recombination rate
    complete_positions_df["recombination_rate_cM_Mb"] = (
        complete_positions_df["delta_cM"]
        / (complete_positions_df["delta_position"] / 1e6)
    )

    return complete_positions_df


def filter_individuals(chromosome_df, max_missing=0.4):
    """
    Filter individuals based on the proportion of missing values in the chromosome data.

    This function filters individuals based on the proportion of missing values in the chromosome data.
    Individuals with a proportion of missing values greater than the threshold are removed.

    Parameters
    ----------
    chromosome_df : pandas.DataFrame
        DataFrame containing chromosome data for multiple individuals.
    max_missing : float, optional
        Maximum proportion of missing values allowed for an individual, by default 0.1.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame containing individuals with missing values below the threshold.

    Notes
    -----
    The function calculates the proportion of missing values for each individual and filters
    out individuals with a proportion greater than the threshold.
    """

    # Calculate the proportion of missing values for each individual
    missing_values = chromosome_df.isna().mean(axis=1)

    # Filter out individuals with missing values above the threshold
    filtered_df = chromosome_df[missing_values <= max_missing]

    print(
        f"Warning: {((len(chromosome_df) - len(filtered_df)) / len(chromosome_df)):.2%}% individuals removed based on unassinged threshold of {max_missing:.2%}%"
    )

    return filtered_df


def create_recombination_map_dataframe_for_chromosome(
    chromosome_pedigree_df,
    model,
    penalty,
    jump=1,
    min_size=1,
    filter_individuals_threshold=0.4,
    min_distance_bp=1_000_000,
    chromosome_X=False,
):

    # Step 1: Filter the dataframe for individuals with too many missing values
    chromosome_pedigree_df_filtered = filter_individuals(
        chromosome_pedigree_df,
        max_missing=filter_individuals_threshold,
    )

    # Step 1: Find all valid recombination events
    recombination_indices = find_all_valid_recombination_events(
        chromosome_pedigree_df_filtered,
        model,
        penalty,
        jump=jump,
        min_size=min_size,
        chromosome_X=chromosome_X,
    )

    # Step 2: Count the recombination events per region
    recombination_counts = recombination_counts_per_region(
        chromosome_pedigree_df_filtered, recombination_indices
    )

    # Step 3: Calculate the genetic map distances
    genetic_map = calculate_genetic_map(
        recombination_counts,
        chromosome_pedigree_df_filtered.shape[0],
        chromosome_X=chromosome_X,
    )

    # Step 4: Interpolate the genetic distances
    interpolated_genetic_map = interpolate_genetic_distances(
        genetic_map, chromosome_pedigree_df_filtered
    )

    return interpolated_genetic_map


# ---------------------------------------------------------------------------
# Imputation (assign each segment to its modal genotype)
# ---------------------------------------------------------------------------
def impute_row(valid_recombination_events, row):
    """Imputes values in a row based on segment-wise mode calculation using recombination events.

    This function modifies the input row directly by setting all values within each segment
    to the mode value of that segment. Segments are defined by the recombination event positions.

    Parameters
    ----------
    valid_recombination_events : list
        List of positions where recombination events occur. These positions define
        the boundaries of segments for imputation.
    row : pandas.Series
        Input row containing values to be imputed. The row is modified in-place.

    Notes
    -----
    - The function operates directly on the input row (in-place modification)
    - Recombination events are assumed to be 0-indexed
    - Each segment's values are replaced with the mode of that segment
    """
    # remove duplicates from valid recombination events
    # this happens because the original code sets double transitions to the same position
    # but we don't care for imputation purposes
    valid_recombination_events = np.unique(
        valid_recombination_events
    ).tolist()
    working_row = row.copy()

    # create arrays for segment iteration
    segment_array_with_start = [0] + valid_recombination_events
    segment_array_with_end = valid_recombination_events + [
        working_row.shape[0] - 1
    ]

    for segment_start, segment_end in zip(
        segment_array_with_start,
        segment_array_with_end,
    ):

        # get mode of the segment
        segment_mode = working_row.iloc[
            segment_start:segment_end
        ].mode()

        if segment_mode.shape[0] == 0:
            print(
                f"Warning: No mode found, segment start: {segment_start}, segment end: {segment_end}"
            )
            continue

        segment_mode = segment_mode.iloc[0]

        # set the whole segment to the mode
        # results from ruptures are 1-indexed, but get_recobination_events is 0-indexed
        # since it subtracts 1
        working_row.iloc[segment_start:segment_end] = segment_mode

    return working_row


def run_imputation_per_chromosome_individual(
    row, na_threshold, penalty, chromosome_X=False
):
    """
    Performs imputation on a single chromosome for an individual based on recombination events.
    This function processes a single chromosome's genotype data, identifies valid recombination
    events, and performs imputation if the missing data doesn't exceed the threshold.
    Parameters
    ----------
    row : pandas.Series
        A series containing genotype data for a single chromosome
    na_threshold : float
        Maximum allowable proportion of missing values (between 0 and 1)
    penalty : float
        Penalty parameter for the recombination event detection algorithm
    Returns
    -------
    None
        The function modifies the input row in-place through imputation
    Notes
    -----
    The function will skip imputation if the proportion of missing values exceeds
    the na_threshold. Imputation is performed using valid recombination events
    detected with an L1 penalty approach.
    """
    working_row = row.copy()

    # skip imputation if chromosomome has too many missing values
    if working_row.isna().sum() / working_row.shape[0] > na_threshold:
        return working_row

    valid_recombination_events = (
        get_valid_recombination_events_chromosome(
            working_row,
            "l1",
            penalty=penalty,
            chromosome_X=chromosome_X,
        )
    )

    # impute the row
    final_row = impute_row(valid_recombination_events, working_row)

    return final_row


def impute_chromosome_df(
    chromosome_df, penalty, na_threshold, chromosome_X=False
):

    working_df = chromosome_df.copy()
    for idx, row in chromosome_df.iterrows():
        imputed_row = run_imputation_per_chromosome_individual(
            row, na_threshold, penalty, chromosome_X=chromosome_X
        )

        working_df.loc[idx] = imputed_row

    return working_df
