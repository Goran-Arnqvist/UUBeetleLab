#!/usr/bin/env python
"""
Author: Athanasios Toros
Email: athanasiostoros@gmail.com
"""

import re
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import polars as pl
import gffutils
from pybedtools import BedTool
from goatools.obo_parser import GODag


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
def _to_namespace(obj):
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: _to_namespace(v) for k, v in obj.items()})
    return obj


def load_config(path=Path(__file__).parent / "config.yaml"):
    """Load config.yaml into a nested SimpleNamespace (e.g. cfg.inputs.interpro2go)."""
    import yaml

    with open(path) as fh:
        raw = yaml.safe_load(fh)
    return _to_namespace(raw)


cfg = load_config()


def _write_csv(df, path, **kwargs):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(path, **kwargs)


def build_obsolete_remap(go_ids, obo_path):
    """Map obsolete GO terms found in `go_ids` to their successor(s).

    Returns (remap, completely_obsolete): `remap` sends each obsolete term to
    its successor(s); `completely_obsolete` lists terms with no successor (to be
    dropped). Run once for the InterPro-only set and once for the merged set.
    """
    godag = GODag(
        obo_path,
        load_obsolete=True,
        optional_attrs={"consider", "replaced_by"},
    )
    remap = {}
    completely_obsolete = []
    stats = {
        "go_terms_with_successors": 0,
        "go_terms_without_successors": 0,
        "completely_obsolete": 0,
    }
    for go_id, term in godag.items():
        if go_id not in go_ids:
            continue
        if term.is_obsolete:
            # 'replace_by' returns a list of successor IDs, from oboInOwl
            successors = term.replaced_by  # GOATOOLS support
            if successors:
                remap[go_id] = successors
                stats["go_terms_with_successors"] += 1
            elif term.consider:
                remap[go_id] = list(term.consider)[0]
                stats["go_terms_without_successors"] += 1
            else:
                completely_obsolete.append(go_id)
                stats["completely_obsolete"] += 1

    print("Total GO terms with successors:", stats["go_terms_with_successors"])
    print("Total GO terms without successors:", stats["go_terms_without_successors"])
    print("Total completely obsolete GO terms:", stats["completely_obsolete"])
    return remap, completely_obsolete


# ---------------------------------------------------------------------------
# InterPro to GO from the NCBI GFF
# ---------------------------------------------------------------------------

db = gffutils.FeatureDB(cfg.inputs.gff_feature_db)


# pull all gene features & their GO tags
def interpro_from_attr(attr):
    attr = attr[0]
    # InterPro IDs are in the form of "InterPro:IPR000000"
    m = re.findall(r"InterPro:IPR\d+", attr)
    return (
        {interpro_id.split(":")[1] for interpro_id in m}
        if m
        else set()
    )


cds_records_go_dict = {}
for gene in db.features_of_type("CDS"):

    gos = interpro_from_attr(gene.attributes.get("inference", [""]))
    protein_id = None
    if "protein_id" in gene.attributes:
        protein_id = gene.attributes["protein_id"][0]
    else:
        protein_id = gene.id

    cds_records_go_dict.setdefault(protein_id, set()).update(gos)

gff_df = pl.DataFrame(
    {
        "gene_id": list(cds_records_go_dict.keys()),
        "interpro_id": list(
            [
                list(set_val)
                for set_val in cds_records_go_dict.values()
            ]
        ),
    }
)
gff_df = gff_df.explode("interpro_id")

# The annotation only carries the InterPro IDs, so the associated GO terms are
# pulled from interpro2go manually.
interpro2go_dict = {}
with open(cfg.inputs.interpro2go, "r") as f:
    for line in f:
        # Split the line into InterPro ID and GO ID
        interpro_id = line.split(" ")[0]
        go_id = line.split(";")[-1].strip()

        interpro2go_dict.setdefault(interpro_id, set()).add(go_id)

interpro2go_dict_lists = {
    key: list(value) for key, value in interpro2go_dict.items()
}
interpro2go_df = pl.DataFrame(
    {
        "interpro_id": list(interpro2go_dict.keys()),
        "go_id": list(interpro2go_dict_lists.values()),
    }
)
interpro2go_df = interpro2go_df.explode("go_id")
interpro2go_df = interpro2go_df.select(
    pl.col("interpro_id").str.replace_all(r"InterPro:", ""),
    pl.col("go_id"),
)

gff_df_with_go = gff_df.join(
    interpro2go_df,
    on="interpro_id",
    how="left",
)

stats = {}
stats["total_genes"] = gff_df_with_go.shape[0]
stats["genes_without_interpro"] = gff_df_with_go.filter(
    pl.col("interpro_id").is_null()
).shape[0]
stats["genes_without_go"] = gff_df_with_go.filter(
    pl.col("go_id").is_null()
).shape[0]

stats["genes_with_interpro"] = (
    stats["total_genes"] - stats["genes_without_interpro"]
)
stats["genes_with_go"] = (
    stats["total_genes"] - stats["genes_without_go"]
)

print("Total genes with InterPro IDs:", stats["genes_with_interpro"])
print("Total genes with GO IDs:", stats["genes_with_go"])
print(
    "Percentage of genes with InterPro IDs, and not GO IDs:",
    (stats["genes_with_interpro"] - stats["genes_with_go"])
    / stats["total_genes"]
    * 100,
)

gff_df_with_go = (
    gff_df_with_go.select(
        pl.col("gene_id"),
        pl.col("go_id"),
    )
    .filter(
        pl.col("go_id").is_not_null(),
    )
    .unique()
)
# Drop the cds-prefixed ids (not protein_ids).
gff_df_with_go = gff_df_with_go.filter(
    ~pl.col("gene_id").str.contains("cds")
)

_write_csv(
    gff_df_with_go,
    cfg.outputs.interpro_only_gene2go,
    include_header=False,
    separator="\t",
)


# ---------------------------------------------------------------------------
# CDS regions BED -> (chrom, start, end, protein_id)
# ---------------------------------------------------------------------------
# the gff file was turned into a BED with only the CDS regions
bed_cds_regions_df = pl.read_csv(
    cfg.inputs.cds_regions_bed,
    separator="\t",
    has_header=False,
    new_columns=[
        "chrom",
        "start",
        "end",
        "col_4",
        "col_5",
        "col_6",
        "EMBL",
        "CDS",
        "no_idea",
        "gff_attributes",
    ],
)

bed_cds_fixed = bed_cds_regions_df.select(
    pl.col("chrom"),
    pl.col("start"),
    pl.col("end"),
    pl.col("gff_attributes")
    .str.extract(r"protein_id=([^;]+)")
    .alias("gene_id"),
)
# The genes with no protein_id are all PSEUDO; drop them.
bed_cds_fixed = bed_cds_fixed.filter(~pl.col("gene_id").is_null())

_write_csv(
    bed_cds_fixed,
    cfg.outputs.cds_regions_fixed_bed,
    separator="\t",
    include_header=False,
)


# ---------------------------------------------------------------------------
# Obsolete-term remap for the InterPro-only GO set
# ---------------------------------------------------------------------------
# Some terms from the first InterProScan are obsolete; remap to the current ones.
go_ids_for_genes = set(gff_df_with_go["go_id"].unique().to_list())

remap, completely_obsolete = build_obsolete_remap(
    go_ids_for_genes, cfg.inputs.go_basic_obo
)

gff_df_with_go = gff_df_with_go.select(
    pl.col("gene_id"), pl.col("go_id").replace(remap)
).filter(~pl.col("go_id").is_in(completely_obsolete))

# Write the updated DataFrame to a new file
_write_csv(
    gff_df_with_go,
    cfg.outputs.interpro_only_gene2go_updated,
    include_header=False,
    separator="\t",
)


# ---------------------------------------------------------------------------
# EggNOG-mapper GO universe (RefSeq CDS)
# ---------------------------------------------------------------------------
eggnog_df = pl.read_csv(
    cfg.inputs.eggnog_refseq_annotations,
    separator="\t",
    has_header=True,
)

eggnog_go_universe = eggnog_df.select(
    pl.col("query"), pl.col("GOs").str.split(",").alias("go_ids")
)
eggnog_go_universe = eggnog_go_universe.explode("go_ids").select(
    pl.col("query"), pl.col("go_ids").str.replace_all(r"-", "")
)
eggnog_go_universe = eggnog_go_universe.filter(
    ~(pl.col("go_ids") == ""),
).select(
    pl.col("go_ids"),
    pl.col("query"),
)

_write_csv(
    eggnog_go_universe,
    cfg.outputs.eggnog_go_universe,
    include_header=True,
    separator="\t",
)


# ---------------------------------------------------------------------------
# Merge InterPro and EggNOG-mapper GO annotations
# ---------------------------------------------------------------------------
# The two GO annotations were done on different genome assemblies; the
# protein_ids should match.
protein_ids_go = set(eggnog_go_universe["query"].to_list())
protein_ids_interpro = set(gff_df_with_go["gene_id"].to_list())

stats = {}
stats["total_proteins"] = len(protein_ids_go | protein_ids_interpro)
stats["common_proteins"] = len(protein_ids_go & protein_ids_interpro)
stats["eggnog_only_proteins"] = len(
    protein_ids_go - protein_ids_interpro
)
stats["interpro_only_proteins"] = len(
    protein_ids_interpro - protein_ids_go
)
print("Total proteins:", stats["total_proteins"])
print(
    "Percentage of common proteins:",
    stats["common_proteins"] / stats["total_proteins"] * 100,
)
print(
    "Percentage of eggnog-only proteins:",
    stats["eggnog_only_proteins"] / stats["total_proteins"] * 100,
)
print(
    "Percentage of interpro-only proteins:",
    stats["interpro_only_proteins"] / stats["total_proteins"] * 100,
)

# Merge the two DataFrames in the exploded format, and deduplicate
merged_df = pl.concat(
    [
        eggnog_go_universe,
        gff_df_with_go.select(
            pl.col("go_id").alias("go_ids"),
            pl.col("gene_id").alias("query"),
        ),
    ],
    how="vertical",
).unique()

merged_df = merged_df.filter(
    (~(pl.col("go_ids") == "")) & (~(pl.col("go_ids").is_null()))
)

merged_df = merged_df.select(
    pl.col("query").alias("GID"),
    pl.col("go_ids").alias("GO"),
)

# Obsolete-term remap for the merged set (same algorithm as above).
total_go_ids = set(merged_df["GO"].to_list())

remap, completely_obsolete = build_obsolete_remap(
    total_go_ids, cfg.inputs.go_basic_obo
)

interpro_eggnog_annotation_filtered_with = merged_df.select(
    pl.col("GID"), pl.col("GO").replace(remap)
).filter(~pl.col("GO").is_in(completely_obsolete))

annotations_out_dir = Path(cfg.outputs.results_dir)
annotations_out_dir.mkdir(parents=True, exist_ok=True)
interpro_eggnog_annotation_filtered_with.write_csv(
    cfg.outputs.go_universe,
    include_header=True,
    separator="\t",
)


# ---------------------------------------------------------------------------
# Chromosome metadata
# ---------------------------------------------------------------------------
annotations_with_chromosome_df = (
    interpro_eggnog_annotation_filtered_with.join(
        bed_cds_fixed,
        left_on="GID",
        right_on="gene_id",
        how="left",
    )
    .select(
        pl.col("GID"),
        pl.col("chrom"),
    )
    .unique()
)

annotations_with_chromosome_df = (
    annotations_with_chromosome_df.select(
        pl.col("GID"),
        pl.col("chrom").alias("CHROM"),
    )
)

annotations_with_chromosome_df.write_csv(
    cfg.outputs.chromosome_annotation,
    include_header=True,
    separator="\t",
)


# ---------------------------------------------------------------------------
# Blacklist / quality metadata
# ---------------------------------------------------------------------------
# Each exclusion BED is intersected with the CDS BED (>= F of the CDS feature
# covered).
bed_cds_fixed_bt = BedTool(cfg.outputs.cds_regions_fixed_bed)


def genes_overlapping_exclusion(exclusion_bed, gid_column, label):
    tmp = BedTool(exclusion_bed).intersect(
        bed_cds_fixed_bt, wb=True, F=cfg.params.blacklist_overlap_fraction
    )
    df = pl.read_csv(tmp.fn, separator="\t", has_header=False)
    return df.select(
        pl.col(gid_column).alias("GID"),
        pl.lit(label).alias("exlusion_reason"),
    ).unique()


low_coverage_genes_df = genes_overlapping_exclusion(
    cfg.inputs.exclusion_beds.low_coverage, "column_13", "low_coverage"
)
high_coverage_genes_df = genes_overlapping_exclusion(
    cfg.inputs.exclusion_beds.high_coverage, "column_13", "high_coverage"
)
low_quality_genes_df = genes_overlapping_exclusion(
    cfg.inputs.exclusion_beds.low_quality, "column_7", "low_quality"
)
near_gap_genes_df = genes_overlapping_exclusion(
    cfg.inputs.exclusion_beds.near_gap, "column_7", "near_gap"
)

# Merge all the exclusion regions
exclusion_regions_df = pl.concat(
    [
        low_quality_genes_df,
        low_coverage_genes_df,
        high_coverage_genes_df,
        near_gap_genes_df,
    ],
    how="vertical",
).unique()

exclusion_regions_df_table = exclusion_regions_df.with_columns(
    pl.lit(True).alias("value")
)
exclusion_regions_df_table = exclusion_regions_df_table.pivot(
    on="exlusion_reason",
    index="GID",
)
exclusion_regions_df_table = exclusion_regions_df_table.fill_null(
    False
)

# Add genes without any exclusion to the table
all_genes = set(annotations_with_chromosome_df["GID"].to_list())

genes_with_exclusion = set(
    exclusion_regions_df_table["GID"].to_list()
)
genes_without_exclusion = all_genes - genes_with_exclusion
genes_without_exclusion = pl.DataFrame(
    {
        "GID": list(genes_without_exclusion),
    }
)
genes_without_exclusion = genes_without_exclusion.with_columns(
    pl.lit(False).alias("low_quality"),
    pl.lit(False).alias("low_coverage"),
    pl.lit(False).alias("high_coverage"),
    pl.lit(False).alias("near_gap"),
)
exclusion_regions_total_table = pl.concat(
    [
        exclusion_regions_df_table,
        genes_without_exclusion,
    ],
    how="vertical",
).unique()
exclusion_regions_total_table = exclusion_regions_total_table.sort(
    "GID"
)

stats = {}
stats["total_genes"] = exclusion_regions_total_table.shape[0]
stats["genes_with_exclusion"] = exclusion_regions_total_table.filter(
    pl.any_horizontal(pl.exclude("GID"))
).shape[0]
stats["low_coverage_genes"] = exclusion_regions_total_table.filter(
    pl.col("low_coverage") == True
).shape[0]
stats["low_quality_genes"] = exclusion_regions_total_table.filter(
    pl.col("low_quality") == True
).shape[0]
stats["high_coverage_genes"] = exclusion_regions_total_table.filter(
    pl.col("high_coverage") == True
).shape[0]
stats["near_gap_genes"] = exclusion_regions_total_table.filter(
    pl.col("near_gap") == True
).shape[0]
print("Total genes:", stats["total_genes"])
print(
    "Percentage of genes with exclusion regions:",
    stats["genes_with_exclusion"] / stats["total_genes"] * 100,
)
print(
    "Percentage of low coverage genes:",
    stats["low_coverage_genes"] / stats["total_genes"] * 100,
)
print(
    "Percentage of low quality genes:",
    stats["low_quality_genes"] / stats["total_genes"] * 100,
)
print(
    "Percentage of high coverage genes:",
    stats["high_coverage_genes"] / stats["total_genes"] * 100,
)
print(
    "Percentage of near gap genes:",
    stats["near_gap_genes"] / stats["total_genes"] * 100,
)

exclusion_regions_total_table.write_csv(
    cfg.outputs.quality_annotations,
    include_header=True,
    separator="\t",
)


# ---------------------------------------------------------------------------
# Per-replicate sweep gene lists (blacklist-filtered)
# ---------------------------------------------------------------------------
interpro_go_bed = BedTool(cfg.outputs.cds_regions_fixed_bed)

bed_files_dir = Path(cfg.inputs.sweep_calls_dir)
bed_files = sorted(list(bed_files_dir.glob("*.bed")))


def get_sample_id(bed_file):
    # Extract the sample ID from the filename
    return bed_file.stem.split("_")[0]


genes_to_exclude = set(
    exclusion_regions_total_table.filter(
        pl.any_horizontal(pl.exclude("GID"))
    )["GID"].to_list()
)

per_replicate_dir = Path(cfg.outputs.per_replicate_genes_dir)
per_replicate_dir.mkdir(parents=True, exist_ok=True)
for bed in bed_files:

    hits = BedTool(bed).intersect(interpro_go_bed, wa=True, wb=True)

    # Isolate last column and write to file
    sample_id = get_sample_id(bed)
    output_file = per_replicate_dir / f"{sample_id}_filtered_for_all.genes.txt"
    tmp_df = pd.read_csv(hits.fn, sep="\t", header=None)
    last_col = tmp_df.iloc[:, -1]
    last_col = last_col.sort_values(ascending=True)
    # Filter out genes to exclude
    last_col = last_col[
        ~last_col.isin(genes_to_exclude)
    ]  # Exclude genes in the list
    last_col.to_csv(
        output_file,
        sep="\t",
        index=False,
        header=False,
    )


# ---------------------------------------------------------------------------
# Consistently-selected genes per regime (hit in > threshold replicate samples)
# ---------------------------------------------------------------------------
outpur_dir = Path(cfg.outputs.consistently_selected_dir)
outpur_dir.mkdir(parents=True, exist_ok=True)

late_population_bed_files = sorted(list(bed_files_dir.glob("L*.bed")))
all_gene_selection_events_Late = []
for bed in late_population_bed_files:

    hits = BedTool(bed).intersect(interpro_go_bed, wa=True, wb=True)

    tmp_df = pd.read_csv(hits.fn, sep="\t", header=None)
    # deduplicate the last column
    # This is neccessary since the GFF --> BED conversion
    # can create multiple lines for the same gene. It does not
    # affect the downstream analysis, since we assume 1 gene in one sample
    last_col = tmp_df.iloc[:, -1]
    last_col = last_col.drop_duplicates()

    # Add intersected genes to the list
    all_gene_selection_events_Late.append(last_col.to_list())

all_gene_selection_events_Late = pl.Series(
    "genes_late",
    all_gene_selection_events_Late,
).explode()
late_value_counts = all_gene_selection_events_Late.value_counts()

all_gene_selection_events_early = []
early_bed_files = sorted(list(bed_files_dir.glob("E*.bed")))
for bed in early_bed_files:

    hits = BedTool(bed).intersect(interpro_go_bed, wa=True, wb=True)

    tmp_df = pd.read_csv(hits.fn, sep="\t", header=None)
    # deduplicate the last column
    # This is neccessary since the GFF --> BED conversion
    # can create multiple lines for the same gene. It does not
    # affect the downstream analysis, since we assume 1 gene in one sample
    last_col = tmp_df.iloc[:, -1]
    last_col = last_col.drop_duplicates()

    # Add intersected genes to the list
    all_gene_selection_events_early.append(last_col.to_list())

all_gene_selection_events_early = pl.Series(
    "genes_early",
    all_gene_selection_events_early,
).explode()

early_value_counts = all_gene_selection_events_early.value_counts()
early_value_counts = early_value_counts.rename(
    {
        "genes_early": "gene_id",
        "count": "early_selection_events",
    }
)

# filter and save the vectors
early_value_counts_filtered = (
    early_value_counts.filter(
        pl.col("early_selection_events") > cfg.params.selection_event_threshold
    )
    .select(pl.col("gene_id"))
    .unique()
)
late_value_counts_filtered = (
    late_value_counts.filter(
        pl.col("count") > cfg.params.selection_event_threshold
    )
    .select(pl.col("genes_late"))
    .unique()
)

# save the vectors
early_value_counts_filtered.write_csv(
    cfg.outputs.consistently_selected_early,
    include_header=False,
    separator="\t",
)
late_value_counts_filtered.write_csv(
    cfg.outputs.consistently_selected_late,
    include_header=False,
    separator="\t",
)
