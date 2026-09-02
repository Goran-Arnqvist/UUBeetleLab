#!/usr/bin/env python
"""
Author: Athanasios Toros
Email: athanasiostoros@gmail.com
"""

import polars as pl
from pathlib import Path
from types import SimpleNamespace


def _to_namespace(obj):
    if isinstance(obj, dict):
        return SimpleNamespace(**{k: _to_namespace(v) for k, v in obj.items()})
    return obj


def load_config(path=Path(__file__).parent / "config.yaml"):
    """Load config.yaml into a nested SimpleNamespace."""
    import yaml

    with open(path) as fh:
        raw = yaml.safe_load(fh)
    return _to_namespace(raw)


cfg = load_config()


def load_gene_set(path: str, sep: str = "\t", has_header: bool = False) -> set[str]:
    df = pl.read_csv(path, separator=sep, has_header=has_header, infer_schema_length=0)
    col = df.columns[0]
    vals = df[col].to_list()
    return {str(x).strip() for x in vals if x is not None and str(x).strip() != ""}


# Inputs
early_cons_path = cfg.outputs.consistently_selected_early
late_cons_path  = cfg.outputs.consistently_selected_late

early_del_path  = cfg.inputs.cnv_gene_lists.early_deletions
late_del_path   = cfg.inputs.cnv_gene_lists.late_deletions

early_gain_path = cfg.inputs.cnv_gene_lists.early_gain
late_gain_path  = cfg.inputs.cnv_gene_lists.late_gain

# Load
consistently_selected_early = load_gene_set(early_cons_path, sep="\t", has_header=False)
consistently_selected_late  = load_gene_set(late_cons_path,  sep="\t", has_header=False)

deletions_early = load_gene_set(early_del_path,  sep="\t", has_header=False)
deletions_late  = load_gene_set(late_del_path,   sep="\t", has_header=False)

gains_early = load_gene_set(early_gain_path, sep="\t", has_header=False)
gains_late  = load_gene_set(late_gain_path,  sep="\t", has_header=False)

def overlap_stats(A: set[str], B: set[str], nameA: str, nameB: str) -> set[str]:
    inter = A & B
    nA, nB, nI = len(A), len(B), len(inter)
    pct_A_in_B = (nI / nA * 100) if nA else 0.0
    pct_B_in_A = (nI / nB * 100) if nB else 0.0
    jaccard = (nI / len(A | B) * 100) if (A or B) and len(A | B) else 0.0
    print(f"{nameA}: {nA}")
    print(f"{nameB}: {nB}")
    print(f"Common {nameA} ∩ {nameB}: {nI}")
    print(f"% of {nameA} in {nameB}: {pct_A_in_B:.1f}%")
    print(f"% of {nameB} in {nameA}: {pct_B_in_A:.1f}%")
    print(f"Jaccard({nameA},{nameB}): {jaccard:.1f}%")
    print("-" * 40)
    return inter

print("===== Counts =====")
print(f"Early Consistently Selected Genes: {len(consistently_selected_early)}")
print(f"Late Consistently Selected Genes:  {len(consistently_selected_late)}")
print(f"Early Deletions:                  {len(deletions_early)}")
print(f"Late Deletions:                   {len(deletions_late)}")
print(f"Early Gains:                      {len(gains_early)}")
print(f"Late Gains:                       {len(gains_late)}")
print("=" * 40)

print("===== Overlaps =====")
inter_early_del  = overlap_stats(consistently_selected_early, deletions_early, "ConsistentlySelected(Early)", "Deletions(Early)")
inter_late_del   = overlap_stats(consistently_selected_late,  deletions_late,  "ConsistentlySelected(Late)",  "Deletions(Late)")
inter_early_gain = overlap_stats(consistently_selected_early, gains_early,     "ConsistentlySelected(Early)", "Gains(Early)")
inter_late_gain  = overlap_stats(consistently_selected_late,  gains_late,      "ConsistentlySelected(Late)",  "Gains(Late)")
