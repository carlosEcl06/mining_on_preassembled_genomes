#!/usr/bin/env python3
"""
01_extract_domains.py

Reads the comBGC summary table and extracts Pfam domain arrays for each BGC.
Applies pre-clustering filters:
  - deepBGC BGCs with BGC_probability < 0.50 are excluded
  - BGCs with empty PFAM_domains are excluded
Outputs a TSV with one BGC per row, ready for pairwise Jaccard computation.

backbone_class is taken directly from Product_class, the classification each
prediction tool already assigns (antiSMASH's rule-based typing, deepBGC's ML
classifier, GECCO's classifier), rather than being re-derived from a
hand-picked Pfam accession list. An earlier heuristic based on such a list
was dropped after cross-tabulation showed it mislabeling most deepBGC
"Saccharide"/"Other" BGCs as "Terpene" (its marker accessions weren't
specific to terpene biosynthesis).

Product_class is missing (NaN) for most deepBGC BGCs — deepBGC only
populates it when its classifier is confident enough. This collapses to
backbone_class = "Unknown" here, same as GECCO's explicit "Unknown" label;
it only affects the descriptive backbone_class column, since backbone
identity computation (00_extract_bgc_fastas.py) only covers antiSMASH BGCs,
whose Product_class is always populated.

bgc_id convention: sample_id + contig_id alone is not unique, since a contig
can carry more than one antiSMASH region or be called by more than one
prediction tool. BGC_start disambiguates these cases, so bgc_id is built here
as sample_id__contig_id__BGC_start and written out as its own column — every
downstream script (00, 02, 03, 04) reads this column directly.

Usage:
    python 01_extract_domains.py \
        --input  ../../R/copied_from_funcscan_results/combgc_complete_summary.tsv \
        --output ../results/bgc_domain_arrays.tsv
"""

import argparse
import pathlib
import re

import pandas as pd


DEEPBGC_PROB_THRESHOLD = 0.50


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",  required=True, help="Path to combgc_complete_summary.tsv")
    parser.add_argument("--output", required=True, help="Path to output TSV")
    return parser.parse_args()


def parse_domain_array(raw: str) -> list[str]:
    """Split a PFAM_domains string into a deduplicated, sorted list of domain IDs.

    Note: this field mixes Pfam accessions (PFxxxxx) and Pfam family names
    (e.g. 'YcaO') inconsistently by tool — antiSMASH reports the name,
    deepBGC/GECCO report the accession, for the same Pfam family. Doesn't
    affect backbone_class (sourced from Product_class), but is a known
    limitation for any future cross-tool Jaccard comparison in 02.
    """
    if pd.isna(raw) or str(raw).strip() == "":
        return []
    tokens = re.split(r"[;,\s]+", str(raw).strip())
    return sorted(set(t.strip() for t in tokens if t.strip()))


def build_bgc_id(df: pd.DataFrame) -> pd.Series:
    """Canonical, unique BGC identifier: sample_id__contig_id__BGC_start.
    Falls back to sample_id__contig_id__NA if BGC_start is missing."""
    start = df["BGC_start"].fillna("NA").astype(str)
    return df["sample_id"] + "__" + df["contig_id"] + "__" + start


def main():
    args = parse_args()

    print(f"[01] Reading {args.input} ...", flush=True)
    df = pd.read_csv(args.input, sep="\t", dtype=str)
    print(f"[01] {len(df)} BGCs loaded.", flush=True)

    # --- filters ---
    # 1. deepBGC probability threshold
    df["BGC_probability"] = pd.to_numeric(df["BGC_probability"], errors="coerce")
    low_prob = (df["Prediction_tool"] == "deepBGC") & (df["BGC_probability"] < DEEPBGC_PROB_THRESHOLD)
    n_dropped_prob = low_prob.sum()
    df = df[~low_prob].copy()
    print(f"[01] Dropped {n_dropped_prob} deepBGC BGCs below probability "
          f"threshold ({DEEPBGC_PROB_THRESHOLD}).", flush=True)

    # 2. Empty domain arrays
    empty_domains = df["PFAM_domains"].apply(lambda x: parse_domain_array(x) == [])
    n_dropped_empty = empty_domains.sum()
    df = df[~empty_domains].copy()
    print(f"[01] Dropped {n_dropped_empty} BGCs with empty PFAM_domains.", flush=True)

    print(f"[01] {len(df)} BGCs retained after filtering.", flush=True)

    # --- domain array (still needed for Jaccard similarity in 02) ---
    df["domain_array"] = df["PFAM_domains"].apply(
        lambda x: "|".join(parse_domain_array(x))
    )

    # --- backbone_class: each tool's own Product_class (see module docstring) ---
    n_missing_product_class = df["Product_class"].isna().sum()
    df["backbone_class"] = df["Product_class"].fillna("Unknown")
    print(f"[01] backbone_class sourced from Product_class "
          f"({n_missing_product_class} BGCs had no Product_class from their "
          f"tool and were set to 'Unknown').", flush=True)
    print("[01] backbone_class distribution:", flush=True)
    for cls, count in df["backbone_class"].value_counts().items():
        print(f"      {count:5d}  {cls}", flush=True)

    # --- canonical bgc_id ---
    df["bgc_id"] = build_bgc_id(df)
    n_dup = df["bgc_id"].duplicated().sum()
    if n_dup:
        print(f"[01] WARNING: {n_dup} duplicate bgc_id values remain even with "
              f"BGC_start included. These are likely exact-duplicate rows in "
              f"the comBGC summary (identical sample_id/contig_id/BGC_start) — "
              f"inspect manually before proceeding.", flush=True)
    else:
        print("[01] bgc_id is unique across all retained BGCs.", flush=True)

    # --- output ---
    cols_out = [
        "bgc_id", "sample_id", "contig_id", "Prediction_tool", "Product_class",
        "BGC_probability", "BGC_complete", "BGC_start", "BGC_end",
        "BGC_length", "CDS_count", "domain_array", "backbone_class",
    ]
    cols_out = [c for c in cols_out if c in df.columns.tolist() + ["bgc_id", "domain_array", "backbone_class"]]

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df[cols_out].to_csv(args.output, sep="\t", index=False)
    print(f"[01] Output written to {args.output}", flush=True)


if __name__ == "__main__":
    main()