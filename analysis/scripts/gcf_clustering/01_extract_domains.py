#!/usr/bin/env python3
"""
01_extract_domains.py

Reads the comBGC summary table and extracts Pfam domain arrays for each BGC.
Applies pre-clustering filters:
  - deepBGC BGCs with BGC_probability < 0.50 are excluded
  - BGCs with empty PFAM_domains are excluded
Outputs a TSV with one BGC per row, ready for pairwise Jaccard computation.

backbone_class convention (revised 2026-07-29): backbone_class is now taken
DIRECTLY from Product_class, the classification already assigned by each
prediction tool (antiSMASH's curated rule-based cluster typing, deepBGC's ML
classifier, GECCO's classifier) — not re-derived from a hand-picked Pfam
accession list (the old BACKBONE_MAP/classify_backbone() approach).

This replaces an approach that was empirically shown to be unreliable: a
cross-tabulation of the old classify_backbone() output against Product_class
found that of BGCs labeled "Saccharide" by deepBGC, 94% (117/125) were
mislabeled "Terpene" by our own domain-matching heuristic — and 100% (83/83)
of "Other" BGCs were mislabeled "Terpene" too — because the deepBGC "Terpene"
marker accessions we'd picked (PF04909, PF03544, PF13088) are not specific to
terpene biosynthesis and appear incidentally across unrelated BGC classes.
Trusting each tool's own Product_class avoids this false-positive problem
entirely, since it comes from purpose-built classification logic rather than
a post-hoc accession list assembled without this kind of verification.

Product_class is missing (NaN) for 1,005 of 1,465 retained deepBGC BGCs
(68.6%) — deepBGC only populates it when its classifier is confident enough,
leaving it blank otherwise (distinct from GECCO's "Unknown", which IS an
explicit label). Both cases collapse to backbone_class = "Unknown" here,
since both mean "the tool did not confidently assign a class" — this has no
effect on backbone identity computation (00_extract_bgc_fastas.py only
extracts antiSMASH BGCs, whose Product_class is always populated), only on
the descriptive backbone_class column in the final output.

bgc_id convention: sample_id + contig_id ALONE is not unique — a single contig
can carry more than one antiSMASH region (region001, region002, ...) or be
called independently by more than one prediction tool, producing multiple
rows with the same (sample_id, contig_id) pair. BGC_start (the genomic
coordinate where the region/cluster begins) disambiguates these cases, so
bgc_id is built here as sample_id__contig_id__BGC_start and written out as
its own column — every downstream script (00, 02, 03, 04) should read this
column directly rather than re-deriving bgc_id from sample_id/contig_id alone.

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

    Note: this field mixes Pfam accessions (PFxxxxx) and Pfam family NAMEs
    (e.g. 'YcaO', 'PP-binding') inconsistently BY TOOL — antiSMASH always
    reports the NAME, deepBGC/GECCO always report the ACC, for the same
    underlying Pfam family (verified 2026-07-29 by checking co-occurrence of
    several NAME/ACC pairs — e.g. antiSMASH's 'PP-binding' vs deepBGC's
    'PF00550' never co-occur within the same tool). This does NOT affect
    backbone_class (now sourced from Product_class, see module docstring),
    but is a known limitation of the Jaccard similarity in 02_jaccard_matrix.py
    for any future CROSS-TOOL domain-content comparison — tracked as a
    low-priority follow-up, not fixed here.
    """
    if pd.isna(raw) or str(raw).strip() == "":
        return []
    tokens = re.split(r"[;,\s]+", str(raw).strip())
    return sorted(set(t.strip() for t in tokens if t.strip()))


def build_bgc_id(df: pd.DataFrame) -> pd.Series:
    """
    Canonical, unique BGC identifier: sample_id__contig_id__BGC_start.
    BGC_start disambiguates multiple regions/calls on the same contig
    (see module docstring). Falls back to sample_id__contig_id__NA if
    BGC_start is missing, with a warning left to the caller to surface.
    """
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

    # --- backbone_class: trust each tool's own Product_class directly ---
    # (see module docstring for why this replaced the old Pfam-matching heuristic)
    n_missing_product_class = df["Product_class"].isna().sum()
    df["backbone_class"] = df["Product_class"].fillna("Unknown")
    print(f"[01] backbone_class sourced from Product_class "
          f"({n_missing_product_class} BGCs had no Product_class from their "
          f"tool and were set to 'Unknown').", flush=True)
    print("[01] backbone_class distribution:", flush=True)
    for cls, count in df["backbone_class"].value_counts().items():
        print(f"      {count:5d}  {cls}", flush=True)

    # --- canonical bgc_id (see module docstring: sample_id__contig_id is NOT unique) ---
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