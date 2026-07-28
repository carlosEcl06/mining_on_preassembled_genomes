#!/usr/bin/env python3
"""
01_extract_domains.py

Reads the comBGC summary table and extracts Pfam domain arrays for each BGC.
Applies pre-clustering filters:
  - deepBGC BGCs with BGC_probability < 0.50 are excluded
  - BGCs with empty PFAM_domains are excluded
Outputs a TSV with one BGC per row, ready for pairwise Jaccard computation.

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

# Backbone domain map built from empirical top-domain analysis of this dataset.
# Classes are mutually inclusive: a BGC may belong to more than one class.
# NRPS_like reflects the absence of canonical C/A/T domains in H. pylori NRP BGCs;
# the dominant domains suggest an NIS (NRPS-independent siderophore) or NRPS-like pathway.
BACKBONE_MAP = {
    "Terpene": {
        "PT_FPPS_like",  # farnesyl-PP synthase-like — antiSMASH Terpene-precursor
        "PF02353",       # terpenoid cyclase
        "PF04909",       # UbiA prenyltransferase — deepBGC Terpene
        "PF03544",       # terpene synthase metal-binding — deepBGC Terpene
        "PF13088",       # Rieske-like, terpene oxidation — deepBGC Terpene
    },
    "NRPS_like": {
        "PF00107",       # alcohol dehydrogenase — dominant in NRP BGCs (n=41)
        "PF01262",       # subtilase — dominant in NRP BGCs (n=41)
        "PF02826",       # D-alanyl carrier — dominant in NRP BGCs (n=41)
        "PF03446",       # NAD-binding — dominant in NRP BGCs (n=41)
        "PF08240",       # alcohol DH-like — dominant in NRP BGCs (n=41)
        "PF01501",       # glycosyl hydrolase — NRP BGCs (n=35)
    },
    "PKS": {
        "PF00109",       # beta-ketoacyl synthase (KS) — canonical PKS backbone
        "PF02801",       # acyltransferase (AT) — canonical PKS backbone
        "PF00550",       # phosphopantetheine attachment site (T/ACP)
    },
    "RiPP": {
        "PF00398",       # thioredoxin-like — most specific for RiPP in this dataset
        "PF02624",       # YcaO-like, thiazole/oxazole biosynthesis — RiPP backbone
        "PF02463",       # RecA-like, associated with sactipeptide RiPPs
    },
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",  required=True, help="Path to combgc_complete_summary.tsv")
    parser.add_argument("--output", required=True, help="Path to output TSV")
    return parser.parse_args()


def parse_domain_array(raw: str) -> list[str]:
    """Split a PFAM_domains string into a deduplicated, sorted list of domain IDs."""
    if pd.isna(raw) or str(raw).strip() == "":
        return []
    tokens = re.split(r"[;,\s]+", str(raw).strip())
    return sorted(set(t.strip() for t in tokens if t.strip()))


def classify_backbone(domains: list[str]) -> str:
    """
    Assign backbone class(es) based on empirically validated Pfam accessions.
    Returns pipe-separated class names, or 'Unknown' if no backbone domain is found.
    """
    domain_set = set(domains)
    classes = [cls for cls, markers in BACKBONE_MAP.items() if domain_set & markers]
    return "|".join(classes) if classes else "Unknown"


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

    # --- domain extraction and backbone classification ---
    df["domain_array"]   = df["PFAM_domains"].apply(
        lambda x: "|".join(parse_domain_array(x))
    )
    df["backbone_class"] = df["domain_array"].apply(
        lambda x: classify_backbone(x.split("|") if x else [])
    )

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