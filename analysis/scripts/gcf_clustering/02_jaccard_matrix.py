#!/usr/bin/env python3
"""
02_jaccard_matrix.py

Computes pairwise Jaccard similarity between all BGC domain arrays.
BGC pairs with Jaccard similarity below MIN_JACCARD are excluded (set to 0).
Outputs a sparse pairwise TSV (bgc_i, bgc_j, jaccard) for use in 03 and 04.

bgc_id is read directly from the bgc_id column written by 01_extract_domains.py,
not re-derived from sample_id/contig_id (see 01's module docstring).

Usage:
    python 02_jaccard_matrix.py \
        --input  ../results/bgc_domain_arrays.tsv \
        --output ../results/jaccard_pairs.tsv
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import pairwise_distances
from tqdm import tqdm


MIN_JACCARD = 0.10


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",  required=True, help="Path to bgc_domain_arrays.tsv")
    parser.add_argument("--output", required=True, help="Path to output sparse pairs TSV")
    parser.add_argument("--min-jaccard", type=float, default=MIN_JACCARD,
                        help=f"Minimum Jaccard similarity threshold (default: {MIN_JACCARD})")
    return parser.parse_args()


def domain_array_to_binary(domain_arrays: pd.Series) -> np.ndarray:
    """
    Convert pipe-separated domain arrays into a binary occurrence matrix.
    Each column is a Pfam domain; each row is a BGC.
    Uses CountVectorizer with binary=True so duplicate domains are collapsed.
    """
    # Replace pipes with spaces so each domain becomes a token
    tokenized = domain_arrays.str.replace("|", " ", regex=False)
    vec = CountVectorizer(binary=True, token_pattern=r"[^\s]+")
    return vec.fit_transform(tokenized).toarray().astype(np.float32)


def jaccard_from_binary(mat: np.ndarray) -> np.ndarray:
    """
    Compute pairwise Jaccard similarity from a binary matrix.
    Jaccard(A, B) = |A ∩ B| / |A ∪ B|
    Uses sklearn pairwise_distances with metric='jaccard' on boolean arrays,
    which returns distance (1 - similarity), so we subtract from 1.
    """
    dist = pairwise_distances(mat.astype(bool), metric="jaccard", n_jobs=-1)
    return (1.0 - dist).astype(np.float32)


def main():
    args = parse_args()

    print(f"[02] Reading {args.input} ...", flush=True)
    df = pd.read_csv(args.input, sep="\t", dtype=str)
    n = len(df)
    print(f"[02] {n} BGCs loaded.", flush=True)

    if "bgc_id" not in df.columns:
        raise SystemExit(
            "[02] ERROR: no 'bgc_id' column in input. Re-run 01_extract_domains.py "
            "with the current version, which writes this column directly."
        )

    if df["bgc_id"].duplicated().any():
        n_dup = df["bgc_id"].duplicated().sum()
        print(f"[02] WARNING: {n_dup} duplicate bgc_id values found even after "
              "including BGC_start. Check for exact-duplicate rows upstream.",
              flush=True)

    # --- binary domain matrix ---
    print("[02] Building binary domain matrix ...", flush=True)
    mat = domain_array_to_binary(df["domain_array"])
    print(f"[02] Matrix shape: {mat.shape} ({mat.shape[0]} BGCs × {mat.shape[1]} domains).",
          flush=True)

    # --- pairwise Jaccard ---
    print("[02] Computing pairwise Jaccard similarities ...", flush=True)
    sim = jaccard_from_binary(mat)

    # --- sparsify: keep only upper triangle above threshold ---
    print(f"[02] Filtering pairs with Jaccard < {args.min_jaccard} ...", flush=True)
    rows_idx, cols_idx = np.triu_indices(n, k=1)
    values = sim[rows_idx, cols_idx]

    mask = values >= args.min_jaccard
    rows_idx = rows_idx[mask]
    cols_idx = cols_idx[mask]
    values   = values[mask]

    n_pairs_total    = len(np.triu_indices(n, k=1)[0])
    n_pairs_retained = mask.sum()
    print(f"[02] {n_pairs_retained:,} / {n_pairs_total:,} pairs retained "
        f"({100 * n_pairs_retained / n_pairs_total:.1f}%).", flush=True)

    # --- output with progress bar ---
    bgc_ids = df["bgc_id"].to_numpy()
    chunk_size = 100_000
    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    with open(args.output, "w") as fh:
        fh.write("bgc_i\tbgc_j\tjaccard\n")
        for start in tqdm(range(0, len(rows_idx), chunk_size),
                        desc="[02] Writing pairs", unit="chunk"):
            end = start + chunk_size
            chunk = pd.DataFrame({
                "bgc_i":   bgc_ids[rows_idx[start:end]],
                "bgc_j":   bgc_ids[cols_idx[start:end]],
                "jaccard": values[start:end].round(6),
            })
            chunk.to_csv(fh, sep="\t", index=False, header=False)

    print(f"[02] Output written to {args.output}", flush=True)

if __name__ == "__main__":
    main()