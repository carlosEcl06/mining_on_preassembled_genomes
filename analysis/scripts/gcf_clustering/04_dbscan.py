#!/usr/bin/env python3
"""
04_dbscan.py

Computes the final similarity score from jaccard_pairs.tsv, builds a pairwise
distance matrix, and runs DBSCAN to delineate Gene Cluster Families (GCFs).

When --identity is provided (backbone_identity.tsv from 03_backbone_identity.py),
pairs with a computed sequence_identity use the full Robey et al. (2021) score:

    Similarity Score = sqrt(0.8 * sequence_identity + 0.2 * Jaccard)

Pairs without a computed sequence_identity (or when --identity is omitted
entirely) fall back to:

    Similarity Score = sqrt(0.2 * Jaccard)

which corresponds to the Robey et al. (2021) formula with sequence_identity = 0.
This is the appropriate fallback for H. pylori, whose biosynthetic repertoire
is dominated by BGC classes lacking canonical backbone enzyme architecture.

BGCs not present in any pair above the Jaccard threshold (singletons by isolation)
are assigned GCF_id = -1 before DBSCAN runs. DBSCAN may also produce singletons
(label = -1) for BGCs that are connected but fail the min_samples density criterion.

Output: gcf_assignments.tsv — one row per BGC with GCF assignment and metadata.

Usage:
    python 04_dbscan.py \
        --domains ../results/bgc_domain_arrays.tsv \
        --pairs   ../results/jaccard_pairs.tsv \
        --identity ../results/backbone_identity.tsv \
        --output  ../results/gcf_assignments.tsv
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix, csr_matrix
from sklearn.cluster import DBSCAN
from tqdm import tqdm

# epsilon = 0.56 corresponds to grouping BGCs with identical Pfam domain content
# under the Jaccard-only fallback (Jaccard = 1.0, distance = 1 - sqrt(0.2) ≈ 0.5528).
# Under the full score (identity=1.0, Jaccard=1.0), distance = 1 - sqrt(1.0) = 0.0,
# so 0.56 remains a valid (more permissive, never more restrictive) upper bound
# when --identity is used.
DBSCAN_EPSILON = 0.56
DBSCAN_MIN_SAMPLES = 2


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--domains",      required=True, help="Path to bgc_domain_arrays.tsv")
    parser.add_argument("--pairs",        required=True, help="Path to jaccard_pairs.tsv")
    parser.add_argument("--output",       required=True, help="Path to output TSV")
    parser.add_argument("--epsilon",      type=float, default=DBSCAN_EPSILON,
                        help=f"DBSCAN epsilon in distance space (default: {DBSCAN_EPSILON})")
    parser.add_argument("--min-samples",  type=int,   default=DBSCAN_MIN_SAMPLES,
                        help=f"DBSCAN min_samples (default: {DBSCAN_MIN_SAMPLES})")
    parser.add_argument("--identity", default=None,
                        help="Path to backbone_identity.tsv from 03_backbone_identity.py "
                             "(optional). When provided, pairs with a computed "
                             "sequence_identity use the full Robey et al. score "
                             "sqrt(0.8*identity + 0.2*Jaccard); pairs without one "
                             "fall back to sqrt(0.2*Jaccard), matching "
                             "experimental_design.Rmd.")
    return parser.parse_args()


def similarity_score_jaccard_only(jaccard: np.ndarray) -> np.ndarray:
    """Fallback for BGCs without identifiable backbone domains: sqrt(0.2 * Jaccard)."""
    return np.sqrt(0.2 * jaccard)


def similarity_score_full(jaccard: np.ndarray, identity: np.ndarray) -> np.ndarray:
    """Robey et al. (2021) full score: sqrt(0.8 * identity + 0.2 * Jaccard)."""
    return np.sqrt(0.8 * identity + 0.2 * jaccard)


def main():
    args = parse_args()

    # --- load BGC metadata ---
    print(f"[04] Reading domain arrays from {args.domains} ...", flush=True)
    df = pd.read_csv(args.domains, sep="\t", dtype=str)
    df["bgc_id"] = df["sample_id"] + "__" + df["contig_id"]
    bgc_ids  = df["bgc_id"].tolist()
    bgc_index = {bid: i for i, bid in enumerate(bgc_ids)}
    n = len(bgc_ids)
    print(f"[04] {n} BGCs loaded.", flush=True)

    # --- load pairs ---
    print(f"[04] Reading pairs from {args.pairs} ...", flush=True)
    pairs = pd.read_csv(args.pairs, sep="\t",
                        dtype={"bgc_i": str, "bgc_j": str, "jaccard": float})
    print(f"[04] {len(pairs):,} pairs loaded.", flush=True)

    # --- merge backbone sequence identity, if provided ---
    if args.identity:
        print(f"[04] Reading backbone identities from {args.identity} ...", flush=True)
        identity_df = pd.read_csv(
            args.identity, sep="\t",
            dtype={"bgc_i": str, "bgc_j": str, "sequence_identity": float},
        )
        pairs = pairs.merge(
            identity_df[["bgc_i", "bgc_j", "sequence_identity"]],
            on=["bgc_i", "bgc_j"], how="left",
        )
        n_with_identity = pairs["sequence_identity"].notna().sum()
        print(f"[04] {n_with_identity:,} / {len(pairs):,} pairs have a computed "
              f"sequence_identity; the rest fall back to Jaccard-only.", flush=True)
    else:
        pairs["sequence_identity"] = np.nan
        print("[04] No --identity provided; using Jaccard-only for all pairs "
              "(sqrt(0.2 * Jaccard)).", flush=True)

    # --- build distance matrix ---
    print("[04] Building sparse distance matrix ...", flush=True)
    # Use lil_matrix for efficient incremental construction, then convert to csr
    # Pairs not in the table have distance = 1.0 (Jaccard = 0, below threshold)
    dist_sparse = lil_matrix((n, n), dtype=np.float32)

    # Fill diagonal with 0 (self-distance)
    for i in range(n):
        dist_sparse[i, i] = 0.0

    skipped = 0
    for row in tqdm(pairs.itertuples(index=False),
                    total=len(pairs),
                    desc="[04] Filling distance matrix",
                    unit="pair"):
        i = bgc_index.get(row.bgc_i)
        j = bgc_index.get(row.bgc_j)
        if i is None or j is None:
            skipped += 1
            continue
        if pd.notna(row.sequence_identity):
            sim = similarity_score_full(
                np.array([row.jaccard]), np.array([row.sequence_identity])
            )[0]
        else:
            sim = similarity_score_jaccard_only(np.array([row.jaccard]))[0]
        d = float(1.0 - sim)
        dist_sparse[i, j] = d
        dist_sparse[j, i] = d

    if skipped:
        print(f"[04] WARNING: {skipped} pairs skipped (bgc_id not found in domain array).",
              flush=True)

    # Pairs absent from the sparse table are implicitly distance = 1.0 in DBSCAN
    # with metric='precomputed'. We set them explicitly to 1.0 for correctness.
    dist_csr = dist_sparse.tocsr()

    print("[04] Running DBSCAN ...", flush=True)
    db = DBSCAN(
        eps=args.epsilon,
        min_samples=args.min_samples,
        metric="precomputed",
        n_jobs=-1,
    )
    # DBSCAN with precomputed requires a dense or explicit sparse matrix.
    # With n=~2000, dense is fine (~16 MB float32).
    dist_dense = dist_csr.toarray()
    # Fill unset off-diagonal entries (still 0 from lil_matrix init) with 1.0
    off_diag_zero = (dist_dense == 0) & (np.eye(n, dtype=bool) == False)
    dist_dense[off_diag_zero] = 1.0

    labels = db.fit_predict(dist_dense)

    n_gcfs      = len(set(labels)) - (1 if -1 in labels else 0)
    n_singletons = (labels == -1).sum()
    print(f"[04] {n_gcfs} GCFs found, {n_singletons} singletons (label = -1).", flush=True)

    # --- format GCF labels ---
    # Singletons → "singleton"; clusters → "GCF_0001", "GCF_0002", ...
    gcf_ids = []
    for label in labels:
        if label == -1:
            gcf_ids.append("singleton")
        else:
            gcf_ids.append(f"GCF_{label + 1:04d}")

    # --- output ---
    df["gcf_id"] = gcf_ids
    df["dbscan_label"] = labels

    cols_out = [
        "bgc_id", "sample_id", "contig_id", "Prediction_tool", "Product_class",
        "BGC_probability", "BGC_complete", "BGC_start", "BGC_end", "BGC_length",
        "CDS_count", "backbone_class", "gcf_id", "dbscan_label",
    ]
    cols_out = [c for c in cols_out if c in df.columns]

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df[cols_out].to_csv(args.output, sep="\t", index=False)
    print(f"[04] Output written to {args.output}", flush=True)

    # --- summary to stdout ---
    print("\n[04] GCF size distribution:", flush=True)
    gcf_sizes = (
        df[df["gcf_id"] != "singleton"]
        .groupby("gcf_id")
        .size()
        .value_counts()
        .sort_index()
    )
    for size, count in gcf_sizes.items():
        print(f"      {count:4d} GCFs with {size:3d} BGC(s)", flush=True)
    print(f"      {n_singletons:4d} singletons", flush=True)


if __name__ == "__main__":
    main()