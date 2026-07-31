#!/usr/bin/env python3
"""
05_epsilon_sweep.py

Diagnostic sweep of the DBSCAN epsilon parameter, to characterize whether GCF
clustering granularity is dominated by a single global epsilon poorly suited
to heterogeneous cluster density — e.g. a near-clonal Terpene-precursor
"family" (554 BGCs differing by a gene or fragment) sitting alongside
genuinely diverse, sparser families (Azole-containing-RiPP, deepBGC/GECCO
Product_class values such as Polyketide, NRP, RiPP, Saccharide).

This does NOT re-run hmmscan/hmmfetch/hmmalign (03_backbone_identity.py) — it
reuses the already-computed jaccard_pairs.tsv and backbone_identity.tsv, builds
the pairwise distance matrix ONCE (same logic as 04_dbscan.py), and re-runs only
DBSCAN.fit_predict() for each candidate epsilon. This is why it's fast (~1 min
for a few dozen epsilon values) compared to the ~90 min step 03 takes.

For each epsilon, records:
  - n_gcfs                : number of non-singleton clusters
  - n_singletons          : BGCs with label = -1
  - largest_gcf_size      : size of the single biggest cluster
  - second_largest_size   : size of the next biggest (0 if only one cluster exists)
  - median_gcf_size_excl_largest : median size of all clusters EXCLUDING the largest
                            (tracks whether the "everything else" clusters are
                            collapsing into singletons as epsilon shrinks, which
                            would indicate artificial fragmentation rather than
                            genuine resolution of the giant cluster)

Usage:
    python 05_epsilon_sweep.py \
        --domains  ../results/bgc_domain_arrays.tsv \
        --pairs    ../results/jaccard_pairs.tsv \
        --identity ../results/backbone_identity.tsv \
        --output   ../results/epsilon_sweep.tsv \
        --plot     ../results/epsilon_sweep.png
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix
from sklearn.cluster import DBSCAN


DEFAULT_MIN_SAMPLES = 2


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--domains",  required=True, help="Path to bgc_domain_arrays.tsv")
    parser.add_argument("--pairs",    required=True, help="Path to jaccard_pairs.tsv")
    parser.add_argument("--identity", default=None,
                        help="Path to backbone_identity.tsv (optional, same as 04_dbscan.py)")
    parser.add_argument("--output",   required=True, help="Path to output sweep summary TSV")
    parser.add_argument("--plot",     default=None,
                        help="Optional path to save a diagnostic PNG plot")
    parser.add_argument("--min-samples", type=int, default=DEFAULT_MIN_SAMPLES,
                        help=f"DBSCAN min_samples, held fixed across the sweep "
                             f"(default: {DEFAULT_MIN_SAMPLES})")
    parser.add_argument("--eps-min",  type=float, default=0.05, help="Sweep start (default: 0.05)")
    parser.add_argument("--eps-max",  type=float, default=0.90, help="Sweep end (default: 0.90)")
    parser.add_argument("--eps-step", type=float, default=0.025, help="Sweep step (default: 0.025)")
    return parser.parse_args()


def similarity_score_jaccard_only(jaccard: np.ndarray) -> np.ndarray:
    return np.sqrt(0.2 * jaccard)


def similarity_score_full(jaccard: np.ndarray, identity: np.ndarray) -> np.ndarray:
    return np.sqrt(0.8 * identity + 0.2 * jaccard)


def build_distance_matrix(domains_path: str, pairs_path: str, identity_path: str | None):
    """Same construction as 04_dbscan.py — kept in sync deliberately, duplicated
    here rather than imported so this script stays a standalone diagnostic tool
    consistent with the rest of gcf_clustering/'s style (each script is
    self-contained)."""
    print(f"[05] Reading domain arrays from {domains_path} ...", flush=True)
    df = pd.read_csv(domains_path, sep="\t", dtype=str)
    if "bgc_id" not in df.columns:
        raise SystemExit("[05] ERROR: no 'bgc_id' column in --domains.")
    bgc_ids = df["bgc_id"].tolist()
    bgc_index = {bid: i for i, bid in enumerate(bgc_ids)}
    n = len(bgc_ids)
    print(f"[05] {n} BGCs loaded.", flush=True)

    print(f"[05] Reading pairs from {pairs_path} ...", flush=True)
    pairs = pd.read_csv(pairs_path, sep="\t",
                        dtype={"bgc_i": str, "bgc_j": str, "jaccard": float})

    if identity_path:
        print(f"[05] Reading backbone identities from {identity_path} ...", flush=True)
        identity_df = pd.read_csv(
            identity_path, sep="\t",
            dtype={"bgc_i": str, "bgc_j": str, "sequence_identity": float},
        )
        pairs = pairs.merge(
            identity_df[["bgc_i", "bgc_j", "sequence_identity"]],
            on=["bgc_i", "bgc_j"], how="left",
        )
    else:
        pairs["sequence_identity"] = np.nan

    print("[05] Building distance matrix (once, reused for every epsilon) ...", flush=True)
    dist_sparse = lil_matrix((n, n), dtype=np.float32)
    for i in range(n):
        dist_sparse[i, i] = 0.0

    for row in pairs.itertuples(index=False):
        i = bgc_index.get(row.bgc_i)
        j = bgc_index.get(row.bgc_j)
        if i is None or j is None:
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

    dist_dense = dist_sparse.tocsr().toarray()
    off_diag_zero = (dist_dense == 0) & (np.eye(n, dtype=bool) == False)
    dist_dense[off_diag_zero] = 1.0

    return dist_dense, n


def sweep(dist_dense: np.ndarray, n: int, eps_values: np.ndarray, min_samples: int) -> pd.DataFrame:
    rows = []
    for eps in eps_values:
        db = DBSCAN(eps=float(eps), min_samples=min_samples, metric="precomputed", n_jobs=-1)
        labels = db.fit_predict(dist_dense)

        n_singletons = int((labels == -1).sum())
        cluster_labels = labels[labels != -1]
        sizes = pd.Series(cluster_labels).value_counts().sort_values(ascending=False).tolist()

        n_gcfs = len(sizes)
        largest = sizes[0] if sizes else 0
        second_largest = sizes[1] if len(sizes) > 1 else 0
        sizes_excl_largest = sizes[1:] if len(sizes) > 1 else []
        median_excl_largest = float(np.median(sizes_excl_largest)) if sizes_excl_largest else 0.0

        rows.append({
            "epsilon": round(float(eps), 4),
            "n_gcfs": n_gcfs,
            "n_singletons": n_singletons,
            "largest_gcf_size": largest,
            "second_largest_size": second_largest,
            "median_gcf_size_excl_largest": median_excl_largest,
        })
        print(f"[05] eps={eps:.4f}  n_gcfs={n_gcfs:4d}  singletons={n_singletons:4d}  "
              f"largest={largest:4d}  2nd_largest={second_largest:4d}  "
              f"median_excl_largest={median_excl_largest:.1f}", flush=True)

    return pd.DataFrame(rows)


def make_plot(df: pd.DataFrame, plot_path: str, current_epsilon: float = 0.56):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axes[0].plot(df["epsilon"], df["n_gcfs"], marker="o", ms=3)
    axes[0].axvline(current_epsilon, color="red", linestyle="--", linewidth=1,
                     label=f"current epsilon = {current_epsilon}")
    axes[0].set_ylabel("Number of GCFs\n(excl. singletons)")
    axes[0].legend(fontsize=8)

    axes[1].plot(df["epsilon"], df["largest_gcf_size"], marker="o", ms=3,
                 label="largest GCF", color="tab:red")
    axes[1].plot(df["epsilon"], df["second_largest_size"], marker="o", ms=3,
                 label="2nd largest GCF", color="tab:orange")
    axes[1].axvline(current_epsilon, color="red", linestyle="--", linewidth=1)
    axes[1].set_ylabel("GCF size (# BGCs)")
    axes[1].legend(fontsize=8)

    axes[2].plot(df["epsilon"], df["median_gcf_size_excl_largest"], marker="o", ms=3,
                 color="tab:green", label="median size, excl. largest GCF")
    axes[2].plot(df["epsilon"], df["n_singletons"], marker="o", ms=3,
                 color="tab:gray", label="singletons")
    axes[2].axvline(current_epsilon, color="red", linestyle="--", linewidth=1)
    axes[2].set_ylabel("Count / median size")
    axes[2].set_xlabel("DBSCAN epsilon")
    axes[2].legend(fontsize=8)

    fig.suptitle("Epsilon sensitivity sweep — GCF clustering granularity")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    print(f"[05] Plot saved to {plot_path}", flush=True)


def main():
    args = parse_args()

    dist_dense, n = build_distance_matrix(args.domains, args.pairs, args.identity)

    eps_values = np.arange(args.eps_min, args.eps_max + args.eps_step / 2, args.eps_step)
    print(f"[05] Sweeping {len(eps_values)} epsilon values from "
          f"{args.eps_min} to {args.eps_max} (step {args.eps_step}) ...", flush=True)

    df = sweep(dist_dense, n, eps_values, args.min_samples)

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[05] Sweep summary written to {args.output}", flush=True)

    if args.plot:
        make_plot(df, args.plot)


if __name__ == "__main__":
    main()