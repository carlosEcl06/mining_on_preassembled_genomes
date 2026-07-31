#!/usr/bin/env python3
"""
discover_deepbgc_backbone_domains.py

Diagnostic tool (run manually, NOT part of run_clustering.sh): empirically
discovers candidate backbone Pfam domains per deepBGC Product_class, by
aggregating each sample's own {sample}.pfam.tsv (deepBGC's internal Pfam
domain calls) and computing per-class domain prevalence — the same
"find the dominant, specific marker gene" methodology used manually via
hmmscan spot-checks for the antiSMASH classes (Terpene-precursor -> PF00348,
Azole-containing-RiPP -> PF02624), but automated here since deepBGC's own
Pfam calls are already clean PFxxxxx accessions (no NAME/ACC inconsistency
like comBGC's PFAM_domains column — see 01_extract_domains.py's docstring).

This does NOT modify BACKBONE_DOMAINS in 03_backbone_identity.py — it only
produces a candidates table for review. Populating BACKBONE_DOMAINS from
these results should still involve checking that the top domain makes
biological sense for that Product_class, the same discipline applied to
every class so far this project (the Terpene "false positive" episode
happened precisely because that check was skipped once).

For each BGC row (Prediction_tool == deepBGC), this script:
  1. Reads that sample's {sample_id}.pfam.tsv (cached per sample).
  2. Selects domain hits whose gene coordinates overlap [BGC_start, BGC_end]
     on the matching contig (sequence_id == contig_id) — a standard interval
     overlap test (gene_start < BGC_end AND gene_end > BGC_start).
  3. Cross-checks this against the file's own in_cluster flag and counts how
     often they disagree, as a sanity check that the coordinate systems
     actually line up (rather than silently trusting a possible mismatch).
  4. Records which Pfam accessions were found for that BGC.

Then, per Product_class, computes what fraction of BGCs in that class
contain each Pfam accession, and writes a candidates TSV. Crucially, this
also computes a SPECIFICITY score per (class, domain): own-class prevalence
minus the highest prevalence that same domain reaches in any OTHER class.
High within-class prevalence alone is not sufficient evidence of a genuine
backbone marker — a first version of this analysis (2026-07-30) found that
several top "Saccharide"/"Other"/"Unknown" candidates were generic,
promiscuous tailoring-enzyme domains (a family of Methyltransferase_N Pfam
entries) shared at similarly high prevalence across those three unrelated
classes, echoing the exact false-positive pattern already caught once for
antiSMASH's Terpene markers (see 01_extract_domains.py's docstring). Results
are sorted by specificity, not raw prevalence, so genuinely class-specific
candidates surface first even if their raw prevalence is lower than a
promiscuous domain's.

Usage:
    python discover_deepbgc_backbone_domains.py \
        --domains     ../results/bgc_domain_arrays.tsv \
        --deepbgc-dir ../../funcscan/results/bgc/deepbgc \
        --output      ../results/deepbgc_backbone_candidates.tsv \
        --min-prevalence 0.3
"""

import argparse
import pathlib
from collections import defaultdict

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--domains",     required=True, help="Path to bgc_domain_arrays.tsv")
    parser.add_argument("--deepbgc-dir", required=True,
                        help="Path to funcscan deepBGC results dir "
                             "(contains {sample_id}/{sample_id}.pfam.tsv)")
    parser.add_argument("--output",      required=True, help="Path to output candidates TSV")
    parser.add_argument("--min-prevalence", type=float, default=0.30,
                        help="Only report accessions present in at least this fraction "
                             "of BGCs within a class (default: 0.30)")
    return parser.parse_args()


def load_sample_pfam(pfam_tsv_path: pathlib.Path):
    if not pfam_tsv_path.exists():
        return None
    try:
        df = pd.read_csv(pfam_tsv_path, sep="\t", dtype=str)
    except pd.errors.EmptyDataError:
        return None
    for col in ("gene_start", "gene_end"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def main():
    args = parse_args()

    print(f"[discover] Reading {args.domains} ...", flush=True)
    df = pd.read_csv(args.domains, sep="\t", dtype=str)
    df = df[df["Prediction_tool"] == "deepBGC"].copy()
    df["BGC_start"] = pd.to_numeric(df["BGC_start"], errors="coerce")
    df["BGC_end"]   = pd.to_numeric(df["BGC_end"], errors="coerce")
    print(f"[discover] {len(df)} deepBGC BGCs to scan.", flush=True)

    deepbgc_dir = pathlib.Path(args.deepbgc_dir)

    pfam_cache: dict = {}
    class_domain_bgc_counts: dict = defaultdict(lambda: defaultdict(int))
    class_total_bgcs: dict = defaultdict(int)

    n_missing_tsv = 0
    n_no_hits = 0
    n_mismatch_warned = 0

    for row in df.itertuples(index=False):
        sample_id = row.sample_id
        cls = row.Product_class if pd.notna(row.Product_class) else "Unknown"
        class_total_bgcs[cls] += 1

        if sample_id not in pfam_cache:
            pfam_cache[sample_id] = load_sample_pfam(
                deepbgc_dir / sample_id / f"{sample_id}.pfam.tsv"
            )
        pfam_df = pfam_cache[sample_id]
        if pfam_df is None:
            n_missing_tsv += 1
            continue

        sub = pfam_df[pfam_df["sequence_id"] == row.contig_id]
        if pd.notna(row.BGC_start) and pd.notna(row.BGC_end):
            in_window = sub[
                (sub["gene_start"] < row.BGC_end) & (sub["gene_end"] > row.BGC_start)
            ]
        elif "in_cluster" in sub.columns:
            in_window = sub[sub["in_cluster"] == "1"]
        else:
            in_window = sub.iloc[0:0]

        if len(in_window) == 0:
            n_no_hits += 1
            continue

        # Sanity cross-check against the file's own in_cluster flag
        if "in_cluster" in in_window.columns:
            disagreement = (in_window["in_cluster"] != "1").mean()
            if disagreement > 0.5:
                n_mismatch_warned += 1

        domains_here = set(in_window["pfam_id"].dropna())
        for dom in domains_here:
            class_domain_bgc_counts[cls][dom] += 1

    print(f"[discover] {n_missing_tsv} BGCs skipped (missing/empty pfam.tsv).", flush=True)
    print(f"[discover] {n_no_hits} BGCs had no domain hits in their coordinate window.",
          flush=True)
    if n_mismatch_warned:
        print(f"[discover] WARNING: {n_mismatch_warned} BGCs had >50% of their "
              f"coordinate-matched hits NOT flagged in_cluster=1 by deepBGC's own "
              f"output — coordinate matching may be unreliable for these; "
              f"inspect before trusting results.", flush=True)

    # --- build full prevalence matrix (all classes x all domains, unfiltered) ---
    # needed so we can compute, for each candidate, how prevalent it is in
    # OTHER classes too — high within-class prevalence alone isn't enough
    # evidence of a genuine backbone marker if the same domain is equally
    # common elsewhere (see module docstring: this is exactly the trap that
    # produced false "Terpene" markers for Other/Saccharide BGCs previously).
    prevalence_matrix: dict[str, dict[str, float]] = {}
    for cls, domain_counts in class_domain_bgc_counts.items():
        total = class_total_bgcs[cls]
        prevalence_matrix[cls] = {
            dom: (n / total if total else 0.0) for dom, n in domain_counts.items()
        }

    rows_out = []
    for cls, domain_counts in class_domain_bgc_counts.items():
        total = class_total_bgcs[cls]
        for dom, n_bgcs in domain_counts.items():
            own_prevalence = n_bgcs / total if total else 0.0
            if own_prevalence < args.min_prevalence:
                continue

            other_prevalences = {
                other_cls: prevalence_matrix[other_cls].get(dom, 0.0)
                for other_cls in prevalence_matrix
                if other_cls != cls
            }
            if other_prevalences:
                other_max_class = max(other_prevalences, key=other_prevalences.get)
                other_max_prevalence = other_prevalences[other_max_class]
            else:
                other_max_class, other_max_prevalence = None, 0.0

            rows_out.append({
                "Product_class": cls,
                "pfam_id": dom,
                "n_bgcs_with_domain": n_bgcs,
                "n_bgcs_total_in_class": total,
                "prevalence": round(own_prevalence, 3),
                "other_max_prevalence": round(other_max_prevalence, 3),
                "other_max_class": other_max_class,
                "specificity": round(own_prevalence - other_max_prevalence, 3),
            })

    out_df = pd.DataFrame(rows_out)
    if not out_df.empty:
        # Sort by specificity within class, not raw prevalence — a domain
        # that's 100% prevalent in its own class but ALSO 90%+ prevalent in
        # some other class is a promiscuous tailoring enzyme, not a backbone
        # marker, and should sort BELOW a domain that's 60% prevalent here
        # and near-0% everywhere else.
        out_df = out_df.sort_values(["Product_class", "specificity"], ascending=[True, False])

    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"[discover] {len(out_df)} candidate (class, domain) rows written to {args.output}",
          flush=True)
    print("\n[discover] Preview (sorted by specificity, not raw prevalence):", flush=True)
    print(out_df.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()