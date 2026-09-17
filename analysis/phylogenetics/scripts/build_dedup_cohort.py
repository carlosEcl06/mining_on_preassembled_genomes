#!/usr/bin/env python3
"""
Defines the deduplicated cohort for the primary GCF-phenotype association
analysis (see plan: near-duplicate pairs found in QC turned out, on NCBI
BioSample follow-up, to be near-entirely non-independent samples -- literal
duplicate deposits of the same isolate, same-patient paired tumor/non-tumor
biopsies, or samples sharing a research biobank across two BioProjects --
not independent biological replicates). Clusters the 40 near-duplicate pairs
into connected components (a pair can chain into a larger near-clonal group,
e.g. 3-4 genomes all within 5 SNPs of each other), keeps one representative
per component (highest Snippy callable-fraction, a QC metric already
computed and a principled, non-arbitrary tie-breaker), and writes the
exclusion list. The full 552-genome cohort remains available as the named
secondary/sensitivity analysis -- this script only produces the exclusion
list, it does not decide which cohort is "correct" on its own.
"""
import csv

BASE = "analysis/phylogenetics/results/grafgen"


def find(parent, x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(parent, a, b):
    ra, rb = find(parent, a), find(parent, b)
    if ra != rb:
        parent[ra] = rb


def main():
    pairs = list(csv.DictReader(open(f"{BASE}/near_duplicates.tsv"), delimiter="\t"))
    parent = {}
    for p in pairs:
        union(parent, p["sample_a"], p["sample_b"])

    clusters = {}
    for node in parent:
        clusters.setdefault(find(parent, node), []).append(node)
    clusters = list(clusters.values())

    callable_frac = {}
    with open(f"{BASE}/qc_callable_fraction.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            callable_frac[row["ID"]] = float(row["callable_fraction"])

    kept, excluded = [], []
    rows = []
    for cluster in clusters:
        cluster_sorted = sorted(cluster, key=lambda g: -callable_frac.get(g, 0))
        representative = cluster_sorted[0]
        kept.append(representative)
        for g in cluster_sorted[1:]:
            excluded.append(g)
        for g in cluster_sorted:
            rows.append({
                "sample_id": g,
                "cluster_size": len(cluster),
                "callable_fraction": callable_frac.get(g, ""),
                "role": "representative (kept)" if g == representative else "excluded (near-duplicate)",
            })

    with open(f"{BASE}/dedup_cohort_exclusions.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample_id", "cluster_size", "callable_fraction", "role"])
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: -r["cluster_size"]))

    print(f"{len(clusters)} clusters from {len(pairs)} pairs, {sum(len(c) for c in clusters)} genomes involved")
    print(f"Kept (representatives): {len(kept)}")
    print(f"Excluded (near-duplicates): {len(excluded)}")
    print(f"Dedup cohort size: 552 - {len(excluded)} = {552 - len(excluded)}")


if __name__ == "__main__":
    main()
