"""
Minimum-independent-cluster check for GCF carrier sets against the 521-genome tree.

For a given GCF, "independent clusters" = the minimum number of monophyletic
groups needed to contain all carrier genomes -- i.e. the maximal subtrees whose
tips are entirely carriers. This tells apart "N independent acquisition/loss
events" from "one clade of close relatives counted N times".
"""
import csv
import sys
from Bio import Phylo

TREE_PATH = "analysis/phylogenetics/results/gubbins/full_hybrid_521/full_hybrid_521.final_tree.tre"
CARRIERS_PATH = "analysis/phylogenetics/results/gcf_roles_followup/significant_gcf_carriers.tsv"
OUT_TSV = "analysis/phylogenetics/results/gcf_roles_followup/independent_cluster_check.tsv"


def pure_clusters(clade, carriers):
    """Return list of maximal clades whose full tip set is a subset of `carriers`."""
    terms = clade.get_terminals()
    names = {t.name for t in terms}
    if names and names.issubset(carriers):
        return [clade]
    if clade.is_terminal():
        return []
    out = []
    for child in clade.clades:
        out.extend(pure_clusters(child, carriers))
    return out


def main():
    tree = Phylo.read(TREE_PATH, "newick")
    all_tips = {t.name for t in tree.get_terminals()}

    gcf_carriers = {}
    with open(CARRIERS_PATH) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            gcf_carriers.setdefault(row["gcf_id"], set()).add(row["sample_id"])

    results = {}
    for gcf, carriers in sorted(gcf_carriers.items()):
        missing = carriers - all_tips
        if missing:
            print(f"WARNING: {gcf} has {len(missing)} carriers missing from tree: {missing}", file=sys.stderr)
        clusters = pure_clusters(tree.root, carriers)
        sizes = sorted((len(c.get_terminals()) for c in clusters), reverse=True)
        results[gcf] = {
            "n_carriers": len(carriers),
            "n_clusters": len(clusters),
            "cluster_sizes": sizes,
            "clusters": clusters,
        }
        print(f"{gcf}: n_carriers={len(carriers)} n_independent_clusters={len(clusters)} cluster_sizes={sizes}")

    with open(OUT_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gcf_id", "n_carriers", "n_clusters", "pct_independent", "cluster_sizes"])
        for gcf, res in results.items():
            pct = 100 * res["n_clusters"] / res["n_carriers"]
            w.writerow([gcf, res["n_carriers"], res["n_clusters"], f"{pct:.1f}",
                        ";".join(str(s) for s in res["cluster_sizes"])])
    print("wrote", OUT_TSV)

    return results, tree, gcf_carriers


if __name__ == "__main__":
    main()
