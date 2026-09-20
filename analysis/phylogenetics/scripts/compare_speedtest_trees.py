#!/usr/bin/env python3
"""
Speed-test verdict: does a faster Gubbins tree builder reproduce the 21h57m
reference run (pilot B, --tree-builder iqtree) closely enough to use for the
full run?

  python3 compare_speedtest_trees.py speed_fasttree speed_iqtreefast

Each argument is a run name; expects, under results/gubbins/speedtest_<suffix>/
(fetched from the server first):
  speed_<suffix>.final_tree.tre and speed_<suffix>.per_branch_statistics.csv
Fetch with:
  rsync -a -e "ssh -p 2222 -i ~/.ssh/id_ed25519_labserver" \
    'clima@143.107.143.27:/data2/projects/LGMB-009/NGS/analysis/phylogenetics/results/gubbins/speedtest_fasttree/speed_fasttree.{final_tree.tre,per_branch_statistics.csv,log}' \
    analysis/phylogenetics/results/gubbins/speedtest_fasttree/
Reference: results/gubbins/pilot_b/pilot_b.{final_tree.tre,per_branch_statistics.csv}.

Decision bar (same as the Panaroo cross-check): patristic Pearson r >= 0.98, and a
recombination fraction of the same order as the reference (0.57% of SNPs, 62.7% of
branches with >=1 block). Wall time is read from the PBS out file, not here.
"""
import csv
import sys
import numpy as np
import dendropy
from dendropy.calculate import treecompare
from scipy.stats import pearsonr, spearmanr

BASE = "analysis/phylogenetics/results/gubbins"
REF_TREE = f"{BASE}/pilot_b/pilot_b.final_tree.tre"
REF_STATS = f"{BASE}/pilot_b/pilot_b.per_branch_statistics.csv"
R_BAR = 0.98


def recomb_summary(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    snps = sum(int(r["Total SNPs"]) for r in rows)
    inrec = sum(int(r["Number of SNPs Inside Recombinations"]) for r in rows)
    blocks = sum(int(r["Number of Recombination Blocks"]) for r in rows)
    withb = sum(1 for r in rows if int(r["Number of Recombination Blocks"]) > 0)
    return dict(branches=len(rows), pct_snps=100 * inrec / snps, blocks=blocks, pct_branches=100 * withb / len(rows))


def patristic(tree, labels):
    pdm = tree.phylogenetic_distance_matrix()
    tx = {t.label: t for t in tree.taxon_namespace if t.label in labels}
    lab = sorted(labels)
    return np.array([pdm.distance(tx[lab[i]], tx[lab[j]]) for i in range(len(lab)) for j in range(i + 1, len(lab))])


def compare(name):
    d = f"{BASE}/speedtest_{name.replace('speed_', '')}"
    tns = dendropy.TaxonNamespace()
    ref = dendropy.Tree.get(path=REF_TREE, schema="newick", taxon_namespace=tns)
    tst = dendropy.Tree.get(path=f"{d}/{name}.final_tree.tre", schema="newick", taxon_namespace=tns)
    shared = {l.taxon.label for l in ref.leaf_node_iter()} & {l.taxon.label for l in tst.leaf_node_iter()}
    ref.retain_taxa_with_labels(shared); tst.retain_taxa_with_labels(shared)
    x, y = patristic(ref, shared), patristic(tst, shared)
    r, rho = pearsonr(x, y)[0], spearmanr(x, y)[0]
    ref.encode_bipartitions(); tst.encode_bipartitions()
    rf = treecompare.symmetric_difference(ref, tst)
    rf_norm = rf / (2 * (len(shared) - 3))
    rec = recomb_summary(f"{d}/{name}.per_branch_statistics.csv")
    print(f"\n=== {name} vs reference (pilot B, iqtree), {len(shared)} shared taxa ===")
    print(f"patristic Pearson r = {r:.4f}   Spearman = {rho:.4f}   normalized RF = {rf_norm:.3f}")
    print(f"recombination: {rec['pct_snps']:.2f}% of SNPs, {rec['blocks']} blocks, "
          f"{rec['pct_branches']:.1f}% of branches (reference: 0.57%, 506 blocks, 62.7%)")
    print("PASS (r >= %.2f)" % R_BAR if r >= R_BAR else "FAIL (r < %.2f)" % R_BAR)


if __name__ == "__main__":
    ref_rec = recomb_summary(REF_STATS)
    print(f"reference: {ref_rec['pct_snps']:.2f}% SNPs, {ref_rec['blocks']} blocks, {ref_rec['pct_branches']:.1f}% branches")
    for n in sys.argv[1:]:
        compare(n)
