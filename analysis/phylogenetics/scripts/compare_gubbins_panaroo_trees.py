#!/usr/bin/env python3
"""
Cross-check (task #8, per the plan): compares the Gubbins pilot-A tree
(recombination-corrected, Snippy reference-mapped SNPs) against the Panaroo
pilot-A tree (pangenome core-gene alignment, no recombination correction,
built independently of any single reference) on the same 50 study genomes.
Agreement between the two supports the reference-mapping approach; material
disagreement would flag a reference-bias or recombination-masking problem.

Metrics: unweighted Robinson-Foulds distance (topology-only) and a
patristic-distance (cophenetic) correlation (structure-preserving, less
sensitive to individual weakly-supported bipartitions than raw RF).
"""
import dendropy
from dendropy.calculate import treecompare
import numpy as np
from scipy.stats import pearsonr, spearmanr

GUBBINS_TREE = "analysis/phylogenetics/results/gubbins/pilot_a/pilot_a.final_tree.tre"
PANAROO_TREE = "analysis/phylogenetics/results/panaroo/panaroo_pilot_a_tree.treefile"

tns = dendropy.TaxonNamespace()
t1 = dendropy.Tree.get(path=GUBBINS_TREE, schema="newick", taxon_namespace=tns)
t2 = dendropy.Tree.get(path=PANAROO_TREE, schema="newick", taxon_namespace=tns)

gubbins_labels = {leaf.taxon.label for leaf in t1.leaf_node_iter()}
panaroo_labels = {leaf.taxon.label for leaf in t2.leaf_node_iter()}

print(f"Gubbins tree: {len(gubbins_labels)} tips (includes 'Reference' if present)")
print(f"Panaroo tree: {len(panaroo_labels)} tips")

shared = gubbins_labels & panaroo_labels
only_gubbins = gubbins_labels - panaroo_labels
only_panaroo = panaroo_labels - gubbins_labels
print(f"Shared taxa: {len(shared)}")
print(f"Only in Gubbins tree (expected: 'Reference'): {only_gubbins}")
print(f"Only in Panaroo tree: {only_panaroo}")

# Prune both trees down to the shared taxon set for a fair comparison
t1.retain_taxa_with_labels(shared)
t2.retain_taxa_with_labels(shared)
t1.encode_bipartitions()
t2.encode_bipartitions()

rf = treecompare.symmetric_difference(t1, t2)
n_internal_edges = len(shared) - 3  # max possible bipartitions for an unrooted binary tree
print(f"\n=== Robinson-Foulds distance ===")
print(f"RF = {rf} (out of a max possible {2 * n_internal_edges} for {len(shared)} taxa)")
print(f"Normalized RF = {rf / (2 * n_internal_edges):.3f} (0 = identical topology, 1 = maximally different)")

# Patristic distance correlation
pdm1 = t1.phylogenetic_distance_matrix()
pdm2 = t2.phylogenetic_distance_matrix()
taxa = sorted(shared)
taxon_objs1 = {t.label: t for t in t1.taxon_namespace if t.label in shared}
taxon_objs2 = {t.label: t for t in t2.taxon_namespace if t.label in shared}

d1, d2 = [], []
for i in range(len(taxa)):
    for j in range(i + 1, len(taxa)):
        a, b = taxa[i], taxa[j]
        d1.append(pdm1.distance(taxon_objs1[a], taxon_objs1[b]))
        d2.append(pdm2.distance(taxon_objs2[a], taxon_objs2[b]))

d1, d2 = np.array(d1), np.array(d2)
pearson_r, pearson_p = pearsonr(d1, d2)
spearman_r, spearman_p = spearmanr(d1, d2)
print(f"\n=== Patristic distance correlation ({len(d1)} pairs) ===")
print(f"Pearson r  = {pearson_r:.3f} (p={pearson_p:.2e})")
print(f"Spearman r = {spearman_r:.3f} (p={spearman_p:.2e})")
