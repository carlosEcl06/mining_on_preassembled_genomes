#!/usr/bin/env python3
"""Render the final Gubbins tree (521-genome cohort) as a tip-colored rectangular
cladogram, phenotype-coded (cancer vs. non-cancer), using the same palette as the
BSB paper's figures (#d7191c cancer, #2c7bb6 non-cancer).

Run with ~/.venvs/phylo/bin/python (see memory: mint-toolchain).
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[3]
TREE = ROOT / "analysis/phylogenetics/results/gubbins/full_hybrid_521/full_hybrid_521.final_tree.tre"
METADATA = ROOT / "rawdata/metadata/gbk_hp_genomes.tsv"
OUT = ROOT / "analysis/phylogenetics/results/gubbins/full_hybrid_521/full_hybrid_521.final_tree.pdf"

CANCER_COLOR = "#d7191c"
NONCANCER_COLOR = "#2c7bb6"
UNKNOWN_COLOR = "#999999"

phenotype = {}
with open(METADATA, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        is_cancer = row["Phenotype"].strip() == "Gastric cancer"
        phenotype[row["Genome_ID"]] = "cancer" if is_cancer else "non-cancer"

tree = Phylo.read(TREE, "newick")
tree.ladderize()
tips = tree.get_terminals()
n = len(tips)
print(f"Tips in tree: {n}")

counts = {"cancer": 0, "non-cancer": 0, "unknown": 0}
label_colors = {}
for tip in tips:
    grp = phenotype.get(tip.name, "unknown")
    counts[grp] += 1
    label_colors[tip.name] = {
        "cancer": CANCER_COLOR,
        "non-cancer": NONCANCER_COLOR,
        "unknown": UNKNOWN_COLOR,
    }[grp]
print("Phenotype counts among tree tips:", counts)

fig_height = max(11, n * 0.115)
fig, ax = plt.subplots(figsize=(14, fig_height))

Phylo.draw(
    tree,
    axes=ax,
    do_show=False,
    label_colors=label_colors,
    label_func=lambda c: c.name if c.is_terminal() else "",
    branch_labels=None,
)

for text in ax.texts:
    if text.get_text() in label_colors:
        text.set_fontsize(3.3)

ax.set_title(
    f"Gubbins recombination-masked core-genome tree -- 521-genome cohort "
    f"({counts['cancer']} cancer, {counts['non-cancer']} non-cancer)",
    fontsize=12,
)
ax.set_xlabel("Branch length (substitutions per site, recombination-masked)")
ax.set_ylabel("")
ax.set_yticks([])

legend_handles = [
    Line2D([0], [0], marker="o", color="none", markerfacecolor=CANCER_COLOR, markersize=8, label="Cancer"),
    Line2D([0], [0], marker="o", color="none", markerfacecolor=NONCANCER_COLOR, markersize=8, label="Non-cancer"),
]
ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize=10)

fig.tight_layout()
fig.savefig(OUT, dpi=150)
print(f"Wrote {OUT}")
