#!/usr/bin/env python3
"""make_clade_report_figs.py -- every useful figure out of the three root-split
comparison tables (genome_clade_assignments.tsv, gcf_association_with_clade.tsv,
clade_defining_snps_annotated.tsv), for the clade-split report.

Run with ~/.venvs/phylo/bin/python (see memory: mint-toolchain).
"""
import math
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path("/home/carlos/lab/projects/bgcs_in_cancer/mt/mining_on_preassembled_genomes")
DATA = ROOT / "analysis/phylogenetics/results/root_split_comparison"
OUT = DATA / "figures"
OUT.mkdir(exist_ok=True)

# Palette: phenotype/continent colors match the camera-ready paper and the
# circular tree already shown to the user, for visual continuity across documents.
CANCER, NONCANCER = "#d7191c", "#2c7bb6"
CONTINENT_COLORS = {
    "Asia": "#E69F00", "Europe": "#009E73",
    "North America": "#8B4CB8", "South America": "#4D4D4D",
}
ACCENT, MUTED, BLUE = "#eb6834", "#999999", "#2a78d6"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "savefig.facecolor": "#fcfcfb", "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "axes.labelcolor": "#222222",
    "text.color": "#222222", "xtick.color": "#444444", "ytick.color": "#444444",
})

members = pd.read_csv(DATA / "genome_clade_assignments.tsv", sep="\t")
gcf = pd.read_csv(DATA / "gcf_association_with_clade.tsv", sep="\t")
snps = pd.read_csv(DATA / "clade_defining_snps_annotated.tsv", sep="\t")

# --- Fig 1: continent composition by clade ---
cont_order = ["Asia", "Europe", "North America", "South America"]
ct = pd.crosstab(members["clade"], members["continent"])[cont_order]
fig, ax = plt.subplots(figsize=(7, 4.5))
x = np.arange(len(ct.index))
width = 0.19
for i, cont in enumerate(cont_order):
    bars = ax.bar(x + (i - 1.5) * width, ct[cont], width, label=cont, color=CONTINENT_COLORS[cont])
    ax.bar_label(bars, padding=2, fontsize=8.5, color="#222222")
ax.set_xticks(x); ax.set_xticklabels([f"Clade {c}" for c in ct.index])
ax.set_ylabel("Genomes")
ax.set_title("Continent of origin by clade")
ax.legend(frameon=False, ncol=2, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.12))
fig.tight_layout()
fig.savefig(OUT / "fig1_continent_by_clade.png", dpi=200)
plt.close(fig)

# --- Fig 2: phenotype composition by clade ---
pt = pd.crosstab(members["clade"], members["phenotype"])[["Cancer", "Non-cancer"]]
_, p_pheno = fisher_exact(pt.values)
fig, ax = plt.subplots(figsize=(5.5, 4.5))
x = np.arange(len(pt.index))
width = 0.32
for i, (ph, col) in enumerate([("Cancer", CANCER), ("Non-cancer", NONCANCER)]):
    bars = ax.bar(x + (i - 0.5) * width, pt[ph], width, label=ph, color=col)
    ax.bar_label(bars, padding=2, fontsize=9, color="#222222")
ax.set_xticks(x); ax.set_xticklabels([f"Clade {c}" for c in pt.index])
ax.set_ylabel("Genomes")
ax.set_title(f"Phenotype by clade (Fisher p = {p_pheno:.2f}, no association)")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "fig2_phenotype_by_clade.png", dpi=200)
plt.close(fig)

# --- Fig 3: GCF volcano (log2 odds ratio vs -log10 p_adj) ---
def log2_or(row):
    o = row["odds_ratio"]
    if o == 0: return -6.0  # cap for plotting (clade-B-restricted)
    if math.isinf(o): return 6.0  # cap (clade-A-restricted)
    return math.log2(o)

gcf = gcf.copy()
gcf["log2_or"] = gcf.apply(log2_or, axis=1)
gcf["neglog10_p"] = -np.log10(gcf["p_adj"].clip(lower=1e-300))
sig = gcf["p_adj"] < 0.05

fig, ax = plt.subplots(figsize=(7, 5.5))
ax.scatter(gcf.loc[~sig, "log2_or"], gcf.loc[~sig, "neglog10_p"], s=28, color=MUTED, alpha=0.7, label="n.s.")
ax.scatter(gcf.loc[sig, "log2_or"], gcf.loc[sig, "neglog10_p"], s=48, color=ACCENT, label="p_adj < 0.05", zorder=3)
# stack labels vertically when points share (near-)identical coordinates
sig_sorted = gcf.loc[sig].sort_values(["log2_or", "neglog10_p"])
seen_xy = {}
for _, r in sig_sorted.iterrows():
    key = (round(r["log2_or"], 1), round(r["neglog10_p"], 1))
    stack = seen_xy.get(key, 0)
    seen_xy[key] = stack + 1
    ax.annotate(r["gcf_id"], (r["log2_or"], r["neglog10_p"]), textcoords="offset points",
                xytext=(10, 6 - 12 * stack), fontsize=8.5, color="#222222")
ax.axhline(-np.log10(0.05), color="#bbbbbb", linestyle="--", linewidth=1)
ax.set_xlabel("log2(odds ratio), Clade A vs Clade B  [capped at ±6 for 0/Inf cases]")
ax.set_ylabel("-log10(BH-adjusted p)")
ax.set_title("GCF presence/absence association with clade (82 qualifying GCFs)")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "fig3_gcf_volcano.png", dpi=200)
plt.close(fig)

# --- Fig 4: carrier counts for significant GCFs ---
top = gcf.loc[sig].sort_values("p_adj")
fig, ax = plt.subplots(figsize=(7, 4.5))
y = np.arange(len(top))
height = 0.35
ax.barh(y + height/2, top["n_A"], height, label="Clade A", color=MUTED)
ax.barh(y - height/2, top["n_B"], height, label="Clade B", color=BLUE)
ax.set_yticks(y); ax.set_yticklabels(top["gcf_id"])
ax.invert_yaxis()
ax.set_xlabel("Carrier genomes")
ax.set_title("Carrier counts for clade-associated GCFs (p_adj < 0.05)")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "fig4_significant_gcf_carriers.png", dpi=200)
plt.close(fig)

# --- Fig 5: genome-wide distribution of clade-defining SNPs ---
GENOME_LEN = 1_667_867
fig, ax = plt.subplots(figsize=(9, 3.2))
ax.hist(snps["real_pos_NC_000915"], bins=80, range=(0, GENOME_LEN), color=BLUE, edgecolor="#fcfcfb", linewidth=0.3)
ax.set_xlabel("Position on NC_000915 (H. pylori 26695 reference, bp)")
ax.set_ylabel("SNPs / 20.8 kb bin")
ax.set_title("Genome-wide distribution of the 455 clade-defining SNPs")
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v/1e6:.1f}M"))
fig.tight_layout()
fig.savefig(OUT / "fig5_snp_genome_distribution.png", dpi=200)
plt.close(fig)

# --- Fig 6: top genes by SNP count ---
gene_counts = (snps.loc[snps["gene"] != "intergenic"]
               .groupby(["locus_tag", "gene", "product"]).size()
               .reset_index(name="n").sort_values("n", ascending=False).head(15))
labels = [f"{g} ({l})" if g else l for g, l in zip(gene_counts["gene"], gene_counts["locus_tag"])]
fig, ax = plt.subplots(figsize=(8, 5.5))
y = np.arange(len(gene_counts))
ax.barh(y, gene_counts["n"], color=BLUE)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("Clade-defining SNPs in this gene")
ax.set_title("Genes carrying the most clade-defining SNPs")
fig.tight_layout()
fig.savefig(OUT / "fig6_top_genes.png", dpi=200)
plt.close(fig)

# --- Fig 7: genic vs intergenic ---
is_intergenic = snps["gene"] == "intergenic"
counts = [(~is_intergenic).sum(), is_intergenic.sum()]
fig, ax = plt.subplots(figsize=(4.5, 4.5))
bars = ax.bar(["Genic", "Intergenic"], counts, color=[BLUE, MUTED])
ax.bar_label(bars, fontsize=10)
ax.set_ylabel("SNPs")
ax.set_title("Clade-defining SNPs: genic vs. intergenic")
fig.tight_layout()
fig.savefig(OUT / "fig7_genic_vs_intergenic.png", dpi=200)
plt.close(fig)

# --- Fig 8: substitution spectrum (ancestral -> derived) ---
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
order = []
for a in "ACGT":
    for d in "ACGT":
        if a != d:
            order.append((a, d))
change_counts = Counter(zip(snps["ancestral"], snps["derived"]))
vals = [change_counts.get(c, 0) for c in order]
colors = [BLUE if c in TRANSITIONS else ACCENT for c in order]
labels = [f"{a}>{d}" for a, d in order]
n_ts = sum(v for c, v in zip(order, vals) if c in TRANSITIONS)
n_tv = sum(vals) - n_ts
fig, ax = plt.subplots(figsize=(8, 4.5))
bars = ax.bar(labels, vals, color=colors)
ax.bar_label(bars, fontsize=8)
ax.set_ylabel("SNPs")
ax.set_title(f"Substitution spectrum (Ts/Tv = {n_ts/max(n_tv,1):.2f}; blue = transition, orange = transversion)")
fig.tight_layout()
fig.savefig(OUT / "fig8_substitution_spectrum.png", dpi=200)
plt.close(fig)

print("Wrote 8 figures to", OUT)
print(f"Fisher phenotype-by-clade p = {p_pheno:.4f}")
print(f"Ts/Tv = {n_ts/max(n_tv,1):.3f} ({n_ts} transitions, {n_tv} transversions)")
print(f"Genic: {(~is_intergenic).sum()}, Intergenic: {is_intergenic.sum()}")
