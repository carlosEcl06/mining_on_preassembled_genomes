# Builds the report figures (pilot timing, tree agreement, GrafGen ancestry, QC). Needs results/ from the pipeline.
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import dendropy
from scipy.stats import pearsonr, spearmanr

ROOT = "/home/clima/lab/projects/bgcs_in_cancer/mt/mining_on_preassembled_genomes"
OUT = ROOT + "/analysis/phylogenetics/results/report_figs"  # results/ is gitignored
os.makedirs(OUT, exist_ok=True)
os.chdir(ROOT)

# validated categorical slots (light) + text/surface tokens
SURF, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8983", "#e4e3df"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GRAY = "#b9b8b1"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 8.5, "text.color": INK,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
    "axes.titlesize": 9, "axes.titleweight": "semibold", "axes.titlelocation": "left",
})

def panel_tag(ax, t):
    ax.text(-0.02, 1.13, t, transform=ax.transAxes, fontsize=9, fontweight="bold", color=INK, va="bottom", ha="right")

def rbar(ax, y, w, h, color, x0=0.0, r=0.0):
    ax.barh(y, w, height=h, left=x0, color=color, linewidth=0, zorder=3)

# ------------------------------------------------------------------ data
anc = pd.read_csv("analysis/phylogenetics/results/grafgen/ancestry.tsv", sep="\t")
meta = pd.read_csv("rawdata/metadata/gbk_hp_genomes.tsv", sep="\t")
meta["Continent"] = meta["Continent"].str.strip(); meta["Phenotype"] = meta["Phenotype"].str.strip()
meta["cancer"] = (meta["Phenotype"] == "Gastric cancer")
d = anc.merge(meta, left_on="Sample", right_on="Genome_ID")
qc = pd.read_csv("analysis/phylogenetics/results/grafgen/qc_callable_fraction.tsv", sep="\t")

# ================================================================ FIG 1: pilots
fig, (a, b) = plt.subplots(1, 2, figsize=(6.9, 2.9), gridspec_kw={"width_ratios": [1.25, 1]}, constrained_layout=True)
n = np.array([50, 175]); t = np.array([2 + 43.6/60 + 0.0, 21 + 57.3/60 + 0.0])  # 2h43m37s, 21h57m17s
t = np.array([2 + 43/60 + 37/3600, 21 + 57/60 + 17/3600])
p = np.log(t[1] / t[0]) / np.log(n[1] / n[0])
xs = np.geomspace(50, 620, 200)
proj = lambda x, e: t[1] * (x / 175) ** e
lo, hi = proj(552, p - 0.15), proj(552, p + 0.15)
a.fill_between(np.geomspace(175, 620, 100), proj(np.geomspace(175, 620, 100), p - 0.15),
               proj(np.geomspace(175, 620, 100), p + 0.15), color=ORANGE, alpha=0.12, linewidth=0, zorder=1)
a.plot(np.geomspace(175, 620, 100), proj(np.geomspace(175, 620, 100), p), color=ORANGE, lw=2, zorder=2)
a.plot(n, t, color=BLUE, lw=2, zorder=2)
for xx, yy in zip(n, t):
    a.plot(xx, yy, "o", ms=8, mfc=BLUE, mec=SURF, mew=2, zorder=4)
a.plot(552, proj(552, p), "o", ms=8, mfc=ORANGE, mec=SURF, mew=2, zorder=4)
a.set_xscale("log"); a.set_yscale("log")
a.set_xticks([50, 100, 175, 300, 552]); a.set_xticklabels(["50", "100", "175", "300", "552"])
a.set_yticks([1, 10, 100, 1000]); a.set_yticklabels(["1 h", "10 h", "100 h", "1000 h"])
a.minorticks_off(); a.set_xlim(40, 700); a.set_ylim(0.8, 600)
a.grid(True, axis="y", zorder=0); a.set_axisbelow(True)
a.set_xlabel("Genomes in the alignment (N)")
a.set_title("Wall time grows faster than N")
a.annotate("Pilot A (N=50):\n2 h 44 min", (50, t[0]), xytext=(9, -4), textcoords="offset points", color=INK, fontsize=8, va="top")
a.annotate("Pilot B (N=175):\n21 h 57 min", (175, t[1]), xytext=(9, -4), textcoords="offset points", color=INK, fontsize=8, va="top")
a.annotate(f"N=552: ~{proj(552, p):.0f} h\n(≈{proj(552, p)/24:.0f} days)", (552, proj(552, p)), xytext=(-6, 12), textcoords="offset points",
           color=INK, fontsize=8, ha="right", va="bottom")
a.text(0.02, 0.97, f"fitted exponent ≈ {p:.2f}\nband: ±0.15 (illustrative)", transform=a.transAxes, va="top", fontsize=7.5, color=INK2)
panel_tag(a, "a")

its = ["1", "2", "3", "4", "5"]
dur = [5.79, 5.78, 1.73, 3.10, 5.55]
bw = 0.55
for i, v in enumerate(dur):
    b.bar(i, v, width=bw, color=BLUE, linewidth=0, zorder=3)
    b.text(i, v + 0.15, f"{v:.1f} h", ha="center", va="bottom", fontsize=8, color=INK)
b.set_xticks(range(5)); b.set_xticklabels(its)
b.set_ylim(0, 7.4); b.set_yticks([0, 2, 4, 6]); b.set_yticklabels(["0", "2 h", "4 h", "6 h"])
b.grid(True, axis="y", zorder=0); b.set_axisbelow(True)
b.set_xlabel("Gubbins iteration (pilot B, N=175)")
b.set_title("Per-iteration time is irregular")
panel_tag(b, "b")
fig.savefig(f"{OUT}/fig1_pilots.png", dpi=220); plt.close(fig)

# ================================================================ FIG 2: tree agreement + recombination
gt = "analysis/phylogenetics/results/gubbins/pilot_a/pilot_a.final_tree.tre"
pt = "analysis/phylogenetics/results/panaroo/panaroo_pilot_a_tree.treefile"
tns = dendropy.TaxonNamespace()
t1 = dendropy.Tree.get(path=gt, schema="newick", taxon_namespace=tns)
t2 = dendropy.Tree.get(path=pt, schema="newick", taxon_namespace=tns)
shared = {l.taxon.label for l in t1.leaf_node_iter()} & {l.taxon.label for l in t2.leaf_node_iter()}
t1.retain_taxa_with_labels(shared); t2.retain_taxa_with_labels(shared)
pd1, pd2 = t1.phylogenetic_distance_matrix(), t2.phylogenetic_distance_matrix()
tx1 = {x.label: x for x in t1.taxon_namespace if x.label in shared}
tx2 = {x.label: x for x in t2.taxon_namespace if x.label in shared}
lab = sorted(shared); x, y = [], []
for i in range(len(lab)):
    for j in range(i + 1, len(lab)):
        x.append(pd1.distance(tx1[lab[i]], tx1[lab[j]])); y.append(pd2.distance(tx2[lab[i]], tx2[lab[j]]))
x, y = np.array(x), np.array(y)
r, rho = pearsonr(x, y)[0], spearmanr(x, y)[0]

fig, (a, b, c) = plt.subplots(1, 3, figsize=(6.9, 2.7), gridspec_kw={"width_ratios": [1.5, 1, 1]}, constrained_layout=True)
a.scatter(x, y, s=5, color=BLUE, alpha=0.35, linewidths=0, zorder=3)
a.set_xlabel("Gubbins tree: pairwise distance (substitutions)"); a.set_ylabel("Panaroo tree: pairwise distance\n(substitutions per site)")
a.grid(True, zorder=0); a.set_axisbelow(True)
a.set_title("Independent trees agree")
a.text(0.04, 0.96, f"r = {r:.3f} (Pearson)\nρ = {rho:.3f} (Spearman)\n{len(x):,} genome pairs", transform=a.transAxes, va="top", fontsize=8, color=INK)
panel_tag(a, "a")
for ax, vals, title, lab_, ymax, fmt in [
    (b, [0.28, 0.57], "SNPs called\nrecombinant", "", 0.75, "{:.2f}%"),
    (c, [54.5, 62.7], "Branches with\n≥1 recombination", "", 80, "{:.0f}%")]:
    for i, v in enumerate(vals):
        ax.bar(i, v, width=0.5, color=BLUE, linewidth=0, zorder=3)
        ax.text(i, v + ymax * 0.02, fmt.format(v), ha="center", va="bottom", fontsize=8, color=INK)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["N=50\n(A)", "N=175\n(B)"])
    ax.set_ylim(0, ymax); ax.grid(True, axis="y", zorder=0); ax.set_axisbelow(True)
    ax.set_title(title)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}%"))
panel_tag(b, "b"); panel_tag(c, "c")
fig.savefig(f"{OUT}/fig2_trees.png", dpi=220); plt.close(fig)
print("scatter r", r, rho, len(x))

# ================================================================ FIG 3: ancestry
fig = plt.figure(figsize=(6.9, 3.5), constrained_layout=True)
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1, 1])
a = fig.add_subplot(gs[:, 0]); b = fig.add_subplot(gs[0, 1]); c = fig.add_subplot(gs[1, 1])

cnt = d["Refpop"].value_counts().sort_values()
names = [s.replace("hpgp", "") for s in cnt.index]
for i, v in enumerate(cnt.values):
    a.barh(i, v, height=0.55, color=BLUE, linewidth=0, zorder=3)
    a.text(v + 3, i, str(v), va="center", fontsize=8, color=INK)
a.set_yticks(range(len(cnt))); a.set_yticklabels(names)
a.set_xlim(0, 240); a.grid(True, axis="x", zorder=0); a.set_axisbelow(True)
a.set_title("Closest reference population\n(GrafGen, n=552)")
a.set_xlabel("Genomes")
panel_tag(a, "a")

order = ["Europe", "North America", "South America", "Asia"]
comp = d.groupby("Continent")[["E_percent", "F_percent", "A_percent"]].mean().loc[order]
ns = d["Continent"].value_counts()
cols = [BLUE, ORANGE, AQUA]; labs = ["European", "African", "Asian"]
for i, cont in enumerate(order[::-1]):
    left = 0
    for k, col in enumerate(["E_percent", "F_percent", "A_percent"]):
        w = comp.loc[cont, col]
        b.barh(i, w - 0.6, height=0.55, left=left + 0.3, color=cols[k], linewidth=0, zorder=3)
        if w >= 9:
            b.text(left + w / 2, i, f"{w:.0f}", ha="center", va="center", fontsize=7.5, color=("#ffffff" if k == 0 else INK))
        left += w
b.set_yticks(range(4)); b.set_yticklabels([f"{c_} (n={ns[c_]})" for c_ in order[::-1]])
b.set_xlim(0, 100); b.set_xticks([0, 50, 100]); b.set_xticklabels(["0", "50", "100%"])
b.set_title("Mean ancestry by reported continent")
for k, (lb, cl) in enumerate(zip(labs, cols)):
    b.plot([], [], "s", color=cl, ms=6, label=lb)
b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=False, fontsize=7.5, handlelength=0.8, columnspacing=1.2)
panel_tag(b, "b")

eu = d[d["Continent"] == "Europe"]
rng = np.random.default_rng(3)
for i, (flag, col, nm) in enumerate([(False, BLUE, "Non-cancer"), (True, ORANGE, "Cancer")]):
    s = eu.loc[eu["cancer"] == flag, "A_percent"].values
    jy = i + rng.uniform(-0.22, 0.22, len(s))
    c.scatter(s, jy, s=9, color=col, alpha=0.6, linewidths=0, zorder=3)
    med = np.median(s)
    c.plot([med, med], [i - 0.32, i + 0.32], color=INK, lw=2, zorder=4)
    c.text(25, i, f"{nm}  n={len(s)}\nmedian {med:.1f}%", va="center", fontsize=7.5, color=INK, clip_on=False)
c.set_yticks([]); c.set_ylim(-0.5, 1.5); c.set_xlim(-1, 24)
c.set_xlabel("Asian-ancestry component (%)")
c.grid(True, axis="x", zorder=0); c.set_axisbelow(True)
c.set_title("European genomes only")
c.spines["left"].set_visible(False)
panel_tag(c, "c")
fig.savefig(f"{OUT}/fig3_ancestry.png", dpi=220); plt.close(fig)

# ================================================================ FIG 4: QC
fig, (a, b) = plt.subplots(1, 2, figsize=(6.9, 2.9), gridspec_kw={"width_ratios": [1.6, 1]}, constrained_layout=True)
q = qc.copy(); q["Refpop"] = q["Refpop"].str.replace("hpgp", "")
med = q.groupby("Refpop")["callable_fraction"].median().sort_values()
rng = np.random.default_rng(5)
for i, rp in enumerate(med.index):
    v = q.loc[q["Refpop"] == rp, "callable_fraction"].values * 100
    a.scatter(v, i + rng.uniform(-0.2, 0.2, len(v)), s=7, color=BLUE, alpha=0.5, linewidths=0, zorder=3)
    a.plot([np.median(v)] * 2, [i - 0.32, i + 0.32], color=INK, lw=2, zorder=4)
a.set_yticks(range(len(med))); a.set_yticklabels([f"{rp} (n={int((q['Refpop']==rp).sum())})" for rp in med.index])
a.set_xlim(66, 88); a.grid(True, axis="x", zorder=0); a.set_axisbelow(True)
a.set_xlabel("Reference bases callable per genome (%)")
a.set_title("Callable genome by ancestry group")
a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}%"))
panel_tag(a, "a")

cats = [("Duplicate deposit\n(same sample ID)", 12), ("Likely same patient", 12), ("Possibly same patient", 11), ("Cross-project, shared\nMexican isolate pool", 5)]
for i, (lb, v) in enumerate(cats[::-1]):
    b.barh(i, v, height=0.55, color=BLUE, linewidth=0, zorder=3)
    b.text(v + 0.3, i, str(v), va="center", fontsize=8, color=INK)
b.set_yticks(range(4)); b.set_yticklabels([c_[0] for c_ in cats[::-1]], fontsize=7.5)
b.set_xlim(0, 15); b.grid(True, axis="x", zorder=0); b.set_axisbelow(True)
b.set_xlabel("Near-duplicate pairs (of 40)")
b.set_title("The 40 near-identical\npairs, by origin")
panel_tag(b, "b")
fig.savefig(f"{OUT}/fig4_qc.png", dpi=220); plt.close(fig)
print("done")
