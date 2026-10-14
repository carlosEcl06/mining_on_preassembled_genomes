# Figures for Progress Report #4 (Gubbins tree-builder speed tests). Needs results/gubbins/speedtest_*/ fetched from the server.
# Run with ~/.venvs/phylo/bin/python. Same palette/typography as make_report_figs.py; one fixed colour per run in every figure.
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import dendropy
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[3]
os.chdir(ROOT)
G = "analysis/phylogenetics/results/gubbins"
OUT = ROOT / "analysis/phylogenetics/results/report4/figs"  # results/ is gitignored
OUT.mkdir(parents=True, exist_ok=True)

SURF, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8983", "#e4e3df"
BLUE, ORANGE, AQUA, YELLOW, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a8983"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "DejaVu Sans", "font.size": 8.5, "text.color": INK,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": GRID, "grid.linewidth": 0.8, "axes.titlesize": 9, "axes.titleweight": "bold",
    "axes.titlelocation": "left",
})

# entity -> (label, colour); the colour follows the run in every figure
RUNS = {
    "ref":   ("Pilot B, IQ-TREE only (reference)", GRAY),
    "hyb":   ("Hybrid: FastTree + RAxML", AQUA),
    "fast":  ("FastTree only", BLUE),
    "hyq":   ("FastTree + IQ-TREE", ORANGE),
    "iqf":   ("IQ-TREE --fast", YELLOW),
}
H = 3600.0
WALL = {"ref": 21 + 57 / 60 + 17 / 3600, "hyb": 5 + 16 / 60 + 45 / 3600, "fast": 7614.74 / H,
        "hyq": 2 + 3 / 60 + 55 / 3600, "iqf": 36 / 60 + 12 / 3600}
CRASHED = {"hyq": "crashed in iteration 3", "iqf": "crashed in iteration 5"}
ITER_H = {  # hours per Gubbins iteration, from each run's PBS .out
    "ref": [5.79, 5.78, 1.73, 3.10, 5.55],
    "fast": [0.434, 0.419, 0.423, 0.418, 0.422],
    "hyb": [0.367, 1.476, 1.312, 1.107, 1.014],
    "hyq": [0.374, 1.678],
    "iqf": [0.133, 0.122, 0.169, 0.156],
}

# ------------------------------------------------ FIG 1: time and per-iteration time
fig, (a, b) = plt.subplots(1, 2, figsize=(6.9, 3.0), gridspec_kw={"width_ratios": [1.15, 1]}, constrained_layout=True)
order = ["ref", "hyb", "fast", "hyq", "iqf"]
ys = np.arange(len(order))[::-1]
for y, k in zip(ys, order):
    lab, col = RUNS[k]
    crashed = k in CRASHED
    a.barh(y, WALL[k], height=0.55, color=SURF if crashed else col, edgecolor=col, linewidth=1.4 if crashed else 0,
           hatch="////" if crashed else None, zorder=3)
    txt = f"{WALL[k]:.1f} h" if WALL[k] >= 1 else f"{WALL[k]*60:.0f} min"
    if crashed:
        txt += f"  ({CRASHED[k]})"
    a.text(WALL[k] + 0.4, y, txt, va="center", ha="left", fontsize=8, color=INK2)
a.set_yticks(ys); a.set_yticklabels([RUNS[k][0] for k in order], fontsize=8)
a.set_xlim(0, 33); a.set_xlabel("Wall time on the 175-genome alignment (h)")
a.grid(axis="x", zorder=0); a.tick_params(axis="y", length=0)
a.set_title("Wall time")
for k in order:
    lab, col = RUNS[k]
    xs = np.arange(1, len(ITER_H[k]) + 1)
    b.plot(xs, ITER_H[k], color=col, lw=2, zorder=3, solid_capstyle="round")
    b.plot(xs, ITER_H[k], "o", ms=5.5, mfc=col, mec=SURF, mew=1.5, zorder=4)
    if k in CRASHED:
        b.plot(xs[-1] + 0.35, ITER_H[k][-1], marker="x", ms=7, color=INK2, mew=1.6, zorder=5)
b.set_yscale("log"); b.set_yticks([0.1, 0.3, 1, 3, 10]); b.set_yticklabels(["0.1", "0.3", "1", "3", "10"])
b.set_xticks(range(1, 6)); b.set_xlim(0.7, 5.3)
b.set_xlabel("Gubbins iteration"); b.set_ylabel("Hours per iteration (log)")
b.grid(axis="y", zorder=0)
lab_at = {"ref": (5, "reference"), "hyb": (5, "FastTree+RAxML"), "fast": (5, "FastTree"),
          "hyq": (2, "FastTree+IQ-TREE"), "iqf": (4, "IQ-TREE --fast")}
off = {"ref": (0, 8), "hyb": (0, -14), "fast": (0, 8), "hyq": (8, 12), "iqf": (0, -13)}
for k, (it, t) in lab_at.items():
    b.annotate(t, (it, ITER_H[k][it - 1]), xytext=off[k], textcoords="offset points", fontsize=7.5, color=INK2,
               ha="left" if k == "hyq" else ("right" if k in ("ref", "hyb", "fast") else "center"))
b.set_title("Time per iteration")
fig.savefig(OUT / "fig1_time.png", dpi=220); plt.close(fig)

# ------------------------------------------------ FIG 2: patristic agreement with the reference
REF = f"{G}/pilot_b/pilot_b.final_tree.tre"
def pairs(path):
    tns = dendropy.TaxonNamespace()
    ref = dendropy.Tree.get(path=REF, schema="newick", taxon_namespace=tns)
    tst = dendropy.Tree.get(path=path, schema="newick", taxon_namespace=tns)
    sh = {l.taxon.label for l in ref.leaf_node_iter()} & {l.taxon.label for l in tst.leaf_node_iter()}
    ref.retain_taxa_with_labels(sh); tst.retain_taxa_with_labels(sh)
    out = []
    for t in (ref, tst):
        pdm = t.phylogenetic_distance_matrix()
        tx = {x.label: x for x in t.taxon_namespace if x.label in sh}
        lab = sorted(sh)
        out.append(np.array([pdm.distance(tx[lab[i]], tx[lab[j]]) for i in range(len(lab)) for j in range(i + 1, len(lab))]))
    return out
panels = [("hyb", f"{G}/speedtest_hybrid_raxml/speed_hybrid_raxml.final_tree.tre", "FastTree + RAxML (hybrid)"),
          ("fast", f"{G}/speedtest_fasttree/speed_fasttree.final_tree.tre", "FastTree only")]
fig, axs = plt.subplots(1, 2, figsize=(6.9, 3.2), constrained_layout=True, sharex=True, sharey=True)
for ax, (k, path, title) in zip(axs, panels):
    x, y = pairs(path)
    r = pearsonr(x, y)[0]
    ax.plot([0, x.max()], [0, x.max()], color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
    ax.scatter(x, y, s=3, color=RUNS[k][1], alpha=0.12, linewidths=0, zorder=2, rasterized=True)
    ax.set_title(title)
    slope = (x @ y) / (x @ x)
    ax.text(0.04, 0.93, f"Pearson r = {r:.3f}\n{len(x):,} genome pairs\nscale = {slope:.2f} x reference (same at all distances)",
            transform=ax.transAxes, va="top", fontsize=8, color=INK)
    ax.set_xlabel("Patristic distance, reference tree (SNPs)")
    ax.grid(zorder=0)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v/1000)}k")); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v/1000)}k"))
axs[0].set_ylabel("Patristic distance, speed-test tree (SNPs)")
fig.savefig(OUT / "fig2_agreement.png", dpi=220); plt.close(fig)

# ------------------------------------------------ FIG 3: the branch-length scale blow-up
def total_len(path):
    t = dendropy.Tree.get(path=path, schema="newick")
    return sum(e.length for e in t.postorder_edge_iter() if e.length is not None)
series = {
    "hyq": [f"{G}/speedtest_hybrid/pilot_b_n175.iteration_{i}.tre" for i in (1, 2)],
    "iqf": [f"{G}/speedtest_iqtreefast/pilot_b_n175.iteration_{i}.tre" for i in (1, 2, 3, 4)],
    "fast": [f"{G}/speedtest_fasttree/pilot_b_n175.iteration_{i}.tre" for i in (1, 2, 3, 4)] + [f"{G}/speedtest_fasttree/speed_fasttree.final_tree.tre"],
    "hyb": [f"{G}/speedtest_hybrid_raxml/pilot_b_n175.iteration_{i}.tre" for i in (1, 2, 3, 4)] + [f"{G}/speedtest_hybrid_raxml/speed_hybrid_raxml.final_tree.tre"],
}
ref_total = total_len(REF)
fig, ax = plt.subplots(figsize=(6.9, 3.1), constrained_layout=True)
ax.axhline(ref_total, color=GRAY, lw=1.2, ls=(0, (4, 3)), zorder=1)
vals = {}
for k, paths in series.items():
    v = [total_len(p) for p in paths]; vals[k] = v
    xs = np.arange(1, len(v) + 1)
    ax.plot(xs, v, color=RUNS[k][1], lw=2, zorder=3, solid_capstyle="round")
    ax.plot(xs, v, "o", ms=5.5, mfc=RUNS[k][1], mec=SURF, mew=1.5, zorder=4)
    if k in CRASHED:
        ax.plot(xs[-1] + 0.3, v[-1] * (1 if k == "iqf" else 1), marker="x", ms=7, color=INK2, mew=1.6, zorder=5)
ax.set_yscale("log"); ax.set_xlim(0.7, 7.0); ax.set_xticks(range(1, 6))
ax.set_ylim(3e5, 2e17)
ax.set_xlabel("Gubbins iteration (tree written at the end of the iteration)")
ax.set_ylabel("Sum of branch lengths (log)")
ax.grid(axis="y", zorder=0)
ax.annotate("FastTree + IQ-TREE", (2, vals["hyq"][1]), xytext=(-6, 8), textcoords="offset points", ha="right", fontsize=7.5, color=INK2)
ax.annotate("IQ-TREE --fast", (4, vals["iqf"][3]), xytext=(-8, -3), textcoords="offset points", ha="right", va="top", fontsize=7.5, color=INK2)
def side_label(txt, y_data, x_from, y_frac):
    ax.annotate(txt, (x_from, y_data), xytext=(5.45, y_frac), textcoords=("data", "axes fraction"), fontsize=7.5, color=INK2,
                va="center", ha="left", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7, shrinkA=1, shrinkB=3))
side_label("reference final tree", ref_total, 5.0, 0.24)
side_label("FastTree + RAxML", vals["hyb"][4], 5.0, 0.14)
side_label("FastTree only", vals["fast"][4], 5.0, 0.04)
ax.set_title("Sum of branch lengths per iteration")
fig.savefig(OUT / "fig3_scale.png", dpi=220); plt.close(fig)
print("wrote", sorted(p.name for p in OUT.glob("*.png")))
print({k: [f"{x:.3g}" for x in v] for k, v in vals.items()}, f"ref {ref_total:.3g}")
