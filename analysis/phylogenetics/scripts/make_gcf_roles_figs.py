"""Figures for the gcf_roles_followup report: independent-cluster-check bar chart."""
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
SURFACE = "#fcfcfb"
INK = "#2b2b2b"

CLUSTER_PATH = "analysis/phylogenetics/results/gcf_roles_followup/independent_cluster_check.tsv"
OUT_PATH = "analysis/phylogenetics/results/gcf_roles_followup/figs/fig_independent_clusters.png"

rows = []
with open(CLUSTER_PATH) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        rows.append(row)

rows.sort(key=lambda r: -int(r["n_carriers"]))

gcfs = [r["gcf_id"] for r in rows]
carriers = [int(r["n_carriers"]) for r in rows]
clusters = [int(r["n_clusters"]) for r in rows]

fig, ax = plt.subplots(figsize=(7.2, 4.2), facecolor=SURFACE)
ax.set_facecolor(SURFACE)

x = range(len(gcfs))
w = 0.38
b1 = ax.bar([i - w/2 for i in x], carriers, width=w, color=BLUE, label="Carrier genomes")
b2 = ax.bar([i + w/2 for i in x], clusters, width=w, color=ORANGE, label="Independent clusters")

for i, (c, k) in enumerate(zip(carriers, clusters)):
    ax.text(i - w/2, c + 0.4, str(c), ha="center", va="bottom", fontsize=9, color=INK)
    ax.text(i + w/2, k + 0.4, str(k), ha="center", va="bottom", fontsize=9, color=INK)
    pct = 100 * k / c
    ax.text(i, max(c, k) + 2.6, f"{pct:.0f}%\nindependent", ha="center", va="bottom",
            fontsize=8, color=AQUA if pct >= 85 else YELLOW, fontweight="bold")

ax.set_xticks(list(x))
ax.set_xticklabels(gcfs)
ax.set_ylabel("Genomes")
ax.set_ylim(0, max(carriers) + 7)
ax.set_title("Carrier count vs. minimum independent clusters on the tree", fontsize=11, color=INK)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
ax.legend(frameon=False, loc="upper right", fontsize=9)
ax.tick_params(colors=INK)
ax.yaxis.label.set_color(INK)
ax.title.set_color(INK)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=200, facecolor=SURFACE)
print("wrote", OUT_PATH)
