"""
Volcano plot of the 20 clade-significant amino-acid substitutions from
make_aa_substitution_table.py: log2(odds ratio, Clade A : Clade B) on x,
-log10(BH p_adj) on y, labeled by their table `id` (matching
aa_substitutions_table.tsv's `id` column 1:1).

Clade direction is already encoded by the x-axis sign (positive = Clade A,
negative = Clade B), so coloring points by clade on top of that would just be
redundant with an axis the plot already has. Points are colored by GENE
instead (metG vs speA) -- using this project's secondary palette colors
(aqua/gold), deliberately NOT blue/orange, since those are already the
established clade colors elsewhere in this project and reusing them here
would misleadingly suggest the color means clade again.
"""
import math
import sys

sys.path.insert(0, "analysis/phylogenetics/scripts")
import matplotlib.pyplot as plt

import make_aa_substitution_table as tbl

OUT = "analysis/phylogenetics/results/gcf_roles_followup/figs/fig_aa_substitution_volcano.png"

GENE_COLORS = {
    "metG": "#1baf7a",  # aqua
    "speA": "#eda100",  # gold/yellow
}
INK = "#2b2b2b"
MUTED = "#8a8a8a"
SURFACE = "#fcfcfb"
ALPHA = 0.05

# marker shape by physicochemical-change category -- a 4-level combination of the
# table's polarity_change/charge_change flags, independent of (and layered on top
# of) the gene color above.
SHAPE_BY_CHANGE = {
    ("no", "no"): ("o", "Conservative (no polarity/charge change)"),
    ("yes", "no"): ("^", "Polarity change only"),
    ("no", "yes"): ("D", "Charge change only"),
    ("yes", "yes"): ("s", "Polarity + charge change"),
}


def main():
    rows = tbl.main()  # recomputes + rewrites the TSV, returns the same rows

    fig, ax = plt.subplots(figsize=(8, 6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    # slight per-shape size bump so triangle/diamond/square read as visually
    # similar in "ink area" to the circle at the same nominal marker size
    SHAPE_SIZE = {"o": 46, "^": 60, "D": 42, "s": 42}

    for r in rows:
        x = r["log2_odds_ratio_A_vs_B"]
        y = r["neg_log10_p_adj"]
        color = GENE_COLORS[r["gene"]]
        shape, _ = SHAPE_BY_CHANGE[(r["polarity_change"], r["charge_change"])]
        ax.scatter(x, y, s=SHAPE_SIZE[shape], marker=shape, color=color, edgecolor="white",
                   linewidth=0.6, zorder=3)
        # offset labels slightly so they don't sit on top of the marker
        ax.annotate(str(r["id"]), (x, y), xytext=(4, 3), textcoords="offset points",
                    fontsize=7.3, color=INK, zorder=4)

    thresh_y = -math.log10(ALPHA)
    ax.axhline(thresh_y, color=MUTED, linewidth=0.9, linestyle="--", zorder=1)
    ax.axvline(0, color=MUTED, linewidth=0.8, linestyle=":", zorder=1)

    xmax = max(abs(r["log2_odds_ratio_A_vs_B"]) for r in rows) * 1.15
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(0, max(r["neg_log10_p_adj"] for r in rows) * 1.2)

    ax.set_xlabel("log2(odds ratio), Clade A : Clade B carriers", fontsize=10, color=INK)
    ax.set_ylabel("-log10(BH-adjusted p-value)", fontsize=10, color=INK)
    ax.set_title("Clade-significant amino-acid substitutions (metG + speA)", fontsize=12,
                 color=INK, pad=12)

    ax.text(0.99, 0.97, "Clade A enriched →", transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, color=INK)
    ax.text(0.01, 0.97, "← Clade B enriched", transform=ax.transAxes, ha="left", va="top",
            fontsize=8.5, color=INK)
    ax.text(0.01, (thresh_y / ax.get_ylim()[1]) + 0.015, f"p_adj = {ALPHA}", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=7.5, color=MUTED)

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=8.5)

    gene_handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=GENE_COLORS["metG"], markersize=7,
                   label="metG"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=GENE_COLORS["speA"], markersize=7,
                   label="speA"),
    ]
    gene_legend = ax.legend(handles=gene_handles, title="Gene", loc="upper center",
                             bbox_to_anchor=(0.27, -0.22), ncol=2, frameon=False, fontsize=9,
                             title_fontsize=9)
    ax.add_artist(gene_legend)

    shape_handles = [
        plt.Line2D([0], [0], marker=shape, color="none", markerfacecolor=MUTED, markeredgecolor="white",
                   markersize=8 if shape != "^" else 9, label=desc)
        for shape, desc in SHAPE_BY_CHANGE.values()
    ]
    ax.legend(handles=shape_handles, title="Residue property change", loc="upper center",
              bbox_to_anchor=(0.75, -0.22), ncol=1, frameon=False, fontsize=8.5, title_fontsize=9)

    plt.tight_layout()
    plt.savefig(OUT, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    print("wrote", OUT, file=sys.stderr)


if __name__ == "__main__":
    main()
