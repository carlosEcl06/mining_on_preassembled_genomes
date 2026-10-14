"""Schematic of the nested clade substructure and where the five GCFs' carriers fall."""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
YELLOW = "#eda100"
PURPLE = "#7B2D8E"
PINK = "#C2185B"
NAVY = "#1A3C8C"
SURFACE = "#fcfcfb"
INK = "#2b2b2b"
GREY = "#9a9a9a"

fig, ax = plt.subplots(figsize=(11, 6.7), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
ax.set_xlim(0, 11)
ax.set_ylim(0, 6.7)
ax.axis("off")

def box(x, y, w, h, title, lines, fc, ec=INK, fontsize=8.3):
    rect = mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                                    linewidth=1.3, edgecolor=ec, facecolor=fc)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h - 0.16, title, ha="center", va="top", fontsize=fontsize+0.8,
             fontweight="bold", color=INK)
    ax.text(x + w/2, y + h/2 - 0.12, "\n".join(lines), ha="center", va="center",
             fontsize=fontsize, color=INK, linespacing=1.35)
    return (x, y, w, h)

def connect(b1, b2, color=GREY):
    x1, y1, w1, h1 = b1
    x2, y2, w2, h2 = b2
    ax.plot([x1 + w1, x2], [y1 + h1/2, y2 + h2/2], color=color, lw=1.2, zorder=0)

# Root
root = box(0.1, 2.7, 1.0, 0.8, "Root", ["521 genomes", "+ Reference"], "white")

# Clade A / Clade B
cladeA = box(1.5, 4.4, 1.7, 1.1, "Clade A", ["260 genomes", "Europe/NAm/SAm", "no Asia"], "#dbe9fb")
cladeB = box(1.5, 0.7, 1.7, 1.1, "Clade B", ["261 genomes", "Europe/Asia", "almost no NAm/SAm"], "#fbe3d6")
connect(root, cladeA)
connect(root, cladeB)

# Clade A children
aOff = box(3.6, 5.15, 1.9, 0.95, "A offshoot", ["n=13", "Europe 11, NAm 2", "GCF_0030: 0  GCF_0031: 0"], "white")
aCore = box(3.6, 4.0, 1.9, 0.95, "A core", ["n=247", "NAm/Europe/SAm mixed", "(continues below)"], "white")
connect(cladeA, aOff)
connect(cladeA, aCore)

aCoreOff = box(5.9, 4.55, 2.3, 0.85, "A-core offshoot", ["n=11, 100% Europe", "(p=3.9e-6 vs rest)",
                                                           "GCF_0030: 0  GCF_0031: 0"], "white", fontsize=7.8)
aCoreCore = box(5.9, 3.55, 2.3, 0.85, "A-core remainder", ["n=236, NAm/Eur/SAm mixed",
                                                             "GCF_0030: 19  GCF_0031: 22"], "#fde9c8", fontsize=7.8)
connect(aCore, aCoreOff)
connect(aCore, aCoreCore)

# Clade B children
bOff = box(3.6, 1.65, 1.9, 0.95, "B offshoot", ["n=10", "100% Asia, 100%", "Non-cancer",
                                                 "GCF_0017:0 GCF_0019:2 GCF_0025:0"], "white", fontsize=7.6)
bCore = box(3.6, 0.5, 1.9, 0.95, "B core", ["n=251", "Europe/Asia/NAm/SAm", "(continues below)"], "white")
connect(cladeB, bOff)
connect(cladeB, bCore)

bCoreOff = box(5.9, 1.05, 2.3, 0.85, "B-core offshoot", ["n=11, Asia 9 / Europe 2", "(p=6.3e-6 vs rest)",
                                                           "GCF_0017:0 GCF_0019:0 GCF_0025:1"], "white", fontsize=7.6)
bCoreCore = box(5.9, 0.05, 2.3, 0.85, "B-core remainder", ["n=240, Europe/NAm/SAm/Asia mixed",
                                                             "GCF_0017:16 GCF_0019:15 GCF_0025:17"], "#cdeee0", fontsize=7.6)
connect(bCore, bCoreOff)
connect(bCore, bCoreCore)

ax.text(8.6, 5.6, "Carriers of the Clade-A GCFs\nland almost entirely in the\nlarge, mixed-continent\n\"core\" -- not in either\nmonophyletic continent-pure\noffshoot.", fontsize=8.2,
        ha="left", va="top", color=INK,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=GREY, lw=0.8))
ax.text(8.6, 2.0, "Same pattern in Clade B:\ncarriers concentrate in the\n240-genome mixed core, not\nthe Asia-dominated basal\noffshoots that create the\n\"visible substructure\".", fontsize=8.2,
        ha="left", va="top", color=INK,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=GREY, lw=0.8))

ax.text(0.1, 6.45, "Nested clade structure and where the five GCFs' carriers fall", fontsize=12,
        fontweight="bold", color=INK)

plt.tight_layout()
plt.savefig("analysis/phylogenetics/results/gcf_roles_followup/figs/fig_substructure_schematic.png",
            dpi=200, facecolor=SURFACE, bbox_inches="tight")
print("wrote fig_substructure_schematic.png")
