"""
Informative-nucleotide-positions alignment figure (classic "MEGA/Sequencher-style"
variable-sites table) for two of the five clade-significant GCF loci, with the
521-genome tree on the left. Rendered with PIL (direct bitmap drawing) rather
than matplotlib, because one ax.text() call per cell (~130,000 calls for a
522-row x 203-col matrix) froze the machine on the first attempt.

Scope (chosen after checking what's actually in the core-genome alignment):
 - GCF_0019 & GCF_0031 share one locus (confirmed by k-mer mapping against the
   reference): cfaS/metG/ferrochelatase/cmoB, NC_000915:429330-434418. Panel 1
   tests all 144 variable positions in the metG gene (430668-432621, present in
   essentially every genome -- metG is a housekeeping gene, so this is real
   core-genome variation, not a sparse accessory region) and displays only the
   ones carrying a clade-significant allele.
 - GCF_0030 = speA (arginine decarboxylase), NC_000915:436051-437899. Panel 2
   tests all 524 variable positions across the full gene and displays only the
   ones carrying a clade-significant allele -- showing every variable position
   regardless of significance would be too dense to read for either gene.
 - GCF_0017 (ketoacyl-ACP synthase III, 206682-207678) and GCF_0025 (hofH outer
   membrane protein, 1232404-1233822) are NOT in this figure: neither locus has
   ANY rows in snippy-core's core521.tab, meaning the region isn't present/
   aligned across enough of the 521 genomes to be part of the core alignment at
   all -- these two are genuinely accessory loci, not just under-called by
   deepBGC, so there is no genome-wide nucleotide alignment to show for them.

Row order = tree leaf order (Bio.Phylo get_terminals(), matching the Newick's own
child order; the "Reference" tip is pruned out of the tree/matrix entirely since
its bases are already shown in the header's REF row -- keeping it as a 522nd leaf
was redundant). Cells are only highlighted when the specific allele at that column
is statistically associated with clade membership: for every (position,
non-reference allele) carried by >=3 genomes, a Fisher's exact test (allele
present/absent x clade A/B) is run, p-values are BH-corrected across all such
tests, and cells carrying an allele significant at p_adj < 0.05 are colored on a
gradient toward the enriched clade's color -- pale near p_adj = 0.05, fully
saturated blue/orange by p_adj <= 1e-12 (log-scaled in between), so how dark a
cell is shows how strong the clade association is, not just whether it cleared
the threshold. Both panels are flagged `filter_to_significant_only=True`:
only columns with at least one significant allele are displayed at all -- this
replaces an earlier version that either grayed out any run of >=3 tree-adjacent
identical calls regardless of significance, or truncated speA to an arbitrary
first-200bp slice, or (one revision ago) still displayed all 144 of metG's
variable positions even though only some were significant.
"""
import csv
import math
import sys
from Bio import Phylo
from PIL import Image, ImageDraw, ImageFont
from scipy.stats import fisher_exact

TREE_PATH = "analysis/phylogenetics/results/gubbins/full_hybrid_521/full_hybrid_521.final_tree.tre"
WINDOWS_TSV = "analysis/phylogenetics/results/gcf_roles_followup/variable_sites_raw.tsv"
CLADE_TSV = "analysis/phylogenetics/results/root_split_comparison/genome_clade_assignments.tsv"
OUT = "analysis/phylogenetics/results/gcf_roles_followup/figs/fig_variable_sites.png"
FONT_PATH = "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"
FONT_BOLD_PATH = "/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf"

GENE_PANELS = [
    # (label, start, end, filter_to_significant_only, description)
    ("metG  (GCF_0019 / GCF_0031 locus), clade-significant positions only", 430668, 432621, True,
     "Methionine--tRNA ligase"),
    ("speA  (GCF_0030), clade-significant positions only", 436051, 437899, True,
     "Arginine decarboxylase"),
]

INK = (43, 43, 43)
SURFACE = (252, 252, 251)
CLADE_A_COLOR = (42, 120, 214)    # blue (tree branches)
CLADE_B_COLOR = (235, 104, 52)    # orange (tree branches)
CLADE_A_PALE = (235, 242, 252)    # cell background gradient, Clade A end: pale (p_adj just under ALPHA)
CLADE_A_DEEP = (90, 150, 225)     # cell background gradient, Clade A end: saturated (p_adj <= P_FLOOR)
CLADE_B_PALE = (253, 240, 231)    # cell background gradient, Clade B end: pale
CLADE_B_DEEP = (240, 150, 95)     # cell background gradient, Clade B end: saturated
GUIDE_COLOR = (200, 200, 200)     # light gray dotted tip-alignment guide
MIN_CARRIERS = 3                  # "repetitive" threshold for testing an allele
ALPHA = 0.05                      # BH-adjusted significance threshold
P_FLOOR = 1e-12                   # p_adj at/below this maps to the fully-saturated color

SCALE = 2               # overall resolution multiplier (crisper text when zoomed); bump
                        # this, not the individual pixel constants below, to resharpen
LINE_W = max(1, round(SCALE))  # tree branches, guide dashes, legend strokes

CELL_W = 7 * SCALE
CELL_H = 9 * SCALE
TREE_W_PX = 650 * SCALE
PANEL_GAP_PX = 28 * SCALE
LEFT_MARGIN = 10 * SCALE
RIGHT_MARGIN = 20 * SCALE
TOP_MARGIN = 10 * SCALE
LABEL_ROW_PX = 14 * SCALE
DESC_ROW_PX = 12 * SCALE   # short protein-name line, between the gene label and the position digits
LEGEND_H = 68 * SCALE   # reserved strip at the very top of the image for the color-scale legend
LEGEND_BAR_W = 260 * SCALE
LEGEND_BAR_H = 14 * SCALE
LEGEND_GAP = 40 * SCALE  # horizontal gap between the two clades' legend bars


def load_tree():
    tree = Phylo.read(TREE_PATH, "newick")
    ref = next(t for t in tree.get_terminals() if t.name == "Reference")
    tree.prune(ref)
    terminals = tree.get_terminals()
    return tree, terminals


def compute_xy(tree, terminals):
    y = {t: i for i, t in enumerate(terminals)}
    x = {}

    def assign_x(clade, parent_x):
        bl = clade.branch_length or 0.0
        cx = parent_x + bl
        x[clade] = cx
        for c in clade.clades:
            assign_x(c, cx)

    assign_x(tree.root, 0.0)

    def assign_y_internal(clade):
        if clade.is_terminal():
            return y[clade]
        ys = [assign_y_internal(c) for c in clade.clades]
        y[clade] = sum(ys) / len(ys)
        return y[clade]

    assign_y_internal(tree.root)
    return x, y


def load_clades():
    label_of = {}
    with open(CLADE_TSV) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            label_of[row["id"]] = row["clade"]
    return label_of


def build_subtree_labels(tree, label_of):
    """Post-order: clade -> frozenset of clade labels ('A'/'B'/None) among its terminals."""
    memo = {}

    def visit(clade):
        if clade in memo:
            return memo[clade]
        if clade.is_terminal():
            s = frozenset([label_of.get(clade.name)])
        else:
            s = frozenset().union(*(visit(c) for c in clade.clades))
        memo[clade] = s
        return s

    visit(tree.root)
    return memo


def color_for(clade, subtree_labels):
    s = subtree_labels[clade]
    if s == {"A"}:
        return CLADE_A_COLOR
    if s == {"B"}:
        return CLADE_B_COLOR
    return INK


def load_matrix():
    with open(WINDOWS_TSV) as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        rows = list(r)
    sample_names = header[3:]
    by_pos = {}
    for row in rows:
        pos = int(row[1])
        ref = row[2]
        calls = row[3:]
        by_pos[pos] = (ref, dict(zip(sample_names, calls)))
    return by_pos


def bh_adjust(pvals):
    """Benjamini-Hochberg adjusted p-values, same order as input."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj_sorted = [0.0] * m
    prev = 1.0
    for k in range(m - 1, -1, -1):
        rank = k + 1
        val = pvals[order[k]] * m / rank
        prev = min(prev, val)
        adj_sorted[k] = prev
    padj = [0.0] * m
    for rank_idx, orig_idx in enumerate(order):
        padj[orig_idx] = adj_sorted[rank_idx]
    return padj


def lerp_color(c0, c1, t):
    return tuple(round(c0[i] + (c1[i] - c0[i]) * t) for i in range(3))


def padj_to_t(padj):
    """0.0 at padj==ALPHA (just significant) to 1.0 at padj<=P_FLOOR (highly significant),
    log-scaled, clamped."""
    padj = max(padj, P_FLOOR)
    lo, hi = math.log10(ALPHA), math.log10(P_FLOOR)
    return min(max((math.log10(padj) - lo) / (hi - lo), 0.0), 1.0)


def gradient_color(direction, padj):
    """Map a BH-adjusted p-value to a color on the pale->saturated gradient for the
    given clade direction. padj == ALPHA -> pale end; padj <= P_FLOOR -> deep end;
    log-scaled in between."""
    t = padj_to_t(padj)
    pale, deep = (CLADE_A_PALE, CLADE_A_DEEP) if direction == "A" else (CLADE_B_PALE, CLADE_B_DEEP)
    return lerp_color(pale, deep, t)


def draw_gradient_legend(draw, x0, y0, direction, title, font_title, font_tick):
    """Draws one clade's pale->saturated gradient bar with tick labels at fixed
    p_adj reference points, starting at (x0, y0)."""
    pale, deep = (CLADE_A_PALE, CLADE_A_DEEP) if direction == "A" else (CLADE_B_PALE, CLADE_B_DEEP)
    draw.text((x0, y0), title, fill=INK, font=font_title, anchor="lt")
    bar_y = y0 + 16 * SCALE
    for px in range(LEGEND_BAR_W):
        t = px / (LEGEND_BAR_W - 1)
        draw.line([(x0 + px, bar_y), (x0 + px, bar_y + LEGEND_BAR_H)], fill=lerp_color(pale, deep, t))
    draw.rectangle([x0, bar_y, x0 + LEGEND_BAR_W, bar_y + LEGEND_BAR_H], outline=INK, width=LINE_W)
    tick_vals = [ALPHA, 1e-3, 1e-6, 1e-9, P_FLOOR]
    tick_labels = ["0.05", "1e-3", "1e-6", "1e-9", "<=1e-12"]
    for val, lab in zip(tick_vals, tick_labels):
        t = padj_to_t(val)
        tx_ = x0 + t * LEGEND_BAR_W
        draw.line([(tx_, bar_y + LEGEND_BAR_H), (tx_, bar_y + LEGEND_BAR_H + 3 * SCALE)], fill=INK, width=LINE_W)
        anchor = "la" if val == ALPHA else ("ra" if val == P_FLOOR else "ma")
        draw.text((tx_, bar_y + LEGEND_BAR_H + 4 * SCALE), lab, fill=INK, font=font_tick, anchor=anchor)


def compute_clade_significant_cells(panels, by_pos, label_of):
    """For every (position, non-ref allele) carried by >=MIN_CARRIERS genomes, test
    allele-presence x clade with Fisher's exact, BH-correct across all tests, and
    return {(panel_label, pos, genome_name): ('A'|'B', p_adj)} for cells whose allele
    is significantly enriched in that clade at p_adj < ALPHA."""
    clade_A = {n for n, c in label_of.items() if c == "A"}
    clade_B = {n for n, c in label_of.items() if c == "B"}
    nA, nB = len(clade_A), len(clade_B)

    tests = []  # [label, pos, allele, pval, direction, carrier_names]
    for label, positions in panels:
        for p in positions:
            ref, calls = by_pos[p]
            allele_carriers = {}
            for name, base in calls.items():
                if base == ref or base in ("N", "-", "", "X"):  # 'X' = ambiguous aa (shared w/ make_variable_aa_fig.py)
                    continue
                allele_carriers.setdefault(base, []).append(name)
            for allele, names in allele_carriers.items():
                if len(names) < MIN_CARRIERS:
                    continue
                a = sum(1 for n in names if n in clade_A)
                b = sum(1 for n in names if n in clade_B)
                if a + b == 0:
                    continue
                c, d = nA - a, nB - b
                _, pval = fisher_exact([[a, b], [c, d]])
                direction = "A" if (a / nA if nA else 0) > (b / nB if nB else 0) else "B"
                tests.append([label, p, allele, pval, direction, names])

    if not tests:
        return {}, 0, 0

    padj = bh_adjust([t[3] for t in tests])
    sig_cells = {}
    n_sig_tests = 0
    for (label, p, allele, pval, direction, names), pa in zip(tests, padj):
        if pa < ALPHA:
            n_sig_tests += 1
            for n in names:
                sig_cells[(label, p, n)] = (direction, pa)
    print(f"clade-association tests: {len(tests)} total, {n_sig_tests} significant "
          f"at BH p_adj<{ALPHA}", file=sys.stderr)
    return sig_cells, len(tests), n_sig_tests


def dotted_hline(draw, y, x0, x1, color, dash=None, gap=None, width=None):
    dash = dash if dash is not None else 3 * SCALE
    gap = gap if gap is not None else 3 * SCALE
    width = width if width is not None else LINE_W
    x = x0
    while x < x1:
        xe = min(x + dash, x1)
        draw.line([(x, y), (xe, y)], fill=color, width=width)
        x += dash + gap


def main():
    tree, terminals = load_tree()
    x, y = compute_xy(tree, terminals)
    max_x = max(x.values())
    label_of = load_clades()
    subtree_labels = build_subtree_labels(tree, label_of)

    by_pos = load_matrix()

    desc_by_label = {label: desc for label, s, e, filt, desc in GENE_PANELS}

    # candidate positions per panel, used for the significance tests
    test_panels = []
    for label, s, e, _, _desc in GENE_PANELS:
        positions = sorted(p for p in by_pos if s <= p <= e)
        test_panels.append((label, positions))

    sig_cells, n_tests, n_sig = compute_clade_significant_cells(test_panels, by_pos, label_of)
    sig_positions_by_label = {}
    for label, p, _name in sig_cells:
        sig_positions_by_label.setdefault(label, set()).add(p)

    # display positions per panel: filtered down to clade-significant columns
    # only for panels flagged filter_to_significant_only
    panels = []
    for (label, s, e, filt, _desc), (_, test_positions) in zip(GENE_PANELS, test_panels):
        if filt:
            positions = sorted(p for p in test_positions if p in sig_positions_by_label.get(label, set()))
        else:
            positions = test_positions
        panels.append((label, positions))
    total_cols = sum(len(p[1]) for p in panels)
    print("panel sizes (displayed):", [(l, len(p)) for l, p in panels], "total:", total_cols, file=sys.stderr)

    max_digits = max(len(str(p)) for _, positions in panels for p in positions)
    n_rows = len(terminals)

    header_px = LABEL_ROW_PX + DESC_ROW_PX + max_digits * CELL_H + CELL_H  # label + desc + digits + REF row
    matrix_w_px = total_cols * CELL_W + PANEL_GAP_PX * (len(panels) - 1)
    img_w = LEFT_MARGIN + TREE_W_PX + matrix_w_px + RIGHT_MARGIN
    img_h = LEGEND_H + TOP_MARGIN + header_px + n_rows * CELL_H + 10 * SCALE
    print("image size:", img_w, "x", img_h, file=sys.stderr)

    img = Image.new("RGB", (img_w, img_h), SURFACE)
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype(FONT_PATH, 7 * SCALE)
    font_bold = ImageFont.truetype(FONT_BOLD_PATH, 7 * SCALE)
    font_label = ImageFont.truetype(FONT_BOLD_PATH, 10 * SCALE)
    font_desc = ImageFont.truetype(FONT_PATH, 8 * SCALE)
    font_legend_title = ImageFont.truetype(FONT_BOLD_PATH, 9 * SCALE)
    font_legend_tick = ImageFont.truetype(FONT_PATH, 7 * SCALE)

    draw_gradient_legend(draw, LEFT_MARGIN, 6 * SCALE, "A", "Clade A association (BH p_adj)",
                          font_legend_title, font_legend_tick)
    draw_gradient_legend(draw, LEFT_MARGIN + LEGEND_BAR_W + LEGEND_GAP, 6 * SCALE, "B",
                          "Clade B association (BH p_adj)", font_legend_title, font_legend_tick)

    mat_x0 = LEFT_MARGIN + TREE_W_PX
    mat_y0 = LEGEND_H + TOP_MARGIN + header_px
    hdr_top = LEGEND_H + TOP_MARGIN

    def row_y(ri):
        return mat_y0 + ri * CELL_H

    # --- tree ---
    def tx(v):
        return LEFT_MARGIN + (v / max_x) * (TREE_W_PX - 20 * SCALE)

    def ty(v):
        return mat_y0 + v * CELL_H + CELL_H / 2

    for clade in x:
        if clade.is_terminal():
            continue
        cx = tx(x[clade])
        child_ys = [ty(y[c]) for c in clade.clades]
        draw.line([(cx, min(child_ys)), (cx, max(child_ys))], fill=color_for(clade, subtree_labels), width=LINE_W)
        for c in clade.clades:
            draw.line([(cx, ty(y[c])), (tx(x[c]), ty(y[c]))], fill=color_for(c, subtree_labels), width=LINE_W)

    # extend each tip's branch with a light dotted guide line out to the matrix's
    # left edge, so every row lines up with its alignment columns regardless of
    # how deep that tip's actual branch length placed it
    for t in terminals:
        tip_x = tx(x[t])
        ty_v = ty(y[t])
        if tip_x < mat_x0:
            dotted_hline(draw, ty_v, tip_x, mat_x0, GUIDE_COLOR)

    # --- column index + panel bounds ---
    col_cursor = 0
    col_index = {}
    panel_bounds = []
    for label, positions in panels:
        start_px = mat_x0 + col_cursor * CELL_W + panels.index((label, positions)) * PANEL_GAP_PX
        for p in positions:
            col_index[(label, p)] = col_cursor
            col_cursor += 1
        panel_bounds.append((label, positions, start_px))

    def col_px(panel_i, local_col):
        return mat_x0 + sum(len(panels[k][1]) for k in range(panel_i)) * CELL_W + panel_i * PANEL_GAP_PX + local_col * CELL_W

    # --- header: label, short protein description, digits, ref row ---
    for pi, (label, positions) in enumerate(panels):
        x0 = col_px(pi, 0)
        x1 = col_px(pi, len(positions) - 1) + CELL_W
        draw.text(((x0 + x1) / 2, hdr_top), label, fill=INK, font=font_label, anchor="mt")
        draw.text(((x0 + x1) / 2, hdr_top + LABEL_ROW_PX), desc_by_label.get(label, ""),
                   fill=INK, font=font_desc, anchor="mt")
        for li, p in enumerate(positions):
            cx = col_px(pi, li)
            digits = str(p)
            for d_i, d in enumerate(digits):
                dy = hdr_top + LABEL_ROW_PX + DESC_ROW_PX + d_i * CELL_H
                draw.text((cx + CELL_W / 2, dy + CELL_H / 2), d, fill=INK, font=font, anchor="mm")
            ref, calls = by_pos[p]
            ref_y = hdr_top + LABEL_ROW_PX + DESC_ROW_PX + max_digits * CELL_H
            draw.text((cx + CELL_W / 2, ref_y + CELL_H / 2), ref, fill=INK, font=font_bold, anchor="mm")

    # --- highlight cells whose allele is significantly clade-associated ---
    # (sig_cells was already computed above, against the full test_panels ranges)
    for pi, (label, positions) in enumerate(panels):
        for li, p in enumerate(positions):
            cx = col_px(pi, li)
            for ri, term in enumerate(terminals):
                name = term.name
                hit = sig_cells.get((label, p, name))
                if hit is None:
                    continue
                direction, padj = hit
                tint = gradient_color(direction, padj)
                y0 = row_y(ri)
                draw.rectangle([cx, y0, cx + CELL_W, y0 + CELL_H], fill=tint)

    # --- per-cell text (one row at a time via multiline-safe single draws) ---
    for ri, term in enumerate(terminals):
        name = term.name
        ry = row_y(ri)
        for pi, (label, positions) in enumerate(panels):
            for li, p in enumerate(positions):
                ref, calls = by_pos[p]
                cx = col_px(pi, li)
                base = calls.get(name, "N")
                show = "." if base == ref else base
                draw.text((cx + CELL_W / 2, ry + CELL_H / 2), show, fill=INK, font=font, anchor="mm")

    img.save(OUT)
    print("wrote", OUT, file=sys.stderr)


if __name__ == "__main__":
    main()
