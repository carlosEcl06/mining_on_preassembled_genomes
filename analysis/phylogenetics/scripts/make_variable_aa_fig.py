"""
Amino-acid analog of make_variable_sites_fig.py: same tree + clade-significance-
gradient alignment figure, but columns are translated amino-acid positions in
metG and speA rather than nucleotide positions. Reuses all the generic
tree/legend/significance-testing machinery from make_variable_sites_fig.py
(itself unmodified) -- only the data layer (nucleotide-SNP-table -> per-genome
sequence -> translation -> variable-codon table) is new.

How each genome's codon sequence is reconstructed: snippy-core's core521.tab
(via variable_sites_raw.tsv, already saved locally) lists every *variable*
core-genome position with the reference base and each genome's call; positions
it doesn't list are -- by construction of a core-genome alignment -- identical
to the reference in every one of the 521 genomes. So for each gene window,
take the reference CDS nucleotide sequence and, per genome, substitute in just
the variable positions that fall inside it. Both genes are on the + strand
with no indels in the core alignment, so this is a straight same-length
substitution, then Bio.Seq.translate(table=11).

**Reference-sequence gotcha (cost an hour to find):** the nucleotide figure's
gene coordinates were originally found by k-mer-matching a representative
GCF member against `26695.gbk` -- which turns out to be accession CP079087,
*not* NC_000915 (same strain, different assembly/submission). That's
close enough for picking out "which ~2kb window is this gene" (which is all
the nucleotide figure needed), but translation needs an exact, frame-correct
sequence, and CP079087 and NC_000915 are offset from each other here -- using
CP079087's sequence for the per-genome reconstruction produced frameshifted
garbage (near-universal fake stop codons partway through an essential gene,
the unmistakable signature of an off-by-some-indels frame error). This is
the *exact* CP079087-vs-NC_000915 mixup Report #5 (Section 3) already flagged
once for the GrafGen reference step -- same two accessions, same project,
different pipeline stage. Fixed by using `NC_000915.fasta` (the actual
sequence `core521.tab` is coordinate-based on) and re-locating each gene's
exact boundaries in *that* sequence by the same k-mer seed-and-vote method,
rather than reusing the CP079087-derived coordinates. Confirmed clean (no
internal stop codons, correct annotated protein length) before use below.
"""
import csv
import sys

sys.path.insert(0, "analysis/phylogenetics/scripts")
from Bio import SeqIO
from Bio.Seq import Seq

import make_variable_sites_fig as base

REF_FASTA = "analysis/phylogenetics/rawdata/reference/NC_000915.fasta"
WINDOWS_TSV = "analysis/phylogenetics/results/gcf_roles_followup/variable_sites_raw.tsv"
OUT = "analysis/phylogenetics/results/gcf_roles_followup/figs/fig_variable_aa.png"

# (label, cds_start, cds_end, filter_to_significant_only, description) -- 0-based
# python-string-index CDS bounds *in NC_000915.fasta*, i.e. nc_seq[cds_start:cds_end]
# is the gene. Located by k-mer seed-matching the CP079087 CDS (whose boundaries
# the nucleotide figure's GENE_PANELS uses) against NC_000915.fasta: both genes
# landed on one unambiguous diagonal, offset +230bp from the CP079087 coordinates
# (430668->430898 for metG, 436051->436281 for speA) -- a single local indel
# between the two assemblies upstream of both genes, not two separate shifts.
AA_GENE_PANELS = [
    ("metG  (GCF_0019 / GCF_0031 locus), clade-significant positions only", 430898, 432851, True,
     "Methionine--tRNA ligase"),
    ("speA  (GCF_0030), clade-significant positions only", 436281, 438129, True,
     "Arginine decarboxylase"),
]


def load_reference_seq():
    rec = next(SeqIO.parse(REF_FASTA, "fasta"))
    return str(rec.seq)


def load_variant_rows():
    with open(WINDOWS_TSV) as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        rows = list(r)
    sample_names = header[3:]
    return sample_names, rows


def build_aa_matrix(label, cds_start, cds_end, ref_genome_seq, sample_names, variant_rows):
    """Returns a by_pos dict shaped like make_variable_sites_fig.load_matrix()'s
    output, but keyed by 1-based amino-acid position: {aa_pos: (ref_aa, {name: aa})},
    restricted to amino-acid positions that are actually variable."""
    ref_nt = list(ref_genome_seq[cds_start:cds_end])
    n_nt = len(ref_nt)
    assert n_nt % 3 == 0, f"{label}: CDS length {n_nt} is not a multiple of 3"

    # per-genome nucleotide sequence, starting as a copy of the reference, with this
    # window's variable positions substituted in per genome
    seqs = {name: ref_nt.copy() for name in sample_names}
    for row in variant_rows:
        pos = int(row[1])
        i = pos - cds_start - 1  # POS is 1-based; cds_start is Biopython 0-based
        if not (0 <= i < n_nt):
            continue  # outside this gene (e.g. the neighboring-gene edge case)
        calls = row[3:]
        for name, base in zip(sample_names, calls):
            # non-ACGT calls (N, -, multi-base ambiguity codes from multi-allelic
            # sites) are written through as 'N' so the codon translates to 'X' and
            # gets skipped downstream, rather than silently assuming it matches the
            # reference.
            seqs[name][i] = base if base in ("A", "C", "G", "T") else "N"

    ref_aa = str(Seq("".join(ref_nt)).translate(table=11))
    aa_calls = {name: str(Seq("".join(seqs[name])).translate(table=11)) for name in sample_names}

    by_pos = {}
    n_aa = len(ref_aa)
    for codon_i in range(n_aa):
        ref_a = ref_aa[codon_i]
        col = {name: aa_calls[name][codon_i] for name in sample_names}
        if any(a != ref_a for a in col.values()):
            by_pos[codon_i + 1] = (ref_a, col)  # 1-based amino-acid position
    print(f"{label}: {n_aa} codons, {len(by_pos)} variable amino-acid positions", file=sys.stderr)
    return by_pos


def main():
    tree, terminals = base.load_tree()
    x, y = base.compute_xy(tree, terminals)
    max_x = max(x.values())
    label_of = base.load_clades()
    metadata = base.load_metadata()
    subtree_labels = base.build_subtree_labels(tree, label_of)

    ref_genome_seq = load_reference_seq()
    sample_names, variant_rows = load_variant_rows()

    desc_by_label = {label: desc for label, s, e, filt, desc in AA_GENE_PANELS}

    # compute_clade_significant_cells() takes one flat by_pos dict per call, keyed
    # by position only -- since metG and speA codon numbering both start at 1, each
    # gene is tested separately against its own by_pos and results merged (the
    # label is preserved in the returned keys, so no collision).
    by_pos_per_label = {}
    test_panels = []
    sig_cells = {}
    for label, s, e, _filt, _desc in AA_GENE_PANELS:
        rows_in_window = [row for row in variant_rows if s <= int(row[1]) <= e]
        by_pos = build_aa_matrix(label, s, e, ref_genome_seq, sample_names, rows_in_window)
        by_pos_per_label[label] = by_pos
        positions = sorted(by_pos.keys())
        test_panels.append((label, positions))
        cells, n_tests, n_sig = base.compute_clade_significant_cells([(label, positions)], by_pos, label_of)
        sig_cells.update(cells)

    sig_positions_by_label = {}
    for label, p, _name in sig_cells:
        sig_positions_by_label.setdefault(label, set()).add(p)

    panels = []
    for (label, s, e, filt, _desc), (_, test_positions) in zip(AA_GENE_PANELS, test_panels):
        if filt:
            positions = sorted(p for p in test_positions if p in sig_positions_by_label.get(label, set()))
        else:
            positions = test_positions
        panels.append((label, positions))
    total_cols = sum(len(p[1]) for p in panels)
    print("panel sizes (displayed):", [(l, len(p)) for l, p in panels], "total:", total_cols, file=sys.stderr)

    by_pos_lookup = {}
    for label, _ in test_panels:
        for p, v in by_pos_per_label[label].items():
            by_pos_lookup[(label, p)] = v

    max_digits = max(len(str(p)) for _, positions in panels for p in positions)
    n_rows = len(terminals)

    from PIL import Image, ImageDraw, ImageFont
    font = ImageFont.truetype(base.FONT_PATH, 7 * base.SCALE)
    font_bold = ImageFont.truetype(base.FONT_BOLD_PATH, 7 * base.SCALE)
    font_label = ImageFont.truetype(base.FONT_BOLD_PATH, 10 * base.SCALE)
    font_desc = ImageFont.truetype(base.FONT_PATH, 8 * base.SCALE)
    font_legend_title = ImageFont.truetype(base.FONT_BOLD_PATH, 9 * base.SCALE)
    font_legend_tick = ImageFont.truetype(base.FONT_PATH, 7 * base.SCALE)
    font_strip_header = ImageFont.truetype(base.FONT_BOLD_PATH, 6 * base.SCALE)

    # With few displayed columns (amino-acid panels filter hard), the gene-name/
    # description text can be wider than the data block itself -- reserve at least
    # enough width per panel for its own label/description, or labels overlap.
    alloc_w = []
    for label, positions in panels:
        data_w = len(positions) * base.CELL_W
        label_w = font_label.getlength(label)
        desc_w = font_desc.getlength(desc_by_label.get(label, ""))
        alloc_w.append(round(max(data_w, label_w, desc_w)) + 6 * base.SCALE)

    header_px = base.LABEL_ROW_PX + base.DESC_ROW_PX + max_digits * base.CELL_H + base.CELL_H
    matrix_w_px = sum(alloc_w) + base.PANEL_GAP_PX * (len(panels) - 1)
    img_w = (base.LEFT_MARGIN + base.TREE_W_PX + 2 * base.ANNOT_STRIP_W + base.ANNOT_GAP_AFTER
             + matrix_w_px + base.RIGHT_MARGIN)
    img_h = base.LEGEND_H + base.TOP_MARGIN + header_px + n_rows * base.CELL_H + 10 * base.SCALE
    print("image size:", img_w, "x", img_h, file=sys.stderr)

    img = Image.new("RGB", (img_w, img_h), base.SURFACE)
    draw = ImageDraw.Draw(img)

    base.draw_gradient_legend(draw, base.LEFT_MARGIN, 6 * base.SCALE, "A", "Clade A association (BH p_adj)",
                               font_legend_title, font_legend_tick)
    base.draw_gradient_legend(draw, base.LEFT_MARGIN + base.LEGEND_BAR_W + base.LEGEND_GAP, 6 * base.SCALE, "B",
                               "Clade B association (BH p_adj)", font_legend_title, font_legend_tick)
    legend3_x0 = base.LEFT_MARGIN + 2 * base.LEGEND_BAR_W + 2 * base.LEGEND_GAP
    legend3_end = base.draw_categorical_legend(draw, legend3_x0, 6 * base.SCALE, "Phenotype", base.PHENOTYPE_ORDER,
                                                base.PHENOTYPE_COLORS, font_legend_title, font_legend_tick)
    base.draw_categorical_legend(draw, legend3_end + base.LEGEND_GAP, 6 * base.SCALE, "Continent",
                                  base.CONTINENT_ORDER, base.CONTINENT_COLORS, font_legend_title, font_legend_tick)

    strip_x0 = base.LEFT_MARGIN + base.TREE_W_PX
    mat_x0 = strip_x0 + 2 * base.ANNOT_STRIP_W + base.ANNOT_GAP_AFTER
    mat_y0 = base.LEGEND_H + base.TOP_MARGIN + header_px
    hdr_top = base.LEGEND_H + base.TOP_MARGIN

    def row_y(ri):
        return mat_y0 + ri * base.CELL_H

    def tx(v):
        return base.LEFT_MARGIN + (v / max_x) * (base.TREE_W_PX - 20 * base.SCALE)

    def ty(v):
        return mat_y0 + v * base.CELL_H + base.CELL_H / 2

    for clade in x:
        if clade.is_terminal():
            continue
        cx = tx(x[clade])
        child_ys = [ty(y[c]) for c in clade.clades]
        draw.line([(cx, min(child_ys)), (cx, max(child_ys))], fill=base.color_for(clade, subtree_labels),
                   width=base.LINE_W)
        for c in clade.clades:
            draw.line([(cx, ty(y[c])), (tx(x[c]), ty(y[c]))], fill=base.color_for(c, subtree_labels),
                       width=base.LINE_W)

    for t in terminals:
        tip_x = tx(x[t])
        ty_v = ty(y[t])
        if tip_x < strip_x0:
            base.dotted_hline(draw, ty_v, tip_x, strip_x0, base.GUIDE_COLOR)

    base.draw_annotation_strips(draw, img, terminals, metadata, strip_x0, row_y, hdr_top, header_px,
                                 font_strip_header)

    def panel_x0(panel_i):
        return mat_x0 + sum(alloc_w[k] for k in range(panel_i)) + panel_i * base.PANEL_GAP_PX

    def col_px(panel_i, local_col):
        data_w = len(panels[panel_i][1]) * base.CELL_W
        data_offset = (alloc_w[panel_i] - data_w) / 2  # center the data block in its slot
        return panel_x0(panel_i) + data_offset + local_col * base.CELL_W

    for pi, (label, positions) in enumerate(panels):
        px0 = panel_x0(pi)
        px1 = px0 + alloc_w[pi]
        draw.text(((px0 + px1) / 2, hdr_top), label, fill=base.INK, font=font_label, anchor="mt")
        draw.text(((px0 + px1) / 2, hdr_top + base.LABEL_ROW_PX), desc_by_label.get(label, ""),
                   fill=base.INK, font=font_desc, anchor="mt")
        for li, p in enumerate(positions):
            cx = col_px(pi, li)
            digits = str(p)
            for d_i, d in enumerate(digits):
                dy = hdr_top + base.LABEL_ROW_PX + base.DESC_ROW_PX + d_i * base.CELL_H
                draw.text((cx + base.CELL_W / 2, dy + base.CELL_H / 2), d, fill=base.INK, font=font, anchor="mm")
            ref_a, _ = by_pos_lookup[(label, p)]
            ref_y = hdr_top + base.LABEL_ROW_PX + base.DESC_ROW_PX + max_digits * base.CELL_H
            draw.text((cx + base.CELL_W / 2, ref_y + base.CELL_H / 2), ref_a, fill=base.INK, font=font_bold,
                       anchor="mm")

    for pi, (label, positions) in enumerate(panels):
        for li, p in enumerate(positions):
            cx = col_px(pi, li)
            for ri, term in enumerate(terminals):
                name = term.name
                hit = sig_cells.get((label, p, name))
                if hit is None:
                    continue
                direction, padj = hit
                tint = base.gradient_color(direction, padj)
                y0 = row_y(ri)
                draw.rectangle([cx, y0, cx + base.CELL_W, y0 + base.CELL_H], fill=tint)

    for ri, term in enumerate(terminals):
        name = term.name
        ry = row_y(ri)
        for pi, (label, positions) in enumerate(panels):
            for li, p in enumerate(positions):
                ref_a, col = by_pos_lookup[(label, p)]
                cx = col_px(pi, li)
                aa = col.get(name, "X")
                show = "." if aa == ref_a else aa
                draw.text((cx + base.CELL_W / 2, ry + base.CELL_H / 2), show, fill=base.INK, font=font, anchor="mm")

    img.save(OUT)
    print("wrote", OUT, file=sys.stderr)


if __name__ == "__main__":
    main()
