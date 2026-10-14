"""
Table of every clade-significant amino-acid substitution found by
make_variable_aa_fig.py (metG/speA), in formal notation (e.g. E256A), with
physicochemical properties for both the reference and substituted residue,
the clade driving the association, and the BH-adjusted significance.

Reuses make_variable_aa_fig.py's/make_variable_sites_fig.py's data-loading and
significance-testing code directly rather than re-implementing it, so this
table is always consistent with whatever the figure itself shows.
"""
import csv
import math
import sys

sys.path.insert(0, "analysis/phylogenetics/scripts")
import make_variable_aa_fig as aa
import make_variable_sites_fig as base

OUT_TSV = "analysis/phylogenetics/results/gcf_roles_followup/aa_substitutions_table.tsv"

# Standard amino-acid properties. Polarity/charge per the common Lehninger-style
# classification (nonpolar-hydrophobic / polar-uncharged / acidic / basic);
# hydrophobicity = Kyte-Doolittle hydropathy index (higher = more hydrophobic).
AA_INFO = {
    "A": ("Ala", "Alanine",       "Nonpolar",        "Neutral",   1.8),
    "R": ("Arg", "Arginine",      "Polar",           "Positive", -4.5),
    "N": ("Asn", "Asparagine",    "Polar",           "Neutral",  -3.5),
    "D": ("Asp", "Aspartate",     "Polar (acidic)",  "Negative", -3.5),
    "C": ("Cys", "Cysteine",      "Polar",           "Neutral",   2.5),
    "Q": ("Gln", "Glutamine",     "Polar",           "Neutral",  -3.5),
    "E": ("Glu", "Glutamate",     "Polar (acidic)",  "Negative", -3.5),
    "G": ("Gly", "Glycine",       "Nonpolar",        "Neutral",  -0.4),
    "H": ("His", "Histidine",     "Polar (basic)",   "Positive", -3.2),
    "I": ("Ile", "Isoleucine",    "Nonpolar",        "Neutral",   4.5),
    "L": ("Leu", "Leucine",       "Nonpolar",        "Neutral",   3.8),
    "K": ("Lys", "Lysine",        "Polar (basic)",   "Positive", -3.9),
    "M": ("Met", "Methionine",    "Nonpolar",        "Neutral",   1.9),
    "F": ("Phe", "Phenylalanine", "Nonpolar",        "Neutral",   2.8),
    "P": ("Pro", "Proline",       "Nonpolar",        "Neutral",  -1.6),
    "S": ("Ser", "Serine",        "Polar",           "Neutral",  -0.8),
    "T": ("Thr", "Threonine",     "Polar",           "Neutral",  -0.7),
    "W": ("Trp", "Tryptophan",    "Nonpolar",        "Neutral",  -0.9),
    "Y": ("Tyr", "Tyrosine",      "Polar",           "Neutral",  -1.3),
    "V": ("Val", "Valine",        "Nonpolar",        "Neutral",   4.2),
}


def main():
    label_of = base.load_clades()
    ref_genome_seq = aa.load_reference_seq()
    sample_names, variant_rows = aa.load_variant_rows()
    n_clade_A = sum(1 for c in label_of.values() if c == "A")
    n_clade_B = sum(1 for c in label_of.values() if c == "B")

    rows_out = []
    for label, s, e, _filt, _desc in aa.AA_GENE_PANELS:
        gene = label.split()[0]
        rows_in_window = [r for r in variant_rows if s <= int(r[1]) <= e]
        by_pos = aa.build_aa_matrix(label, s, e, ref_genome_seq, sample_names, rows_in_window)
        positions = sorted(by_pos.keys())
        sig_cells, n_tests, n_sig = base.compute_clade_significant_cells([(label, positions)], by_pos, label_of)

        # group sig_cells by (position, direction, padj) to recover one row per
        # distinct significant allele (usually 1, occasionally 2, per position)
        groups = {}
        for (lbl, pos, name), (direction, padj) in sig_cells.items():
            groups.setdefault((pos, direction, round(padj, 300)), []).append(name)

        for (pos, direction, padj), names in groups.items():
            ref_a, col = by_pos[pos]
            alt_a = col[names[0]]
            ref3, ref_full, ref_pol, ref_chg, ref_kd = AA_INFO[ref_a]
            alt3, alt_full, alt_pol, alt_chg, alt_kd = AA_INFO[alt_a]
            n_a = sum(1 for n in names if label_of.get(n) == "A")
            n_b = sum(1 for n in names if label_of.get(n) == "B")
            # Haldane-Anscombe-corrected (+0.5/cell) odds ratio, carrier vs non-carrier,
            # Clade A : Clade B -- so OR > 1 means Clade-A-enriched, OR < 1 Clade-B-enriched;
            # the correction keeps this finite for the rows with 0 carriers on one side.
            a, b = n_a + 0.5, n_b + 0.5
            c, d = (n_clade_A - n_a) + 0.5, (n_clade_B - n_b) + 0.5
            odds_ratio = (a / c) / (b / d)
            log2_or = math.log2(odds_ratio)
            rows_out.append({
                "gene": gene,
                "substitution": f"{ref_a}{pos}{alt_a}",
                "substitution_3letter": f"{ref3}{pos}{alt3}",
                "ref_aa_1letter": ref_a, "ref_aa_3letter": ref3, "ref_aa_name": ref_full,
                "ref_polarity": ref_pol, "ref_charge": ref_chg, "ref_hydrophobicity_KD": ref_kd,
                "alt_aa_1letter": alt_a, "alt_aa_3letter": alt3, "alt_aa_name": alt_full,
                "alt_polarity": alt_pol, "alt_charge": alt_chg, "alt_hydrophobicity_KD": alt_kd,
                "delta_hydrophobicity_KD": round(alt_kd - ref_kd, 1),
                "polarity_change": "yes" if ref_pol.split()[0] != alt_pol.split()[0] else "no",
                "charge_change": "yes" if ref_chg != alt_chg else "no",
                "clade_driving_significance": direction,
                "n_carriers": len(names),
                "n_carriers_cladeA": n_a,
                "n_carriers_cladeB": n_b,
                "odds_ratio_A_vs_B": round(odds_ratio, 4),
                "log2_odds_ratio_A_vs_B": round(log2_or, 3),
                "neg_log10_p_adj": round(-math.log10(max(padj, 1e-300)), 3),
                "p_adj": padj,
            })

    rows_out.sort(key=lambda r: r["p_adj"])
    for i, r in enumerate(rows_out, start=1):
        r["id"] = i
    # move id to the front
    rows_out = [{"id": r["id"], **{k: v for k, v in r.items() if k != "id"}} for r in rows_out]

    fieldnames = list(rows_out[0].keys())
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    print(f"wrote {OUT_TSV}: {len(rows_out)} substitutions", file=sys.stderr)
    return rows_out


if __name__ == "__main__":
    main()
