"""
Generates the sortable/searchable HTML <table> partial embedded in
visual_report.qmd, from the same data make_aa_substitution_table.py produces
(imported directly, not re-derived) -- so the table in the Quarto report is
always consistent with aa_substitutions_table.tsv.

Output is a bare <table id="subTable"> (no <html>/<head>), meant to be pulled
into the .qmd via a Quarto `{{< include >}}` shortcode; DataTables (loaded via
CDN in the .qmd itself) turns it into a sortable/searchable widget at render
time in the browser -- no R/Python runtime needed to view the report.
"""
import sys

sys.path.insert(0, "analysis/phylogenetics/scripts")
import make_aa_substitution_table as tbl

OUT = "analysis/phylogenetics/results/gcf_roles_followup/visual_report/table_partial.html"

COLUMNS = [
    ("id", "#"),
    ("gene", "Gene"),
    ("substitution", "Substitution"),
    ("ref_aa_3letter", "Ref"),
    ("ref_polarity", "Ref polarity"),
    ("ref_charge", "Ref charge"),
    ("alt_aa_3letter", "Alt"),
    ("alt_polarity", "Alt polarity"),
    ("alt_charge", "Alt charge"),
    ("clade_driving_significance", "Clade"),
    ("n_carriers_cladeA", "Carriers A"),
    ("n_carriers_cladeB", "Carriers B"),
    ("odds_ratio_A_vs_B", "OR (A:B)"),
    ("p_adj", "p_adj"),
]


def fmt(key, val):
    if key == "p_adj":
        return f"{val:.2e}"
    if key == "odds_ratio_A_vs_B":
        return f"{val:.2f}"
    if key == "clade_driving_significance":
        return f"Clade {val}"
    return str(val)


def main():
    rows = tbl.main()

    lines = ['<table id="subTable" class="display" style="width:100%">', "  <thead><tr>"]
    for _, label in COLUMNS:
        lines.append(f"    <th>{label}</th>")
    lines.append("  </tr></thead>")
    lines.append("  <tbody>")
    for r in rows:
        lines.append(f'    <tr data-id="{r["id"]}">')
        for key, _ in COLUMNS:
            lines.append(f"      <td>{fmt(key, r[key])}</td>")
        lines.append("    </tr>")
    lines.append("  </tbody>")
    lines.append("</table>")

    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {OUT}: {len(rows)} rows", file=sys.stderr)


if __name__ == "__main__":
    main()
