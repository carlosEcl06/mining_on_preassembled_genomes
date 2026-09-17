#!/usr/bin/env python3
"""
Converts a whole-genome antiSMASH .gbk (one per genome, already confirmed to
be the full assembly, not just a BGC region -- see plan's "What was
independently re-verified" section) into a Prokka-style GFF3 file with an
appended ##FASTA section, which is what Panaroo expects as input. Used for
the Panaroo cross-check on the Gubbins pilot subsample (task #8): Panaroo
needs per-genome gene annotations, not just raw contigs, so this is the
minimal conversion needed -- no re-annotation, just a format change of
antiSMASH's own CDS calls (which are themselves inherited from the genome's
original NCBI annotation, since antiSMASH annotates on top of existing
assemblies rather than calling genes itself for annotated inputs).
"""
import argparse
import sys
from pathlib import Path

from Bio import SeqIO


def gbk_to_gff3(gbk_path: Path, out_path: Path, sample_id: str):
    records = list(SeqIO.parse(gbk_path, "genbank"))
    lines = ["##gff-version 3"]
    cds_count = 0
    for rec in records:
        lines.append(f"##sequence-region {rec.id} 1 {len(rec.seq)}")
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            cds_count += 1
            start = int(feat.location.start) + 1  # GFF3 is 1-based
            end = int(feat.location.end)
            strand = "+" if feat.location.strand == 1 else "-"
            locus_tag = feat.qualifiers.get("locus_tag", [f"{sample_id}_{cds_count:05d}"])[0]
            product = feat.qualifiers.get("product", ["hypothetical protein"])[0].replace(";", ",").replace("=", "")
            attrs = f"ID={locus_tag};locus_tag={locus_tag};product={product}"
            lines.append(
                "\t".join([rec.id, "antiSMASH", "CDS", str(start), str(end), ".", strand, "0", attrs])
            )
    lines.append("##FASTA")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
        for rec in records:
            SeqIO.write(rec, f, "fasta")
    return cds_count, len(records)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gbk-dir", required=True, help="dir of <sample_id>.gbk files")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ids-file", required=True, help="one sample_id per line")
    args = ap.parse_args()

    gbk_dir = Path(args.gbk_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ids = [l.strip() for l in open(args.ids_file) if l.strip()]

    ok, failed = 0, []
    for sid in ids:
        gbk_path = gbk_dir / f"{sid}.gbk"
        if not gbk_path.exists():
            failed.append(sid)
            continue
        try:
            n_cds, n_contigs = gbk_to_gff3(gbk_path, out_dir / f"{sid}.gff", sid)
            print(f"{sid}: {n_contigs} contigs, {n_cds} CDS -> {sid}.gff")
            ok += 1
        except Exception as e:
            print(f"{sid}: FAILED ({e})", file=sys.stderr)
            failed.append(sid)

    print(f"\n{ok}/{len(ids)} converted successfully.")
    if failed:
        print(f"Failed: {failed}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
