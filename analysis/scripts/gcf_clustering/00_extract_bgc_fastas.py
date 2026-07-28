#!/usr/bin/env python3
"""
00_extract_bgc_fastas.py

Extracts per-BGC protein FASTA files from antiSMASH region GenBank files, for use
as --fasta-dir input to 03_backbone_identity.py.

antiSMASH output layout (confirmed on this project's server):
    {antismash_dir}/{sample_id}/{contig_id}.region001.gbk

Each region GBK contains CDS features with /translation qualifiers. This script
writes one FASTA file per BGC, named {sample_id}__{contig_id}.faa, matching the
bgc_id convention used throughout scripts/gcf_clustering (see 02_jaccard_matrix.py:
build_bgc_id).

Only BGCs present in bgc_domain_arrays.tsv (i.e. that survived the 01 filters) are
processed, and only BGCs predicted by antiSMASH have a region GBK to extract from
(deepBGC/GECCO BGCs have no antiSMASH region file and are skipped with a warning —
their sequence_identity will fall back to 0.0 in 03/04, which is already the
documented behavior for BGCs without backbone domains).

Usage:
    python 00_extract_bgc_fastas.py \
        --domains      ../results/bgc_domain_arrays.tsv \
        --antismash-dir ../../funcscan/results/bgc/antismash \
        --output-dir   ../../rawdata/fasta/proteins
"""

import argparse
import pathlib
import sys

from Bio import SeqIO  # requires: biopython (see requirements.txt)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--domains", required=True, help="Path to bgc_domain_arrays.tsv")
    parser.add_argument(
        "--antismash-dir", required=True,
        help="Path to funcscan antiSMASH results dir "
             "(contains {sample_id}/{contig_id}.region*.gbk)",
    )
    parser.add_argument("--output-dir", required=True, help="Directory to write {bgc_id}.faa files")
    parser.add_argument(
        "--only-tool", default="antiSMASH",
        help="Only extract BGCs whose Prediction_tool matches this value "
             "(default: antiSMASH; the only tool with region GBKs)",
    )
    return parser.parse_args()


def find_region_gbks(sample_dir: pathlib.Path, contig_id: str) -> list[pathlib.Path]:
    """
    A contig can have multiple antiSMASH regions (region001, region002, ...).
    Return all matching region GBKs for this contig, sorted.
    """
    return sorted(sample_dir.glob(f"{contig_id}.region*.gbk"))


def extract_translations(gbk_path: pathlib.Path) -> dict[str, str]:
    """Parse a region GBK and return {protein_id: translation} for all CDS features."""
    seqs = {}
    for record in SeqIO.parse(str(gbk_path), "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            translation = feature.qualifiers.get("translation", [None])[0]
            if not translation:
                continue
            protein_id = (
                feature.qualifiers.get("locus_tag", [None])[0]
                or feature.qualifiers.get("protein_id", [None])[0]
            )
            if not protein_id:
                continue
            seqs[protein_id] = translation
    return seqs


def main():
    args = parse_args()

    import pandas as pd

    print(f"[00] Reading {args.domains} ...", flush=True)
    df = pd.read_csv(args.domains, sep="\t", dtype=str)

    if args.only_tool:
        n_before = len(df)
        df = df[df["Prediction_tool"] == args.only_tool].copy()
        print(f"[00] Filtered to Prediction_tool == '{args.only_tool}': "
              f"{len(df)} / {n_before} BGCs.", flush=True)

    antismash_dir = pathlib.Path(args.antismash_dir)
    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_written = 0
    n_missing_dir = 0
    n_missing_gbk = 0
    n_empty = 0

    for row in df.itertuples(index=False):
        sample_id = row.sample_id
        contig_id = row.contig_id
        bgc_id = f"{sample_id}__{contig_id}"

        sample_dir = antismash_dir / sample_id
        if not sample_dir.is_dir():
            print(f"[00] WARNING: no antiSMASH dir for sample {sample_id}", flush=True)
            n_missing_dir += 1
            continue

        region_gbks = find_region_gbks(sample_dir, contig_id)
        if not region_gbks:
            print(f"[00] WARNING: no region GBK for {bgc_id} in {sample_dir}", flush=True)
            n_missing_gbk += 1
            continue

        seqs = {}
        for gbk in region_gbks:
            seqs.update(extract_translations(gbk))

        if not seqs:
            print(f"[00] WARNING: no CDS translations found for {bgc_id}", flush=True)
            n_empty += 1
            continue

        out_path = output_dir / f"{bgc_id}.faa"
        with open(out_path, "w") as fh:
            for protein_id, seq in seqs.items():
                fh.write(f">{protein_id}\n{seq}\n")
        n_written += 1

    print(f"\n[00] Done. {n_written} FASTA files written to {output_dir}", flush=True)
    if n_missing_dir or n_missing_gbk or n_empty:
        print(f"[00] Skipped: {n_missing_dir} missing sample dir, "
              f"{n_missing_gbk} missing region GBK, {n_empty} empty translations.",
              flush=True)
    if n_written == 0:
        print("[00] ERROR: no FASTA files were written. Check --antismash-dir path.",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()