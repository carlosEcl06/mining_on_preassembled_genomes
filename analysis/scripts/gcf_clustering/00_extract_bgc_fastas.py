#!/usr/bin/env python3
"""
00_extract_bgc_fastas.py

Extracts per-BGC protein FASTA files from antiSMASH region GenBank files, for use
as --fasta-dir input to 03_backbone_identity.py.

antiSMASH output layout (confirmed on this project's server):
    {antismash_dir}/{sample_id}/{contig_id}.region001.gbk

A single contig can have MORE THAN ONE antiSMASH region (region001, region002, ...).
Each region corresponds to a distinct row in bgc_domain_arrays.tsv, distinguished
by BGC_start/BGC_end (see 01_extract_domains.py's bgc_id convention:
sample_id__contig_id__BGC_start). Earlier versions of this script globbed ALL
region*.gbk files for a contig and merged their CDS into a single FASTA regardless
of which row was being processed — this silently mixed proteins from unrelated
regions together whenever a contig had more than one region. This version instead
matches each row to its own region file by comparing BGC_start/BGC_end against the
region's genomic span (its 'source' feature location, which antiSMASH GBKs preserve
in original contig coordinates), and falls back to the single candidate file when
there is only one.

Only BGCs present in bgc_domain_arrays.tsv (i.e. that survived the 01 filters) are
processed, and only BGCs predicted by antiSMASH have a region GBK to extract from
(deepBGC/GECCO BGCs have no antiSMASH region file and are skipped — their
sequence_identity falls back to 0.0 in 03/04, which is the documented behavior
for BGCs without backbone domains).

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


def region_span(gbk_path: pathlib.Path):
    """
    Return (start, end) genomic coordinates for a region GBK, taken from its
    'source' feature location (antiSMASH preserves original contig coordinates
    here). Returns None if unavailable.
    """
    try:
        record = next(SeqIO.parse(str(gbk_path), "genbank"))
    except StopIteration:
        return None
    for feature in record.features:
        if feature.type == "source":
            return int(feature.location.start), int(feature.location.end)
    return None


def pick_region_gbk(candidates: list[pathlib.Path], bgc_start, bgc_end):
    """
    Choose the region GBK whose genomic span best overlaps [bgc_start, bgc_end].
    With a single candidate, return it directly (no coordinates needed).
    """
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    if bgc_start is None or bgc_end is None:
        # Can't disambiguate without coordinates; bail rather than guess wrong.
        return None

    best, best_overlap = None, -1
    for gbk in candidates:
        span = region_span(gbk)
        if span is None:
            continue
        start, end = span
        overlap = min(end, bgc_end) - max(start, bgc_start)
        if overlap > best_overlap:
            best, best_overlap = gbk, overlap

    return best if best_overlap is not None and best_overlap > 0 else None


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

    if "bgc_id" not in df.columns:
        print("[00] ERROR: no 'bgc_id' column in input. Re-run 01_extract_domains.py "
              "with the current version.", file=sys.stderr)
        sys.exit(1)

    df["BGC_start"] = pd.to_numeric(df["BGC_start"], errors="coerce")
    df["BGC_end"]   = pd.to_numeric(df["BGC_end"], errors="coerce")

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
    n_ambiguous = 0
    n_empty = 0

    for row in df.itertuples(index=False):
        bgc_id = row.bgc_id
        sample_id = row.sample_id
        contig_id = row.contig_id

        sample_dir = antismash_dir / sample_id
        if not sample_dir.is_dir():
            print(f"[00] WARNING: no antiSMASH dir for sample {sample_id}", flush=True)
            n_missing_dir += 1
            continue

        candidates = find_region_gbks(sample_dir, contig_id)
        if not candidates:
            print(f"[00] WARNING: no region GBK for {bgc_id} in {sample_dir}", flush=True)
            n_missing_gbk += 1
            continue

        gbk = pick_region_gbk(candidates, row.BGC_start, row.BGC_end)
        if gbk is None:
            print(f"[00] WARNING: could not disambiguate region GBK for {bgc_id} "
                  f"among {len(candidates)} candidates in {sample_dir} "
                  f"(BGC_start={row.BGC_start}, BGC_end={row.BGC_end})", flush=True)
            n_ambiguous += 1
            continue

        seqs = extract_translations(gbk)
        if not seqs:
            print(f"[00] WARNING: no CDS translations found for {bgc_id} in {gbk.name}",
                  flush=True)
            n_empty += 1
            continue

        out_path = output_dir / f"{bgc_id}.faa"
        with open(out_path, "w") as fh:
            for protein_id, seq in seqs.items():
                fh.write(f">{protein_id}\n{seq}\n")
        n_written += 1

    print(f"\n[00] Done. {n_written} FASTA files written to {output_dir}", flush=True)
    if n_missing_dir or n_missing_gbk or n_ambiguous or n_empty:
        print(f"[00] Skipped: {n_missing_dir} missing sample dir, "
              f"{n_missing_gbk} missing region GBK, "
              f"{n_ambiguous} ambiguous multi-region matches, "
              f"{n_empty} empty translations.", flush=True)
    if n_written == 0:
        print("[00] ERROR: no FASTA files were written. Check --antismash-dir path.",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()