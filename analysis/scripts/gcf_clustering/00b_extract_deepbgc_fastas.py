#!/usr/bin/env python3
"""
00b_extract_deepbgc_fastas.py

Extracts per-BGC protein FASTA files from deepBGC's own pre-cut BGC GenBank
output, for use as --fasta-dir input to 03_backbone_identity.py. Writes into
the SAME directory 00_extract_bgc_fastas.py uses for antiSMASH — bgc_id is
globally unique across tools (sample_id__contig_id__BGC_start), so sharing
one directory carries no collision risk.

deepBGC output layout:
    {deepbgc_dir}/{sample_id}/{sample_id}.bgc.gbk

Unlike antiSMASH (one file per region), deepBGC packs all of a sample's
predicted BGC regions into a single multi-record GenBank file — one
LOCUS/record per region, already cut down to that region's DNA. If a
contig_id appears as only one LOCUS in the sample's .bgc.gbk, that record IS
the BGC; coordinate matching (via the record's ACCESSION/VERSION field, set
by deepBGC to a synthetic "{contig_id}_{start}-{end}[.version]" string, e.g.
'contig_21_6-9057.1') is only needed when two or more predicted regions land
on the same contig_id within one sample.

Note: records are grouped by Biopython's record.name (from the LOCUS line,
equal to the real contig_id), not record.id (which holds that synthetic
ACCESSION string).

Only BGCs present in bgc_domain_arrays.tsv, Prediction_tool == deepBGC, are
processed here (antiSMASH is 00_extract_bgc_fastas.py; GECCO has no embedded
protein sequence in its own output, so it's not covered).

Usage:
    python 00b_extract_deepbgc_fastas.py \
        --domains     ../results/bgc_domain_arrays.tsv \
        --deepbgc-dir ../../funcscan/results/bgc/deepbgc \
        --output-dir  ../../../rawdata/fastas/proteins
"""

import argparse
import pathlib
import re
import sys

from Bio import SeqIO  # requires: biopython (see requirements.txt)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--domains", required=True, help="Path to bgc_domain_arrays.tsv")
    parser.add_argument(
        "--deepbgc-dir", required=True,
        help="Path to funcscan deepBGC results dir "
             "(contains {sample_id}/{sample_id}.bgc.gbk)",
    )
    parser.add_argument("--output-dir", required=True,
                        help="Directory to write {bgc_id}.faa files "
                             "(same directory 00_extract_bgc_fastas.py writes to)")
    parser.add_argument(
        "--only-tool", default="deepBGC",
        help="Only extract BGCs whose Prediction_tool matches this value "
             "(default: deepBGC)",
    )
    return parser.parse_args()


def load_sample_records(gbk_path: pathlib.Path) -> dict[str, list]:
    """Parse a sample's multi-record .bgc.gbk once, grouped by contig_id.

    Groups by record.name (parsed from the LOCUS line, the real contig_id),
    NOT record.id (the synthetic ACCESSION/VERSION string deepBGC assigns).
    """
    by_contig: dict[str, list] = {}
    for record in SeqIO.parse(str(gbk_path), "genbank"):
        by_contig.setdefault(record.name, []).append(record)
    return by_contig


ACCESSION_COORDS_RE = re.compile(r"^.+_(\d+)-(\d+)(?:\.\d+)?$")


def record_span(record):
    """
    Return (start, end) original genomic coordinates for a record.

    Primary source: the ACCESSION/VERSION field, which deepBGC sets to
    "{contig_id}_{start}-{end}[.version]" (e.g. 'contig_21_6-9057.1').
    Falls back to the 'source' feature location if that doesn't parse.
    """
    for accession in record.annotations.get("accessions", []):
        m = ACCESSION_COORDS_RE.match(accession)
        if m:
            return int(m.group(1)), int(m.group(2))
    m = ACCESSION_COORDS_RE.match(record.id)
    if m:
        return int(m.group(1)), int(m.group(2))

    for feature in record.features:
        if feature.type == "source":
            return int(feature.location.start), int(feature.location.end)
    return None


def pick_record(candidates: list, bgc_start, bgc_end, bgc_length):
    """
    Choose the correct record for this row among same-contig_id candidates.
    - Single candidate: return it directly, no ambiguity possible.
    - Multiple candidates: try genomic-coordinate overlap via the 'source'
      feature first; if that's unavailable or inconclusive, fall back to the
      closest record-length match against BGC_length (tolerance 5%). Returns
      None (rather than guessing) if neither signal resolves the ambiguity.
    """
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    if bgc_start is not None and bgc_end is not None:
        best, best_overlap = None, 0
        for rec in candidates:
            span = record_span(rec)
            if span is None:
                continue
            start, end = span
            overlap = min(end, bgc_end) - max(start, bgc_start)
            if overlap > best_overlap:
                best, best_overlap = rec, overlap
        if best is not None:
            return best

    if bgc_length is not None:
        best, best_diff = None, None
        for rec in candidates:
            diff = abs(len(rec.seq) - bgc_length)
            if best_diff is None or diff < best_diff:
                best, best_diff = rec, diff
        if best is not None and best_diff <= 0.05 * bgc_length:
            return best

    return None


def extract_translations(record) -> dict[str, str]:
    """Return {protein_id: translation} for all CDS features in a record."""
    seqs = {}
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

    print(f"[00b] Reading {args.domains} ...", flush=True)
    df = pd.read_csv(args.domains, sep="\t", dtype=str)

    if "bgc_id" not in df.columns:
        print("[00b] ERROR: no 'bgc_id' column in input. Re-run 01_extract_domains.py "
              "with the current version.", file=sys.stderr)
        sys.exit(1)

    df["BGC_start"]  = pd.to_numeric(df["BGC_start"], errors="coerce")
    df["BGC_end"]    = pd.to_numeric(df["BGC_end"], errors="coerce")
    df["BGC_length"] = pd.to_numeric(df["BGC_length"], errors="coerce")

    if args.only_tool:
        n_before = len(df)
        df = df[df["Prediction_tool"] == args.only_tool].copy()
        print(f"[00b] Filtered to Prediction_tool == '{args.only_tool}': "
              f"{len(df)} / {n_before} BGCs.", flush=True)

    deepbgc_dir = pathlib.Path(args.deepbgc_dir)
    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_written = 0
    n_missing_gbk = 0
    n_missing_contig = 0
    n_ambiguous = 0
    n_empty = 0

    sample_cache: dict[str, dict] = {}  # parsed .bgc.gbk per sample, grouped by contig_id

    for row in df.itertuples(index=False):
        bgc_id = row.bgc_id
        sample_id = row.sample_id
        contig_id = row.contig_id

        if sample_id not in sample_cache:
            gbk_path = deepbgc_dir / sample_id / f"{sample_id}.bgc.gbk"
            sample_cache[sample_id] = load_sample_records(gbk_path) if gbk_path.exists() else {}

        by_contig = sample_cache[sample_id]
        if not by_contig:
            print(f"[00b] WARNING: no .bgc.gbk (or empty) for sample {sample_id}", flush=True)
            n_missing_gbk += 1
            continue

        candidates = by_contig.get(contig_id, [])
        if not candidates:
            print(f"[00b] WARNING: no record for contig {contig_id} in "
                  f"{sample_id}.bgc.gbk ({bgc_id})", flush=True)
            n_missing_contig += 1
            continue

        record = pick_record(candidates, row.BGC_start, row.BGC_end, row.BGC_length)
        if record is None:
            print(f"[00b] WARNING: could not disambiguate {len(candidates)} records "
                  f"for contig {contig_id} in {sample_id}.bgc.gbk ({bgc_id})", flush=True)
            n_ambiguous += 1
            continue

        seqs = extract_translations(record)
        if not seqs:
            print(f"[00b] WARNING: no CDS translations found for {bgc_id}", flush=True)
            n_empty += 1
            continue

        out_path = output_dir / f"{bgc_id}.faa"
        with open(out_path, "w") as fh:
            for protein_id, seq in seqs.items():
                fh.write(f">{protein_id}\n{seq}\n")
        n_written += 1

    print(f"\n[00b] Done. {n_written} FASTA files written to {output_dir}", flush=True)
    if n_missing_gbk or n_missing_contig or n_ambiguous or n_empty:
        print(f"[00b] Skipped: {n_missing_gbk} missing .bgc.gbk, "
              f"{n_missing_contig} missing contig record, "
              f"{n_ambiguous} ambiguous multi-record matches, "
              f"{n_empty} empty translations.", flush=True)
    if n_written == 0:
        print("[00b] ERROR: no FASTA files were written. Check --deepbgc-dir path.",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()