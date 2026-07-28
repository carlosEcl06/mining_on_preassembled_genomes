#!/usr/bin/env python3
"""
03_backbone_identity.py

Estimates pairwise sequence identity between backbone enzyme domains for BGC pairs
that passed the Jaccard threshold (from 02_jaccard_matrix.py).

For each pair of BGCs sharing a backbone class, hmmalign is used to align backbone
domains against their respective Pfam HMM profiles. Sequence identity is then
calculated as the mean pairwise identity between matched domain pairs. When multiple
copies of a domain are present, the Hungarian algorithm is used to find the
highest-similarity pairing configuration.

Only pairs present in jaccard_pairs.tsv are processed. BGC pairs without any
backbone domain in common are assigned sequence_identity = 0.0.

Requirements:
    - HMMER (hmmalign) must be in PATH
    - Pfam HMM database (Pfam-A.hmm) must be available (see --pfam-hmm)
    - Input FASTA files per BGC (protein sequences) must be available (see --fasta-dir)

Usage:
    python 03_backbone_identity.py \
        --domains   ../results/bgc_domain_arrays.tsv \
        --pairs     ../results/jaccard_pairs.tsv \
        --fasta-dir ../../rawdata/fasta/proteins \
        --pfam-hmm  ../../rawdata/Pfam-A.hmm \
        --output    ../results/backbone_identity.tsv
"""

import argparse
import pathlib
import subprocess
import tempfile
from collections import defaultdict
from io import StringIO

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


# Backbone domains to use for hmmalign, per class.
# Only domains present in BACKBONE_MAP in 01_extract_domains.py are considered.
BACKBONE_DOMAINS = {
    "Terpene":   ["PT_FPPS_like", "PF02353", "PF04909", "PF03544", "PF13088"],
    "NRPS_like": ["PF00107", "PF01262", "PF02826", "PF03446", "PF08240", "PF01501"],
    "PKS":       ["PF00109", "PF02801", "PF00550"],
    "RiPP":      ["PF00398", "PF02624", "PF02463"],
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--domains",   required=True, help="Path to bgc_domain_arrays.tsv")
    parser.add_argument("--pairs",     required=True, help="Path to jaccard_pairs.tsv")
    parser.add_argument("--fasta-dir", required=True,
                        help="Directory with per-BGC protein FASTA files "
                             "(named as {bgc_id}.faa)")
    parser.add_argument("--pfam-hmm",  required=True, help="Path to Pfam-A.hmm")
    parser.add_argument("--output",    required=True, help="Path to output TSV")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# FASTA / alignment helpers
# ---------------------------------------------------------------------------

def read_fasta(path: pathlib.Path) -> dict[str, str]:
    """Parse a FASTA file into {header: sequence} dict."""
    seqs = {}
    header, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(buf)
                header, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if header:
        seqs[header] = "".join(buf)
    return seqs


def extract_domain_sequences(fasta: dict[str, str],
                              domain: str,
                              pfam_hmm: str) -> dict[str, str]:
    """
    Run hmmscan to identify which proteins in `fasta` contain `domain`,
    then return those sequences keyed by protein ID.
    """
    if not fasta:
        return {}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as tmp_fa:
        for hdr, seq in fasta.items():
            tmp_fa.write(f">{hdr}\n{seq}\n")
        tmp_fa_path = tmp_fa.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".tbl", delete=False) as tmp_tbl:
        tmp_tbl_path = tmp_tbl.name

    try:
        subprocess.run(
            ["hmmscan", "--domtblout", tmp_tbl_path, "--noali", "-E", "1e-5",
             pfam_hmm, tmp_fa_path],
            capture_output=True, check=True
        )
        hits = {}
        with open(tmp_tbl_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                query_domain = parts[0]   # Pfam accession
                protein_id   = parts[3]   # query protein
                if query_domain == domain and protein_id in fasta:
                    hits[protein_id] = fasta[protein_id]
        return hits
    except subprocess.CalledProcessError:
        return {}
    finally:
        pathlib.Path(tmp_fa_path).unlink(missing_ok=True)
        pathlib.Path(tmp_tbl_path).unlink(missing_ok=True)


def align_sequences(seqs: dict[str, str], domain: str, pfam_hmm: str) -> dict[str, str]:
    """
    Run hmmalign to align sequences against a Pfam HMM profile.
    Returns {seq_id: aligned_sequence} with gap characters intact.
    """
    if not seqs:
        return {}

    # Extract single HMM profile for this domain
    with tempfile.NamedTemporaryFile(mode="w", suffix=".hmm", delete=False) as tmp_hmm:
        tmp_hmm_path = tmp_hmm.name

    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa", delete=False) as tmp_fa:
        for hdr, seq in seqs.items():
            tmp_fa.write(f">{hdr}\n{seq}\n")
        tmp_fa_path = tmp_fa.name

    try:
        # hmmfetch extracts one profile from Pfam-A.hmm by accession
        fetch = subprocess.run(
            ["hmmfetch", pfam_hmm, domain],
            capture_output=True, check=True
        )
        with open(tmp_hmm_path, "wb") as fh:
            fh.write(fetch.stdout)

        result = subprocess.run(
            ["hmmalign", "--trim", "--outformat", "afa", tmp_hmm_path, tmp_fa_path],
            capture_output=True, check=True
        )
        aligned = {}
        header, buf = None, []
        for line in result.stdout.decode().splitlines():
            if line.startswith(">"):
                if header:
                    aligned[header] = "".join(buf)
                header, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if header:
            aligned[header] = "".join(buf)
        return aligned

    except subprocess.CalledProcessError:
        return {}
    finally:
        pathlib.Path(tmp_hmm_path).unlink(missing_ok=True)
        pathlib.Path(tmp_fa_path).unlink(missing_ok=True)


def pairwise_identity(seq_a: str, seq_b: str) -> float:
    """
    Compute percent identity between two aligned sequences.
    Only positions where neither sequence has a gap are considered.
    """
    matches = total = 0
    for a, b in zip(seq_a, seq_b):
        if a != "-" and b != "-":
            total += 1
            if a == b:
                matches += 1
    return matches / total if total > 0 else 0.0


def hungarian_best_identity(seqs_i: list[str], seqs_j: list[str]) -> float:
    """
    Given two lists of aligned domain sequences, use the Hungarian algorithm
    to find the pairing configuration that maximises total sequence identity,
    then return the mean identity of the best pairs.
    """
    n, m = len(seqs_i), len(seqs_j)
    size  = max(n, m)
    cost  = np.zeros((size, size), dtype=np.float64)
    for r, sa in enumerate(seqs_i):
        for c, sb in enumerate(seqs_j):
            cost[r, c] = 1.0 - pairwise_identity(sa, sb)   # cost = 1 - identity

    row_ind, col_ind = linear_sum_assignment(cost)
    valid = [(r, c) for r, c in zip(row_ind, col_ind) if r < n and c < m]
    if not valid:
        return 0.0
    identities = [1.0 - cost[r, c] for r, c in valid]
    return float(np.mean(identities))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    print(f"[03] Reading domain arrays from {args.domains} ...", flush=True)
    df_domains = pd.read_csv(args.domains, sep="\t", dtype=str)
    df_domains["bgc_id"] = df_domains["sample_id"] + "__" + df_domains["contig_id"]
    domain_lookup  = df_domains.set_index("bgc_id")["domain_array"].to_dict()
    backbone_lookup = df_domains.set_index("bgc_id")["backbone_class"].to_dict()

    print(f"[03] Reading pairs from {args.pairs} ...", flush=True)
    df_pairs = pd.read_csv(args.pairs, sep="\t", dtype={"bgc_i": str, "bgc_j": str,
                                                          "jaccard": float})
    print(f"[03] {len(df_pairs):,} pairs to process.", flush=True)

    fasta_dir = pathlib.Path(args.fasta_dir)
    pfam_hmm  = args.pfam_hmm

    # Cache FASTA files and domain sequences to avoid redundant hmmscan calls
    fasta_cache:  dict[str, dict[str, str]] = {}
    domain_cache: dict[tuple, dict[str, str]] = {}   # (bgc_id, domain) → aligned seqs

    def get_fasta(bgc_id: str) -> dict[str, str]:
        if bgc_id not in fasta_cache:
            fa_path = fasta_dir / f"{bgc_id}.faa"
            fasta_cache[bgc_id] = read_fasta(fa_path) if fa_path.exists() else {}
        return fasta_cache[bgc_id]

    def get_aligned_domain_seqs(bgc_id: str, domain: str) -> list[str]:
        key = (bgc_id, domain)
        if key not in domain_cache:
            fasta    = get_fasta(bgc_id)
            hit_seqs = extract_domain_sequences(fasta, domain, pfam_hmm)
            aligned  = align_sequences(hit_seqs, domain, pfam_hmm) if hit_seqs else {}
            domain_cache[key] = list(aligned.values())
        return domain_cache[key]

    # ---------------------------------------------------------------------------
    # Compute sequence identity per pair
    # ---------------------------------------------------------------------------
    seq_identities = []
    n_with_backbone = 0

    for idx, row in df_pairs.iterrows():
        if idx % 10_000 == 0:
            print(f"[03] Processing pair {idx:,} / {len(df_pairs):,} ...", flush=True)

        bgc_i, bgc_j = row["bgc_i"], row["bgc_j"]
        classes_i = set(backbone_lookup.get(bgc_i, "Unknown").split("|"))
        classes_j = set(backbone_lookup.get(bgc_j, "Unknown").split("|"))
        shared_classes = (classes_i & classes_j) - {"Unknown"}

        if not shared_classes:
            seq_identities.append(0.0)
            continue

        # For each shared backbone class, compute mean identity across its domains
        class_identities = []
        for cls in shared_classes:
            domains_in_class = BACKBONE_DOMAINS.get(cls, [])
            domain_identities = []
            for domain in domains_in_class:
                seqs_i = get_aligned_domain_seqs(bgc_i, domain)
                seqs_j = get_aligned_domain_seqs(bgc_j, domain)
                if seqs_i and seqs_j:
                    identity = hungarian_best_identity(seqs_i, seqs_j)
                    domain_identities.append(identity)
            if domain_identities:
                class_identities.append(float(np.mean(domain_identities)))

        if class_identities:
            seq_identities.append(float(np.mean(class_identities)))
            n_with_backbone += 1
        else:
            seq_identities.append(0.0)

    print(f"[03] {n_with_backbone:,} pairs had backbone domain sequences for alignment.",
          flush=True)

    # ---------------------------------------------------------------------------
    # Output
    # ---------------------------------------------------------------------------
    df_pairs["sequence_identity"] = seq_identities
    pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df_pairs.to_csv(args.output, sep="\t", index=False)
    print(f"[03] Output written to {args.output}", flush=True)


if __name__ == "__main__":
    main()