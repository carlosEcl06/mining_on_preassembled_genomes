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
    - HMMER (hmmalign, hmmfetch, hmmscan) must be in PATH
    - Pfam HMM database (Pfam-A.hmm) must be available (see --pfam-hmm)
    - Input FASTA files per BGC (protein sequences) must be available (see --fasta-dir)

Usage:
    python 03_backbone_identity.py \
        --domains   ../results/bgc_domain_arrays.tsv \
        --pairs     ../results/jaccard_pairs.tsv \
        --fasta-dir ../../rawdata/fasta/proteins \
        --pfam-hmm  ../../rawdata/pfam/Pfam-A.hmm \
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


# Backbone domains to use for hmmalign, per class. These MUST be real Pfam
# accessions fetchable via `hmmfetch Pfam-A.hmm <accession>` (versioned form
# resolved automatically at runtime — see build_accession_index()).
#
# Keys are now the ACTUAL Product_class strings from combgc/antiSMASH
# (backbone_class is sourced directly from Product_class as of 2026-07-29 —
# see 01_extract_domains.py's module docstring for why the old Pfam-matching
# heuristic was dropped: it produced false-positive rates as high as 94%
# ("Saccharide" BGCs mislabeled "Terpene"), because the marker accessions it
# used were never independently verified against this dataset).
#
# Only antiSMASH BGCs reach this script at all (00_extract_bgc_fastas.py only
# extracts FASTA for Prediction_tool == antiSMASH), so only antiSMASH's own
# three Product_class values are listed below — deepBGC/GECCO Product_class
# values (Polyketide, Saccharide, Other, RiPP, NRP, Terpene,
# Polyketide-Terpene, Unknown) are structurally unreachable here regardless
# of what accessions might describe them, since there's no FASTA to align.
#
# Provenance per entry:
#   "Terpene-precursor": VERIFIED 2026-07-28 — ran hmmscan directly on 4
#       different antiSMASH Terpene-precursor BGC FASTAs and confirmed
#       PF00348 (Polyprenyl synthetase — the enzyme comBGC's own
#       "PT_FPPS_like" marker refers to) present in all 4/4 spot-checks.
#   "Azole-containing-RiPP": NOT YET independently verified via hmmscan.
#       PF02624 (YcaO domain) is a well-motivated candidate — comBGC's own
#       PFAM_domains already reports the literal name "YcaO" in 19/19 of
#       these BGCs, and the YcaO cyclodehydratase is the textbook
#       backbone-defining enzyme for azole-containing RiPPs — but this is
#       inference from nomenclature, not the same direct verification done
#       for Terpene-precursor. Recommend the same spot-check before trusting
#       results for this class: run hmmscan on 2-3 of these BGCs' .faa files
#       and confirm PF02624 appears.
#   "Arylpolyene": left EMPTY deliberately. n=1 BGC in this dataset, so there
#       can never be a within-class pair to compute identity for regardless
#       of which accession is chosen — not worth guessing.
BACKBONE_DOMAINS = {
    "Terpene-precursor":     ["PF00348"],
    "Azole-containing-RiPP": ["PF02624"],   # unverified — see provenance note above
    "Arylpolyene":           [],            # n=1, no pairs possible; intentionally empty
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


# Module-level counters so we surface the first couple of HMMER failures with
# their actual stderr, instead of silently returning {} for all 382k pairs like
# before — that silence is exactly what hid the real cause of "0 pairs had
# backbone domain sequences" across two previous runs.
_HMMSCAN_FAILURES_SHOWN = 0
_HMMFETCH_FAILURES_SHOWN = 0
_HMMALIGN_FAILURES_SHOWN = 0
_MAX_FAILURES_SHOWN = 3


def run_hmmscan_domtbl(fasta: dict[str, str], pfam_hmm: str) -> list[tuple[str, str]]:
    """
    Run hmmscan ONCE for this BGC's full protein set against the full Pfam-A.hmm
    database, returning (bare_accession, protein_id) for every domain hit.

    IMPORTANT: hmmscan is never told which domain we're looking for — it always
    scans the query proteins against the entire Pfam-A.hmm database regardless.
    The previous version called this once per (bgc_id, domain) combination, but
    since the underlying hmmscan command and its output were byte-for-byte
    identical across those calls (only the post-hoc filter differed), this
    re-ran the same multi-second scan 4-6x per BGC for no benefit. This version
    runs hmmscan once per bgc_id and returns all hit rows unfiltered; filtering
    for a specific domain is then a cheap in-memory operation (see
    get_aligned_domain_seqs) — the set of hits returned is identical to before,
    just computed once instead of once per domain.

    hmmscan --domtblout column layout (whitespace-separated):
        [0] target name        — the PROFILE's human-readable name (e.g. "Alcohol_dh"),
                                  NOT its Pfam accession
        [1] target accession   — the Pfam accession, e.g. "PF00107.30" (versioned)
        [2] tlen
        [3] query name         — the query PROTEIN id
        ...
    """
    if not fasta:
        return []

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
        hits = []
        with open(tmp_tbl_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                accession  = parts[1].split(".")[0]   # Pfam accession, version stripped
                protein_id = parts[3]                  # query protein id
                hits.append((accession, protein_id))
        return hits
    except subprocess.CalledProcessError as e:
        global _HMMSCAN_FAILURES_SHOWN
        if _HMMSCAN_FAILURES_SHOWN < _MAX_FAILURES_SHOWN:
            _HMMSCAN_FAILURES_SHOWN += 1
            print(f"[03] hmmscan FAILED (showing first {_MAX_FAILURES_SHOWN} failures only):\n"
                  f"     cmd: hmmscan --domtblout {tmp_tbl_path} --noali -E 1e-5 "
                  f"{pfam_hmm} {tmp_fa_path}\n"
                  f"     stderr: {e.stderr.decode(errors='replace').strip()}", flush=True)
        return []
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

    except subprocess.CalledProcessError as e:
        global _HMMFETCH_FAILURES_SHOWN, _HMMALIGN_FAILURES_SHOWN
        cmd = e.cmd[0] if e.cmd else "?"
        counter_name = "_HMMFETCH_FAILURES_SHOWN" if cmd == "hmmfetch" else "_HMMALIGN_FAILURES_SHOWN"
        shown = _HMMFETCH_FAILURES_SHOWN if cmd == "hmmfetch" else _HMMALIGN_FAILURES_SHOWN
        if shown < _MAX_FAILURES_SHOWN:
            if cmd == "hmmfetch":
                _HMMFETCH_FAILURES_SHOWN += 1
            else:
                _HMMALIGN_FAILURES_SHOWN += 1
            print(f"[03] {cmd} FAILED for domain {domain} "
                  f"(showing first {_MAX_FAILURES_SHOWN} failures only):\n"
                  f"     stderr: {e.stderr.decode(errors='replace').strip()}", flush=True)
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


def build_accession_index(pfam_hmm_path: str) -> dict[str, str]:
    """
    Map bare Pfam accession (no version) -> full versioned accession, by
    scanning 'ACC' lines in the flat-text Pfam-A.hmm file.

    hmmfetch requires an EXACT match on the indexed key, which is the
    versioned accession (e.g. 'PF00109.28'). Calling hmmfetch with a bare
    accession like 'PF00109' fails with "not found in SSI index" even
    though hmmpress ran correctly and the domain is present — this is what
    caused every backbone domain lookup to fail in earlier runs.
    """
    accession_map = {}
    with open(pfam_hmm_path, errors="replace") as fh:
        for line in fh:
            if line.startswith("ACC "):
                full = line.split()[1].strip()
                bare = full.split(".")[0]
                accession_map[bare] = full
    return accession_map


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

    # --- fail fast: verify HMMER binaries and a pressed Pfam-A.hmm exist ---
    for binary in ("hmmscan", "hmmfetch", "hmmalign"):
        check = subprocess.run(["which", binary], capture_output=True)
        if check.returncode != 0:
            raise SystemExit(f"[03] ERROR: '{binary}' not found in PATH. "
                              f"Activate the gcf_clustering conda env first.")

    pressed_exts = (".h3f", ".h3i", ".h3m", ".h3p")
    missing = [ext for ext in pressed_exts if not pathlib.Path(args.pfam_hmm + ext).exists()]
    if missing:
        raise SystemExit(
            f"[03] ERROR: {args.pfam_hmm} is missing pressed index files "
            f"({', '.join(missing)}). Run: hmmpress {args.pfam_hmm}"
        )

    # Build bare→versioned accession map once, since hmmfetch needs the exact
    # versioned key (see build_accession_index docstring).
    print(f"[03] Indexing Pfam accessions in {args.pfam_hmm} ...", flush=True)
    accession_map = build_accession_index(args.pfam_hmm)
    print(f"[03] Found {len(accession_map):,} Pfam accessions.", flush=True)

    # Resolve BACKBONE_DOMAINS to versioned accessions; warn (once) about any
    # that aren't present in this Pfam-A.hmm release rather than failing later.
    resolved_backbone_domains: dict[str, list[str]] = {}
    for cls, doms in BACKBONE_DOMAINS.items():
        resolved = []
        for d in doms:
            versioned = accession_map.get(d)
            if versioned is None:
                print(f"[03] WARNING: accession {d} (class {cls}) not found in "
                      f"{args.pfam_hmm}; skipping.", flush=True)
                continue
            resolved.append(versioned)
        resolved_backbone_domains[cls] = resolved

    # Smoke-test hmmfetch against a real, correctly versioned accession —
    # just take the first resolved accession from any class (no need to prefer
    # a specific class now that BACKBONE_DOMAINS is keyed by Product_class).
    smoke_domain = next((v for vs in resolved_backbone_domains.values() for v in vs), None)
    if smoke_domain is None:
        raise SystemExit("[03] ERROR: none of the BACKBONE_DOMAINS accessions were "
                          f"found in {args.pfam_hmm}. Check the Pfam release in use.")
    smoke = subprocess.run(["hmmfetch", args.pfam_hmm, smoke_domain], capture_output=True)
    if smoke.returncode != 0:
        raise SystemExit(
            f"[03] ERROR: hmmfetch smoke test failed against {args.pfam_hmm} "
            f"(accession {smoke_domain}).\n"
            f"stderr: {smoke.stderr.decode(errors='replace').strip()}\n"
            f"Check that Pfam-A.hmm was downloaded/pressed correctly."
        )
    print(f"[03] HMMER + Pfam-A.hmm smoke test passed ({smoke_domain}).", flush=True)

    print(f"[03] Reading domain arrays from {args.domains} ...", flush=True)
    df_domains = pd.read_csv(args.domains, sep="\t", dtype=str)
    if "bgc_id" not in df_domains.columns:
        raise SystemExit(
            "[03] ERROR: no 'bgc_id' column in --domains. Re-run 01_extract_domains.py "
            "with the current version, which writes this column directly."
        )
    domain_lookup   = df_domains.set_index("bgc_id")["domain_array"].to_dict()
    backbone_lookup = df_domains.set_index("bgc_id")["backbone_class"].to_dict()

    print(f"[03] Reading pairs from {args.pairs} ...", flush=True)
    df_pairs = pd.read_csv(args.pairs, sep="\t", dtype={"bgc_i": str, "bgc_j": str,
                                                          "jaccard": float})
    print(f"[03] {len(df_pairs):,} pairs to process.", flush=True)

    fasta_dir = pathlib.Path(args.fasta_dir)
    pfam_hmm  = args.pfam_hmm

    # Cache FASTA files and hmmscan hits per bgc_id (ONE hmmscan call per BGC,
    # not per domain — see run_hmmscan_domtbl docstring), plus per-(bgc_id,
    # domain) aligned sequences to avoid redundant hmmfetch/hmmalign calls.
    fasta_cache:   dict[str, dict[str, str]] = {}
    hmmscan_cache: dict[str, list[tuple[str, str]]] = {}
    domain_cache:  dict[tuple, list[str]] = {}   # (bgc_id, domain) → aligned seqs

    def get_fasta(bgc_id: str) -> dict[str, str]:
        if bgc_id not in fasta_cache:
            fa_path = fasta_dir / f"{bgc_id}.faa"
            fasta_cache[bgc_id] = read_fasta(fa_path) if fa_path.exists() else {}
        return fasta_cache[bgc_id]

    def get_hmmscan_hits(bgc_id: str) -> list[tuple[str, str]]:
        if bgc_id not in hmmscan_cache:
            fasta = get_fasta(bgc_id)
            hmmscan_cache[bgc_id] = run_hmmscan_domtbl(fasta, pfam_hmm)
        return hmmscan_cache[bgc_id]

    def get_aligned_domain_seqs(bgc_id: str, domain: str) -> list[str]:
        key = (bgc_id, domain)
        if key not in domain_cache:
            fasta       = get_fasta(bgc_id)
            domain_bare = domain.split(".")[0]
            all_hits    = get_hmmscan_hits(bgc_id)
            hit_seqs    = {pid: fasta[pid] for acc, pid in all_hits
                           if acc == domain_bare and pid in fasta}
            aligned     = align_sequences(hit_seqs, domain, pfam_hmm) if hit_seqs else {}
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
            domains_in_class = resolved_backbone_domains.get(cls, [])
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