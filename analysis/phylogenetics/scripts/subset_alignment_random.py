#!/usr/bin/env python3
"""
subset_alignment_random.py

Randomly subsamples N sequences from a FASTA alignment, for the Gubbins pilot
runs (--n 50 for pilot A, --n 150-200 for pilot B) — used to estimate
full-scale (552-genome) Gubbins runtime before committing a long-walltime
full-node job to it. See the project plan's "Risks" section: Gubbins'
iteration count to convergence isn't simply linear in N, so two pilots (not
one) are used to fit an actual scaling curve.

The Reference sequence (present in every snippy-core alignment) is always
kept, on top of the N randomly-sampled study genomes, so the pilot alignment
has the same structure (N+1 sequences) as a real subset of the full cohort.

Usage:
    python subset_alignment_random.py \
        --input core.aln \
        --output pilot_a_subset.aln \
        --n 50 \
        --seed 42
"""

import argparse
import random

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n", type=int, default=50, help="number of study genomes to sample (Reference is always kept on top of this)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--reference-id", default="Reference", help="sequence ID to always keep (default: 'Reference', snippy-core's convention)")
    return parser.parse_args()


def main():
    args = parse_args()
    records = list(SeqIO.parse(args.input, "fasta"))

    ref_records = [r for r in records if r.id == args.reference_id]
    study_records = [r for r in records if r.id != args.reference_id]

    if args.n >= len(study_records):
        print(f"[subset] --n ({args.n}) >= study genomes available ({len(study_records)}); using all of them.")
        chosen = study_records
    else:
        random.seed(args.seed)
        chosen = random.sample(study_records, args.n)

    out_records = ref_records + chosen
    SeqIO.write(out_records, args.output, "fasta")
    print(f"[subset] wrote {len(chosen)} study genome(s) + {len(ref_records)} reference "
          f"= {len(out_records)} total sequences (of {len(records)} available) to {args.output}")


if __name__ == "__main__":
    main()
