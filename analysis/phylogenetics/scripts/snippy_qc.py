#!/usr/bin/env python3
"""
Required Snippy QC diagnostics (per the plan's "Verification" and "Risks"
sections): (a) per-sample callable-reference-fraction cross-tabbed against
GrafGen's discrete Refpop and against self-reported Continent, to check for
ancestry/geography-correlated coverage bias from mapping against a single
reference; (b) flag genomes below 50% callable fraction as contamination/
misassembly candidates the length filter alone can't catch; (c) a pairwise
core-SNP-distance scan on core.aln for near-duplicate genomes; (d) confirm
the reference isn't itself one of the 554 study genomes.
"""
import numpy as np
import pandas as pd
from Bio import SeqIO

BASE = "analysis/phylogenetics/results/snippy/core"
core_txt = pd.read_csv(f"{BASE}/core.txt", sep="\t")
core_txt = core_txt[core_txt["ID"] != "Reference"].copy()
core_txt["callable_fraction"] = core_txt["ALIGNED"] / core_txt["LENGTH"]

grafgen = pd.read_csv("analysis/phylogenetics/results/grafgen/ancestry.tsv", sep="\t")
meta = pd.read_csv("rawdata/metadata/gbk_hp_genomes.tsv", sep="\t")
meta["Continent"] = meta["Continent"].str.strip()

df = core_txt.merge(grafgen[["Sample", "Refpop"]], left_on="ID", right_on="Sample", how="left")
df = df.merge(meta[["Genome_ID", "Continent"]], left_on="ID", right_on="Genome_ID", how="left")

print("=== (d) Reference-in-cohort check ===")
print("'NC_000915'/'26695' in study accession list:",
      any(df["ID"].str.contains("000915|26695", case=False, na=False)))

print("\n=== (a) Callable fraction summary ===")
print(df["callable_fraction"].describe())

print("\n=== (a) Callable fraction by Refpop (GrafGen) ===")
print(df.groupby("Refpop")["callable_fraction"].agg(["mean", "std", "count"]).sort_values("mean"))

print("\n=== (a) Callable fraction by Continent (self-reported) ===")
print(df.groupby("Continent")["callable_fraction"].agg(["mean", "std", "count"]).sort_values("mean"))

print("\n=== (b) Genomes below 50% callable fraction ===")
low = df[df["callable_fraction"] < 0.50]
print(f"{len(low)} genome(s) below 50% callable:")
if len(low):
    print(low[["ID", "callable_fraction", "Refpop", "Continent"]].to_string(index=False))

df.to_csv(f"{BASE}/../../grafgen/qc_callable_fraction.tsv", sep="\t", index=False)

print("\n=== (c) Near-duplicate pairwise-SNP-distance scan ===")
records = {r.id: str(r.seq).upper() for r in SeqIO.parse(f"{BASE}/core.aln", "fasta") if r.id != "Reference"}
ids = list(records.keys())
n = len(ids)
L = len(next(iter(records.values())))
mat = np.frombuffer("".join(records[i] for i in ids).encode(), dtype=np.uint8).reshape(n, L)

# Only compare at variable (non-constant) columns -- core.aln is already the
# variant-site alignment (snippy-core's *core SNP* alignment), so all columns
# are informative; no further column filtering needed.
near_dup_pairs = []
CHUNK = 50
for i0 in range(0, n, CHUNK):
    i1 = min(i0 + CHUNK, n)
    block = mat[i0:i1]  # (chunk, L)
    # Hamming distance of this block against all samples j > i (upper triangle only, done globally below)
    diffs = (block[:, None, :] != mat[None, :, :]).sum(axis=2)  # (chunk, n)
    for bi, i in enumerate(range(i0, i1)):
        for j in range(i + 1, n):
            d = diffs[bi, j]
            if d <= 5:  # near-zero pairwise SNP distance threshold
                near_dup_pairs.append((ids[i], ids[j], int(d)))

print(f"Core alignment: {n} samples x {L} sites")
print(f"Pairs with <=5 SNP core-genome distance: {len(near_dup_pairs)}")
for a, b, d in sorted(near_dup_pairs, key=lambda x: x[2])[:30]:
    print(f"  {a}  <->  {b}   dist={d}")

pd.DataFrame(near_dup_pairs, columns=["sample_a", "sample_b", "snp_distance"]).to_csv(
    f"{BASE}/../../grafgen/near_duplicates.tsv", sep="\t", index=False
)
