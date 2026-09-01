# Mining Biosynthetic Gene Clusters in *H. pylori* Genomes

Genome mining of biosynthetic gene clusters (BGCs) across 554 preassembled
*Helicobacter pylori* genomes (179 from gastric-cancer samples, 375 from
precancerous cases — chronic gastritis and intestinal metaplasia), followed by
tool-agnostic Gene Cluster Family (GCF) clustering and statistical testing for
association with clinical phenotype.

This is a Master's project at the Laboratory of Molecular Genetics and Omics
(LGMO), Ribeirão Preto Medical School – University of São Paulo (USP), building
on an earlier undergraduate pilot study that mined BGCs from two
metagenome-assembled *H. pylori* genomes.

**Authors:** Carlos Eugênio Costa de Lima, Wilson Araújo da Silva Junior

**Notice:** This README has been mostly AI-generated, thus it may contain some inconsistencies.

## Pipeline overview

1. **Genome retrieval** (`rawdata/`) — 554 preassembled genomes downloaded from
   NCBI via `datasets download genome`.
2. **BGC mining** (`analysis/funcscan/`) — [nf-core/funcscan](https://nf-co.re/funcscan)
   (v3.0.0), running antiSMASH 8.0.1, DeepBGC 0.1.31, and GECCO 0.9.10, with
   Bakta/InterProScan functional annotation. Results are aggregated by comBGC
   into a single harmonized table of Pfam domains per BGC
   (`analysis/R/copied_from_funcscan_results/combgc_complete_summary.tsv`).
3. **GCF clustering** (`analysis/scripts/gcf_clustering/`) — a tool-agnostic
   clustering workflow (Python), since BiG-SCAPE requires antiSMASH-formatted
   input and is incompatible with DeepBGC/GECCO output. Pairwise BGC similarity
   is adapted from Robey et al. (2021):

   - Pfam domain arrays are compared via Jaccard similarity (pairs below 0.10
     excluded).
   - For BGCs with a verified backbone domain (antiSMASH: Terpene-precursor,
     Azole-containing-RiPP; DeepBGC: NRP, Saccharide, Polyketide, Terpene,
     Polyketide-Terpene), sequence identity between backbone domains is
     estimated via `hmmalign` against Pfam-A, and the full score is used:
     `sqrt(0.8 * identity + 0.2 * Jaccard)`.
   - All other BGCs fall back to `sqrt(0.2 * Jaccard)`.
   - Final GCFs are delineated with DBSCAN (`epsilon = 0.56`,
     `min_samples = 2`, chosen via an epsilon-sensitivity sweep) on the
     resulting distance matrix.

   See `analysis/scripts/gcf_clustering/run_clustering.sh` for the full
   orchestration and `requirements.txt` for dependencies.
4. **Statistical analysis and reporting** (`analysis/R/RMarkdown/`) — Quarto
   notebooks characterizing the BGC/GCF inventory and testing GCF–phenotype
   association (Fisher's exact test, Cochran–Mantel–Haenszel, and multivariable
   logistic regression, all corrected for multiple testing across every
   qualifying GCF — not a pre-selected subset).

## Repository structure

```
rawdata/                          genome download script and sample metadata
analysis/
├── funcscan/                     nf-core/funcscan run config (params.yaml, custom.config, samplesheet)
├── scripts/gcf_clustering/       GCF clustering pipeline (Python) + results/
├── R/
│   ├── copied_from_funcscan_results/   comBGC summary table + MultiQC report
│   └── RMarkdown/                      analysis notebooks (.qmd, Quarto) and rendered output
```

Large intermediate files (genome FASTAs, the original MultiQC/Excel exports,
Pfam-A.hmm, funcscan work/results directories) are excluded via `.gitignore`
and are not part of this repository.

## Reproducing the analysis

```bash
# 1. Genome mining (requires Nextflow + Apptainer)
cd analysis/funcscan && bash nf-core_funcscan.sh

# 2. GCF clustering (requires conda)
cd analysis/scripts/gcf_clustering
bash setup_env.sh && conda activate gcf_clustering
bash run_clustering.sh

# 3. Statistical analysis (requires Quarto CLI)
cd analysis/R/RMarkdown && quarto render
```

## Status

Preliminary results were submitted as a short paper to the Brazilian
Symposium on Bioinformatics (BSB). BGC/GCF richness was modestly but
significantly higher in cancer genomes overall; no individual GCF reached
significance after correction for the full set of GCFs tested.
