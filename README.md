# Mining Biosynthetic Gene Clusters in *H. pylori* Genomes

Code and data for the short paper and poster presented at the Brazilian
Symposium on Bioinformatics (BSB 2026): genome mining of biosynthetic gene
clusters (BGCs) in preassembled *Helicobacter pylori* genomes from
gastric-cancer and precancerous samples, followed by tool-agnostic Gene Cluster
Family (GCF) clustering and testing for association with clinical phenotype.

Of the 554 genomes retrieved from NCBI, 33 were excluded (2 assemblies under
1.4 Mb and 31 near-duplicates, which violate the independence assumption of the
tests), leaving **521 genomes** (173 cancer, 348 non-cancer). The exclusion list
is in [`analysis/cohort/excluded_33_genomes.tsv`](analysis/cohort/excluded_33_genomes.tsv).

This is a Master's project at the Laboratory of Molecular Genetics and Omics
(LGMO), Ribeirão Preto Medical School – University of São Paulo (USP).

**Authors:** Carlos Eugênio Costa de Lima, Wilson Araújo da Silva Junior

**Notice:** This README has been mostly AI-generated, thus it may contain some inconsistencies.

## Pipeline overview

1. **Genome retrieval** (`rawdata/`) — preassembled genomes downloaded from NCBI
   via `datasets download genome`.
2. **BGC mining** (`analysis/funcscan/`) — [nf-core/funcscan](https://nf-co.re/funcscan)
   (v3.0.0), running antiSMASH 8.0.1, DeepBGC 0.1.31, and GECCO 0.9.10, with
   Bakta/InterProScan functional annotation. comBGC aggregates the results into
   a table of Pfam domains per BGC.
3. **Cohort definition** (`analysis/cohort/`) — near-duplicate pairs are grouped
   into connected components and one representative per component is kept
   (highest Snippy callable fraction).
4. **GCF clustering** (`analysis/scripts/gcf_clustering/`) — a tool-agnostic
   workflow, since BiG-SCAPE requires antiSMASH-formatted input. Pairwise BGC
   similarity is adapted from Robey et al. (2021):

   - Pfam domain arrays are compared via Jaccard similarity (pairs below 0.10
     excluded).
   - For BGCs with a verified backbone domain (antiSMASH: Terpene-precursor,
     Azole-containing-RiPP; DeepBGC: NRP, Saccharide, Polyketide, Terpene,
     Polyketide-Terpene), sequence identity between backbone domains is
     estimated via `hmmalign` against Pfam-A, and the full score is used:
     `sqrt(0.8 * identity + 0.2 * Jaccard)`.
   - All other BGCs fall back to `sqrt(0.2 * Jaccard)`.
   - Final GCFs are delineated with DBSCAN (`epsilon = 0.56`,
     `min_samples = 2`, chosen via an epsilon-sensitivity sweep).

5. **Statistical analysis** (`analysis/R/`) — Wilcoxon tests of BGC and GCF
   richness, Fisher's exact test, Cochran–Mantel–Haenszel, multivariable
   logistic regression, and quasi-Poisson regression of richness on phenotype
   and age, all corrected for multiple testing across every qualifying GCF.

## Repository structure

```
rawdata/                              genome download script and sample metadata
analysis/
├── funcscan/                         nf-core/funcscan run configuration
├── cohort/                           exclusion list, 521-genome ID list, and the script behind them
├── scripts/
│   ├── gcf_clustering/               GCF clustering pipeline
│   └── results_521/                  clustering outputs for the 521-genome cohort
└── R/
    ├── copied_from_funcscan_results/ comBGC summary table and software versions
    ├── gcf_results_521_comparison.R  statistical tests
    └── figures/                      paper and poster figures
```

Genome FASTAs, Pfam-A.hmm and funcscan work/results directories are not part of
this repository.

## Reproducing the analysis

The statistics and figures run directly from the tables in this repository
(requires R with `tidyverse` and `broom`):

```bash
# from the repository root
Rscript analysis/R/gcf_results_521_comparison.R

# from analysis/R/figures/
Rscript regen_figures_521.R
Rscript regen_poster_figures.R
```

Reproducing the full pipeline requires downloading the genomes and databases:

```bash
# 1. Genomes (requires NCBI datasets)
cd rawdata && bash download_genomes.sh

# 2. Genome mining (requires Nextflow + Apptainer)
cd analysis/funcscan && bash nf-core_funcscan.sh

# 3. GCF clustering on the 521-genome cohort (requires conda and HMMER)
cd analysis/scripts/gcf_clustering
bash setup_env.sh && conda activate gcf_clustering
PFAM_HMM=/path/to/Pfam-A.hmm bash recluster_521.sh
```

The sample metadata in `rawdata/metadata/` comes from public NCBI records.

## License

Code is released under the [MIT License](LICENSE).
