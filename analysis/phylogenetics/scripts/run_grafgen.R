#!/usr/bin/env Rscript
# Runs GrafGen on the snippy-core VCF, producing continuous ancestry
# proportions (F/E/A_percent) per genome for the GCF-phenotype confounder
# correction. AVX-limited server hardware can't run this (see plan's
# "Where AVX-dependent steps actually run") so it runs locally.
suppressMessages(library(GrafGen))

args <- commandArgs(trailingOnly = TRUE)
vcf_gz <- if (length(args) >= 1) args[[1]] else "analysis/phylogenetics/results/snippy/core/core.vcf.gz"
out_tsv <- if (length(args) >= 2) args[[2]] else "analysis/phylogenetics/results/grafgen/ancestry.tsv"

res <- grafGen(vcf_gz, print = 0)
write.table(res$table, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote %d genomes' ancestry calls to %s\n", nrow(res$table), out_tsv))
