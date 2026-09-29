#!/usr/bin/env Rscript
# gcf_results_521_comparison.R
#
# Re-runs the exact statistical-testing sequence from gcf_results_full.qmd (Clinical
# context section) on the deduplicated 521-genome cohort, for the BSB short paper's
# camera-ready reconsideration (2026-09-28 discussion). Prints comparable summary
# numbers to what's already written in main.tex's Section 4 (554-genome cohort), so
# the two can be diffed directly.
#
# IMPORTANT correctness point (why this isn't just "point the qmd at results_521/"):
# gcf_results_full.qmd computes n_cancer/n_noncancer and builds its CMH/GLM covariate
# tables from the FULL `meta` (554 genomes), not from whichever genomes happen to be
# in gcf_assignments.tsv. Swapping only the GCF input files would leave those
# denominators and contingency/stratification tables still counting the 33 excluded
# genomes (2 <1.4Mb + 31 near-duplicate/same-patient pairs) -- exactly the
# independence-assumption problem this whole exercise exists to fix. So `meta` itself
# is filtered to the 521-genome cohort here, FIRST, before any test is run.
#
# Run with quarto-env R (matches gcf_results_full.qmd's toolchain):
#   ~/miniconda3/envs/quarto-env/bin/Rscript analysis/R/RMarkdown/gcf_results_521_comparison.R

suppressMessages({
  library(tidyverse)
  library(broom)
})

# --- paths (relative to repo root; run from there) ---
RESULTS_DIR_521 <- file.path("analysis", "scripts", "results_521")
METADATA        <- file.path("rawdata", "metadata", "gbk_hp_genomes.tsv")
DEDUP_EXCLUSIONS <- file.path("analysis", "phylogenetics", "results", "grafgen",
                               "dedup_cohort_exclusions.tsv")
SHORT_GENOME_EXCLUSIONS <- c("GCA_024409305.1", "GCA_015904695.1")  # <1.4Mb, see nextflow/main.nf

required_files <- c(
  file.path(RESULTS_DIR_521, "gcf_assignments.tsv"),
  file.path(RESULTS_DIR_521, "bgc_domain_arrays.tsv"),
  METADATA, DEDUP_EXCLUSIONS
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(sprintf("Required data file(s) not found: %s", paste(missing_files, collapse = ", ")))
}

# --- build the 521-genome cohort list (same logic as the phylogenomics side) ---
near_dup_excluded <- read_tsv(DEDUP_EXCLUSIONS, show_col_types = FALSE) |>
  filter(role == "excluded (near-duplicate)") |>
  pull(sample_id)
excluded_521 <- c(SHORT_GENOME_EXCLUSIONS, near_dup_excluded)
stopifnot(length(excluded_521) == 33)

# --- load data, restricted to the 521-genome cohort from the start ---
gcf     <- read_tsv(file.path(RESULTS_DIR_521, "gcf_assignments.tsv"),  show_col_types = FALSE)
domains <- read_tsv(file.path(RESULTS_DIR_521, "bgc_domain_arrays.tsv"), show_col_types = FALSE)
meta_full <- read_tsv(METADATA, show_col_types = FALSE) |>
  rename(sample_id = Genome_ID) |>
  mutate(Phenotype = str_trim(Phenotype),
         group     = if_else(Phenotype == "Gastric cancer", "Cancer", "Non-cancer"))

meta <- meta_full |> filter(!sample_id %in% excluded_521)
stopifnot(nrow(meta) == 521)
stopifnot(all(!gcf$sample_id %in% excluded_521))  # sanity: the re-cluster shouldn't contain any excluded genome

cat(sprintf("Cohort: %d genomes (%d Cancer, %d Non-cancer) -- excluded %d (2 short + 31 near-duplicate)\n\n",
            nrow(meta), sum(meta$group == "Cancer"), sum(meta$group == "Non-cancer"), length(excluded_521)))

gcf_meta <- gcf |> left_join(meta, by = "sample_id")
is_singleton <- gcf_meta$gcf_id == "singleton"

# =============================================================================
# BGC richness per genome (Wilcoxon) -- paper's Section 4, para 1
# =============================================================================
richness <- gcf_meta |> filter(!is.na(group)) |> count(sample_id, group, name = "n_bgcs")
wt <- wilcox.test(n_bgcs ~ group, data = richness)
cat("=== BGC richness per genome (Wilcoxon) ===\n")
cat(sprintf("Median: Cancer = %g, Non-cancer = %g\n",
            median(richness$n_bgcs[richness$group == "Cancer"]),
            median(richness$n_bgcs[richness$group == "Non-cancer"])))
cat(sprintf("W = %.1f, p = %.3e\n", wt$statistic, wt$p.value))
cat(sprintf("[554-genome paper value for comparison: median 4 vs 3, W=40459.5, p=5.96e-05]\n\n"))

# --- by prediction tool ---
cat("=== BGC richness by prediction tool (Wilcoxon) ===\n")
richness_tool <- gcf_meta |> filter(!is.na(group)) |> count(sample_id, group, Prediction_tool, name = "n_bgcs")
tool_results <- map_dfr(unique(gcf_meta$Prediction_tool), function(tool) {
  sub <- richness_tool |> filter(Prediction_tool == tool)
  if (n_distinct(sub$group) < 2) return(tibble(Tool = tool, W = NA_real_, p_value = NA_real_))
  wt_t <- wilcox.test(n_bgcs ~ group, data = sub)
  tibble(Tool = tool, W = wt_t$statistic, p_value = wt_t$p.value)
}) |> mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
print(tool_results, n = Inf)
cat("[554-genome paper: antiSMASH W=33898.5 p_adj=0.387 ns; deepBGC W=37563.0 p_adj=2.61e-05 ***]\n\n")

# =============================================================================
# GCF richness per genome (Wilcoxon) -- paper's Section 4, para 3
# =============================================================================
gcf_richness <- gcf_meta |> filter(!is_singleton & !is.na(group)) |>
  distinct(sample_id, gcf_id, group) |> count(sample_id, group, name = "n_gcfs")
wt2 <- wilcox.test(n_gcfs ~ group, data = gcf_richness)
cat("=== GCF richness per genome (Wilcoxon) ===\n")
cat(sprintf("Median: Cancer = %g, Non-cancer = %g\n",
            median(gcf_richness$n_gcfs[gcf_richness$group == "Cancer"]),
            median(gcf_richness$n_gcfs[gcf_richness$group == "Non-cancer"])))
cat(sprintf("W = %.1f, p = %.3e\n", wt2$statistic, wt2$p.value))
cat("[554-genome paper value: both medians 3, W=39892, p=2.10e-04]\n\n")

# =============================================================================
# GCF-level Fisher / CMH / logistic regression -- paper's Section 4, para 4
# =============================================================================
MIN_GENOMES_TOTAL <- 3
n_cancer    <- sum(meta$group == "Cancer",     na.rm = TRUE)
n_noncancer <- sum(meta$group == "Non-cancer", na.rm = TRUE)

presence <- gcf_meta |> filter(!is_singleton & !is.na(group)) |> distinct(sample_id, gcf_id, group)
gcf_carrier_counts <- presence |> count(gcf_id, name = "n_carrier_genomes")
qualifying_gcfs <- gcf_carrier_counts |> filter(n_carrier_genomes >= MIN_GENOMES_TOTAL) |> pull(gcf_id)
n_gcf_total      <- n_distinct(presence$gcf_id)
n_gcf_qualifying <- length(qualifying_gcfs)

fisher_per_gcf <- function(gcf_id_val) {
  sub <- presence |> filter(gcf_id == gcf_id_val)
  a <- sum(sub$group == "Cancer"); b <- sum(sub$group == "Non-cancer")
  c <- n_cancer - a; d <- n_noncancer - b
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
  tibble(gcf_id = gcf_id_val, n_cancer = a, n_noncancer = b, odds_ratio = ft$estimate, p_value = ft$p.value)
}
gcf_stats <- map_dfr(qualifying_gcfs, fisher_per_gcf) |>
  mutate(p_adj = p.adjust(p_value, method = "BH")) |> arrange(p_adj)

cat("=== GCF-level testing ===\n")
cat(sprintf("Non-singleton GCFs total: %d; qualifying (>= %d carriers): %d\n",
            n_gcf_total, MIN_GENOMES_TOTAL, n_gcf_qualifying))
cat(sprintf("Fisher: %d with p_adj<0.05, %d with p_adj<0.10, min p_adj = %.3e\n",
            sum(gcf_stats$p_adj < 0.05), sum(gcf_stats$p_adj < 0.10), min(gcf_stats$p_adj)))
cat("[554-genome paper: 85 qualifying GCFs, 0 with p_adj<0.05 or <0.10, min p_adj=0.189]\n\n")

# --- CMH + GLM ---
covars <- meta |> select(sample_id, Continent, Age) |> mutate(Age = as.numeric(Age))
group_lookup <- meta |> select(sample_id, group)
richness_per_genome <- gcf_meta |> count(sample_id, name = "n_bgcs_total")
presence_wide <- presence |> left_join(covars, by = "sample_id") |>
  left_join(richness_per_genome, by = "sample_id") |> filter(!is.na(group), !is.na(Continent))
all_genomes <- group_lookup |> left_join(covars, by = "sample_id") |>
  left_join(richness_per_genome, by = "sample_id") |> filter(!is.na(group), !is.na(Continent))

cmh_per_gcf <- function(gcf_id_val) {
  carriers <- presence_wide |> filter(gcf_id == gcf_id_val) |> pull(sample_id)
  tab_df <- all_genomes |> mutate(carrier = sample_id %in% carriers) |> count(Continent, group, carrier)
  strata <- unique(tab_df$Continent); arr_list <- list()
  for (s in strata) {
    sub <- tab_df |> filter(Continent == s)
    m <- matrix(0, nrow = 2, ncol = 2, dimnames = list(carrier = c("TRUE","FALSE"), group = c("Cancer","Non-cancer")))
    for (i in seq_len(nrow(sub))) m[as.character(sub$carrier[i]), sub$group[i]] <- sub$n[i]
    if (all(rowSums(m) > 0) && all(colSums(m) > 0)) arr_list[[s]] <- m
  }
  if (length(arr_list) < 1) return(tibble(gcf_id = gcf_id_val, cmh_p = NA_real_))
  arr <- array(unlist(arr_list), dim = c(2, 2, length(arr_list)))
  res <- tryCatch(mantelhaen.test(arr, correct = TRUE), error = function(e) NULL)
  if (is.null(res)) return(tibble(gcf_id = gcf_id_val, cmh_p = NA_real_))
  tibble(gcf_id = gcf_id_val, cmh_p = res$p.value)
}
cmh_results <- map_dfr(gcf_stats$gcf_id, cmh_per_gcf)

glm_per_gcf <- function(gcf_id_val) {
  carriers <- presence_wide |> filter(gcf_id == gcf_id_val) |> pull(sample_id)
  dat <- all_genomes |> mutate(carrier = as.integer(sample_id %in% carriers),
                                group = factor(group, levels = c("Non-cancer", "Cancer")),
                                Continent = factor(Continent)) |> filter(!is.na(Age))
  if (n_distinct(dat$carrier) < 2 || nrow(dat) < 20) return(tibble(gcf_id = gcf_id_val, glm_p = NA_real_))
  fit <- tryCatch(glm(carrier ~ group + Continent + Age + n_bgcs_total, data = dat, family = binomial()),
                   error = function(e) NULL,
                   warning = function(w) suppressWarnings(glm(carrier ~ group + Continent + Age + n_bgcs_total, data = dat, family = binomial())))
  if (is.null(fit) || !"groupCancer" %in% names(coef(fit))) return(tibble(gcf_id = gcf_id_val, glm_p = NA_real_))
  tibble(gcf_id = gcf_id_val, glm_p = summary(fit)$coefficients["groupCancer", "Pr(>|z|)"])
}
glm_results <- map_dfr(gcf_stats$gcf_id, glm_per_gcf)

comparison <- gcf_stats |> left_join(cmh_results, by = "gcf_id") |> left_join(glm_results, by = "gcf_id") |>
  mutate(cmh_p_adj = p.adjust(cmh_p, method = "BH"), glm_p_adj = p.adjust(glm_p, method = "BH"))

cat(sprintf("CMH: %d testable, %d significant (p_adj<0.05)\n",
            sum(!is.na(comparison$cmh_p)), sum(comparison$cmh_p_adj < 0.05, na.rm = TRUE)))
cat(sprintf("GLM: %d testable, %d significant (p_adj<0.05)\n",
            sum(!is.na(comparison$glm_p)), sum(comparison$glm_p_adj < 0.05, na.rm = TRUE)))
cat("[554-genome paper: CMH 70 testable 0 sig; GLM 85 testable 0 sig]\n\n")

# =============================================================================
# Quasi-Poisson: richness ~ phenotype + age -- paper's Section 4, last paragraph
# =============================================================================
age_richness <- richness |> left_join(meta |> select(sample_id, Age), by = "sample_id") |>
  filter(!is.na(Age)) |> mutate(Age = as.numeric(Age))
mod_dat <- age_richness |> mutate(group = factor(group, levels = c("Non-cancer", "Cancer")))
fit_main <- glm(n_bgcs ~ group + Age, data = mod_dat, family = quasipoisson())
fit_int  <- glm(n_bgcs ~ group * Age, data = mod_dat, family = quasipoisson())
lrt <- anova(fit_main, fit_int, test = "Chisq")
summary_main <- broom::tidy(fit_main, conf.int = TRUE, exponentiate = TRUE)

irr_group <- summary_main |> filter(term == "groupCancer")
irr_age   <- summary_main |> filter(term == "Age")
cor_overall <- cor.test(age_richness$Age, age_richness$n_bgcs, method = "spearman")
cor_cancer  <- cor.test(age_richness$Age[age_richness$group == "Cancer"], age_richness$n_bgcs[age_richness$group == "Cancer"], method = "spearman")
cor_noncancer <- cor.test(age_richness$Age[age_richness$group == "Non-cancer"], age_richness$n_bgcs[age_richness$group == "Non-cancer"], method = "spearman")

cat("=== Quasi-Poisson: BGC richness ~ group + age ===\n")
cat(sprintf("IRR (Cancer vs Non-cancer) = %.3f (95%% CI %.3f-%.3f), p = %.4f\n",
            irr_group$estimate, irr_group$conf.low, irr_group$conf.high, irr_group$p.value))
cat(sprintf("IRR (per year age) = %.4f (95%% CI %.4f-%.4f), p = %.4e\n",
            irr_age$estimate, irr_age$conf.low, irr_age$conf.high, irr_age$p.value))
cat(sprintf("Interaction LRT p = %.3f\n", lrt$`Pr(>Chi)`[2]))
cat(sprintf("Spearman rho: overall = %.2f (p=%.3e), Cancer = %.2f (p=%.3f), Non-cancer = %.2f (p=%.4f)\n",
            cor_overall$estimate, cor_overall$p.value, cor_cancer$estimate, cor_cancer$p.value,
            cor_noncancer$estimate, cor_noncancer$p.value))
cat("[554-genome paper: IRR(age)=1.005 p=0.0047; IRR(cancer)=1.119 (CI 1.000-1.250) p=0.0487; interaction p=0.864;\n")
cat(" Spearman overall rho=0.16 p=1.14e-04, non-cancer rho=0.15 p=0.0048, cancer rho=0.08 p=0.305]\n")
