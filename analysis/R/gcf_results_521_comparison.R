#!/usr/bin/env Rscript
# gcf_results_521_comparison.R
#
# Statistical tests reported in the BSB 2026 short paper, on the deduplicated
# 521-genome cohort. The metadata is filtered to the cohort before any test is run,
# so denominators and contingency/stratification tables do not count the 33 excluded
# genomes (2 <1.4Mb + 31 near-duplicate/same-patient pairs).
#
# Run from the repository root:
#   Rscript analysis/R/gcf_results_521_comparison.R

suppressMessages({
  library(tidyverse)
  library(broom)
})

# --- paths (relative to repo root; run from there) ---
RESULTS_DIR_521 <- file.path("analysis", "scripts", "results_521")
METADATA        <- file.path("rawdata", "metadata", "gbk_hp_genomes.tsv")
DEDUP_EXCLUSIONS <- file.path("analysis", "cohort", "dedup_cohort_exclusions.tsv")
SHORT_GENOME_EXCLUSIONS <- c("GCA_024409305.1", "GCA_015904695.1")  # <1.4Mb

required_files <- c(
  file.path(RESULTS_DIR_521, "gcf_assignments.tsv"),
  file.path(RESULTS_DIR_521, "bgc_domain_arrays.tsv"),
  METADATA, DEDUP_EXCLUSIONS
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(sprintf("Required data file(s) not found: %s", paste(missing_files, collapse = ", ")))
}

# --- build the 521-genome cohort list ---
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
# BGC richness per genome (Wilcoxon)
# =============================================================================
richness <- gcf_meta |> filter(!is.na(group)) |> count(sample_id, group, name = "n_bgcs")
wt <- wilcox.test(n_bgcs ~ group, data = richness)
cat("=== BGC richness per genome (Wilcoxon) ===\n")
cat(sprintf("Median: Cancer = %g, Non-cancer = %g\n",
            median(richness$n_bgcs[richness$group == "Cancer"]),
            median(richness$n_bgcs[richness$group == "Non-cancer"])))
cat(sprintf("W = %.1f, p = %.3e\n", wt$statistic, wt$p.value))
cat("
")

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
cat("
")

# =============================================================================
# GCF richness per genome (Wilcoxon)
# =============================================================================
gcf_richness <- gcf_meta |> filter(!is_singleton & !is.na(group)) |>
  distinct(sample_id, gcf_id, group) |> count(sample_id, group, name = "n_gcfs")
wt2 <- wilcox.test(n_gcfs ~ group, data = gcf_richness)
cat("=== GCF richness per genome (Wilcoxon) ===\n")
cat(sprintf("Median: Cancer = %g, Non-cancer = %g\n",
            median(gcf_richness$n_gcfs[gcf_richness$group == "Cancer"]),
            median(gcf_richness$n_gcfs[gcf_richness$group == "Non-cancer"])))
cat(sprintf("W = %.1f, p = %.3e\n", wt2$statistic, wt2$p.value))
cat("
")

# =============================================================================
# GCF-level Fisher / CMH / logistic regression
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
cat("
")

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
cat("
")

# =============================================================================
# Quasi-Poisson: richness ~ phenotype + age
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
