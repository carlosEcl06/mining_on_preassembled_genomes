#!/usr/bin/env Rscript
# regen_figures_521.R -- regenerates the three short-paper figures (written to figures/) on the
# 521-genome cohort.
#
# Run from analysis/R/figures/:
#   Rscript regen_figures_521.R

suppressMessages(library(tidyverse))

ROOT <- file.path("..", "..", "..")
RESULTS_DIR_521 <- file.path(ROOT, "analysis", "scripts", "results_521")
METADATA         <- file.path(ROOT, "rawdata", "metadata", "gbk_hp_genomes.tsv")
DEDUP_EXCLUSIONS <- file.path(ROOT, "analysis", "cohort", "dedup_cohort_exclusions.tsv")
SHORT_GENOME_EXCLUSIONS <- c("GCA_024409305.1", "GCA_015904695.1")

near_dup_excluded <- read_tsv(DEDUP_EXCLUSIONS, show_col_types = FALSE) |>
  filter(role == "excluded (near-duplicate)") |> pull(sample_id)
excluded_521 <- c(SHORT_GENOME_EXCLUSIONS, near_dup_excluded)

gcf  <- read_tsv(file.path(RESULTS_DIR_521, "gcf_assignments.tsv"), show_col_types = FALSE)
meta <- read_tsv(METADATA, show_col_types = FALSE) |>
  rename(sample_id = Genome_ID) |>
  mutate(Phenotype = str_trim(Phenotype),
         group     = if_else(Phenotype == "Gastric cancer", "Cancer", "Non-cancer")) |>
  filter(!sample_id %in% excluded_521)
stopifnot(nrow(meta) == 521)

gcf_meta <- gcf |> left_join(meta, by = "sample_id")
is_singleton <- gcf_meta$gcf_id == "singleton"

dir.create("figures", showWarnings = FALSE)

# --- fig_gcf_size_dist.png ---
gcf_sizes <- gcf_meta |> filter(!is_singleton) |> count(gcf_id, name = "size")

p1 <- gcf_sizes |>
  filter(size <= 50) |>
  ggplot(aes(x = size)) +
  geom_histogram(binwidth = 1, fill = "#2c7bb6", color = "white") +
  labs(x = "GCF size (BGCs)", y = "Number of GCFs") +
  theme_minimal(base_size = 18)
ggsave("figures/fig_gcf_size_dist.png", p1, width = 8, height = 5, dpi = 192)

# --- fig_bgc_richness.png ---
richness <- gcf_meta |> filter(!is.na(group)) |> count(sample_id, group, name = "n_bgcs")
wt <- wilcox.test(n_bgcs ~ group, data = richness)

richness_ymax <- max(richness$n_bgcs)
bracket_y <- richness_ymax * 1.12
tick <- richness_ymax * 0.035

p2 <- richness |>
  ggplot(aes(x = group, y = n_bgcs, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.2) +
  scale_fill_manual(values = c("Cancer" = "#d7191c", "Non-cancer" = "#2c7bb6")) +
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y) +
  annotate("segment", x = 1, xend = 1, y = bracket_y - tick, yend = bracket_y) +
  annotate("segment", x = 2, xend = 2, y = bracket_y - tick, yend = bracket_y) +
  annotate("text", x = 1.5, y = bracket_y + richness_ymax * 0.05,
           label = "***", size = 6) +
  annotate("text", x = 0.62, y = bracket_y, hjust = 0,
           label = sprintf("W = %.0f", wt$statistic), size = 4.5) +
  labs(x = NULL, y = "BGCs per genome", fill = NULL) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")
ggsave("figures/fig_bgc_richness.png", p2, width = 8, height = 5, dpi = 192)

# --- fig_bgc_richness_by_tool.png ---
richness_tool <- gcf_meta |>
  filter(!is.na(group)) |>
  count(sample_id, group, Prediction_tool, name = "n_bgcs")

tool_ann <- richness_tool |>
  group_by(Prediction_tool) |>
  summarize(ymax = max(n_bgcs), .groups = "drop") |>
  mutate(
    bracket_y = ymax * 1.15,
    tick = ymax * 0.05,
    label = if_else(grepl("deep", tolower(Prediction_tool)), "***", "ns")
  )

p3 <- richness_tool |>
  ggplot(aes(x = group, y = n_bgcs, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.0) +
  scale_fill_manual(values = c("Cancer" = "#d7191c", "Non-cancer" = "#2c7bb6")) +
  geom_segment(data = tool_ann, aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = tool_ann, aes(x = 1, xend = 1, y = bracket_y - tick, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_segment(data = tool_ann, aes(x = 2, xend = 2, y = bracket_y - tick, yend = bracket_y),
               inherit.aes = FALSE) +
  geom_text(data = tool_ann, aes(x = 1.5, y = bracket_y + ymax * 0.07, label = label),
            inherit.aes = FALSE, size = 5) +
  facet_wrap(~ Prediction_tool, scales = "free_y") +
  labs(x = NULL, y = "BGCs per genome", fill = NULL) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none", strip.text = element_text(size = 18, face = "plain"))
ggsave("figures/fig_bgc_richness_by_tool.png", p3, width = 8, height = 5, dpi = 192)

cat(sprintf("Wilcoxon (for cross-check against gcf_results_521_comparison.R): W=%.1f p=%.3e\n",
            wt$statistic, wt$p.value))
cat("Wrote figures/fig_gcf_size_dist.png, figures/fig_bgc_richness.png, figures/fig_bgc_richness_by_tool.png\n")
