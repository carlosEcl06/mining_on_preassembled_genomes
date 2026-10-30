#!/usr/bin/env Rscript
# regen_poster_figures.R -- the BSB 2026 poster's figures (written to figures/), 521-genome
# cohort. Same data and ggplot code as regen_figures_521.R, at poster physical and font sizes.
#
# Run from analysis/R/figures/:
#   Rscript regen_poster_figures.R

suppressMessages(library(tidyverse))

ROOT <- file.path("..", "..", "..")
RESULTS_DIR_521 <- file.path(ROOT, "analysis", "scripts", "results_521")
METADATA         <- file.path(ROOT, "rawdata", "metadata", "gbk_hp_genomes.tsv")
DEDUP_EXCLUSIONS <- file.path(ROOT, "analysis", "cohort", "dedup_cohort_exclusions.tsv")
SHORT_GENOME_EXCLUSIONS <- c("GCA_024409305.1", "GCA_015904695.1")

BASE <- 24  # pt at print size; the poster body text is 28 pt
# Slot sizes (cm) of the figures added after the first draft
POSTER_H_ROW      <- 11.0
POSTER_W_GCF_RICH <- 9.4
POSTER_W_AGE      <- 15.0
POSTER_W_VOLCANO  <- 24.4
POSTER_H_VOLCANO  <- 10.6

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

# --- poster_gcf_size_dist.png ---
gcf_sizes <- gcf_meta |> filter(!is_singleton) |> count(gcf_id, name = "size")

p1 <- gcf_sizes |>
  filter(size <= 50) |>
  ggplot(aes(x = size)) +
  geom_histogram(binwidth = 1, fill = "#2c7bb6", color = "white") +
  labs(x = "GCF size (BGCs)", y = "Number of GCFs") +
  theme_minimal(base_size = BASE)
ggsave("figures/poster_gcf_size_dist.png", p1, width = 25, height = 9, units = "cm", dpi = 300,
       bg = "white")

# --- poster_bgc_richness.png ---
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
           label = "***", size = 9) +
  # W sits over the Cancer box (its whiskers top out well below the bracket), clear of the stars
  annotate("text", x = 1, y = richness_ymax * 0.85,
           label = sprintf("W = %.0f", wt$statistic), size = 6.5) +
  labs(x = NULL, y = "BGCs per genome", fill = NULL) +
  theme_minimal(base_size = BASE) +
  theme(legend.position = "none")
ggsave("figures/poster_bgc_richness.png", p2, width = 10.4, height = 11.5, units = "cm",
       dpi = 300, bg = "white")

# --- poster_bgc_richness_by_tool.png ---
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

tool_names <- c(antiSMASH = "antiSMASH", deepBGC = "DeepBGC")

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
            inherit.aes = FALSE, size = 7.5) +
  scale_x_discrete(labels = c("Cancer" = "Cancer", "Non-cancer" = "Non-\ncancer")) +
  facet_wrap(~ Prediction_tool, scales = "free_y",
             labeller = as_labeller(\(x) coalesce(tool_names[x], x))) +
  labs(x = NULL, y = "BGCs per genome", fill = NULL) +
  theme_minimal(base_size = BASE) +
  theme(legend.position = "none", strip.text = element_text(size = BASE, face = "plain"))
ggsave("figures/poster_bgc_richness_by_tool.png", p3, width = 14.8, height = 11.5, units = "cm",
       dpi = 300, bg = "white")

cat(sprintf("Wilcoxon cross-check (paper: W = 36204.0, p = 1.11e-04): W=%.1f p=%.3e\n",
            wt$statistic, wt$p.value))
cat("Wrote figures/poster_gcf_size_dist.png, figures/poster_bgc_richness.png,",
    "figures/poster_bgc_richness_by_tool.png\n")

# =============================================================================
# Figures beyond the paper's three, drawn from the same tests the paper reports
# (logic copied from analysis/R/gcf_results_521_comparison.R; the printed numbers are
# cross-checks against the paper).
# =============================================================================
GROUP_COLS <- c("Cancer" = "#d7191c", "Non-cancer" = "#2c7bb6")

# --- poster_gcf_richness.png: GCFs per genome (paper: W = 35636.5, p = 4.13e-04) ---
gcf_richness <- gcf_meta |> filter(!is_singleton & !is.na(group)) |>
  distinct(sample_id, gcf_id, group) |> count(sample_id, group, name = "n_gcfs")
wt2 <- wilcox.test(n_gcfs ~ group, data = gcf_richness)
cat(sprintf("GCF richness: W=%.1f p=%.3e (paper: W = 35636.5, p = 4.13e-04)\n",
            wt2$statistic, wt2$p.value))

g_ymax <- max(gcf_richness$n_gcfs)
p4 <- gcf_richness |>
  ggplot(aes(x = group, y = n_gcfs, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.15, height = 0.15, alpha = 0.3, size = 1.0) +
  scale_fill_manual(values = GROUP_COLS) +
  annotate("segment", x = 1, xend = 2, y = g_ymax * 1.1, yend = g_ymax * 1.1) +
  annotate("text", x = 1.5, y = g_ymax * 1.17, label = "***", size = 9) +
  scale_x_discrete(labels = c("Cancer" = "Cancer", "Non-cancer" = "Non-\ncancer")) +
  labs(x = NULL, y = "GCFs per genome") +
  theme_minimal(base_size = BASE) +
  theme(legend.position = "none")
ggsave("figures/poster_gcf_richness.png", p4, width = POSTER_W_GCF_RICH, height = POSTER_H_ROW,
       units = "cm", dpi = 300, bg = "white")

# --- poster_gcf_volcano.png: Fisher test per GCF (paper: 82 tested, min p_adj = 0.130) ---
n_cancer    <- sum(meta$group == "Cancer")
n_noncancer <- sum(meta$group == "Non-cancer")
presence <- gcf_meta |> filter(!is_singleton & !is.na(group)) |> distinct(sample_id, gcf_id, group)
qualifying_gcfs <- presence |> count(gcf_id) |> filter(n >= 3) |> pull(gcf_id)
gcf_stats <- map_dfr(qualifying_gcfs, function(g) {
  sub <- presence |> filter(gcf_id == g)
  a <- sum(sub$group == "Cancer"); b <- sum(sub$group == "Non-cancer")
  ft <- fisher.test(matrix(c(a, b, n_cancer - a, n_noncancer - b), nrow = 2))
  tibble(gcf_id = g, carriers = a + b, odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}) |> mutate(p_adj = p.adjust(p_value, method = "BH"))
cat(sprintf("Volcano: %d GCFs tested, %d with p_adj < 0.05, min p_adj = %.3f (paper: 82, 0, 0.130)\n",
            nrow(gcf_stats), sum(gcf_stats$p_adj < 0.05), min(gcf_stats$p_adj)))

# Fisher's conditional MLE is 0 or Inf when a GCF occurs in only one group; those points are
# drawn at the plot edge (open triangles) rather than dropped.
LOG2_CAP <- 4
volc <- gcf_stats |>
  mutate(log2_or = log2(odds_ratio),
         capped = !is.finite(log2_or) | abs(log2_or) > LOG2_CAP,
         log2_or = pmax(pmin(log2_or, LOG2_CAP), -LOG2_CAP),
         side = if_else(log2_or > 0, "Cancer", "Non-cancer"))
p5 <- volc |>
  ggplot(aes(x = log2_or, y = -log10(p_adj))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(size = carriers, colour = side, shape = capped), alpha = 0.75, stroke = 1.2) +
  scale_colour_manual(values = GROUP_COLS, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2), guide = "none") +
  scale_size_area(max_size = 9, breaks = c(5, 50, 200), name = "Carrier\ngenomes") +
  annotate("text", x = -LOG2_CAP, y = -log10(0.05) + 0.07, hjust = 0, vjust = 0,
           label = "p[adj] == 0.05", parse = TRUE, size = 6.5, colour = "grey30") +
  scale_y_continuous(limits = c(0, 1.6)) +
  labs(x = expression(log[2]~"odds ratio (cancer vs. non-cancer)"),
       y = expression(-log[10]~p[adj])) +
  theme_minimal(base_size = BASE) +
  theme(legend.position = "right")
ggsave("figures/poster_gcf_volcano.png", p5, width = POSTER_W_VOLCANO, height = POSTER_H_VOLCANO,
       units = "cm", dpi = 300, bg = "white")

# --- poster_age_richness.png: quasi-Poisson richness ~ phenotype + age
#     (paper: IRR cancer 1.146, IRR age 1.004/yr, interaction LRT p = 0.640) ---
age_richness <- richness |> left_join(meta |> select(sample_id, Age), by = "sample_id") |>
  filter(!is.na(Age)) |> mutate(Age = as.numeric(Age),
                                 group = factor(group, levels = c("Non-cancer", "Cancer")))
fit_main <- glm(n_bgcs ~ group + Age, data = age_richness, family = quasipoisson())
fit_int  <- glm(n_bgcs ~ group * Age, data = age_richness, family = quasipoisson())
irr <- exp(coef(fit_main))
cat(sprintf("Quasi-Poisson: IRR cancer = %.3f, IRR age = %.4f, LRT p = %.3f (n = %d with age)\n",
            irr["groupCancer"], irr["Age"], anova(fit_main, fit_int, test = "Chisq")$`Pr(>Chi)`[2],
            nrow(age_richness)))

grid <- expand_grid(Age = seq(min(age_richness$Age), max(age_richness$Age), length.out = 100),
                    group = factor(c("Non-cancer", "Cancer"), levels = c("Non-cancer", "Cancer")))
pr <- predict(fit_main, grid, type = "link", se.fit = TRUE)
grid <- grid |> mutate(fit = exp(pr$fit), lo = exp(pr$fit - 1.96 * pr$se.fit),
                       hi = exp(pr$fit + 1.96 * pr$se.fit))
p6 <- ggplot(age_richness, aes(x = Age, y = n_bgcs, colour = group)) +
  geom_jitter(width = 0.3, height = 0.2, alpha = 0.35, size = 1.3) +
  geom_ribbon(data = grid, aes(x = Age, ymin = lo, ymax = hi, fill = group), inherit.aes = FALSE,
              alpha = 0.2) +
  geom_line(data = grid, aes(x = Age, y = fit, colour = group), linewidth = 1.6) +
  scale_colour_manual(values = GROUP_COLS, name = NULL) +
  scale_fill_manual(values = GROUP_COLS, name = NULL) +
  coord_cartesian(ylim = c(0, 12)) +
  labs(x = "Age (years)", y = "BGCs per genome") +
  theme_minimal(base_size = BASE) +
  theme(legend.position = "top")
ggsave("figures/poster_age_richness.png", p6, width = POSTER_W_AGE, height = POSTER_H_ROW,
       units = "cm", dpi = 300, bg = "white")
cat("Wrote figures/poster_gcf_richness.png, figures/poster_gcf_volcano.png,",
    "figures/poster_age_richness.png\n")
