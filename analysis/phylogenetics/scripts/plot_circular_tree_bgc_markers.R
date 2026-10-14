#!/usr/bin/env Rscript
# plot_circular_tree_bgc_markers.R -- variant of plot_circular_tree.R that adds
# two extra rings marking presence/absence of GCF_0017 and GCF_0019, the two
# BGCs found exclusive to one clade in the root-split comparison. Writes to a
# SEPARATE output file; does not touch full_hybrid_521.circular_tree.pdf.
#
# Run with ~/miniconda3/envs/quarto-env/bin/Rscript (see memory: mint-toolchain).

suppressMessages({
  library(ggtree)
  library(treeio)
  library(ggnewscale)
  library(tidyverse)
  library(ape)
})

ROOT <- file.path("..", "..", "..")
TREE <- file.path(ROOT, "analysis/phylogenetics/results/gubbins/full_hybrid_521/full_hybrid_521.final_tree.tre")
METADATA <- file.path(ROOT, "rawdata/metadata/gbk_hp_genomes.tsv")
GCF_ASSIGNMENTS <- file.path(ROOT, "analysis/scripts/results_521/gcf_assignments.tsv")
OUT <- file.path(ROOT, "analysis/phylogenetics/results/gubbins/full_hybrid_521/full_hybrid_521.circular_tree_bgc_markers.pdf")

CANCER_COLOR <- "#d7191c"
NONCANCER_COLOR <- "#2c7bb6"
REFERENCE_COLOR <- "#2CA02C"
CONTINENT_COLORS <- c(
  "Asia" = "#E69F00",
  "Europe" = "#009E73",
  "North America" = "#8B4CB8",
  "South America" = "#4D4D4D"
)
GCF0017_COLOR <- "#1A3C8C"
GCF0019_COLOR <- "#C2185B"
ABSENT_COLOR <- "grey92"

tree <- ape::read.tree(TREE)
cat(sprintf("Tips (including snippy-core Reference): %d\n", length(tree$tip.label)))

meta <- read_tsv(METADATA, show_col_types = FALSE) |>
  rename(id = Genome_ID) |>
  mutate(phenotype = if_else(str_trim(Phenotype) == "Gastric cancer", "Cancer", "Non-cancer"),
         continent = str_trim(Continent)) |>
  select(id, phenotype, continent)

gcf <- read_tsv(GCF_ASSIGNMENTS, show_col_types = FALSE)
richness <- gcf |> count(sample_id, name = "n_bgcs") |> rename(id = sample_id)

carriers_0017 <- gcf |> filter(gcf_id == "GCF_0017") |> distinct(sample_id) |> pull(sample_id)
carriers_0019 <- gcf |> filter(gcf_id == "GCF_0019") |> distinct(sample_id) |> pull(sample_id)
cat(sprintf("GCF_0017 carriers: %d, GCF_0019 carriers: %d\n", length(carriers_0017), length(carriers_0019)))

dat <- tibble(id = tree$tip.label) |>
  left_join(meta, by = "id") |>
  left_join(richness, by = "id") |>
  mutate(
    is_reference = id == "Reference",
    phenotype = if_else(is_reference, "Reference", phenotype),
    n_bgcs = if_else(is_reference, NA_real_, replace_na(n_bgcs, 0)),
    has_0017 = if_else(is_reference, NA_character_,
                        if_else(id %in% carriers_0017, "Present", "Absent")),
    has_0019 = if_else(is_reference, NA_character_,
                        if_else(id %in% carriers_0019, "Present", "Absent"))
  ) |>
  select(-is_reference)

stopifnot(sum(is.na(dat$phenotype)) == 0)
stopifnot(sum(is.na(dat$continent) & dat$id != "Reference") == 0)
cat("Phenotype counts:\n"); print(table(dat$phenotype))
cat("Continent counts:\n"); print(table(dat$continent))
cat("GCF_0017 presence by clade-irrelevant tally:\n"); print(table(dat$has_0017, useNA = "ifany"))
cat("GCF_0019 presence by clade-irrelevant tally:\n"); print(table(dat$has_0019, useNA = "ifany"))
cat(sprintf("BGC richness range: %d - %d (median %d)\n",
            min(dat$n_bgcs, na.rm = TRUE), max(dat$n_bgcs, na.rm = TRUE), median(dat$n_bgcs, na.rm = TRUE)))

# geom_fruit's proportional ring-width math (pwidth) breaks silently on this
# tree -- likely due to its very unevenly distributed branch lengths (two deep
# sublineages) -- producing an invisible, zero-height ring. gheatmap() (older,
# part of core ggtree) uses absolute offsets/proportional widths computed
# directly off the tree's own depth and renders reliably here; verified against
# this exact tree file before use. Each new ring needs its own
# ggnewscale::new_scale_fill() before the next gheatmap() call.
max_depth <- max(node.depth.edgelength(tree))

# Phenotype goes on the tree itself (see plot_circular_tree.R for the
# na.translate=FALSE pitfall this avoids).
p <- ggtree(tree, layout = "circular", size = 0.3, open.angle = 8) %<+% dat +
  aes(color = phenotype) +
  scale_color_manual(values = c(Cancer = CANCER_COLOR, `Non-cancer` = NONCANCER_COLOR,
                                 Reference = REFERENCE_COLOR),
                      breaks = c("Cancer", "Non-cancer", "Reference"),
                      na.value = "grey40",
                      name = "Phenotype", guide = guide_legend(order = 1))

continent_df <- dat |> select(id, continent) |> column_to_rownames("id")
phenotype_df <- dat |> select(id, phenotype) |> column_to_rownames("id")
gcf0017_df <- dat |> select(id, has_0017) |> column_to_rownames("id")
gcf0019_df <- dat |> select(id, has_0019) |> column_to_rownames("id")
richness_df <- dat |> select(id, n_bgcs) |> rename(`BGCs per genome` = n_bgcs) |> column_to_rownames("id")

gap <- max_depth * 0.006
ring_width <- 0.08
marker_width <- ring_width * 2 / 3
pheno_ring_width <- marker_width * 0.5

ring1a_offset <- max_depth * 0.03                                           # phenotype (thin, innermost)
ring1_offset <- ring1a_offset + pheno_ring_width * max_depth + gap          # continent
ring2_offset <- ring1_offset + marker_width * max_depth + gap               # GCF_0017
ring3_offset <- ring2_offset + marker_width * max_depth + gap               # GCF_0019
ring4_offset <- ring3_offset + marker_width * max_depth + gap               # BGC richness

# Thin phenotype ring -- same colors/categories already shown via branch color,
# so its own legend is suppressed (guide="none") to avoid a duplicate
# "Phenotype" box; the existing branch-color legend already explains it.
p1a <- gheatmap(p, phenotype_df, offset = ring1a_offset, width = pheno_ring_width, colnames = FALSE) +
  scale_fill_manual(values = c(Cancer = CANCER_COLOR, `Non-cancer` = NONCANCER_COLOR,
                                Reference = REFERENCE_COLOR),
                     breaks = c("Cancer", "Non-cancer", "Reference"), guide = "none")

p1b <- p1a + new_scale_fill()

p1c <- gheatmap(p1b, continent_df, offset = ring1_offset, width = marker_width, colnames = FALSE) +
  scale_fill_manual(values = CONTINENT_COLORS, breaks = names(CONTINENT_COLORS),
                     name = "Continent", guide = guide_legend(order = 2))

p2 <- p1c + new_scale_fill()

p3 <- gheatmap(p2, gcf0017_df, offset = ring2_offset, width = marker_width, colnames = FALSE) +
  scale_fill_manual(values = c(Present = GCF0017_COLOR, Absent = ABSENT_COLOR),
                     breaks = "Present", name = "GCF_0017", guide = guide_legend(order = 3))

p4 <- p3 + new_scale_fill()

p5 <- gheatmap(p4, gcf0019_df, offset = ring3_offset, width = marker_width, colnames = FALSE) +
  scale_fill_manual(values = c(Present = GCF0019_COLOR, Absent = ABSENT_COLOR),
                     breaks = "Present", name = "GCF_0019", guide = guide_legend(order = 4))

p6 <- p5 + new_scale_fill()

p7 <- gheatmap(p6, richness_df, offset = ring4_offset, width = ring_width, colnames = FALSE) +
  scale_fill_gradient(low = "#fee5d9", high = "#99000d", name = "BGCs per genome",
                       guide = guide_colorbar(order = 5)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 12)
  ) +
  ggtitle(sprintf(
    "Gubbins core-genome tree -- 521-genome cohort (%d cancer, %d non-cancer)\nwith GCF_0017 / GCF_0019 carrier rings",
    sum(dat$phenotype == "Cancer"), sum(dat$phenotype == "Non-cancer")
  ))

p <- p7

ggsave(OUT, p, width = 14, height = 14, dpi = 200, limitsize = FALSE)
cat(sprintf("Wrote %s\n", OUT))
