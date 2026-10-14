#!/usr/bin/env Rscript
# plot_circular_tree.R -- circular ggtree of the final Gubbins tree (521-genome
# cohort). Phenotype (cancer/non-cancer) is painted directly on the tree's
# terminal branches; two outer rings add geographic origin (continent) and
# BGC repertoire (BGC count per genome, from the 521-cohort GCF clustering
# results).
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
OUT <- file.path(ROOT, "analysis/phylogenetics/results/gubbins/full_hybrid_521/full_hybrid_521.circular_tree.pdf")

CANCER_COLOR <- "#d7191c"
NONCANCER_COLOR <- "#2c7bb6"
REFERENCE_COLOR <- "#2CA02C"
CONTINENT_COLORS <- c(
  "Asia" = "#E69F00",
  "Europe" = "#009E73",
  "North America" = "#8B4CB8",
  "South America" = "#4D4D4D"
)

tree <- ape::read.tree(TREE)
cat(sprintf("Tips (including snippy-core Reference): %d\n", length(tree$tip.label)))

meta <- read_tsv(METADATA, show_col_types = FALSE) |>
  rename(id = Genome_ID) |>
  mutate(phenotype = if_else(str_trim(Phenotype) == "Gastric cancer", "Cancer", "Non-cancer"),
         continent = str_trim(Continent)) |>
  select(id, phenotype, continent)

gcf <- read_tsv(GCF_ASSIGNMENTS, show_col_types = FALSE)
richness <- gcf |> count(sample_id, name = "n_bgcs") |> rename(id = sample_id)

dat <- tibble(id = tree$tip.label) |>
  left_join(meta, by = "id") |>
  left_join(richness, by = "id") |>
  mutate(
    is_reference = id == "Reference",
    phenotype = if_else(is_reference, "Reference", phenotype),
    n_bgcs = if_else(is_reference, NA_real_, replace_na(n_bgcs, 0))
  ) |>
  select(-is_reference)

stopifnot(sum(is.na(dat$phenotype)) == 0)
stopifnot(sum(is.na(dat$continent) & dat$id != "Reference") == 0)
cat("Phenotype counts:\n"); print(table(dat$phenotype))
cat("Continent counts:\n"); print(table(dat$continent))
cat(sprintf("BGC richness range: %d - %d (median %d)\n",
            min(dat$n_bgcs, na.rm = TRUE), max(dat$n_bgcs, na.rm = TRUE), median(dat$n_bgcs, na.rm = TRUE)))

# geom_fruit's proportional ring-width math (pwidth) breaks silently on this
# tree -- likely due to its very unevenly distributed branch lengths (two deep
# sublineages) -- producing an invisible, zero-height ring. gheatmap() (older,
# part of core ggtree) uses absolute offsets/proportional widths computed
# directly off the tree's own depth and renders reliably here; verified against
# this exact tree file before use. Two independent fill legends need
# ggnewscale::new_scale_fill() between the two gheatmap() calls.
max_depth <- max(node.depth.edgelength(tree))

# Phenotype goes on the tree itself: %<+% attaches tip metadata, then aes(color=)
# remaps the tree geom's color to it. Internal branches have no phenotype (only
# tips do), so they fall back to na.value -- a neutral grey -- leaving color on
# just the terminal (tip) branches. IMPORTANT: na.translate=FALSE (tried first)
# does not just hide the NA legend key -- it drops every NA-colored geom from
# rendering, silently deleting all internal branches. Use breaks= instead, which
# only restricts the legend, not what's drawn. The snippy-core Reference tip
# gets its own category/color instead of being dropped from the tree.
p <- ggtree(tree, layout = "circular", size = 0.3, open.angle = 8) %<+% dat +
  aes(color = phenotype) +
  scale_color_manual(values = c(Cancer = CANCER_COLOR, `Non-cancer` = NONCANCER_COLOR,
                                 Reference = REFERENCE_COLOR),
                      breaks = c("Cancer", "Non-cancer", "Reference"),
                      na.value = "grey40",
                      name = "Phenotype", guide = guide_legend(order = 1))

continent_df <- dat |> select(id, continent) |> column_to_rownames("id")
richness_df <- dat |> select(id, n_bgcs) |> rename(`BGCs per genome` = n_bgcs) |> column_to_rownames("id")

gap <- max_depth * 0.006
ring_width <- 0.08
continent_width <- ring_width * 2 / 3

ring1_offset <- max_depth * 0.03
ring2_offset <- ring1_offset + continent_width * max_depth + gap

p1 <- gheatmap(p, continent_df, offset = ring1_offset, width = continent_width, colnames = FALSE) +
  scale_fill_manual(values = CONTINENT_COLORS, breaks = names(CONTINENT_COLORS),
                     name = "Continent", guide = guide_legend(order = 2))

p2 <- p1 + new_scale_fill()

p3 <- gheatmap(p2, richness_df, offset = ring2_offset, width = ring_width, colnames = FALSE) +
  scale_fill_gradient(low = "#fee5d9", high = "#99000d", name = "BGCs per genome",
                       guide = guide_colorbar(order = 3)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 12)
  ) +
  ggtitle(sprintf(
    "Gubbins core-genome tree -- 521-genome cohort (%d cancer, %d non-cancer)",
    sum(dat$phenotype == "Cancer"), sum(dat$phenotype == "Non-cancer")
  ))

p <- p3

ggsave(OUT, p, width = 14, height = 14, dpi = 200, limitsize = FALSE)
cat(sprintf("Wrote %s\n", OUT))
