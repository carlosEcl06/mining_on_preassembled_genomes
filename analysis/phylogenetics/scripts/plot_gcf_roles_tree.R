#!/usr/bin/env Rscript
# plot_gcf_roles_tree.R -- circular tree colored by clade (A/B) with five
# presence/absence rings, one per clade-significant GCF (0017,0019,0025,0030,0031).
# Built for the gcf_roles_followup report. Run with:
#   PATH="$HOME/miniconda3/envs/quarto-env/bin:$PATH" Rscript plot_gcf_roles_tree.R

suppressMessages({
  library(ggtree)
  library(treeio)
  library(ggnewscale)
  library(tidyverse)
  library(ape)
})

ROOT <- file.path("..", "..", "..")
TREE <- file.path(ROOT, "analysis/phylogenetics/results/gubbins/full_hybrid_521/full_hybrid_521.final_tree.tre")
CLADE_ASSIGN <- file.path(ROOT, "analysis/phylogenetics/results/root_split_comparison/genome_clade_assignments.tsv")
GCF_ASSIGNMENTS <- file.path(ROOT, "analysis/scripts/results_521/gcf_assignments.tsv")
OUT <- file.path(ROOT, "analysis/phylogenetics/results/gcf_roles_followup/figs/fig_clade_gcf_tree.pdf")

CLADE_A_COLOR <- "#2a78d6"
CLADE_B_COLOR <- "#eb6834"
REFERENCE_COLOR <- "grey50"
ABSENT_COLOR <- "grey92"

GCF_COLORS <- c(
  GCF_0017 = "#1A3C8C",
  GCF_0019 = "#C2185B",
  GCF_0025 = "#1baf7a",
  GCF_0030 = "#eda100",
  GCF_0031 = "#7B2D8E"
)
TARGET_GCFS <- names(GCF_COLORS)

tree <- ape::read.tree(TREE)
cat(sprintf("Tips (including snippy-core Reference): %d\n", length(tree$tip.label)))

clade <- read_tsv(CLADE_ASSIGN, show_col_types = FALSE) |> rename(gid = id)
gcf <- read_tsv(GCF_ASSIGNMENTS, show_col_types = FALSE)

carrier_lists <- lapply(TARGET_GCFS, function(g) {
  gcf |> filter(gcf_id == g) |> distinct(sample_id) |> pull(sample_id)
})
names(carrier_lists) <- TARGET_GCFS
for (g in TARGET_GCFS) cat(sprintf("%s carriers: %d\n", g, length(carrier_lists[[g]])))

dat <- tibble(gid = tree$tip.label) |>
  left_join(clade, by = "gid") |>
  mutate(
    is_reference = gid == "Reference",
    clade = if_else(is_reference, "Reference", clade)
  )

for (g in TARGET_GCFS) {
  dat[[g]] <- if_else(dat$is_reference, NA_character_,
                       if_else(dat$gid %in% carrier_lists[[g]], "Present", "Absent"))
}
dat <- dat |> select(-is_reference)

stopifnot(sum(is.na(dat$clade)) == 0)
cat("Clade counts:\n"); print(table(dat$clade))

max_depth <- max(node.depth.edgelength(tree))

p <- ggtree(tree, layout = "circular", size = 0.3, open.angle = 8) %<+% dat +
  aes(color = clade) +
  scale_color_manual(values = c(A = CLADE_A_COLOR, B = CLADE_B_COLOR, Reference = REFERENCE_COLOR),
                      breaks = c("A", "B", "Reference"),
                      na.value = "grey40",
                      name = "Clade", guide = guide_legend(order = 1))

gap <- max_depth * 0.006
ring_width <- 0.08

offsets <- list()
prev_offset <- max_depth * 0.03
for (g in TARGET_GCFS) {
  offsets[[g]] <- prev_offset
  prev_offset <- prev_offset + ring_width * max_depth + gap
}

p_cur <- p
for (i in seq_along(TARGET_GCFS)) {
  g <- TARGET_GCFS[i]
  df <- dat |> select(gid, !!g) |> column_to_rownames("gid")
  p_cur <- p_cur + new_scale_fill()
  p_cur <- gheatmap(p_cur, df, offset = offsets[[g]], width = ring_width, colnames = FALSE) +
    scale_fill_manual(values = setNames(c(GCF_COLORS[[g]], ABSENT_COLOR), c("Present", "Absent")),
                       breaks = "Present", name = g, guide = guide_legend(order = 1 + i))
}

p_final <- p_cur +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 12)
  ) +
  ggtitle(sprintf(
    "521-genome tree, clade A/B split, with carriers of the five clade-significant GCFs\n(Clade A: %d genomes, Clade B: %d genomes)",
    sum(dat$clade == "A", na.rm = TRUE), sum(dat$clade == "B", na.rm = TRUE)
  ))

ggsave(OUT, p_final, width = 13, height = 13, dpi = 200, limitsize = FALSE)
cat(sprintf("Wrote %s\n", OUT))
