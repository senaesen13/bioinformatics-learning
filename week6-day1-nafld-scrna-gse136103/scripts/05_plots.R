## Week 6 Day 1 — GSE136103: key visualisation plots
## Loads the annotated Seurat object; no QC/clustering/DE re-run.

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(RColorBrewer)

out <- "../plots"

message("Loading annotated Seurat object...")
obj <- readRDS("../results/nafld_seurat_annotated.rds")
message("Loaded: ", ncol(obj), " cells, ", length(unique(obj$cell_type)), " cell types")

# ── Colour palette: 20 distinct cell types ────────────────────────────────────
ct_levels <- c(
  "CD4+ T cells", "CD8+ T cells", "CD8+ T cells (exhausted)", "Naive T cells",
  "NK/NKT cells", "NK cells (cytotoxic)", "NK cells (liver-resident)", "gdT/NK cells",
  "Monocytes", "Dendritic cells", "Kupffer cells",
  "B cells", "Plasma cells",
  "Endothelial cells", "LSEC",
  "Hepatic stellate cells", "Hepatocytes",
  "Proliferating cells", "Cholangiocytes", "Mast cells"
)
obj$cell_type <- factor(obj$cell_type, levels = ct_levels)
Idents(obj) <- "cell_type"

# 20-colour palette: group immune warm, stromal cool, parenchymal neutral
pal <- c(
  # T cells — blues
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  # NK — greens
  "#1A9641", "#4DAC26", "#B8E186", "#F1B6DA",
  # Myeloid — oranges/reds
  "#E66101", "#D73027", "#A50026",
  # B lineage — purples
  "#762A83", "#C51B7D",
  # Endothelial — teals
  "#35978F", "#80CDC1",
  # Stromal/parenchymal — browns/golds
  "#8C510A", "#DFC27D",
  # Other
  "#878787", "#4D4D4D", "#000000"
)
names(pal) <- ct_levels

# ── 1. UMAP: annotated cell types ─────────────────────────────────────────────
message("Plot 1: annotated UMAP")

p1 <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
              cols = pal, pt.size = 0.25, label = TRUE,
              label.size = 3, repel = TRUE) +
  ggtitle("Cell types — GSE136103 NAFLD subset (35,050 cells)") +
  theme_void(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 13, face = "bold"),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave(file.path(out, "umap_annotated.png"), p1,
       width = 13, height = 7, dpi = 180)
message("  Saved: umap_annotated.png")

# ── 2. UMAP: healthy vs NAFLD side-by-side ────────────────────────────────────
message("Plot 2: healthy vs NAFLD split UMAP")

obj$disease_label <- ifelse(obj$disease_status == "healthy",
                            "Healthy (n=5 donors)",
                            "NAFLD cirrhosis (n=2 donors)")

p2 <- DimPlot(obj, reduction = "umap", group.by = "cell_type",
              split.by = "disease_label", cols = pal,
              pt.size = 0.2, label = FALSE) +
  ggtitle("Cell distribution: healthy vs NAFLD cirrhosis") +
  theme_void(base_size = 11) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 13, face = "bold"),
    strip.text   = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 7.5),
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave(file.path(out, "umap_healthy_vs_nafld.png"), p2,
       width = 16, height = 7, dpi = 180)
message("  Saved: umap_healthy_vs_nafld.png")

# ── 3. FeaturePlot: TREM2 + CD9 ───────────────────────────────────────────────
message("Plot 3: FeaturePlot TREM2 + CD9")

fp3 <- FeaturePlot(obj, features = c("TREM2", "CD9"),
                   reduction = "umap", ncol = 2,
                   pt.size = 0.15, order = TRUE,
                   max.cutoff = "q95",
                   cols = c("lightgrey", "#D73027")) &
  theme_void(base_size = 11) &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

fp3 <- fp3 + plot_annotation(
  title = "TREM2 and CD9 expression on UMAP",
  theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
)

ggsave(file.path(out, "featureplot_trem2_cd9.png"), fp3,
       width = 13, height = 6, dpi = 180)
message("  Saved: featureplot_trem2_cd9.png")

# ── 4. FeaturePlot: COL1A1 + GPNMB ───────────────────────────────────────────
message("Plot 4: FeaturePlot COL1A1 + GPNMB")

fp4 <- FeaturePlot(obj, features = c("COL1A1", "GPNMB"),
                   reduction = "umap", ncol = 2,
                   pt.size = 0.15, order = TRUE,
                   max.cutoff = "q95",
                   cols = c("lightgrey", "#1A9641")) &
  theme_void(base_size = 11) &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

fp4 <- fp4 + plot_annotation(
  title = "COL1A1 and GPNMB expression on UMAP",
  theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
)

ggsave(file.path(out, "featureplot_col1a1_gpnmb.png"), fp4,
       width = 13, height = 6, dpi = 180)
message("  Saved: featureplot_col1a1_gpnmb.png")

# ── 5. Composition bar plot ───────────────────────────────────────────────────
message("Plot 5: composition bar plot")

comp <- read.csv("../results/cell_composition_major.csv")

# Pivot to long for ggplot
comp_long <- comp %>%
  select(major_cell_type, healthy, NAFLD_cirrhosis) %>%
  pivot_longer(c(healthy, NAFLD_cirrhosis),
               names_to = "disease", values_to = "mean_pct") %>%
  mutate(
    disease = recode(disease,
                     "healthy"         = "Healthy (n=5)",
                     "NAFLD_cirrhosis" = "NAFLD cirrhosis (n=2)"),
    major_cell_type = factor(major_cell_type,
                             levels = comp %>%
                               arrange(healthy) %>%
                               pull(major_cell_type))
  )

p5 <- ggplot(comp_long,
             aes(x = major_cell_type, y = mean_pct, fill = disease)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c("Healthy (n=5)"         = "#4393C3",
                               "NAFLD cirrhosis (n=2)" = "#D73027")) +
  coord_flip() +
  labs(x = NULL, y = "Mean % of total cells per donor",
       fill = NULL,
       title = "Cell type proportions: healthy vs NAFLD cirrhosis",
       caption = "Mean across donors within each group. Error bars omitted at n=2.") +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "top",
    plot.title       = element_text(face = "bold", size = 13),
    plot.caption     = element_text(size = 8, colour = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out, "composition_barplot.png"), p5,
       width = 10, height = 7, dpi = 180)
message("  Saved: composition_barplot.png")

# ── 6. Dot plot: TREM2, SPP1, GPNMB, COL1A1 across cell types ────────────────
message("Plot 6: dot plot key NAFLD genes")

# Cell type order: parenchymal → stromal → myeloid → lymphoid (bottom to top)
ct_dot_order <- c(
  "Mast cells", "Cholangiocytes", "Proliferating cells",
  "Hepatocytes", "Hepatic stellate cells", "LSEC", "Endothelial cells",
  "Plasma cells", "B cells",
  "Kupffer cells", "Dendritic cells", "Monocytes",
  "gdT/NK cells", "NK cells (liver-resident)", "NK cells (cytotoxic)", "NK/NKT cells",
  "Naive T cells", "CD8+ T cells (exhausted)", "CD8+ T cells", "CD4+ T cells"
)

Idents(obj) <- factor(obj$cell_type, levels = ct_dot_order)

p6 <- DotPlot(obj,
              features  = c("TREM2", "SPP1", "GPNMB", "COL1A1"),
              cols      = c("lightgrey", "#A50026"),
              dot.scale = 7,
              col.min   = 0) +
  coord_flip() +
  labs(x = NULL, y = NULL,
       title = "Cell-type specificity of key NAFLD genes",
       subtitle = "Dot size = % cells expressing; colour = average expression") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y  = element_text(size = 11, face = "bold"),
    plot.title   = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    legend.position = "right"
  )

ggsave(file.path(out, "dotplot_key_genes.png"), p6,
       width = 12, height = 5, dpi = 180)
message("  Saved: dotplot_key_genes.png")

# ── Done ──────────────────────────────────────────────────────────────────────
message("\n=== All 6 plots saved to plots/ ===")
message("  umap_annotated.png")
message("  umap_healthy_vs_nafld.png")
message("  featureplot_trem2_cd9.png")
message("  featureplot_col1a1_gpnmb.png")
message("  composition_barplot.png")
message("  dotplot_key_genes.png")
