## Week 6 Day 1 — GSE136103: cell-type annotation, composition, within-cell-type DE
## Healthy (5 donors) vs NAFLD cirrhosis (Cirrhotic1 + Cirrhotic4, 2 donors)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ── Load object ───────────────────────────────────────────────────────────────
message("Loading Seurat object...")
obj <- readRDS("../results/nafld_seurat_clustered.rds")
message("Loaded: ", ncol(obj), " cells, clusters 0–20")

# ── 1. Cell-type annotation ───────────────────────────────────────────────────
# Based on top marker genes from FindAllMarkers (script 03)
cluster_to_celltype <- c(
  "0"  = "CD4+ T cells",
  "1"  = "NK/NKT cells",
  "2"  = "CD8+ T cells",
  "3"  = "CD4+ T cells",
  "4"  = "CD8+ T cells (exhausted)",
  "5"  = "NK cells (cytotoxic)",
  "6"  = "Monocytes",
  "7"  = "Endothelial cells",
  "8"  = "B cells",
  "9"  = "LSEC",
  "10" = "Dendritic cells",
  "11" = "Kupffer cells",
  "12" = "Hepatic stellate cells",
  "13" = "Hepatocytes",
  "14" = "gdT/NK cells",
  "15" = "Plasma cells",
  "16" = "Naive T cells",
  "17" = "Proliferating cells",
  "18" = "NK cells (liver-resident)",
  "19" = "Cholangiocytes",
  "20" = "Mast cells"
)

obj$cell_type <- unname(cluster_to_celltype[as.character(obj$seurat_clusters)])

# Broader grouping for composition plots
obj$major_cell_type <- dplyr::recode(obj$cell_type,
  "CD4+ T cells"             = "T cells",
  "CD8+ T cells"             = "T cells",
  "CD8+ T cells (exhausted)" = "T cells",
  "Naive T cells"            = "T cells",
  "NK/NKT cells"             = "NK/NKT cells",
  "NK cells (cytotoxic)"     = "NK/NKT cells",
  "NK cells (liver-resident)"= "NK/NKT cells",
  "gdT/NK cells"             = "NK/NKT cells",
  "B cells"                  = "B/Plasma cells",
  "Plasma cells"             = "B/Plasma cells",
  "Monocytes"                = "Monocytes",
  "Dendritic cells"          = "Dendritic cells",
  "Kupffer cells"            = "Kupffer cells",
  "LSEC"                     = "Endothelial/LSEC",
  "Endothelial cells"        = "Endothelial/LSEC",
  "Hepatic stellate cells"   = "Hepatic stellate cells",
  "Hepatocytes"              = "Hepatocytes",
  "Proliferating cells"      = "Proliferating cells",
  "Cholangiocytes"           = "Cholangiocytes",
  "Mast cells"               = "Mast cells"
)

Idents(obj) <- "cell_type"

# Labelled UMAP
p_umap_ct <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 3,
                     repel = TRUE, pt.size = 0.3) +
  ggtitle("UMAP — annotated cell types") +
  theme(legend.position = "right")
ggsave("../plots/umap_cell_types.png", p_umap_ct, width = 14, height = 8, dpi = 150)

saveRDS(obj, "../results/nafld_seurat_annotated.rds")
message("Annotated Seurat object saved.\n")

# ── 2. Cell composition ───────────────────────────────────────────────────────
message("=== Cell composition analysis ===\n")

# Per-donor cell counts per major cell type
comp <- obj@meta.data %>%
  count(major_cell_type, disease_status, donor) %>%
  group_by(donor) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

# Mean ± SD across donors, per disease group
comp_summary <- comp %>%
  group_by(major_cell_type, disease_status) %>%
  summarise(mean_pct = mean(pct), sd_pct = sd(pct), n_donors = n(),
            .groups = "drop")

# Wide table for reporting
comp_wide <- comp_summary %>%
  select(major_cell_type, disease_status, mean_pct) %>%
  pivot_wider(names_from = disease_status,
              values_from = mean_pct, values_fill = 0) %>%
  mutate(
    diff = NAFLD_cirrhosis - healthy,
    fold = ifelse(healthy > 0.1, NAFLD_cirrhosis / healthy, NA)
  ) %>%
  arrange(desc(abs(diff)))

write.csv(comp_wide, "../results/cell_composition_major.csv", row.names = FALSE)

message("Cell composition table (mean % per donor group):\n")
print(as.data.frame(comp_wide %>%
  mutate(across(where(is.numeric), ~ round(., 2)))), row.names = FALSE)

# Plot A: stacked bars per donor (fine-grained cell type)
ct_order <- obj@meta.data %>%
  count(cell_type, disease_status) %>%
  filter(disease_status == "healthy") %>%
  arrange(desc(n)) %>%
  pull(cell_type)

comp_fine <- obj@meta.data %>%
  count(cell_type, disease_status, donor) %>%
  group_by(donor) %>%
  mutate(pct = n / sum(n) * 100,
         disease_label = ifelse(disease_status == "healthy",
                                "Healthy", "NAFLD Cirrhosis")) %>%
  ungroup()
comp_fine$cell_type <- factor(comp_fine$cell_type, levels = rev(ct_order))

p_stack <- ggplot(comp_fine, aes(x = donor, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ disease_label, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "% of cells", fill = "Cell type",
       title = "Cell type composition per donor") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8))
ggsave("../plots/composition_stacked_by_donor.png", p_stack,
       width = 14, height = 7, dpi = 150)

# Plot B: mean % per major cell type, dodge bars with error
disease_colours <- c("healthy" = "#4DAF4A", "NAFLD_cirrhosis" = "#E41A1C")
major_order <- comp_summary %>%
  filter(disease_status == "healthy") %>%
  arrange(mean_pct) %>%
  pull(major_cell_type)

comp_summary$major_cell_type <- factor(comp_summary$major_cell_type,
                                        levels = major_order)

p_dodge <- ggplot(comp_summary,
                  aes(x = major_cell_type, y = mean_pct, fill = disease_status)) +
  geom_bar(stat = "identity", position = position_dodge(0.85), width = 0.8) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_pct - sd_pct), ymax = mean_pct + sd_pct),
    position = position_dodge(0.85), width = 0.3, linewidth = 0.5
  ) +
  coord_flip() +
  scale_fill_manual(values = disease_colours,
                    labels = c("healthy" = "Healthy (n=5)",
                               "NAFLD_cirrhosis" = "NAFLD cirrhosis (n=2)")) +
  labs(x = NULL, y = "Mean % of cells (± SD across donors)",
       fill = NULL,
       title = "Cell type proportions: healthy vs NAFLD cirrhosis",
       caption = "Error bars = SD across donors. NAFLD: only 2 donors — interpret cautiously.") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")
ggsave("../plots/composition_healthy_vs_nafld.png", p_dodge,
       width = 11, height = 8, dpi = 150)

# ── 3. Within-cell-type DE ────────────────────────────────────────────────────
message("\n=== Within-cell-type DE: NAFLD vs healthy ===\n")
message("NOTE: Only 2 NAFLD donors — results indicate direction/magnitude,")
message("but are not powered for formal discovery. Interpret as exploratory.\n")

MIN_CELLS <- 20   # minimum cells per group to attempt DE

# Cell types to test (those with adequate cells in both groups)
Idents(obj) <- "cell_type"
cell_type_counts <- obj@meta.data %>%
  count(cell_type, disease_status) %>%
  pivot_wider(names_from = disease_status, values_from = n, values_fill = 0)

message("Cells per cell type per disease group:\n")
print(as.data.frame(cell_type_counts), row.names = FALSE)

run_de <- cell_type_counts %>%
  filter(healthy >= MIN_CELLS, NAFLD_cirrhosis >= MIN_CELLS) %>%
  pull(cell_type)

message("\nCell types with >= ", MIN_CELLS, " cells per group: ",
        paste(run_de, collapse = ", "), "\n")

all_de <- list()

for (ct in run_de) {
  n_h <- cell_type_counts$healthy[cell_type_counts$cell_type == ct]
  n_n <- cell_type_counts$NAFLD_cirrhosis[cell_type_counts$cell_type == ct]
  message("  DE: ", ct, "  (healthy=", n_h, ", NAFLD=", n_n, ")")

  ct_obj <- subset(obj, idents = ct)
  Idents(ct_obj) <- "disease_status"

  de <- tryCatch(
    FindMarkers(ct_obj, ident.1 = "NAFLD_cirrhosis", ident.2 = "healthy",
                min.pct = 0.10, logfc.threshold = 0.10,
                test.use = "wilcox", verbose = FALSE),
    error = function(e) { message("    ERROR: ", e$message); NULL }
  )

  if (!is.null(de) && nrow(de) > 0) {
    de$gene      <- rownames(de)
    de$cell_type <- ct
    all_de[[ct]] <- de
  }
}

de_all <- bind_rows(all_de)
de_sig  <- de_all %>% filter(p_val_adj < 0.05)

write.csv(de_all,  "../results/within_celltype_de_all.csv",  row.names = FALSE)
write.csv(de_sig,  "../results/within_celltype_de_sig.csv",  row.names = FALSE)
message("\nTotal tested genes across cell types  : ", nrow(de_all))
message("Significant (padj < 0.05)             : ", nrow(de_sig))
message("Cell types with ≥1 sig gene            : ",
        length(unique(de_sig$cell_type)))

# Top 10 up / down per cell type (by avg_log2FC)
top_de <- de_sig %>%
  group_by(cell_type) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  bind_rows(
    de_sig %>%
      group_by(cell_type) %>%
      slice_min(order_by = avg_log2FC, n = 10)
  ) %>%
  distinct() %>%
  arrange(cell_type, desc(avg_log2FC))

write.csv(top_de, "../results/top10_de_per_celltype.csv", row.names = FALSE)

message("\nTop 5 upregulated in NAFLD per cell type:\n")
top5_up <- de_sig %>%
  group_by(cell_type) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  select(cell_type, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
print(as.data.frame(top5_up), row.names = FALSE)

message("\nTop 5 downregulated in NAFLD (i.e., higher in healthy) per cell type:\n")
top5_dn <- de_sig %>%
  group_by(cell_type) %>%
  slice_min(order_by = avg_log2FC, n = 5) %>%
  select(cell_type, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
print(as.data.frame(top5_dn), row.names = FALSE)

# ── 4. Dot plot: key NAFLD genes across cell types ───────────────────────────
nafld_genes <- c(
  # Kupffer cell / macrophage activation
  "TREM2", "GPNMB", "SPP1", "CD9", "LGALS3",
  # Fibrosis / HSC activation
  "ACTA2", "COL1A1", "COL1A2", "LOXL2", "PDGFRB",
  # Hepatocyte stress
  "HMGB1", "TXN", "MT1X", "CYP2E1",
  # Inflammation
  "IL1B", "TNF", "CCL2", "CXCL8"
)
# Keep only genes actually in the object
nafld_genes <- nafld_genes[nafld_genes %in% rownames(obj)]

Idents(obj) <- "cell_type"
p_dot <- DotPlot(obj, features = nafld_genes,
                 group.by = "cell_type", split.by = "disease_status",
                 cols = c("lightgrey", "#E41A1C"),
                 dot.scale = 5) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9)) +
  ggtitle("Key NAFLD-related genes across cell types (split: healthy | NAFLD)")
ggsave("../plots/nafld_genes_dotplot.png", p_dot,
       width = 16, height = 8, dpi = 150)

# ── Summary ───────────────────────────────────────────────────────────────────
message("\n========================================")
message("SUMMARY")
message("========================================")
message("Cell types annotated  : ", length(unique(obj$cell_type)))
message("Major categories      : ", length(unique(obj$major_cell_type)))
message("Healthy cells         : 27,632  (5 donors)")
message("NAFLD cells           : 7,418   (2 donors: Cirrhotic1, Cirrhotic4)")
message("Cell types DE-tested  : ", length(run_de))
message("Significant DE genes  : ", nrow(de_sig), " across ",
        length(unique(de_sig$cell_type)), " cell types")
message("========================================\n")
message("Output files:")
message("  results/nafld_seurat_annotated.rds")
message("  results/cell_composition_major.csv")
message("  results/within_celltype_de_all.csv")
message("  results/within_celltype_de_sig.csv")
message("  results/top10_de_per_celltype.csv")
message("  plots/umap_cell_types.png")
message("  plots/composition_stacked_by_donor.png")
message("  plots/composition_healthy_vs_nafld.png")
message("  plots/nafld_genes_dotplot.png")
