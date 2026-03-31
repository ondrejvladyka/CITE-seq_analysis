# ==============================================================================
# SCRIPT: 01c_Alternative_Visualizations.R
# PURPOSE: Compare UMAP and t-SNE using GLM-PCA embeddings
# ==============================================================================
library(Seurat)
library(ggplot2)
library(patchwork)

message("\n>>> Running Alternative Visualizations (UMAP vs t-SNE)...")

if(!exists("sobj")) stop(">>> ERROR: 'sobj' not found.")

reduc_use <- ifelse("glmpca" %in% names(sobj@reductions), "glmpca", "pca")

message(">>> Calculating t-SNE (This may take a minute)...")
sobj <- RunTSNE(sobj, reduction = reduc_use, dims = 1:20, check_duplicates = FALSE)

p_umap <- DimPlot(sobj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("UMAP (Focus on Global Structure)") + 
  theme_minimal() + NoLegend()

p_tsne <- DimPlot(sobj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("t-SNE (Focus on Local Density)") + 
  theme_minimal() + NoLegend()

umap_vs_tsne <- p_umap | p_tsne
output_file <- paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_UMAP_vs_tSNE.png")

ggsave(output_file, plot = umap_vs_tsne, width = 12, height = 6, dpi = PUB_DPI)
message(">>> Alternative visualizations saved to ", output_file)