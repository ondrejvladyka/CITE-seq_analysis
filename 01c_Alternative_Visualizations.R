# ==============================================================================
# SCRIPT: 01c_Alternative_Visualizations.R
# PURPOSE: Compare UMAP, t-SNE, and FlowSOM using GLM-PCA embeddings
# ==============================================================================
library(Seurat)
library(ggplot2)
library(patchwork)
library(qs2)
library(FlowSOM)

message("\n>>> Running Alternative Visualizations...")

# 1. Load the object (Assuming 00_Master_Control handled the alias, but we will play it safe)
if(!exists("sobj")) { sobj <- qs_read("data/CD4_TargetCohorts_GLMPCA.qs2") }

# Identify which math to use (Prefers glmpca if it exists, otherwise falls back to pca)
reduc_use <- ifelse("glmpca" %in% names(sobj@reductions), "glmpca", "pca")
message(">>> Using reduction: ", reduc_use)

# ---------------------------------------------------------
# VISUALIZATION 1 & 2: UMAP vs t-SNE
# ---------------------------------------------------------
message(">>> Calculating t-SNE (This may take a minute)...")
sobj <- RunTSNE(sobj, reduction = reduc_use, dims = 1:20, check_duplicates = FALSE)

# Build the comparison plot
p_umap <- DimPlot(sobj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("UMAP (Focus on Global Structure)") + 
  theme_minimal() + NoLegend()

p_tsne <- DimPlot(sobj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("t-SNE (Focus on Local Density)") + 
  theme_minimal() + NoLegend()

# Save the side-by-side comparison
umap_vs_tsne <- p_umap | p_tsne
ggsave("CD4_UMAP_vs_tSNE.png", plot = umap_vs_tsne, width = 12, height = 6, dpi = 300)
print(umap_vs_tsne)

# ---------------------------------------------------------
# VISUALIZATION 3: FlowSOM Minimum Spanning Tree
# ---------------------------------------------------------
message(">>> Building FlowSOM Minimum Spanning Tree...")

# 1. Extract the raw GLM-PCA/PCA matrix (Cells = Rows, PCs = Columns)
pca_matrix <- Embeddings(sobj, reduc_use)[, 1:20]

# 2. Build the FlowSOM object step-by-step (The Safest Method)
set.seed(42) # Lock the seed for reproducible trees
fSOM <- ReadInput(pca_matrix)
fSOM <- BuildSOM(fSOM, colsToUse = 1:20)
fSOM <- BuildMST(fSOM)

# 3. Plot the Tree with Pie Charts
message(">>> Saving FlowSOM Tree to PDF...")
pdf("CD4_FlowSOM_Tree.pdf", width = 12, height = 12)

# Extract Seurat clusters as a clean factor for the pie charts
seurat_labels <- as.factor(sobj$seurat_clusters)

# THE FIX: Changed 'main' to 'title', and saved it to an object
p_tree <- PlotPies(fSOM, 
                   cellTypes = seurat_labels, 
                   title = "FlowSOM Tree: CD4 GLM-PCA Clusters")

# Explicitly print the ggplot object into the PDF device
print(p_tree)

dev.off()

message(">>> All alternative visualizations successfully generated!")
