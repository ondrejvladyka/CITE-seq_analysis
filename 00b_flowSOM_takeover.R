# ==============================================================================
# SCRIPT: 00b_flowSOM_takeover.R
# PURPOSE: Recluster using FlowSOM, evaluate K, and hijack Seurat identities
# ==============================================================================
library(FlowSOM)
library(ggplot2)

message("\n>>> HOSTILE TAKEOVER: Running FlowSOM over GLM-PCA...")

if(!exists("sobj")) {
  stop(">>> ERROR: 'sobj' not found. Ensure 00_Master_Control.R loaded the object.")
}

# 1. Extract the math
reduc_use <- ifelse("glmpca" %in% names(sobj@reductions), "glmpca", "pca")
pca_matrix <- Embeddings(sobj, reduc_use)[, 1:20]

# 2. Build FlowSOM Object
set.seed(42)
fSOM <- ReadInput(pca_matrix)
fSOM <- BuildSOM(fSOM, colsToUse = 1:20)
fSOM <- BuildMST(fSOM)

# ---------------------------------------------------------
# OUTPUT 1: ELBOW PLOT GENERATION
# ---------------------------------------------------------
message(">>> Calculating WCSS for FlowSOM Elbow Plot...")
som_codes <- fSOM$map$codes
max_k <- 61
wcss <- sapply(1:max_k, function(k) { kmeans(som_codes, centers = k, nstart = 61)$tot.withinss })

elbow_data <- data.frame(k = 1:max_k, WCSS = wcss)
p_elbow <- ggplot(elbow_data, aes(x = k, y = WCSS)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkred", size = 3) +
  scale_x_continuous(breaks = seq(1, max_k, by = 1)) +
  labs(title = paste(TARGET_CELL_TYPE, "- FlowSOM Metacluster Optimization"),
       subtitle = "Look for the 'Elbow' where variance levels off",
       x = "Number of Metaclusters (k)", y = "Total Within-Cluster Sum of Squares (WCSS)") +
  theme_minimal()

ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_FlowSOM_Elbow_Plot.png"), plot = p_elbow, width = 8, height = 6, dpi = PUB_DPI)

# ---------------------------------------------------------
# OUTPUT 2: ANNOTATED TREE PLOT (Biological Cell Types)
# ---------------------------------------------------------
message(">>> Saving Annotated FlowSOM Tree to PDF...")
pdf(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_FlowSOM_Tree_Annotated.pdf"), width = 15, height = 12)

biological_labels <- as.factor(sobj$starCAT_label)
p_tree_bio <- PlotPies(fSOM, 
                       cellTypes = biological_labels, 
                       title = paste("FlowSOM Tree:", TARGET_CELL_TYPE, "Biological Makeup (starCAT)"))
print(p_tree_bio)

dev.off()

# ---------------------------------------------------------
# HOSTILE TAKEOVER (In RAM)
# ---------------------------------------------------------
# Determine optimal 'k' (Uses Seurat's count unless overridden in Master Script)
if(!exists("FLOWSOM_K")) {
  num_clusters <- length(unique(sobj$seurat_clusters))
  message(">>> Defaulting to Seurat's optimal cluster count: k = ", num_clusters)
} else {
  num_clusters <- FLOWSOM_K
  message(">>> Using manually specified FlowSOM metacluster count: k = ", num_clusters)
}

# Generate metaclusters
meta_clusters <- metaClustering_consensus(fSOM$map$codes, k = num_clusters)
cell_assignments <- meta_clusters[fSOM$map$mapping[, 1]]

# Hijack the Seurat identity in the global environment!
sobj$FlowSOM_clusters <- factor(cell_assignments)
sobj$seurat_clusters <- factor(cell_assignments)
Idents(sobj) <- sobj$seurat_clusters

message(">>> Takeover complete! Downstream scripts will now use FlowSOM clusters.\n")