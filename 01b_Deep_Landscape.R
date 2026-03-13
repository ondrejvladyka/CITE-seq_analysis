# ==============================================================================
# SCRIPT: 01b_Deep_Landscape.R
# PURPOSE: Extended exploratory visualizations for cell states and clonality
# ==============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

message(">>> Generating Deep Landscape Plots...")

# (Assumes 'sobj' is already loaded from 01_Landscape_Analysis.R)

# --- ANGLE 1: Reverse Composition (Who owns each cluster?) ---
rev_comp_data <- as.data.frame(table(Cluster = sobj$seurat_clusters, Group = sobj$Clinical_Group)) %>%
  group_by(Cluster) %>% mutate(Proportion = Freq / sum(Freq)) %>% ungroup()

p_rev_comp <- ggplot(rev_comp_data, aes(x = Cluster, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(title = paste(TARGET_CELL_TYPE, "- Clinical Group Makeup per Cluster"), 
       subtitle = "Which cohorts dominate which biological states?",
       x = "Seurat Cluster", y = "Percentage of Cluster") +
  theme_minimal() + theme(axis.text.x = element_text(size = 10, face = "bold"))

ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Reverse_Composition.png"), plot = p_rev_comp, width = 8, height = 5, dpi = 300)


# --- ANGLE 2: Canonical Marker DotPlot ---
# (Assumes CANONICAL_GENES is already loaded from Master Script)

# Only plot genes that actually exist in the current object
genes_to_plot <- intersect(CANONICAL_GENES, rownames(sobj))

p_dot <- DotPlot(sobj, features = genes_to_plot, group.by = "seurat_clusters") + 
  RotatedAxis() + 
  scale_color_gradient2(low = HEATMAP_LOW, 
                        mid = HEATMAP_MID, 
                        high = HEATMAP_HIGH, 
                        midpoint = 0) +
  
  ggtitle(paste(TARGET_CELL_TYPE, "- Canonical State Markers")) +
  theme(axis.text.x = element_text(size = 10, face = "bold")) 

# Save using the Master Script DPI settings
ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Canonical_DotPlot.png"), 
       plot = p_dot, width = 12, height = 6, dpi = PUB_DPI)


# --- ANGLE 3: Clonal Homeostasis (Stacked Barplot) ---
if("cloneSize" %in% colnames(sobj@meta.data)) {
  # Calculate proportion of clone sizes per clinical group
  clone_comp <- as.data.frame(table(Size = sobj$cloneSize, Group = sobj$Clinical_Group)) %>%
    group_by(Group) %>% mutate(Proportion = Freq / sum(Freq)) %>% ungroup()
  
  # Ensure logical ordering (Single -> Small -> Medium -> Large)
  clone_levels <- c("Single (1)", "Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)", "Hyperexpanded (0.1 < X <= 1)")
  clone_comp$Size <- factor(clone_comp$Size, levels = intersect(clone_levels, unique(clone_comp$Size)))
  
  p_clone_bar <- ggplot(clone_comp, aes(x = Group, y = Proportion, fill = Size)) +
    geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.2) +
    scale_fill_brewer(palette = "Reds") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Clonal Expansion Space", y = "Proportion of Repertoire") +
    theme_minimal()
  
  ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Clonal_Homeostasis.png"), plot = p_clone_bar, width = 6, height = 6, dpi = 300)
}


# --- ANGLE 4: Key Feature UMAPs ---
# Let's map Naive (SELL), Cytotoxic (GNLY), and Stress/Migration (CXCR4)
features_umap <- intersect((KEY_FEATURES), rownames(sobj))
if(length(features_umap) > 0) {
  p_feat <- FeaturePlot(sobj, features = features_umap, ncol = 3, pt.size = 0.5) & 
    scale_color_viridis_c(option = "magma", direction = -1)
  
  ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Key_Features_UMAP.png"), plot = p_feat, width = 12, height = 4, dpi = 300)
}

message(">>> Deep Landscape plots saved successfully to ", OUTPUT_DIR)

# --- ANGLE 5: Cluster Phylogenetic Tree (Phenotypic Similarity) ---

# Calculate hierarchical clustering based on average expression in PCA space
reduc_for_tree <- ifelse("glmpca" %in% names(sobj@reductions), "glmpca", "pca")

message(">>> Generating Cluster Dendrogram using ", reduc_for_tree, "...")
sobj <- BuildClusterTree(sobj, dims = 1:20, reduction = reduc_for_tree)

# Save the tree plot
png(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Cluster_Dendrogram.png"), width = 1200, height = 800, res = PUB_DPI/2)
PlotClusterTree(sobj, edge.width = 2, font.size = 14)
title(main = paste(TARGET_CELL_TYPE, "- Cluster Phenotypic Similarity Tree"))
dev.off()

# --- ANGLE 6: Clinical Density Topography ---
message(">>> Generating Clinical Density Contour Maps...")

# --- BULLETPROOF UMAP EXTRACTION ---
# 1. Pull the raw matrix directly from the dimensional reduction slot
umap_raw <- Embeddings(sobj, reduction = "umap")
umap_coords <- as.data.frame(umap_raw)

# 2. Force the columns to be named exactly UMAP_1 and UMAP_2 so ggplot never fails
colnames(umap_coords)[1:2] <- c("UMAP_1", "UMAP_2")

# 3. Attach the clinical metadata
umap_coords$Clinical_Group <- sobj$Clinical_Group

p_density <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2)) +
  # Draw topographic density lines
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", color = "white", linewidth = 0.1, bins = 12) +
  facet_wrap(~Clinical_Group) +
  scale_fill_viridis_c(option = "magma", name = "Cell Density") +
  theme_minimal() +
  labs(title = paste(TARGET_CELL_TYPE, "- Topographical Density of Cohorts"),
       subtitle = "Visualizing the 'center of gravity' for each clinical group") +
  theme(panel.background = element_rect(fill = "grey10"), 
        panel.grid = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text = element_blank()) 

ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Clinical_Density_Topography.png"), plot = p_density, width = 12, height = 4, dpi = PUB_DPI)