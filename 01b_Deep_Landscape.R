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
  scale_fill_manual(values = c("PRE" = "#66c2a5", "POST" = "#fc8d62", "22q11" = "#8da0cb")) +
  labs(title = paste(TARGET_CELL_TYPE, "- Clinical Group Makeup per Cluster"), 
       subtitle = "Which cohorts dominate which biological states?",
       x = "Seurat Cluster", y = "Percentage of Cluster") +
  theme_minimal() + theme(axis.text.x = element_text(size = 10, face = "bold"))

ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Reverse_Composition.png"), plot = p_rev_comp, width = 8, height = 5, dpi = 300)


# --- ANGLE 2: Canonical Marker DotPlot ---
# Grouped by biological function for a clean visual narrative
canonical_genes <- c(
  # Naive / Memory / Survival
  "CCR7", "SELL", "TCF7", "LEF1", "IL7R", "CD27", "CD28", 
  
  # Effector / Cytotoxicity
  "GZMB", "PRF1", "NKG7", "GNLY", 
  
  # Migration / Virtual Memory / Stress
  "CXCR4", "S100A4", "ANXA1", 
  
  # Exhaustion / Dysfunction / Senescence
  "PDCD1",   # PD-1
  "HAVCR2",  # TIM-3
  "TIGIT",   # TIGIT
  "TOX",     # TOX
  "LAG3",    # LAG-3
  "ENTPD1",  # CD39
  "KLRG1",   # Senescence
  "B3GAT1",  # CD57 proxy
  
  # Proliferation
  "MKI67", "TOP2A"
)

# Only plot genes that actually exist in the current object
genes_to_plot <- intersect(canonical_genes, rownames(sobj))

p_dot <- DotPlot(sobj, features = genes_to_plot, group.by = "seurat_clusters") + 
  RotatedAxis() + 
  scale_color_distiller(palette = "RdYlBu") +
  ggtitle(paste(TARGET_CELL_TYPE, "- Canonical State Markers")) +
  theme(axis.text.x = element_text(size = 10, face = "bold")) # Makes the gene names pop

ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Canonical_DotPlot.png"), plot = p_dot, width = 12, height = 6, dpi = 300)

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
features_umap <- intersect(c("SELL", "GNLY", "CXCR4"), rownames(sobj))
if(length(features_umap) > 0) {
  p_feat <- FeaturePlot(sobj, features = features_umap, ncol = 3, pt.size = 0.5) & 
    scale_color_viridis_c(option = "magma", direction = -1)
  
  ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Key_Features_UMAP.png"), plot = p_feat, width = 12, height = 4, dpi = 300)
}

message(">>> Deep Landscape plots saved successfully to ", OUTPUT_DIR)