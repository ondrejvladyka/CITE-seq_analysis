# ==============================================================================
# SCRIPT: 01_Landscape_Analysis.R
# ==============================================================================
library(Seurat)
library(qs2)
library(ggplot2)
library(patchwork)
library(dplyr)

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

if("Clinical_Group" %in% colnames(sobj@meta.data)) {
  sobj$Clinical_Group <- factor(sobj$Clinical_Group, levels = c("PRE", "POST", "MEM", "22q11"))
}

# --- Create the Universal Cluster Color Dictionary ---
# Figure out exactly how many clusters exist in this specific object
num_clusters <- length(levels(sobj$seurat_clusters))
# Stretch the Wes Anderson palette to fit the exact number of clusters
CLUSTER_COLORS <- wes_palette(WES_PALETTE_NAME, num_clusters, type = "continuous")
names(CLUSTER_COLORS) <- levels(sobj$seurat_clusters)

# --- 1. UMAPs ---
# Add scale_color_manual to force the universal colors
p1 <- DimPlot(sobj, group.by = "seurat_clusters", label = TRUE) + 
  scale_color_manual(values = CLUSTER_COLORS) +
  ggtitle(paste(TARGET_CELL_TYPE, "Clusters")) + NoLegend()

p2 <- DimPlot(sobj, group.by = "starCAT_label") + 
  ggtitle("starCAT Labels")

p3 <- DimPlot(sobj, group.by = "seurat_clusters", split.by = "Clinical_Group", ncol = 3) +
  scale_color_manual(values = CLUSTER_COLORS) +
  theme(legend.position = "bottom")

final_overview <- (p1 | p2) / p3
# Update DPI:
ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Landscape_Overview.png"), plot = final_overview, width = 14, height = 10, dpi = PUB_DPI)

# --- 2. VDJ Expansion ---
# --- 2. VDJ Expansion (If it exists) ---
if("clonalFrequency" %in% colnames(sobj@meta.data)) {
  p4 <- FeaturePlot(sobj, features = "clonalFrequency", split.by = "Clinical_Group", order = TRUE, keep.scale = "all", pt.size = 0.8) & 
    
    # NEW DARJEELING GRADIENT:
    scale_color_gradientn(colors = c(HEATMAP_LOW, HEATMAP_MID, HEATMAP_HIGH), 
                          na.value = "grey90", name = "Clone Size") &
    
    theme(legend.position = "right")
  
  ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_VDJ_Expansion.png"), plot = p4, width = 15, height = 5, dpi = PUB_DPI)
}

# --- 3. Composition Barplot ---
comp_data <- as.data.frame(table(Cluster = sobj$seurat_clusters, Group = sobj$Clinical_Group)) %>%
  group_by(Group) %>% mutate(Proportion = Freq / sum(Freq)) %>% ungroup()

p_comp <- ggplot(comp_data, aes(x = Group, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = CLUSTER_COLORS) +   # Replace scale_fill_brewer with your universal colors!
  labs(title = paste(TARGET_CELL_TYPE, "Cluster Distribution"), x = "Clinical Group", y = "Percentage of Cells") +
  theme_minimal() + theme(axis.text.x = element_text(size = 12, face = "bold"))

ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Composition_Barplot.png"), plot = p_comp, width = 7, height = 7, dpi = PUB_DPI)
message(">>> Landscape plots saved successfully to ", OUTPUT_DIR)