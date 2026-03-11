# ==============================================================================
# SCRIPT: 01_Landscape_Analysis.R
# ==============================================================================
library(Seurat)
library(qs2)
library(ggplot2)
library(patchwork)
library(dplyr)

message(">>> Loading Object: ", INPUT_OBJECT_PATH)
sobj <- qs_read(INPUT_OBJECT_PATH)

if("Clinical_Group" %in% colnames(sobj@meta.data)) {
  sobj$Clinical_Group <- factor(sobj$Clinical_Group, levels = c("PRE", "POST", "22q11"))
}

# --- 1. UMAPs ---
p1 <- DimPlot(sobj, group.by = "seurat_clusters", label = TRUE) + 
  ggtitle(paste(TARGET_CELL_TYPE, "Clusters")) + NoLegend()
p2 <- DimPlot(sobj, group.by = "starCAT_label") + 
  ggtitle("starCAT Labels")
p3 <- DimPlot(sobj, group.by = "seurat_clusters", split.by = "Clinical_Group", ncol = 3) +
  theme(legend.position = "bottom")

final_overview <- (p1 | p2) / p3
# Added OUTPUT_DIR here:
ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Landscape_Overview.png"), plot = final_overview, width = 14, height = 10, dpi = 300)

# --- 2. VDJ Expansion (If it exists) ---
if("clonalFrequency" %in% colnames(sobj@meta.data)) {
  p4 <- FeaturePlot(sobj, features = "clonalFrequency", split.by = "Clinical_Group", order = TRUE, keep.scale = "all", pt.size = 0.8) & 
    scale_color_viridis_c(option = "plasma", na.value = "grey90", direction = -1, name = "Clone Size") &
    theme(legend.position = "right")
  # Added OUTPUT_DIR here:
  ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_VDJ_Expansion.png"), plot = p4, width = 15, height = 5, dpi = 300)
}

# --- 3. Composition Barplot ---
comp_data <- as.data.frame(table(Cluster = sobj$seurat_clusters, Group = sobj$Clinical_Group)) %>%
  group_by(Group) %>% mutate(Proportion = Freq / sum(Freq)) %>% ungroup()

p_comp <- ggplot(comp_data, aes(x = Group, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set3") + 
  labs(title = paste(TARGET_CELL_TYPE, "Cluster Distribution"), x = "Clinical Group", y = "Percentage of Cells") +
  theme_minimal() + theme(axis.text.x = element_text(size = 12, face = "bold"))

# Added OUTPUT_DIR here:
ggsave(paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Composition_Barplot.png"), plot = p_comp, width = 7, height = 7, dpi = 300)
message(">>> Landscape plots saved successfully to ", OUTPUT_DIR)