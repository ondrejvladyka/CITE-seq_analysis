# ==============================================================================
# SCRIPT: 06b_Progeny_Analysis.R
# PURPOSE: Estimate signaling pathway activity using PROGENy (Seurat v5 Safe)
# ==============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(progeny)
library(pheatmap)

message("\n>>> Running PROGENy Pathway Analysis...")

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

# ---------------------------------------------------------
# THE SEURAT V5 FIX & PROGENY MATH
# ---------------------------------------------------------
message(">>> Extracting RNA matrix natively for Seurat v5...")
# Extract the data and forcefully uncompress it to a dense matrix for PROGENy
expr_matrix <- as.matrix(GetAssayData(sobj, assay = "RNA", layer = "data"))

message(">>> Calculating PROGENy signaling activities (Top 500 footprint genes)...")
progeny_matrix <- progeny(expr_matrix, scale = FALSE, organism = "Human", top = 1000)

# Manually inject the PROGENy math back into the Seurat object
sobj[["progeny"]] <- CreateAssayObject(data = t(progeny_matrix))

# Scale the progeny scores so we can compare them visually across cells
DefaultAssay(sobj) <- "progeny"
sobj <- ScaleData(sobj, assay = "progeny", verbose = FALSE)

# Extract the 14 pathway names
progeny_pathways <- rownames(sobj[["progeny"]])

# ---------------------------------------------------------
# PLOT 1: Global Pathway Heatmap (All Clusters)
# ---------------------------------------------------------
message(">>> Generating PROGENy Cluster Heatmap...")

progeny_scores <- as.data.frame(t(GetAssayData(sobj, assay = "progeny", layer = "scale.data")))
progeny_scores$Cluster <- Idents(sobj)

avg_scores <- progeny_scores %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean)) %>%
  tibble::column_to_rownames("Cluster")

heatmap_file <- paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_PROGENy_Global_Heatmap.png")
png(heatmap_file, width = 800, height = 600, res = 150)
pheatmap(t(avg_scores), 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = paste("PROGENy Signaling Activity:", TARGET_CELL_TYPE, "Clusters"),
         angle_col = 0,
         treeheight_row = 15, 
         treeheight_col = 15)
dev.off()

# ---------------------------------------------------------
# PLOT 1.5: Global UMAP Feature Overlays (Cell-by-Cell)
# ---------------------------------------------------------
message(">>> Generating PROGENy UMAP Overlays (This may take a moment)...")

# We map the scaled PROGENy activity directly onto the UMAP coordinates
p_umap <- FeaturePlot(sobj, features = progeny_pathways, ncol = 4, pt.size = 0.5) &
  scale_color_gradientn(colors = c("navy", "grey90", "firebrick3"), name = "Activity") &
  theme(plot.title = element_text(face = "bold", size = 12))

umap_file <- paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_PROGENy_UMAP_Overlay.png")
ggsave(umap_file, plot = p_umap, width = 18, height = 14, dpi = PUB_DPI)

# ---------------------------------------------------------
# PLOT 2: Targeted Signaling Comparison (Separate Violins)
# ---------------------------------------------------------
message(">>> Generating Individual Targeted Pathway Violins for: ", COMPARISON_PREFIX)

DefaultAssay(sobj) <- "RNA" # Revert to RNA for safe metadata operations

group_A_label <- paste0("C", NAME_A)
group_B_label <- ifelse("Rest" %in% CLUSTER_B, "Rest", paste0("C", NAME_B))

if ("Rest" %in% CLUSTER_B) {
  sobj_sub <- sobj
  sobj_sub$Target_Group <- group_B_label
  sobj_sub$Target_Group[sobj_sub$seurat_clusters %in% CLUSTER_A] <- group_A_label
} else {
  sobj_sub <- subset(sobj, seurat_clusters %in% c(CLUSTER_A, CLUSTER_B))
  sobj_sub$Target_Group <- "Unknown"
  sobj_sub$Target_Group[sobj_sub$seurat_clusters %in% CLUSTER_A] <- group_A_label
  sobj_sub$Target_Group[sobj_sub$seurat_clusters %in% CLUSTER_B] <- group_B_label
}

sobj_sub$Target_Group <- factor(sobj_sub$Target_Group, levels = c(group_A_label, group_B_label))

DefaultAssay(sobj_sub) <- "progeny"
expr_data <- FetchData(sobj_sub, vars = c(progeny_pathways, "Target_Group", "Clinical_Group"))
expr_long <- pivot_longer(expr_data, cols = all_of(progeny_pathways), names_to = "Pathway", values_to = "Activity")

target_pal <- setNames(c(HEATMAP_HIGH, HEATMAP_LOW), c(group_A_label, group_B_label))

# --- THE NEW LOOP: Generate a distinct image for every pathway ---
for(pathway in progeny_pathways) {
  
  # Filter data to only the current pathway
  plot_data <- expr_long %>% filter(Pathway == pathway)
  
  p_prog_viol <- ggviolin(plot_data, 
                          x = "Target_Group", 
                          y = "Activity", 
                          fill = "Target_Group",
                          palette = target_pal, 
                          add = "boxplot", 
                          add.params = list(fill = "white", width = 0.1, size = 0.2),
                          trim = TRUE) +
    stat_compare_means(comparisons = list(c(group_A_label, group_B_label)), 
                       label = "p.signif", 
                       method = "wilcox.test") +
    labs(title = paste("PROGENy:", pathway), 
         subtitle = paste(group_A_label, "vs", group_B_label),
         x = "", y = "Pathway Activity Score") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(size = 12, face = "bold"))
  
  # Clean up pathway name just in case it has weird characters for the filesystem
  clean_pathway <- gsub("[^A-Za-z0-9]", "_", pathway)
  viol_file <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_PROGENy_Violin_", clean_pathway, ".png")
  
  ggsave(viol_file, plot = p_prog_viol, width = 5, height = 5, dpi = PUB_DPI)
}

message(">>> PROGENy Analysis complete! 14 individual violins and UMAP overlay saved to: ", OUTPUT_DIR)

# Important: Reset Default Assay back to RNA so downstream scripts don't break
DefaultAssay(sobj) <- "RNA"