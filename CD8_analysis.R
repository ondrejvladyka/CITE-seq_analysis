# 1. Clear the environment and free up RAM
rm(list = ls())
gc()

# 2. Load necessary libraries
library(Seurat)
library(qs2)
library(ggplot2)
library(patchwork)
library(dplyr)





install.packages('devtools')
devtools::install_github('immunogenomics/presto')




# 1. Load the un-processed CD8 Archive
message(">>> Loading CD8 Archive...")
sobj_cd8_raw <- qs_read("CD8_Archive_with_Spike.qs2")

# 2. Assign Clinical Groups
sobj_cd8_raw$Clinical_Group <- "Exclude"
sobj_cd8_raw$Clinical_Group[sobj_cd8_raw$Final_Combined_ID %in% c("BrisudaPRE", "LetakPRE", "SatnikPRE")] <- "PRE"
sobj_cd8_raw$Clinical_Group[sobj_cd8_raw$Final_Combined_ID %in% c("BrisudaPOST", "LetakPOST", "SatnikPOST")] <- "POST"
sobj_cd8_raw$Clinical_Group[sobj_cd8_raw$Final_Combined_ID %in% c("Jelinek", "Prazak", "Stringer")] <- "22q11"

# 3. SUBSET FIRST
message(">>> Subsetting to Target Cohorts...")
sobj_cd8_target <- subset(sobj_cd8_raw, subset = Clinical_Group %in% c("PRE", "POST", "22q11"))

# 4. Streamlined Analysis (No MT regression)
reprocess_target_subset <- function(obj, name_prefix) {
  message(paste0("\n>>> Re-processing ", name_prefix, " (Pure Variance Mode)..."))
  
  DefaultAssay(obj) <- "RNA"
  
  # Normalize and find features based ONLY on these 3 groups
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = 2000, verbose = FALSE)
  
  # Scale without regression
  obj <- ScaleData(obj, verbose = FALSE)
  
  # Dimensionality Reduction
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  # Save
  file_name <- paste0("data/", name_prefix, "_TargetCohorts_Analyzed.qs2")
  qs_save(obj, file_name)
  message(">>> Successfully saved: ", file_name)
  
  return(obj)
}

# 5. Execute
sobj_cd8_final <- reprocess_target_subset(sobj_cd8_target, "CD8")







# 1. Load the fresh, hyper-clean object
sobj <- qs_read("data/CD8_TargetCohorts_Analyzed.qs2")

# Set the factor levels for consistent plotting order
sobj$Clinical_Group <- factor(sobj$Clinical_Group, levels = c("PRE", "POST", "22q11"))

# ---------------------------------------------------------
# FIGURE 1: The Landscape (Clusters & original Atlas Labels)
# ---------------------------------------------------------
p1 <- DimPlot(sobj, group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("New Granular Clusters") + NoLegend()

p2 <- DimPlot(sobj, group.by = "starCAT_label") + 
  ggtitle("Original starCAT Labels")

# ---------------------------------------------------------
# FIGURE 2: Clinical Group Distribution (Split UMAP)
# ---------------------------------------------------------
p3 <- DimPlot(sobj, group.by = "seurat_clusters", split.by = "Clinical_Group", ncol = 3) +
  ggtitle("CD8 Density Shifts: PRE vs POST vs 22q11") +
  theme(legend.position = "bottom")

# ---------------------------------------------------------
# FIGURE 3: VDJ Clonal Expansion Overlay
# ---------------------------------------------------------
# We use order = TRUE to ensure expanded clones are drawn on top of grey singletons
p4 <- FeaturePlot(sobj, 
                  features = "clonalFrequency", 
                  split.by = "Clinical_Group", 
                  order = TRUE,         # Expanded cells are drawn on top
                  keep.scale = "all",    # Ensures 22q11 and POST use the same color logic
                  pt.size = 0.8) & 
  scale_color_viridis_c(option = "plasma", 
                        na.value = "grey90", 
                        direction = -1,
                        name = "Clone Size") &
  theme(legend.position = "right")

# Display the plot
print(p4)

# Save the high-res version
ggsave("CD8_VDJ_Expansion_Numeric.png", plot = p4, width = 15, height = 5, dpi = 300)

comp_data <- as.data.frame(table(Cluster = sobj$seurat_clusters, Group = sobj$Clinical_Group))

# 2. Calculate proportions (Ensuring the pipe works)
comp_data <- comp_data %>%
  group_by(Group) %>%
  mutate(Proportion = Freq / sum(Freq)) %>%
  ungroup()

# 3. Check the result
print(head(comp_data))

p5 <- ggplot(comp_data, aes(x = Group, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Cluster Proportions by Cohort", y = "Percentage of CD8 T cells")

# ---------------------------------------------------------
# SAVE ALL FIGURES
# ---------------------------------------------------------
# Combine the first two for an overview
final_overview <- (p1 | p2) / p3
ggsave("CD8_Final_Overview.png", plot = final_overview, width = 14, height = 10)

# Save VDJ and Composition separately
ggsave("CD8_Final_VDJ_Split.png", plot = p4, width = 15, height = 5)
ggsave("CD8_Final_Composition.png", plot = p5, width = 6, height = 7)

# Display the split UMAP now
print(p3)

p_comp <- ggplot(comp_data, aes(x = Group, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set3") + # Clean, distinct colors
  labs(title = "CD8 Cluster Distribution", 
       subtitle = "Comparing PRE, POST, and 22q11 Cohorts",
       x = "Clinical Group", 
       y = "Percentage of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

print(p_comp)
ggsave("CD8_Final_Composition_Barplot.png", plot = p_comp, width = 7, height = 7)











































# 1. Run the direct comparison without strict logFC filters
# only.pos = FALSE allows us to see what is downregulated as well
markers_full_spectrum <- FindMarkers(sobj_0_1, 
                                     ident.1 = 0, 
                                     ident.2 = 1, 
                                     only.pos = FALSE, 
                                     min.pct = 0.1, 
                                     logfc.threshold = 0.0) 

# 2. Move gene names to a column for easier sorting/saving
markers_full_spectrum$gene <- rownames(markers_full_spectrum)

# 3. Sort from most positive (high in 0) to most negative (high in 1)
markers_sorted <- markers_full_spectrum %>%
  arrange(desc(avg_log2FC))

# 4. Save the entire sorted list to a CSV so you can scroll through it in Excel
write.csv(markers_sorted, "CD8_Cluster0_vs_Cluster1_FullSpectrum.csv", row.names = FALSE)
message(">>> Saved full ranked list to CD8_Cluster0_vs_Cluster1_FullSpectrum.csv")

# 5. Preview the extremes in the console
message("\n>>> TOP 15 GENES HIGHER IN CLUSTER 0 (PRE/POST dominated):")
print(head(markers_sorted, 15))

message("\n>>> TOP 15 GENES HIGHER IN CLUSTER 1 (22q11/POST dominated):")
print(tail(markers_sorted, 15))






# Create a simple category for coloring
markers_sorted$Significance <- "Not Significant"
markers_sorted$Significance[markers_sorted$avg_log2FC > 0.25 & markers_sorted$p_val_adj < 0.05] <- "Up in Cluster 0"
markers_sorted$Significance[markers_sorted$avg_log2FC < -0.25 & markers_sorted$p_val_adj < 0.05] <- "Up in Cluster 1"

p_volcano <- ggplot(markers_sorted, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey80", 
                                "Up in Cluster 0" = "blue", 
                                "Up in Cluster 1" = "red")) +
  theme_minimal() +
  labs(title = "Cluster 0 vs Cluster 1 Transcriptomic Shift",
       x = "Average Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black")

print(p_volcano)
ggsave("CD8_Volcano_0_vs_1.png", plot = p_volcano, width = 8, height = 6)







