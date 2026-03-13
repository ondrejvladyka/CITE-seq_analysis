# ==============================================================================
# SCRIPT: Pre-Processing and Global Landscape Mapping
# ==============================================================================

# 1. Clear the environment and free up RAM
rm(list = ls())
gc()

# 2. Load necessary libraries
library(Seurat)
library(qs2)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratWrappers) # <--- Required for RunGLMPCA

# ==============================================================================
# MASTER CONTROL TWEAK
# ==============================================================================
TARGET_CELL <- "CD4"  

# Set to TRUE to use the advanced GLM-PCA, or FALSE for classical PCA
USE_GLMPCA <- TRUE    

INPUT_ARCHIVE <- paste0(TARGET_CELL, "_Archive_with_Spike.qs2")

# Dynamically set the output name so they don't overwrite each other
if(USE_GLMPCA){
  PROCESSED_FILE <- paste0("data/", TARGET_CELL, "_TargetCohorts_GLMPCA.qs2")
} else {
  PROCESSED_FILE <- paste0("data/", TARGET_CELL, "_TargetCohorts_Analyzed.qs2")
}

# ==============================================================================
# PART 1: Data Loading & Reprocessing
# ==============================================================================
message(">>> Loading ", TARGET_CELL, " Archive...")
sobj_raw <- qs_read(INPUT_ARCHIVE)

# Assign Clinical Groups
sobj_raw$Clinical_Group <- "Exclude"
sobj_raw$Clinical_Group[sobj_raw$Final_Combined_ID %in% c("BrisudaPRE", "LetakPRE", "SatnikPRE")] <- "PRE"
sobj_raw$Clinical_Group[sobj_raw$Final_Combined_ID %in% c("BrisudaPOST", "LetakPOST", "SatnikPOST")] <- "POST"
sobj_raw$Clinical_Group[sobj_raw$Final_Combined_ID %in% c("Jelinek", "Prazak", "Stringer")] <- "22q11"

# SUBSET FIRST
message(">>> Subsetting to Target Cohorts...")
sobj_target <- subset(sobj_raw, subset = Clinical_Group %in% c("PRE", "POST", "22q11"))

# Streamlined Analysis Function
reprocess_target_subset <- function(obj, name_prefix) {
  
  DefaultAssay(obj) <- "RNA"
  
  # 1. Normalize and find features based ONLY on these 3 groups
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = 2000, verbose = FALSE)
  
  # Clear out old dimensional reductions to avoid ghost warnings
  obj[["pca"]] <- NULL
  obj[["umap"]] <- NULL
  
  if(USE_GLMPCA) {
    # ---------------------------------------------------------
    # ROUTE A: GLM-PCA
    # ---------------------------------------------------------
    message(paste0("\n>>> Re-processing ", name_prefix, " (GLM-PCA Mode)..."))
    message(">>> Calculating GLM-PCA (Note: This is computationally heavy)...")
    
    obj <- RunGLMPCA(obj, 
                     features = VariableFeatures(obj), 
                     L = 20, 
                     minibatch = "memoized")
    
    message(">>> Running UMAP and Clustering on GLM-PCA space...")
    obj <- RunUMAP(obj, reduction = "glmpca", dims = 1:20, verbose = FALSE)
    obj <- FindNeighbors(obj, reduction = "glmpca", dims = 1:20, verbose = FALSE)
    
  } else {
    # ---------------------------------------------------------
    # ROUTE B: CLASSICAL PCA
    # ---------------------------------------------------------
    message(paste0("\n>>> Re-processing ", name_prefix, " (Classical PCA Mode)..."))
    
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 20, verbose = FALSE)
    
    message(">>> Running UMAP and Clustering on PCA space...")
    obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  }
  
  # Cluster using whichever graph was just generated
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  # Save using the dynamic filename
  qs_save(obj, PROCESSED_FILE)
  message(">>> Successfully saved: ", PROCESSED_FILE)
  
  return(obj)
}

# Execute
sobj_final <- reprocess_target_subset(sobj_target, TARGET_CELL)

# ==============================================================================
# PART 2: Visualizing the Landscape
# ==============================================================================
# Clear raw objects to save RAM before plotting
rm(sobj_raw, sobj_target, sobj_final)
gc()

# Load the fresh, hyper-clean object
sobj <- qs_read(PROCESSED_FILE)

# Set the factor levels for consistent plotting order
sobj$Clinical_Group <- factor(sobj$Clinical_Group, levels = c("PRE", "POST", "22q11"))

# ---------------------------------------------------------
# FIGURE 1: The Landscape (Clusters & original Atlas Labels)
# ---------------------------------------------------------
p1 <- DimPlot(sobj, group.by = "seurat_clusters", label = TRUE) + 
  ggtitle(paste(TARGET_CELL, "- New Granular Clusters")) + NoLegend()

p2 <- DimPlot(sobj, group.by = "starCAT_label") + 
  ggtitle(paste(TARGET_CELL, "- Original starCAT Labels"))

# ---------------------------------------------------------
# FIGURE 2: Clinical Group Distribution (Split UMAP)
# ---------------------------------------------------------
p3 <- DimPlot(sobj, group.by = "seurat_clusters", split.by = "Clinical_Group", ncol = 3) +
  ggtitle(paste(TARGET_CELL, "Density Shifts: PRE vs POST vs 22q11")) +
  theme(legend.position = "bottom")

# ---------------------------------------------------------
# FIGURE 3: VDJ Clonal Expansion Overlay
# ---------------------------------------------------------
if("clonalFrequency" %in% colnames(sobj@meta.data)) {
  p4 <- FeaturePlot(sobj, 
                    features = "clonalFrequency", 
                    split.by = "Clinical_Group", 
                    order = TRUE,          # Expanded cells are drawn on top
                    keep.scale = "all",    # Ensures 22q11 and POST use the same color logic
                    pt.size = 0.8) & 
    scale_color_viridis_c(option = "plasma", 
                          na.value = "grey90", 
                          direction = -1,
                          name = "Clone Size") &
    theme(legend.position = "right")
  
  print(p4)
  ggsave(paste0(TARGET_CELL, "_VDJ_Expansion_Numeric.png"), plot = p4, width = 15, height = 5, dpi = 300)
} else {
  message(">>> WARNING: 'clonalFrequency' not found. Skipping VDJ overlay plot.")
  p4 <- NULL
}

# ---------------------------------------------------------
# FIGURE 4: Composition Barplot
# ---------------------------------------------------------
comp_data <- as.data.frame(table(Cluster = sobj$seurat_clusters, Group = sobj$Clinical_Group))

# Calculate proportions
comp_data <- comp_data %>%
  group_by(Group) %>%
  mutate(Proportion = Freq / sum(Freq)) %>%
  ungroup()

p_comp <- ggplot(comp_data, aes(x = Group, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set3") + # Clean, distinct colors
  labs(title = paste(TARGET_CELL, "Cluster Distribution"), 
       subtitle = "Comparing PRE, POST, and 22q11 Cohorts",
       x = "Clinical Group", 
       y = "Percentage of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

print(p_comp)

# ---------------------------------------------------------
# SAVE ALL FIGURES
# ---------------------------------------------------------
final_overview <- (p1 | p2) / p3
ggsave(paste0(TARGET_CELL, "_Final_Overview.png"), plot = final_overview, width = 14, height = 10, dpi = 300)
ggsave(paste0(TARGET_CELL, "_Final_Composition_Barplot.png"), plot = p_comp, width = 7, height = 7, dpi = 300)

message(">>> Script Complete!")