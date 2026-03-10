library(Seurat)
library(qs2)

# 1. Define the exact labels based on your output
cd8_labels <- c("CD8_CM", "CD8_EM", "CD8_Naive", "CD8_TEMRA")
cd4_labels <- c("CD4_CM", "CD4_EM", "CD4_Naive", "Treg")
thymus_ids <- c("libT1", "libT2", "libT3")

# ---------------------------------------------------------
# A. THYMUS SUBSETS (Using Final_Combined_ID instead of starCAT)
# ---------------------------------------------------------
message("\n>>> Subsetting Thymus cells...")
# Grab only the cells from the Thymus libraries
sobj_thymus_all <- subset(sobj, subset = Final_Combined_ID %in% thymus_ids)
# Thymus libraries don't contain Spike, but we filter just to be 100% consistent
sobj_thymus_clean <- subset(sobj_thymus_all, subset = Final_Combined_ID != "Spike")

qs_save(sobj_thymus_all, "Thymus_Archive_with_Spike.qs2")
message("Thymus Clean cells: ", ncol(sobj_thymus_clean))

# ---------------------------------------------------------
# B. CD8 SUBSETS (Blood only)
# ---------------------------------------------------------
message("\n>>> Subsetting CD8 cells...")
# Grab CD8 labels, BUT explicitly exclude the Thymus libraries
sobj_cd8_all <- subset(sobj, subset = starCAT_label %in% cd8_labels & !(Final_Combined_ID %in% thymus_ids))
sobj_cd8_clean <- subset(sobj_cd8_all, subset = Final_Combined_ID != "Spike")

qs_save(sobj_cd8_all, "CD8_Archive_with_Spike.qs2")
message("CD8 Clean (No Spike) cells: ", ncol(sobj_cd8_clean))

# ---------------------------------------------------------
# C. CD4 SUBSETS (Blood only)
# ---------------------------------------------------------
message("\n>>> Subsetting CD4 cells...")
# Grab CD4 labels, explicitly exclude Thymus
sobj_cd4_all <- subset(sobj, subset = starCAT_label %in% cd4_labels & !(Final_Combined_ID %in% thymus_ids))
sobj_cd4_clean <- subset(sobj_cd4_all, subset = Final_Combined_ID != "Spike")

qs_save(sobj_cd4_all, "CD4_Archive_with_Spike.qs2")
message("CD4 Clean (No Spike) cells: ", ncol(sobj_cd4_clean))


reprocess_clean_subset <- function(obj, name_prefix) {
  message(paste0("\n>>> Re-processing ", name_prefix, " using available data slots..."))
  
  DefaultAssay(obj) <- "RNA"
  
  # 1. Use "dispersion" because the raw counts were removed to save file size
  obj <- FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = 2000, verbose = FALSE)
  
  # 2. Scale Data (Check if percent.mt exists before regressing to be safe)
  if ("percent.mt" %in% colnames(obj@meta.data)) {
    obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  } else {
    obj <- ScaleData(obj, verbose = FALSE)
  }
  
  # 3. Dimensionality Reduction & Clustering
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  # 4. Save
  file_name <- paste0("data/", name_prefix, "_Clean_Analyzed.qs2")
  qs_save(obj, file_name)
  message(">>> Successfully saved: ", file_name)
  
  return(obj)
}

# Run it on your clean subsets!
sobj_cd8_clean <- reprocess_clean_subset(sobj_cd8_clean, "CD8")
sobj_cd4_clean <- reprocess_clean_subset(sobj_cd4_clean, "CD4")
sobj_thymus_clean <- reprocess_clean_subset(sobj_thymus_clean, "Thymus")
