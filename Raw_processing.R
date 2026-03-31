# ==============================================================================
# PIPELINE CONFIGURATION: FILL IN YOUR PROJECT DETAILS BELOW
# ==============================================================================

# 1. Project Name (used for naming output files)
PROJECT_NAME <- "Study_2026"

# 2. Output Directory (where results and plots will be saved)
OUTPUT_DIR <- "./analysis_results"

# 3. Input Samples (Sample Name = Path to its filtered_feature_bc_matrix folder)
# ADVICE: Ensure paths are absolute or relative to this script's location.
SAMPLE_INPUTS <- list(
  "Patient_A_Baseline" = "/data/run1/outs/filtered_feature_bc_matrix",
  "Patient_A_Post"     = "/data/run2/outs/filtered_feature_bc_matrix",
  "Healthy_Control"    = "/data/control/outs/filtered_feature_bc_matrix"
)

# 4. starCAT Reference Catalog
# Options: "TCAT.V1" (T-cells), "MYELOID.GLIOMA.V1", "BONEMARROW.CD34POS.HSPC.V1"
STARCAT_REF <- "TCAT.V1"

# 5. QC Thresholds (Adjust based on your tissue quality)
MIN_FEATURES <- 500    # Filter cells with very few genes detected
MAX_FEATURES <- 12000   # Filter potential doublets with too many genes
MAX_MT_PCT   <- 12     # Filter dying cells (high mitochondrial content)

# ==============================================================================
# PIPELINE EXECUTION ENGINE (Do not edit below unless modifying the workflow)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(scDblFinder)
  library(reticulate)
  library(dplyr)
  library(ggplot2)
})

# Create output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ------------------------------------------------------------------------------
# STEP 1: LOAD DATA & REMOVE DEAD CELLS (QC)
# ------------------------------------------------------------------------------
message("\n[STATUS] Step 1: Loading samples and performing QC filtering...")

obs_list <- lapply(names(SAMPLE_INPUTS), function(sn) {
  counts <- Read10X(data.dir = SAMPLE_INPUTS[[sn]])
  obj <- CreateSeuratObject(counts = counts, project = sn)
  obj$sample_id <- sn
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  return(obj)
})

# Merge samples into a single object
merged_obj <- merge(obs_list[[1]], y = obs_list[-1], add.cell.ids = names(SAMPLE_INPUTS))
rm(obs_list) # Free up memory

# Filter dead cells and poor quality barcodes
merged_obj <- subset(merged_obj, subset = nFeature_RNA > MIN_FEATURES & 
                       nFeature_RNA < MAX_FEATURES & 
                       percent.mt < MAX_MT_PCT)

message(sprintf("[INFO] Cells remaining after QC: %d", ncol(merged_obj)))

# ------------------------------------------------------------------------------
# STEP 2: DOUBLET DELETION (scDblFinder)
# ------------------------------------------------------------------------------
message("\n[STATUS] Step 2: Running scDblFinder to detect and remove doublets...")

# Convert to SingleCellExperiment object for scDblFinder
sce <- as.SingleCellExperiment(merged_obj)
# Run doublet detection per sample
sce <- scDblFinder(sce, samples = "sample_id")

# Add doublet calls back to Seurat metadata and filter
merged_obj$scDblFinder.class <- sce$scDblFinder.class
doublet_count <- sum(merged_obj$scDblFinder.class == "doublet")
message(sprintf("[INFO] Removing %d identified doublets.", doublet_count))

merged_obj <- subset(merged_obj, subset = scDblFinder.class == "singlet")

# ------------------------------------------------------------------------------
# STEP 3: NORMALIZATION & CELL CYCLE CORRECTION
# ------------------------------------------------------------------------------
message("\n[STATUS] Step 3: Normalizing data and correcting for cell cycle...")

# Initial normalization to score cell cycle phases
merged_obj <- SCTransform(merged_obj, vars.to.regress = "percent.mt", verbose = FALSE)

# Score Cell Cycle phases
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
merged_obj <- CellCycleScoring(merged_obj, s.features = s.genes, g2m.features = g2m.genes, assay = "SCT")

# Re-run SCTransform, this time regressing out cell cycle scores
message("[INFO] Regressing out S.Score and G2M.Score...")
merged_obj <- SCTransform(merged_obj, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE)

# ------------------------------------------------------------------------------
# STEP 4: DIMENSIONALITY REDUCTION & HARMONY INTEGRATION
# ------------------------------------------------------------------------------
message("\n[STATUS] Step 4: Running PCA and Harmony Integration...")

merged_obj <- RunPCA(merged_obj, npcs = 30, verbose = FALSE)

# Integrate across the 'sample_id' variable to harmonize batches
merged_obj <- RunHarmony(merged_obj, group.by.vars = "sample_id", 
                         dims.use = 1:30, assay.use = "SCT", 
                         reduction.save = "harmony")

# Generate UMAP and standard clusters based on the harmonized space
merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:30)
merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.5, verbose = FALSE)

# ------------------------------------------------------------------------------
# STEP 5: STARCAT ANNOTATION
# ------------------------------------------------------------------------------
message("\n[STATUS] Step 5: Running starCAT annotation (Reference: ", STARCAT_REF, ")...")

# Define the starCAT wrapper function
run_starCAT_annotation <- function(seurat_obj, reference, assay = "RNA") {
  counts <- GetAssayData(seurat_obj, assay = assay, layer = "counts")
  
  # Import starcat via reticulate
  # Ensure your active conda/virtual environment has starcat installed (`pip install starcat`)
  sc <- import("starcat")
  estimator <- sc$StarCAT(reference = reference)
  
  # Fit the model and extract usage scores
  results <- estimator$fit_transform(t(as.matrix(counts)))
  
  # Add continuous scores to metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = results)
  
  # Define discrete clusters based on the highest scoring program per cell
  seurat_obj$starCAT_clusters <- colnames(results)[apply(results, 1, which.max)]
  return(seurat_obj)
}

# Run the annotation
merged_obj <- run_starCAT_annotation(merged_obj, reference = STARCAT_REF)

# ------------------------------------------------------------------------------
# STEP 6: SAVE FINAL OUTPUTS
# ------------------------------------------------------------------------------
message("\n[STATUS] Step 6: Saving outputs...")

# Save a UMAP plot comparing standard clustering to starCAT annotation
pdf(file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_UMAP_Comparison.pdf")), width = 12, height = 6)
p1 <- DimPlot(merged_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Standard Clusters")
p2 <- DimPlot(merged_obj, reduction = "umap", group.by = "starCAT_clusters", label = TRUE, repel = TRUE) + ggtitle("starCAT Programs")
print(p1 + p2)
dev.off()

# Save the fully processed Seurat object
save_path <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_final_annotated.rds"))
saveRDS(merged_obj, save_path)

message("\n[SUCCESS] Pipeline complete! Object saved to: ", save_path)