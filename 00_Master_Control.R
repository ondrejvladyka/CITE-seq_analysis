# ==============================================================================
# SCRIPT: 00_Master_Control.R
# PURPOSE: Set target variables, create output directories, and run pipeline
# ==============================================================================

# 1. DEFINE YOUR TARGETS HERE
TARGET_CELL_TYPE <- "CD8"           # Options: "CD8", "CD4", "Thymus"
CLUSTER_A <- 0                      # First cluster (e.g., PRE/POST dominant)
CLUSTER_B <- 1                      # Second cluster (e.g., 22q11 dominant)

# 2. AUTOMATED DIRECTORY & SETUP (Do not touch)
# Generate timestamp (Format: YYYYMMDD_HHMMSS)
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create the specific output folder name
OUTPUT_DIR <- paste0("results_", TARGET_CELL_TYPE, "_", TIMESTAMP, "/")

# Tell your computer to actually create the folder if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Define file paths
INPUT_OBJECT_PATH <- paste0("data/", TARGET_CELL_TYPE, "_TargetCohorts_Analyzed.qs2")
COMPARISON_PREFIX <- paste0(TARGET_CELL_TYPE, "_C", CLUSTER_A, "_vs_C", CLUSTER_B)
DEG_CSV_FILE      <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_FullSpectrum.csv")

message("\n==================================================")
message(">>> PIPELINE READY")
message(">>> TARGET OBJECT: ", TARGET_CELL_TYPE)
message(">>> COMPARISON: Cluster ", CLUSTER_A, " vs Cluster ", CLUSTER_B)
message(">>> ALL OUTPUTS GOING TO: ", OUTPUT_DIR)
message("==================================================\n")

# 3. RUN THE PIPELINE 
source("01_Landscape_Analysis.R")
source("01b_Deep_Landscape.R")
source("02_Run_Differential_Genes.R")
# source("03_GSEA_GOBP.R")

