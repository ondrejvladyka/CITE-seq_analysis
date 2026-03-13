# ==============================================================================
# SCRIPT: 00_Master_Control.R
# PURPOSE: Set target variables, create output directories, and run pipeline
# ==============================================================================

#BiocManager::install("FlowSOM")
#install.packages("ggalluvial")
library(FlowSOM)

# 1. DEFINE YOUR TARGETS HERE
TARGET_CELL_TYPE <- "CD4"           
CLUSTER_A <- c(11)                # Can be a single number (0) or a vector (c(0, 1))
#CLUSTER_B <- c(1) 
CLUSTER_B <- "Rest"




# --- NEW: Collapse vectors into text for clean file naming ---
NAME_A <- paste(CLUSTER_A, collapse = "_")
#NAME_B <- paste(CLUSTER_B, collapse = "_")
NAME_B <- "Rest"

# ==============================================================================
# 2. SCRIPT-SPECIFIC VARIABLES & PUBLICATION SETTINGS
# ==============================================================================
library(wesanderson)

# --- Define Global Colors ---
WES_PALETTE_NAME <- "Darjeeling1" 

# Extract the base 5 colors and drop the alpha to 70% (0.7) to make them muted/pastel
BASE_PAL <- adjustcolor(wes_palette(WES_PALETTE_NAME, 5), alpha.f = 0.7)

# 1. Clinical Group Colors (Teal for PRE, Gold for POST, Red for 22q11)
GROUP_COLORS <- c("PRE" = BASE_PAL[5], 
                  "POST" = BASE_PAL[4], 
                  "22q11" = BASE_PAL[1])

# 2. Heatmap Gradient (First color -> White -> Last color)
# Usually, Blue (5th color) is low expression, and Red (1st color) is high expression
HEATMAP_LOW  <- BASE_PAL[5]  
HEATMAP_MID  <- BASE_PAL[2]
HEATMAP_HIGH <- BASE_PAL[1]  

# --- Universal Image Settings ---
PUB_DPI <- 600           
VOLCANO_LABEL_SIZE <- 3  

# --- Script 01b (Deep Landscape) Variables ---
CANONICAL_GENES <- c("CCR7", "SELL", "TCF7", "LEF1", "IL7R", "CD27", "CD28", 
                     "GZMB", "PRF1", "NKG7", "GNLY", "CXCR4", "S100A4", "ANXA1", 
                     "PDCD1", "HAVCR2", "TIGIT", "TOX", "LAG3", "ENTPD1", "KLRG1", 
                     "B3GAT1", "MKI67", "TOP2A")

KEY_FEATURES <- c("SELL", "GNLY", "CXCR4")

# --- Script 02 (DGE & Volcano) Variables ---
EXTRA_GENES_TO_LABEL <- c("CXCR4", "TOX", "HAVCR2", "S100A4")


# Script 06 (Signature Scoring) Variables ---
# Define lists of genes that represent broad biological states
SIGNATURES <- list(
  "Cytotoxicity" = c("PRF1", "GZMB", "GZMA", "GNLY", "NKG7", "FGFBP2"),
  "Exhaustion"   = c("PDCD1", "HAVCR2", "TIGIT", "LAG3", "TOX", "CTLA4", "ENTPD1")
)

# ==============================================================================
# 3. AUTOMATED DIRECTORY & SETUP (Do not touch)
# ==============================================================================
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- paste0("results_", TARGET_CELL_TYPE, "_", TIMESTAMP, "/")
if (!dir.exists(OUTPUT_DIR)) { dir.create(OUTPUT_DIR, recursive = TRUE) }

INPUT_OBJECT_PATH <- paste0("data/", TARGET_CELL_TYPE, "_TargetCohorts_GLMPCA.qs2")

# Uses the new collapsed names!
COMPARISON_PREFIX <- paste0(TARGET_CELL_TYPE, "_C", NAME_A, "_vs_C", NAME_B)
DEG_CSV_FILE      <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_FullSpectrum.csv")

# ==============================================================================
# GLOBAL OBJECT LOADER & ALIAS CHECK
# ==============================================================================

if(!exists("sobj")) { 
  message(">>> Loading primary Seurat object from: ", INPUT_OBJECT_PATH)
  sobj <- qs_read(INPUT_OBJECT_PATH) 
  
  # --- THE DOWNSTREAM RAM ALIAS ---
  # If the object was processed with GLM-PCA, alias it to "pca" for standard Seurat functions
  if("glmpca" %in% names(sobj@reductions) && !("pca" %in% names(sobj@reductions))) {
    message(">>> GLM-PCA detected. Creating a temporary 'pca' alias in RAM for downstream compatibility...")
    sobj[["pca"]] <- CreateDimReducObject(embeddings = Embeddings(sobj, "glmpca"), 
                                          key = "PC_", 
                                          assay = DefaultAssay(sobj))
  }
}

message("\n==================================================")
message(">>> PIPELINE READY: ", COMPARISON_PREFIX)
message("==================================================\n")

# 4. RUN THE PIPELINE 
source("00b_flowSOM_takeover.R")
source("01_Landscape_Analysis.R")
source("01b_Deep_Landscape.R")
source("01c_Alternative_Visualizations.R")
source("02_Run_Differential_Genes.R")
source("03_Run_GSEA.R")
source("04_Targeted_Clonality.R")
source("05_Targeted_Expression.R")
source("06_Signature_Scoring.R")
source("07_Clonal_Tracking.R") 
source("08_Clonal_Purity.R")













