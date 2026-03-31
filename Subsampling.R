# ==============================================================================
# SCRIPT: Subpopulation_Generation.R
# PURPOSE: Load base data, integrate VDJ, anonymize, subset, and run PCA
# ==============================================================================

# 1. Clear environment and load libraries
rm(list = ls())
gc()

setwd("C:/Users/ondre/OneDrive - Univerzita Karlova/BigData/CITE-seq_analysis")

library(Seurat)
library(qs2)
library(dplyr)
library(scRepertoire) # <-- Required for VDJ

# ==============================================================================
# 1. LOAD BASE DATA
# ==============================================================================
message(">>> Loading the global raw object...")
sobj <- qs_read("data/Klocperk_analysed.qs2") 

# ==============================================================================
# 2. VDJ INTEGRATION (The Nuclear Cleanup Version)
# ==============================================================================
message(">>> Loading VDJ annotation files...")
vdj_lib1 <- read.csv("data/filtered_contig_annotations_lib1.csv")
vdj_lib2 <- read.csv("data/filtered_contig_annotations_lib2.csv")
vdj_lib3 <- read.csv("data/filtered_contig_annotations_lib3.csv")
vdj_lib4 <- read.csv("data/filtered_contig_annotations_lib4.csv")

# scRepertoire v2 strictly requires a named list
vdj_list <- list("lib1" = vdj_lib1, "lib2" = vdj_lib2, "lib3" = vdj_lib3, "lib4" = vdj_lib4)

message(">>> Scrubbing formatting bugs from CellRanger output...")
for (name in names(vdj_list)) {
  # A. Delete the confusing sample column entirely
  vdj_list[[name]]$sample <- NULL 
  
  # B. Bulletproof boolean coercion (Looks for 't' or 'T', ignoring text casing)
  vdj_list[[name]]$productive      <- grepl("[Tt]", as.character(vdj_list[[name]]$productive))
  vdj_list[[name]]$high_confidence <- grepl("[Tt]", as.character(vdj_list[[name]]$high_confidence))
  vdj_list[[name]]$is_cell         <- grepl("[Tt]", as.character(vdj_list[[name]]$is_cell))
  
  # C. Strip trailing white spaces from chains just in case
  vdj_list[[name]]$chain <- trimws(as.character(vdj_list[[name]]$chain))
}

message(">>> Combining TCR data...")
combined_tcr <- combineTCR(vdj_list, samples = c("lib1", "lib2", "lib3", "lib4"))

# Fix barcodes (Strip the -1 suffix so it perfectly matches Seurat)
for (i in 1:length(combined_tcr)) {
  if(nrow(combined_tcr[[i]]) > 0) {
    combined_tcr[[i]]$barcode <- sub("-1$", "", combined_tcr[[i]]$barcode)
  }
}

message(">>> Integrating VDJ into Seurat Object...")
sobj <- combineExpression(combined_tcr, 
                          sobj, 
                          cloneCall = "gene", 
                          proportion = TRUE)

message(">>> Cells with matched TCR data globally: ", sum(!is.na(sobj$CTgene)))

# Clean up raw VDJ files from RAM to save memory
rm(vdj_lib1, vdj_lib2, vdj_lib3, vdj_lib4, vdj_list, combined_tcr)
gc()

# ==============================================================================
# 3. PATIENT ANONYMIZATION (HIPAA / Publication Prep)
# ==============================================================================
message("\n>>> Anonymizing Patient IDs...")

sobj$Anonymous_ID <- as.character(sobj$Final_Combined_ID)

# Dictionary for PRE/POST Controls
# ==============================================================================
# 3. PATIENT ANONYMIZATION (HIPAA / Publication Prep)
# ==============================================================================
message("\n>>> Anonymizing Patient IDs (Including MEM cells)...")

sobj$Anonymous_ID <- as.character(sobj$Final_Combined_ID)

# Dictionary for PRE/POST Controls (Merging MEM seamlessly into POST)
sobj$Anonymous_ID[sobj$Final_Combined_ID == "BrisudaPRE"]  <- "C1_PRE"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "BrisudaPOST"] <- "C1_POST"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "BrisudaMEM"]  <- "C1_POST" 

sobj$Anonymous_ID[sobj$Final_Combined_ID == "LetakPRE"]    <- "C2_PRE"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "LetakPOST"]   <- "C2_POST"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "LetakMEM"]    <- "C2_POST" 

sobj$Anonymous_ID[sobj$Final_Combined_ID == "SatnikPRE"]   <- "C3_PRE"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "SatnikPOST"]  <- "C3_POST"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "SatnikMEM"]   <- "C3_POST" 

# Dictionary for 22q11 Patients
sobj$Anonymous_ID[sobj$Final_Combined_ID == "Jelinek"]  <- "DG_1"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "Prazak"]   <- "DG_2"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "Stringer"] <- "DG_3"

# Overwrite the original column so all downstream plots use the anonymous names
sobj$Final_Combined_ID <- factor(sobj$Anonymous_ID)

# Re-assign Clinical Groups
sobj$Clinical_Group <- "Exclude"
sobj$Clinical_Group[grepl("PRE", sobj$Final_Combined_ID)] <- "PRE"
sobj$Clinical_Group[grepl("POST", sobj$Final_Combined_ID)] <- "POST"
sobj$Clinical_Group[grepl("DG_", sobj$Final_Combined_ID)] <- "22q11"

# Dictionary for 22q11 Patients
sobj$Anonymous_ID[sobj$Final_Combined_ID == "Jelinek"]  <- "DG_1"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "Prazak"]   <- "DG_2"
sobj$Anonymous_ID[sobj$Final_Combined_ID == "Stringer"] <- "DG_3"

# Overwrite the original column so all downstream plots use the anonymous names
sobj$Final_Combined_ID <- factor(sobj$Anonymous_ID)

# Re-assign Clinical Groups
sobj$Clinical_Group <- "Exclude"
sobj$Clinical_Group[grepl("PRE", sobj$Final_Combined_ID)] <- "PRE"
sobj$Clinical_Group[grepl("POST", sobj$Final_Combined_ID)] <- "POST"
sobj$Clinical_Group[grepl("DG_", sobj$Final_Combined_ID)] <- "22q11"

# ==============================================================================
# 4. DEFINE TARGET SUBSETS & REMOVE SPIKE
# ==============================================================================
cd8_labels <- c("CD8_CM", "CD8_EM", "CD8_Naive", "CD8_TEMRA")
cd4_labels <- c("CD4_CM", "CD4_EM", "CD4_Naive", "Treg")
thymus_ids <- c("libT1", "libT2", "libT3")

# Filter out Spike instantly
message(">>> Removing Spike cells...")
sobj_clean <- subset(sobj, subset = Final_Combined_ID != "Spike")

# A. THYMUS
message("\n>>> Subsetting Thymus cells...")
sobj_thymus <- subset(sobj_clean, subset = Final_Combined_ID %in% thymus_ids)
message("Thymus cells: ", ncol(sobj_thymus))

# B. CD8 (Blood only)
message("\n>>> Subsetting CD8 cells...")
sobj_cd8 <- subset(sobj_clean, subset = starCAT_label %in% cd8_labels & !(Final_Combined_ID %in% thymus_ids))
message("CD8 cells: ", ncol(sobj_cd8))

# C. CD4 (Blood only)
message("\n>>> Subsetting CD4 cells...")
sobj_cd4 <- subset(sobj_clean, subset = starCAT_label %in% cd4_labels & !(Final_Combined_ID %in% thymus_ids))
message("CD4 cells: ", ncol(sobj_cd4))

# Free up memory by deleting the massive starting object
rm(sobj, sobj_clean)
gc()

# ==============================================================================
# 5. STANDARD REPROCESSING FUNCTION (Classical PCA)
# ==============================================================================
reprocess_clean_subset <- function(obj, name_prefix) {
  message(paste0("\n>>> Re-processing ", name_prefix, " (Classical PCA Mode)..."))
  
  DefaultAssay(obj) <- "RNA"
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "dispersion", nfeatures = 2000, verbose = FALSE)
  
  if ("percent.mt" %in% colnames(obj@meta.data)) {
    obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  } else {
    obj <- ScaleData(obj, verbose = FALSE)
  }
  
  obj <- RunPCA(obj, npcs = 25, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:25, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  
  file_name <- paste0("data/", name_prefix, "_TargetCohorts_Analyzed.qs2")
  qs_save(obj, file_name)
  message(">>> Successfully saved: ", file_name)
  
  return(obj)
}

# ==============================================================================
# 6. EXECUTE REPROCESSING
# ==============================================================================
sobj_cd8_analyzed <- reprocess_clean_subset(sobj_cd8, "CD8")
sobj_cd4_analyzed <- reprocess_clean_subset(sobj_cd4, "CD4")
sobj_thymus_analyzed <- reprocess_clean_subset(sobj_thymus, "Thymus")

message("\n>>> All subpopulations successfully generated, integrated with VDJ, anonymized, and analyzed!")