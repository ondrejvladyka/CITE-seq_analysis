# ==============================================================================
# SCRIPT: 06_Signature_Scoring.R
# PURPOSE: Calculate and visualize composite biological signature scores
# ==============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

message("\n>>> Running Signature Scoring Analysis...")

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

# 1. Calculate Scores for all defined signatures
valid_signatures <- c() # Keep track of which signatures actually worked

for (sig_name in names(SIGNATURES)) {
  genes <- SIGNATURES[[sig_name]]
  # Only use genes that actually exist in your data matrix
  valid_genes <- intersect(genes, rownames(sobj))
  
  if (length(valid_genes) > 0) {
    message(">>> Calculating score for: ", sig_name, " (", length(valid_genes), " valid genes)")
    
    # Seurat automatically appends a "1" to the end of the name (e.g., "Exhaustion1")
    sobj <- AddModuleScore(sobj, features = list(valid_genes), name = sig_name)
    
    # Rename the column back to the clean name
    colnames(sobj@meta.data)[colnames(sobj@meta.data) == paste0(sig_name, "1")] <- sig_name
    valid_signatures <- c(valid_signatures, sig_name)
    
  } else {
    message(">>> WARNING: No valid genes found for signature: ", sig_name)
  }
}

if(length(valid_signatures) > 0) {
  
  # ---------------------------------------------------------
  # PLOT 1: Global Landscape (Feature UMAPs)
  # ---------------------------------------------------------
  message(">>> Generating UMAP Signature Maps...")
  
  # Map the scores onto the UMAP using your Darjeeling gradient
  p_sig_umap <- FeaturePlot(sobj, features = valid_signatures, ncol = 2, pt.size = 0.5) &
    scale_color_gradientn(colors = c(HEATMAP_LOW, HEATMAP_MID, HEATMAP_HIGH), name = "Score") &
    theme(plot.title = element_text(face = "bold"))
  
  umap_file <- paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Signatures_UMAP.png")
  ggsave(umap_file, plot = p_sig_umap, width = 10, height = 5 * ceiling(length(valid_signatures)/2), dpi = PUB_DPI)
  
  # ---------------------------------------------------------
  # Setup for Targeted Violins (SMART SUBSETTING)
  # ---------------------------------------------------------
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
  
  expr_data <- FetchData(sobj_sub, vars = c(valid_signatures, "Target_Group", "Clinical_Group"))
  expr_long <- pivot_longer(expr_data, cols = all_of(valid_signatures), names_to = "Signature", values_to = "Score")
  
  target_pal <- setNames(c(HEATMAP_HIGH, HEATMAP_LOW), c(group_A_label, group_B_label))
  

  # ---------------------------------------------------------
  # PLOT 2: Targeted State Comparison (Cluster A vs Cluster B)
  # ---------------------------------------------------------
  
  target_pal <- setNames(c(HEATMAP_HIGH, HEATMAP_LOW), c(paste0("C", NAME_A), paste0("C", NAME_B)))
  
  message(">>> Generating Targeted Signature Violins...")
  
  p_sig_viol_state <- ggviolin(expr_long, 
                               x = "Target_Group", 
                               y = "Score", 
                               fill = "Target_Group",
                               palette = target_pal, # Red vs Blue
                               add = "boxplot", 
                               add.params = list(fill = "white", width = 0.1, size = 0.2),
                               trim = TRUE) +
    facet_wrap(~Signature, scales = "free_y", ncol = 2) +
    
    # --- THE FIX: Comma instead of parenthesis, and closed at the very end! ---
    stat_compare_means(comparisons = list(c(group_A_label, group_B_label)), 
                       label = "p.signif", 
                       method = "wilcox.test") +
    
    labs(title = "Composite Signature Scores: State Comparison", 
         x = "", y = "Module Score") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(size = 12, face = "bold", color = "black"),
          strip.background = element_rect(fill = "grey90", color = "black"))
  
  state_file <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Signatures_State_Comparison.png")
  ggsave(state_file, plot = p_sig_viol_state, width = 8, height = 5, dpi = PUB_DPI)
  
  # ---------------------------------------------------------
  # PLOT 3: Clinical Shift (Inside EACH targeted cluster)
  # ---------------------------------------------------------
  message(">>> Generating Clinical Shift Violins for both groups...")
  
  clin_comparisons <- list(c("PRE", "POST"), c("POST", "22q11"), c("PRE", "22q11"))
  
  # Loop through both NAME_A and NAME_B to generate plots for both sides
  for (target_name in c(NAME_A, NAME_B)) {
    
    # Filter only to the current target group
    expr_long_target <- expr_long %>% filter(Target_Group == paste0("C", target_name))
    
    if(nrow(expr_long_target) > 0) {
      p_sig_viol_clin <- ggviolin(expr_long_target, 
                                  x = "Clinical_Group", 
                                  y = "Score", 
                                  fill = "Clinical_Group",
                                  palette = GROUP_COLORS, # Darjeeling Teal, Gold, Red
                                  add = "boxplot", 
                                  add.params = list(fill = "white", width = 0.1, size = 0.2),
                                  trim = TRUE) +
        facet_wrap(~Signature, scales = "free_y", ncol = 2) +
        stat_compare_means(comparisons = clin_comparisons, label = "p.signif", step.increase = 0.1) + 
        labs(title = paste0("Signature Shifts within C", target_name), 
             x = "", y = "Module Score") +
        theme_minimal() +
        theme(legend.position = "none",
              strip.text = element_text(size = 12, face = "bold", color = "black"),
              strip.background = element_rect(fill = "grey90", color = "black"))
      
      # Save with the specific cluster name in the file title
      clin_file <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Signatures_Clinical_Shift_C", target_name, ".png")
      ggsave(clin_file, plot = p_sig_viol_clin, width = 8, height = 5, dpi = PUB_DPI)
    }
  }
  
  message(">>> Signature plots saved successfully to ", OUTPUT_DIR)
}