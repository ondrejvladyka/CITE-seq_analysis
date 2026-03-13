# ==============================================================================
# SCRIPT: 07_Clonal_Tracking.R
# PURPOSE: Track specific TCR clonotypes globally and across target states
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggalluvial) 

message("\n>>> Running Clonal Tracking Analysis...")

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

if("CTgene" %in% colnames(sobj@meta.data)) {
  
  # Extract metadata and filter ONLY for cells that have a valid TCR
  tcr_meta_full <- sobj@meta.data %>% filter(!is.na(CTgene))
  
  if(nrow(tcr_meta_full) > 0) {
    
    # ---------------------------------------------------------
    # PLOT 1: GLOBAL CLONAL FLOW (All Clusters)
    # ---------------------------------------------------------
    message(">>> Generating Global Clonal Alluvial Plot...")
    
    # Find Top 50 clones across the ENTIRE CD8 compartment
    top_clones_global <- tcr_meta_full %>% 
      count(CTgene) %>% 
      top_n(50, n) %>% 
      pull(CTgene)
    
    # Format the data for the global river plot
    alluvial_data_global <- tcr_meta_full %>%
      filter(CTgene %in% top_clones_global) %>%
      group_by(CTgene, Clinical_Group, seurat_clusters) %>%
      summarise(Freq = n(), .groups = "drop")
    
    p_alluvial_global <- ggplot(alluvial_data_global,
                                aes(y = Freq, axis1 = Clinical_Group, axis2 = seurat_clusters)) +
      geom_alluvium(aes(fill = CTgene), width = 1/12, alpha = 0.8, color = "white", linewidth = 0.2) +
      geom_stratum(width = 1/4, fill = "grey20", color = "white") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "white", fontface = "bold", size = 3) +
      scale_x_discrete(limits = c("Clinical Cohort", "All Seurat Clusters"), expand = c(.05, .05)) +
      labs(title = paste(TARGET_CELL_TYPE, "- Global Clonal Architecture"),
           subtitle = "Tracking the Top 50 overall expanded clonotypes across all states",
           y = "Number of Cells") +
      theme_minimal() +
      theme(legend.position = "none", 
            axis.text.x = element_text(size = 12, face = "bold", color = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    
    global_file <- paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Global_Clonal_Alluvial_Top50.png")
    ggsave(global_file, plot = p_alluvial_global, width = 10, height = 12, dpi = PUB_DPI)
    
    # Save a CSV map for the global clones
    write.csv(alluvial_data_global %>% arrange(desc(Freq)), 
              paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Global_Top50_Clones_Map.csv"), 
              row.names = FALSE)
    
    # ---------------------------------------------------------
    # PLOT 2: TARGETED CLONAL FLOW (Group A vs Group B)
    # ---------------------------------------------------------
    message(">>> Generating Targeted Clonal Alluvial Plot (", COMPARISON_PREFIX, ")...")
    
    # Subset to the exact targeted clusters
    tcr_meta_sub <- tcr_meta_full %>% filter(seurat_clusters %in% c(CLUSTER_A, CLUSTER_B))
    
    if(nrow(tcr_meta_sub) > 0) {
      tcr_meta_sub$Target_Group <- "Unknown"
      tcr_meta_sub$Target_Group[tcr_meta_sub$seurat_clusters %in% CLUSTER_A] <- paste0("C", NAME_A)
      tcr_meta_sub$Target_Group[tcr_meta_sub$seurat_clusters %in% CLUSTER_B] <- paste0("C", NAME_B)
      tcr_meta_sub$Target_Group <- factor(tcr_meta_sub$Target_Group, levels = c(paste0("C", NAME_A), paste0("C", NAME_B)))
      
      # Find the Top 50 most massive clones IN THESE SPECIFIC CLUSTERS
      top_clones_target <- tcr_meta_sub %>% 
        count(CTgene) %>% 
        top_n(50, n) %>% 
        pull(CTgene)
      
      alluvial_data_target <- tcr_meta_sub %>%
        filter(CTgene %in% top_clones_target) %>%
        group_by(CTgene, Clinical_Group, Target_Group) %>%
        summarise(Freq = n(), .groups = "drop")
      
      p_alluvial_target <- ggplot(alluvial_data_target,
                                  aes(y = Freq, axis1 = Clinical_Group, axis2 = Target_Group)) +
        geom_alluvium(aes(fill = CTgene), width = 1/12, alpha = 0.8, color = "white", linewidth = 0.2) +
        geom_stratum(width = 1/4, fill = "grey20", color = "white") +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "white", fontface = "bold", size = 4) +
        scale_x_discrete(limits = c("Clinical Cohort", "Biological State"), expand = c(.05, .05)) +
        labs(title = "Targeted Clonal Architecture & Phenotypic Flow",
             subtitle = paste("Tracking the Top 50 expanded clonotypes in", paste0("C", NAME_A), "and", paste0("C", NAME_B)),
             y = "Number of Cells") +
        theme_minimal() +
        theme(legend.position = "none",
              axis.text.x = element_text(size = 12, face = "bold", color = "black"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank())
      
      target_file <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Targeted_Clonal_Alluvial_Top50.png")
      ggsave(target_file, plot = p_alluvial_target, width = 8, height = 8, dpi = PUB_DPI)
      
      write.csv(alluvial_data_target %>% arrange(desc(Freq)), 
                paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Targeted_Top50_Clones_Map.csv"), 
                row.names = FALSE)
      
      message(">>> Clonal tracking saved successfully to ", OUTPUT_DIR)
      
    } else {
      message(">>> WARNING: Target clusters do not contain any mapped TCRs.")
    }
  } else {
    message(">>> WARNING: No TCRs mapped in the entire dataset.")
  }
} else {
  message(">>> WARNING: 'CTgene' not found. Ensure scRepertoire integration ran successfully.")
}
