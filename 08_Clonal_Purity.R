# ==============================================================================
# SCRIPT: 08_Clonal_Purity.R
# PURPOSE: Quantify the cluster-restriction of expanded TCRs (Purity Heatmap)
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

message("\n>>> Running Clonal Purity Analysis...")

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

if("CTgene" %in% colnames(sobj@meta.data)) {
  
  tcr_meta <- sobj@meta.data %>% filter(!is.na(CTgene))
  
  if(nrow(tcr_meta) > 0) {
    
    # 1. Find the Top 20 absolute largest clones in the entire CD8 compartment
    top_clones <- tcr_meta %>% 
      count(CTgene) %>% 
      top_n(20, n) %>% 
      pull(CTgene)
    
    # --- NEW: Calculate the Dominant Clinical Group for each clone ---
    # (Which cohort does this clone primarily belong to?)
    dominant_clinical <- tcr_meta %>%
      filter(CTgene %in% top_clones) %>%
      group_by(CTgene, Clinical_Group) %>%
      summarise(Clin_Count = n(), .groups = "drop") %>%
      group_by(CTgene) %>%
      slice_max(order_by = Clin_Count, n = 1, with_ties = FALSE) %>%
      select(CTgene, Dominant_Clinical = Clinical_Group)
    
    # 2. Calculate the "Purity" (What % of each clone lives in each cluster?)
    purity_data <- tcr_meta %>%
      filter(CTgene %in% top_clones) %>%
      group_by(CTgene, seurat_clusters) %>%
      summarise(Cell_Count = n(), .groups = "drop") %>%
      group_by(CTgene) %>%
      mutate(Total_Clone_Size = sum(Cell_Count),
             Percent_in_Cluster = (Cell_Count / Total_Clone_Size) * 100) %>%
      ungroup()
    
    # 3. Merge the dominant clinical info into our purity data
    purity_data <- left_join(purity_data, dominant_clinical, by = "CTgene")
    
    # Fill in the missing 0% combinations so the heatmap is a complete grid
    purity_matrix <- purity_data %>%
      complete(CTgene, seurat_clusters, fill = list(Percent_in_Cluster = 0, Cell_Count = 0)) %>%
      select(-Dominant_Clinical) %>% # Remove to avoid NAs during complete
      left_join(dominant_clinical, by = "CTgene") # Re-attach so every row has the clinical label
    
    # 4. Order the clones nicely for the heatmap
    # We sort them FIRST by Clinical Group, then by Cluster, then by Size
    dominant_clusters <- purity_data %>%
      group_by(CTgene) %>%
      slice_max(order_by = Percent_in_Cluster, n = 1, with_ties = FALSE) %>%
      arrange(Dominant_Clinical, seurat_clusters, desc(Total_Clone_Size))
    
    purity_matrix$CTgene <- factor(purity_matrix$CTgene, levels = dominant_clusters$CTgene)
    
    # ---------------------------------------------------------
    # PLOT: The Purity Heatmap (With Clinical Secondary Axis)
    # ---------------------------------------------------------
    message(">>> Plotting Clonal Purity Heatmap with Clinical Axis...")
    
    p_purity <- ggplot(purity_matrix, aes(x = seurat_clusters, y = CTgene, fill = Percent_in_Cluster)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = ifelse(Percent_in_Cluster > 0, paste0(round(Percent_in_Cluster, 0), "%"), "")), 
                color = "black", size = 3, fontface = "bold") +
      
      scale_fill_gradient(low = "grey95", high = HEATMAP_HIGH, name = "% of Clone\nin Cluster") +
      
      # --- NEW: Create the Secondary Y-Axis using Facets ---
      # This groups the Y-axis by clinical cohort and puts a label block on the right side
      facet_grid(Dominant_Clinical ~ ., scales = "free_y", space = "free_y") +
      
      labs(title = "Clonal Purity & Origin",
           subtitle = "Cluster distribution and dominant clinical origin of Top 20 clones",
           x = "Seurat Cluster", y = "Top 20 Expanded Clonotypes") +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 12, face = "bold"),
            axis.text.y = element_text(size = 8, face = "bold"),
            
            # Format the secondary axis labels (the facet strips) to be horizontal and readable
            strip.text.y = element_text(size = 10, face = "bold", angle = 0), 
            strip.background = element_rect(fill = "grey90", color = "black"),
            
            # Add a tiny gap between the clinical blocks
            panel.spacing = unit(0.2, "lines"), 
            panel.grid = element_blank())
    
    purity_file <- paste0(OUTPUT_DIR, TARGET_CELL_TYPE, "_Clonal_Purity_Heatmap.png")
    ggsave(purity_file, plot = p_purity, width = 11, height = 8, dpi = PUB_DPI)
    
    message(">>> Clonal Purity Heatmap saved successfully to ", OUTPUT_DIR)
    
  } else {
    message(">>> WARNING: No TCRs mapped in the dataset.")
  }
} else {
  message(">>> WARNING: 'CTgene' not found. Ensure scRepertoire integration ran successfully.")
}