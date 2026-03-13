# ==============================================================================
# SCRIPT: 05_Targeted_Expression.R
# PURPOSE: Statistical Violin Plots for Key Driver Genes
# ==============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

# You may need to run: install.packages("ggpubr")
library(ggpubr) 

message("\n>>> Generating Targeted Expression Violins for: ", COMPARISON_PREFIX)

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

# ---------------------------------------------------------
# SMART SUBSETTING
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

genes_to_plot <- intersect(EXTRA_GENES_TO_LABEL, rownames(sobj_sub))

if(length(genes_to_plot) > 0) {
  
  expr_data <- FetchData(sobj_sub, vars = c(genes_to_plot, "Target_Group", "Clinical_Group"))
  expr_long <- pivot_longer(expr_data, cols = all_of(genes_to_plot), names_to = "Gene", values_to = "Expression")
  
  # Set the dynamic colors for the Violins
  target_pal <- setNames(c(HEATMAP_HIGH, HEATMAP_LOW), c(group_A_label, group_B_label))
  
  # ... (Inside your PLOT 1 ggviolin code, change palette to: palette = target_pal) ...
  # ... (And change stat_compare_means to: comparisons = list(c(group_A_label, group_B_label))) ...
  
  
  message(">>> Plotting Cluster Comparison Violins...")
  
  p_viol_state <- ggviolin(expr_long, 
                           x = "Target_Group", 
                           y = "Expression", 
                           fill = "Target_Group",
                           palette = target_pal, 
                           add = "boxplot", 
                           add.params = list(fill = "white", width = 0.1, size = 0.2),
                           trim = TRUE) +
    facet_wrap(~Gene, scales = "free_y", ncol = 4) +
    
    # Automatically calculates Wilcoxon p-values and draws the significance brackets
    stat_compare_means(comparisons = list(c(group_A_label, group_B_label)), 
                       label = "p.signif", # Uses stars (***) instead of messy numbers
                       method = "wilcox.test") +
    
    labs(title = "Key Gene Drivers: State Comparison", 
         subtitle = paste("Expression distribution across targeted clusters"),
         x = "", y = "Normalized Expression") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(size = 12, face = "bold", color = "black"),
          strip.background = element_rect(fill = "grey90", color = "black"),
          axis.text.x = element_text(size = 12, face = "bold"))
  
  file_viol_state <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Violin_State_Comparison.png")
  ggsave(file_viol_state, plot = p_viol_state, width = 10, height = 5, dpi = PUB_DPI)
  
  # ---------------------------------------------------------
  # PLOT 2: Clinical Shift (Inside EACH targeted cluster)
  # ---------------------------------------------------------
  message(">>> Plotting Clinical Shift Violins for both groups...")
  
  # Define the exact comparisons we want brackets for
  clin_comparisons <- list(c("PRE", "POST"), c("POST", "22q11"), c("PRE", "22q11"))
  
  # Loop through both NAME_A and NAME_B to generate plots for both sides
  for (target_name in c(NAME_A, NAME_B)) {
    
    # Filter only to the current target group
    expr_long_target <- expr_long %>% filter(Target_Group == paste0("C", target_name))
    
    if(nrow(expr_long_target) > 0) {
      p_viol_clin <- ggviolin(expr_long_target, 
                              x = "Clinical_Group", 
                              y = "Expression", 
                              fill = "Clinical_Group",
                              palette = GROUP_COLORS, # Uses your Darjeeling Teal/Gold/Red!
                              add = "boxplot", 
                              add.params = list(fill = "white", width = 0.1, size = 0.2),
                              trim = TRUE) +
        facet_wrap(~Gene, scales = "free_y", ncol = 4) +
        
        # Draw the three statistical brackets
        stat_compare_means(comparisons = clin_comparisons, 
                           label = "p.signif", 
                           step.increase = 0.1) + # Adds vertical spacing so brackets don't overlap
        
        labs(title = paste0("Clinical Shifts within C", target_name), 
             subtitle = paste("How does the cohort impact C", target_name, "?"),
             x = "", y = "Normalized Expression") +
        theme_minimal() +
        theme(legend.position = "none",
              strip.text = element_text(size = 12, face = "bold", color = "black"),
              strip.background = element_rect(fill = "grey90", color = "black"),
              axis.text.x = element_text(size = 10, face = "bold"))
      
      # Save with the specific cluster name in the file title
      file_viol_clin <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Violin_Clinical_Shift_C", target_name, ".png")
      ggsave(file_viol_clin, plot = p_viol_clin, width = 12, height = 5, dpi = PUB_DPI)
    }
  }
}
    