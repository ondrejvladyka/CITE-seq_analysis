# ==============================================================================
# SCRIPT: 04_Targeted_Clonality.R
# PURPOSE: Clonal Homeostasis strictly for the targeted clusters
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)

message("\n>>> Running Targeted Clonal Analysis for: ", COMPARISON_PREFIX)

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

# ---------------------------------------------------------
# SMART SUBSETTING & DYNAMIC LABELS
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

# Ensure VDJ data actually exists before plotting
if("cloneSize" %in% colnames(sobj_sub@meta.data)) { 
  
  clone_levels <- c("Single (1)", "Small (1e-04 < X <= 0.001)", "Medium (0.001 < X <= 0.01)", "Large (0.01 < X <= 0.1)", "Hyperexpanded (0.1 < X <= 1)")
  
  # --- THE SIMPLE FIX ---
  # 1. Convert to text
  sobj_sub$cloneSize <- as.character(sobj_sub$cloneSize)
  
  # 2. If a cell HAS a TCR, but its size label doesn't match your clone_levels list (e.g. it's NA or "Rare"), force it to "Single (1)"
  sobj_sub$cloneSize[!is.na(sobj_sub$CTgene) & !(sobj_sub$cloneSize %in% clone_levels)] <- "Single (1)"
  
  # 3. Factor it strictly using your exact levels list
  sobj_sub$cloneSize <- factor(sobj_sub$cloneSize, levels = clone_levels)
  
  # ---------------------------------------------------------
  # PLOT 1: Overall Clonality (Group A vs Group B)
  # ---------------------------------------------------------
  clone_comp_main <- as.data.frame(table(Size = sobj_sub$cloneSize, Group = sobj_sub$Target_Group)) %>%
    group_by(Group) %>% 
    filter(sum(Freq) > 0) %>% # Silences the "Removed rows" warning
    mutate(Proportion = Freq / sum(Freq)) %>% 
    ungroup()
  
  p_clone_main <- ggplot(clone_comp_main, aes(x = Group, y = Proportion, fill = Size)) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
    scale_fill_brewer(palette = "Reds", drop = FALSE) + # Ensures colors stay consistent even if a bin is empty
    scale_y_continuous(labels = scales::percent) +
    labs(title = paste("Clonal Expansion:", paste0("C", NAME_A), "vs", paste0("C", NAME_B)),
         y = "Proportion of Cells", x = "Targeted Clusters") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 14, face = "bold"))
  
  ggsave(paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Targeted_Clonality_Main.png"), plot = p_clone_main, width = 6, height = 6, dpi = PUB_DPI)
  
  # ---------------------------------------------------------
  # PLOT 2: Faceted by Clinical Group (The Deep Dive)
  # ---------------------------------------------------------
  clone_comp_clin <- as.data.frame(table(Size = sobj_sub$cloneSize, Target = sobj_sub$Target_Group, Clinical = sobj_sub$Clinical_Group)) %>%
    group_by(Target, Clinical) %>% 
    filter(sum(Freq) > 0) %>% # Silences the "Removed rows" warning
    mutate(Proportion = Freq / sum(Freq)) %>% 
    ungroup()
  
  p_clone_clin <- ggplot(clone_comp_clin, aes(x = Clinical, y = Proportion, fill = Size)) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
    facet_wrap(~Target) + 
    scale_fill_brewer(palette = "Reds", drop = FALSE) + 
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Targeted Clonal Expansion by Clinical Group",
         subtitle = "Comparing expansion profiles within the target clusters",
         y = "Proportion of Cells", x = "Clinical Cohort") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
          strip.text = element_text(size = 12, face = "bold"),
          strip.background = element_rect(fill = "grey90", color = "black"))
  
  ggsave(paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Targeted_Clonality_Clinical_Split.png"), plot = p_clone_clin, width = 8, height = 6, dpi = PUB_DPI)
  
  message(">>> Targeted Clonal plots saved to ", OUTPUT_DIR)
  
} else {
  message(">>> WARNING: 'cloneSize' not found in metadata. Skipping targeted clonal plots.")
}