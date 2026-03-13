# ==============================================================================
# SCRIPT: 02_Run_Differential_Genes.R
# PURPOSE: Full-spectrum DGE and Labeled Volcano Plot (Optimized & Cleaned)
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel) 

message(">>> Running Differential Expression: ", COMPARISON_PREFIX)

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

# ---------------------------------------------------------
# 1. SMART SUBSETTING: Pairwise vs One-vs-All
# ---------------------------------------------------------
if ("Rest" %in% CLUSTER_B) {
  message(">>> MODE: Target vs Rest of Population")
  sobj_sub <- sobj # DO NOT subset! Keep the whole object.
  sobj_sub$DGE_Group <- "Target_B" # Default everyone to "The Rest"
  sobj_sub$DGE_Group[sobj_sub$seurat_clusters %in% CLUSTER_A] <- "Target_A" # Overwrite our specific target
} else {
  message(">>> MODE: Pairwise Comparison")
  sobj_sub <- subset(sobj, seurat_clusters %in% c(CLUSTER_A, CLUSTER_B))
  sobj_sub$DGE_Group <- "Unknown"
  sobj_sub$DGE_Group[sobj_sub$seurat_clusters %in% CLUSTER_A] <- "Target_A"
  sobj_sub$DGE_Group[sobj_sub$seurat_clusters %in% CLUSTER_B] <- "Target_B"
}

# Apply the NEW Super Groups
Idents(sobj_sub) <- "DGE_Group"

# Count the cells
n_A <- sum(Idents(sobj_sub) == "Target_A")
n_B <- sum(Idents(sobj_sub) == "Target_B")

# Calculate a perfectly balanced cap based on the smaller group (Max 500)
balanced_cap <- min(n_A, n_B, 2000)

message(">>> Group A Cells: ", n_A, " | Group B Cells (Background): ", n_B)
message(">>> Balancing statistical power! Downsampling to exactly ", balanced_cap, " cells per group.")

# --- THE FIX: Use the NEW names in FindMarkers! ---
markers_full <- FindMarkers(sobj_sub, 
                            ident.1 = "Target_A", 
                            ident.2 = "Target_B",
                            test.use = "t",
                            only.pos = FALSE, 
                            min.pct = 0.1, 
                            logfc.threshold = 0.1,  # Cut out the 0.0 logFC noise!
                            max.cells.per.ident = balanced_cap)

markers_full$gene <- rownames(markers_full)
markers_sorted <- markers_full %>% arrange(desc(avg_log2FC))

# Save the full CSV
write.csv(markers_sorted, DEG_CSV_FILE, row.names = FALSE)
message(">>> Saved full ranked list to: ", DEG_CSV_FILE)


# ---------------------------------------------------------
# 2. Volcano Plot Configuration (Clean & Bug-Free)
# ---------------------------------------------------------

# --- BUG FIX 1: Prevent the "Infinity" Decapitation ---
min_valid_p <- min(markers_sorted$p_val_adj[markers_sorted$p_val_adj > 0], na.rm = TRUE)
markers_sorted$p_val_adj[markers_sorted$p_val_adj == 0] <- min_valid_p * 1e-10

# Assign Significance
markers_sorted$Significance <- "Not Significant"
markers_sorted$Significance[markers_sorted$avg_log2FC > 0.25 & markers_sorted$p_val_adj < 0.05] <- paste0("Up in C", NAME_A)
markers_sorted$Significance[markers_sorted$avg_log2FC < -0.25 & markers_sorted$p_val_adj < 0.05] <- paste0("Up in C", NAME_B)

# Grab the top 12 genes
top_A <- markers_sorted %>% filter(Significance == paste0("Up in C", NAME_A)) %>% slice_max(order_by = avg_log2FC, n = 12)
top_B <- markers_sorted %>% filter(Significance == paste0("Up in C", NAME_B)) %>% slice_min(order_by = avg_log2FC, n = 12)

# --- BUG FIX 2: Clean Labels using NA ---
markers_sorted$Label <- NA 
genes_to_label <- unique(c(top_A$gene, top_B$gene, EXTRA_GENES_TO_LABEL))
markers_sorted$Label[markers_sorted$gene %in% genes_to_label] <- markers_sorted$gene[markers_sorted$gene %in% genes_to_label]

# Set colors (Using your beautiful Darjeeling Heatmap Colors)
color_mapping <- setNames(
  c("grey85", HEATMAP_HIGH, HEATMAP_LOW), 
  c("Not Significant", paste0("Up in C", NAME_A), paste0("Up in C", NAME_B))
)

p_volcano <- ggplot(markers_sorted, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.2) + 
  scale_color_manual(values = color_mapping) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", linewidth = 0.3) +
  
  # --- NEW: THE Y-AXIS SQUEEZE ---
  # This physically squishes the massive values at the top, giving the bottom room to breathe!
  scale_y_sqrt(breaks = c(0, 10, 25, 75, 200)) + 
  
  geom_text_repel(aes(label = Label),
                  size = VOLCANO_LABEL_SIZE,  
                  fontface = "bold",          
                  color = "black",            
                  box.padding = 0.8,          
                  point.padding = 0.3,        
                  min.segment.length = 0,     
                  max.overlaps = 20,          
                  na.rm = TRUE,               
                  show.legend = FALSE) +      
  
  theme_minimal() +
  labs(title = paste(COMPARISON_PREFIX, "Transcriptomic Shift"),
       x = "Average Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme(legend.position = "bottom",           
        legend.title = element_blank()) +
  
  coord_cartesian(xlim = c(-3, 4))

volcano_file <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Volcano.png")
ggsave(volcano_file, plot = p_volcano, width = 10, height = 8, dpi = PUB_DPI)
message(">>> Saved Labeled Volcano Plot to: ", volcano_file)
