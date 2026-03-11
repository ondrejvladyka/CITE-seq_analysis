# ==============================================================================
# SCRIPT: 02_Run_Differential_Genes.R
# PURPOSE: Full-spectrum DGE and Labeled Volcano Plot
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)

# IMPORTANT: You may need to run install.packages("ggrepel") once in your console
library(ggrepel) 

message(">>> Running Differential Expression: ", COMPARISON_PREFIX)

if(!exists("sobj")) { sobj <- qs_read(INPUT_OBJECT_PATH) }

# 1. Subset and Run DGE
sobj_sub <- subset(sobj, idents = c(CLUSTER_A, CLUSTER_B))

markers_full <- FindMarkers(sobj_sub, ident.1 = CLUSTER_A, ident.2 = CLUSTER_B, 
                            only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.0) 

markers_full$gene <- rownames(markers_full)
markers_sorted <- markers_full %>% arrange(desc(avg_log2FC))

# Save the full CSV
write.csv(markers_sorted, DEG_CSV_FILE, row.names = FALSE)
message(">>> Saved full ranked list to: ", DEG_CSV_FILE)


# ---------------------------------------------------------
# 2. Volcano Plot Configuration with Labels
# ---------------------------------------------------------
markers_sorted$Significance <- "Not Significant"
markers_sorted$Significance[markers_sorted$avg_log2FC > 0.25 & markers_sorted$p_val_adj < 0.05] <- paste0("Up in C", CLUSTER_A)
markers_sorted$Significance[markers_sorted$avg_log2FC < -0.25 & markers_sorted$p_val_adj < 0.05] <- paste0("Up in C", CLUSTER_B)

# --- NEW LABELING LOGIC ---
# Grab the top 15 genes pushing Right (Cluster A) and Left (Cluster B)
top_A <- markers_sorted %>% 
  filter(Significance == paste0("Up in C", CLUSTER_A)) %>% 
  slice_max(order_by = avg_log2FC, n = 25)

top_B <- markers_sorted %>% 
  filter(Significance == paste0("Up in C", CLUSTER_B)) %>% 
  slice_min(order_by = avg_log2FC, n = 25)

# Create an empty label column, then fill it ONLY for the top genes
markers_sorted$Label <- ""
genes_to_label <- c(top_A$gene, top_B$gene)
markers_sorted$Label[markers_sorted$gene %in% genes_to_label] <- markers_sorted$gene[markers_sorted$gene %in% genes_to_label]

# Set colors dynamically
color_mapping <- setNames(c("grey80", "blue", "red"), 
                          c("Not Significant", paste0("Up in C", CLUSTER_A), paste0("Up in C", CLUSTER_B)))

# Build the labeled plot
p_volcano <- ggplot(markers_sorted, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = color_mapping) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
  
  # --- ADD THE LABELS ---
  geom_text_repel(aes(label = Label),
                  size = 4,                   # Font size
                  fontface = "bold",          # Make text bold
                  color = "black",            # Keep text black regardless of dot color
                  box.padding = 0.5,          # Space around text
                  max.overlaps = Inf,         # Force it to draw every label
                  show.legend = FALSE) +      # Don't put "a" in the legend
  
  theme_minimal() +
  labs(title = paste(COMPARISON_PREFIX, "Transcriptomic Shift"),
       x = "Average Log2 Fold Change", y = "-Log10 Adjusted P-value") 

# Save slightly wider so the labels don't get chopped off on the edges
volcano_file <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Volcano.png")
ggsave(volcano_file, plot = p_volcano, width = 10, height = 8, dpi = 300)
message(">>> Saved Labeled Volcano Plot to: ", volcano_file)