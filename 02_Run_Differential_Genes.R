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
# 2. Volcano Plot Configuration with Universal Colors
# ---------------------------------------------------------
num_clusters <- length(levels(sobj$seurat_clusters))
CLUSTER_COLORS <- wes_palette(WES_PALETTE_NAME, num_clusters, type = "continuous")
names(CLUSTER_COLORS) <- levels(sobj$seurat_clusters)

markers_sorted$Significance <- "Not Significant"

# USE NAME_A and NAME_B for the text labels
markers_sorted$Significance[markers_sorted$avg_log2FC > 0.25 & markers_sorted$p_val_adj < 0.05] <- paste0("Up in C", NAME_A)
markers_sorted$Significance[markers_sorted$avg_log2FC < -0.25 & markers_sorted$p_val_adj < 0.05] <- paste0("Up in C", NAME_B)

top_A <- markers_sorted %>% filter(Significance == paste0("Up in C", NAME_A)) %>% slice_max(order_by = avg_log2FC, n = 25)
top_B <- markers_sorted %>% filter(Significance == paste0("Up in C", NAME_B)) %>% slice_min(order_by = avg_log2FC, n = 25)

markers_sorted$Label <- ""
genes_to_label <- unique(c(top_A$gene, top_B$gene, EXTRA_GENES_TO_LABEL))
markers_sorted$Label[markers_sorted$gene %in% genes_to_label] <- markers_sorted$gene[markers_sorted$gene %in% genes_to_label]

# USE CLUSTER_A[1] so it always successfully picks exactly one color, even if it's a vector!
color_mapping <- setNames(
  c("grey80", CLUSTER_COLORS[as.character(CLUSTER_A[1])], CLUSTER_COLORS[as.character(CLUSTER_B[1])]), 
  c("Not Significant", paste0("Up in C", NAME_A), paste0("Up in C", NAME_B))
)

# ---------------------------------------------------------
# 2. Volcano Plot Configuration (Optimized & Bug-Free)
# ---------------------------------------------------------

# --- BUG FIX 1: Prevent the "Infinity" Decapitation ---
# Find the smallest p-value that isn't exactly zero
min_valid_p <- min(markers_sorted$p_val_adj[markers_sorted$p_val_adj > 0], na.rm = TRUE)
# Replace exact 0s with a number slightly smaller than the minimum so they sit at the very top of the peak
markers_sorted$p_val_adj[markers_sorted$p_val_adj == 0] <- min_valid_p * 1e-10

markers_sorted$Significance <- "Not Significant"
markers_sorted$Significance[markers_sorted$avg_log2FC > 0.25 & markers_sorted$p_val_adj < 0.05] <- paste0("Up in C", NAME_A)
markers_sorted$Significance[markers_sorted$avg_log2FC < -0.25 & markers_sorted$p_val_adj < 0.05] <- paste0("Up in C", NAME_B)

# Grab the top 12 genes so it's not a cluttered spiderweb
top_A <- markers_sorted %>% filter(Significance == paste0("Up in C", NAME_A)) %>% slice_max(order_by = avg_log2FC, n = 12)
top_B <- markers_sorted %>% filter(Significance == paste0("Up in C", NAME_B)) %>% slice_min(order_by = avg_log2FC, n = 12)

# --- BUG FIX 2: Use NA instead of "" so ggrepel ignores non-labeled dots ---
markers_sorted$Label <- NA 
genes_to_label <- unique(c(top_A$gene, top_B$gene, EXTRA_GENES_TO_LABEL))
markers_sorted$Label[markers_sorted$gene %in% genes_to_label] <- markers_sorted$gene[markers_sorted$gene %in% genes_to_label]

# Set colors (Using your Muted Red/Blue Heatmap Colors!)
color_mapping <- setNames(
  c("grey85", HEATMAP_HIGH, HEATMAP_LOW), 
  c("Not Significant", paste0("Up in C", NAME_A), paste0("Up in C", NAME_B))
)

p_volcano <- ggplot(markers_sorted, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.2) + 
  scale_color_manual(values = color_mapping) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", linewidth = 0.3) +
  
  geom_text_repel(aes(label = Label),
                  size = VOLCANO_LABEL_SIZE,  
                  fontface = "bold",          
                  color = "black",            
                  box.padding = 0.8,          
                  point.padding = 0.3,        
                  min.segment.length = 0,     
                  max.overlaps = 20,          
                  na.rm = TRUE,               # Tells ggplot to ignore the NA labels entirely
                  show.legend = FALSE) +      
  
  theme_minimal() +
  labs(title = paste(COMPARISON_PREFIX, "Transcriptomic Shift"),
       x = "Average Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme(legend.position = "bottom",           
        legend.title = element_blank())

volcano_file <- paste0(OUTPUT_DIR, COMPARISON_PREFIX, "_Volcano.png")


ggsave(volcano_file, plot = p_volcano, width = 10, height = 8, dpi = PUB_DPI)
message(">>> Saved Labeled Volcano Plot to: ", volcano_file)