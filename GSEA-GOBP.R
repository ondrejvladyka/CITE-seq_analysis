#
#
#GSEA - GO: BIological process


library(dplyr)
library(ggplot2)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)

# 2. Load and Prepare the Ranked List
message(">>> Loading Ranked Gene List...")
deg_data <- read.csv("CD8_Cluster0_vs_Cluster1_FullSpectrum.csv")

# Create the strictly sorted named vector
ranked_genes <- deg_data$avg_log2FC
names(ranked_genes) <- deg_data$gene
ranked_genes <- na.omit(ranked_genes)
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# 3. Fetch the Robust GO:BP Dataset (Gene Ontology: Biological Process)
message(">>> Fetching GO:BP Pathways...")
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, gene_symbol)

# Clean up the pathway names so they look nice on the plot (Remove "GOBP_" and underscores)
m_t2g$gs_name <- gsub("GOBP_", "", m_t2g$gs_name)
m_t2g$gs_name <- gsub("_", " ", m_t2g$gs_name)

# 4. Run clusterProfiler GSEA
message(">>> Running Advanced GSEA...")
gsea_res <- GSEA(geneList = ranked_genes, 
                 TERM2GENE = m_t2g, 
                 minGSSize = 15, 
                 maxGSSize = 500, 
                 pvalueCutoff = 0.05, 
                 verbose = FALSE)

message(">>> Found ", nrow(gsea_res@result %>% filter(p.adjust < 0.05)), " significant robust pathways!")

# 5. Calculate Similarity between pathways for the Network Plot
# This mathematically determines which pathways share the same genes
message(">>> Calculating Pathway Network Topology...")
gsea_res <- pairwise_termsim(gsea_res)

# ---------------------------------------------------------
# VISUALIZATION 1: The Network Plot (Enrichment Map)
# ---------------------------------------------------------
# Nodes are pathways. Lines connect pathways that share many genes.
# Color is NES (Red = Up in Cluster 0/PRE, Blue = Up in Cluster 1/22q11)
p_network <- emapplot(gsea_res, 
                      showCategory = 40,  # Show the top 40 pathways
                      color = "NES",      # Color by direction of fold-change
                      layout = "nicely",
                      node_label = "category") +
  ggtitle("GO:BP Pathway Network (Enrichment Map)",
          subtitle = "Red = Enriched in Cluster 0 (PRE/POST) | Blue = Enriched in Cluster 1 (22q11/POST)")

print(p_network)
ggsave("CD8_GSEA_NetworkPlot.png", plot = p_network, width = 28, height = 24, dpi = 300)

# ---------------------------------------------------------
# VISUALIZATION 2: The Classic Dot Plot
# ---------------------------------------------------------
# 1. Define the custom labels we want to use instead of the defaults
custom_labels <- c("activated" = "Upregulated in C0 (PRE/POST)", 
                   "suppressed" = "Upregulated in C1 (22q11/POST)")

# 2. Generate the plot and overwrite the facet labels
p_dot <- dotplot(gsea_res, 
                 showCategory = 30, 
                 split = ".sign", 
                 font.size = 10) + 
  facet_grid(.~.sign, labeller = as_labeller(custom_labels)) + # This line applies the new names
  ggtitle("Top Up/Down Regulated GO:BP Pathways") +
  theme(strip.text = element_text(size = 8, face = "bold", color = "black"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        axis.text.y = element_text(size = 8))



# 3. Print and Save
print(p_dot)
ggsave("CD8_GSEA_DotPlot_ClearLabels.png", plot = p_dot, width = 12, height = 8, dpi = 300)

# 1. Rebuild the plot with a much wider label format (e.g., 80 characters)
p_dot_wide_text <- dotplot(gsea_res, 
                           showCategory = 30, 
                           split = ".sign", 
                           label_format = 70, 
                           font.size = 8) + 
  facet_grid(.~.sign, labeller = as_labeller(custom_labels)) +
  ggtitle("Top Up/Down Regulated GO:BP Pathways") +
  theme(strip.text = element_text(size = 8, face = "bold", color = "black"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        axis.text.y = element_text(size = 8)) 

# 2. Save it. 
# Because the text is now physically wider, keeping the width at 9 or 10 
# will naturally squish the dot panels, giving you that "slimmer graph" look.
ggsave("CD8_GSEA_DotPlot_Widetext.png", plot = p_dot_wide_text, width = 10, height = 8, dpi = 300)

