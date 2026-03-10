# ==============================================================================
# SCRIPT: 02_GSEA_Analysis.R
# PURPOSE: Run Gene Set Enrichment Analysis on Cluster 0 vs 1
# ==============================================================================

# 1. Install required packages if you don't have them
#install.packages("BiocManager")
#BiocManager::install(c("fgsea", "msigdbr"))

# 2. Load Libraries
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)

# 3. Load your full-spectrum ranked list
message(">>> Loading Ranked Gene List...")
deg_data <- read.csv("CD8_Cluster0_vs_Cluster1_FullSpectrum.csv")

# 4. Create the Ranked Vector
# fgsea requires a named numeric vector sorted from highest to lowest.
# Positive values = enriched in Cluster 0 (PRE/POST)
# Negative values = enriched in Cluster 1 (22q11/POST)
ranked_genes <- deg_data$avg_log2FC
names(ranked_genes) <- deg_data$gene

# Sort strictly descending, handle any potential duplicates or NAs safely
ranked_genes <- na.omit(ranked_genes)
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# 5. Fetch the Hallmark Pathways for Humans
message(">>> Fetching Hallmark Pathways...")
m_df <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(x = m_df$gene_symbol, f = m_df$gs_name)

# 6. Run the GSEA Algorithm
message(">>> Running Fast GSEA...")
fgsea_res <- fgsea(pathways = pathways, 
                   stats = ranked_genes,
                   minSize = 15,   # Ignore tiny pathways
                   maxSize = 500)  # Ignore massively broad pathways

# 7. Clean up the results for plotting
# We filter for significant pathways (adjusted p-value < 0.05)
significant_pathways <- fgsea_res %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES)) # NES = Normalized Enrichment Score

# If the names are too long, remove the "HALLMARK_" prefix for cleaner plots
significant_pathways$pathway <- gsub("HALLMARK_", "", significant_pathways$pathway)

message(">>> Found ", nrow(significant_pathways), " significant pathways!")

# 8. Plot the Results
p_gsea <- ggplot(significant_pathways, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = NES > 0), color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"), 
                    labels = c("TRUE" = "Up in Cluster 0 (PRE/POST)", 
                               "FALSE" = "Up in Cluster 1 (22q11/POST)")) +
  labs(title = "GSEA: Hallmark Pathways",
       subtitle = "Cluster 0 vs Cluster 1",
       x = "Pathway",
       y = "Normalized Enrichment Score (NES)") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")



print(p_gsea)
ggsave("CD8_GSEA_Hallmark_0_vs_1.png", plot = p_gsea, width = 10, height = 8, dpi = 300)

# 9. Save the data to look at the exact genes driving the pathways
write.csv(apply(significant_pathways,2,as.character), "CD8_GSEA_Significant_Pathways.csv", row.names = FALSE)









