library(Seurat)
library(qs2)

# Load your GitHub data
sobj <- qs_read("data/Klocperk_analysed.qs2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scRepertoire")

library(scRepertoire)
library(Seurat)

# 1. Load VDJ annotation files
vdj_lib1 <- read.csv("data/filtered_contig_annotations_lib1.csv")
vdj_lib2 <- read.csv("data/filtered_contig_annotations_lib2.csv")
vdj_lib3 <- read.csv("data/filtered_contig_annotations_lib3.csv")
vdj_lib4 <- read.csv("data/filtered_contig_annotations_lib4.csv")






# 1. Manually format the raw barcodes to perfectly match Seurat
# Adds "libX_" to the front, and strips "-1" from the back
vdj_lib1$barcode <- paste0("lib1_", gsub("-1$", "", vdj_lib1$barcode))
vdj_lib2$barcode <- paste0("lib2_", gsub("-1$", "", vdj_lib2$barcode))
vdj_lib3$barcode <- paste0("lib3_", gsub("-1$", "", vdj_lib3$barcode))
vdj_lib4$barcode <- paste0("lib4_", gsub("-1$", "", vdj_lib4$barcode))

# 1. Define the prefixes in the exact order of your vdj_list (lib1, lib2, lib3, lib4)
library_prefixes <- c("lib1_", "lib2_", "lib3_", "lib4_")

# 2. Surgically fix the barcodes in the combined_tcr list
for (i in 1:length(combined_tcr)) {
  # Strip the "-1" from the end
  clean_barcodes <- sub("-1$", "", combined_tcr[[i]]$barcode)
  
  # Add the "libX_" prefix to the beginning
  combined_tcr[[i]]$barcode <- paste0(library_prefixes[i], clean_barcodes)
}

# 3. Verify the fix worked (This should now print "lib1_AAACCAAAGATGAGGT")
message(">>> FIXED VDJ BARCODE EXAMPLE:")
print(head(combined_tcr[[1]]$barcode, 1))

# 4. Integrate!
sobj <- combineExpression(combined_tcr, 
                          sobj, 
                          cloneCall = "gene", 
                          proportion = TRUE)

# 5. The moment of truth:
message(">>> Cells with matched TCR data: ", sum(!is.na(sobj$CTgene)))






# Create a clean plotting column
sobj$Has_TCR <- ifelse(!is.na(sobj$CTgene), "Has VDJ", "No VDJ")

# Generate the UMAP
p_sanity <- DimPlot(sobj, group.by = "Has_TCR", reduction = "umap") +
  scale_color_manual(values = c("Has VDJ" = "#D55E00", "No VDJ" = "grey85")) +
  ggtitle("VDJ Integration Success", subtitle = "Cells with paired TCR Alpha/Beta chains") +
  theme_minimal()

# Save and print
ggsave("VDJ_SanityCheck_UMAP.png", plot = p_sanity, width = 8, height = 6, dpi = 300)
print(p_sanity)

























library(dplyr)
library(tidyr)

# 1. Global Mapping Summary
mapped_cells <- sum(!is.na(sobj$CTgene))
total_cells <- ncol(sobj)
message(sprintf(">>> GLOBAL SUMMARY: %d out of %d cells (%.1f%%) have TCR data.", 
                mapped_cells, total_cells, (mapped_cells/total_cells)*100))

# 2. Breakdown by Library
message("\n>>> MAPPING BY LIBRARY:")
print(table(Library = sobj$library_id, Has_TCR = !is.na(sobj$CTgene)))

message("\n>>> MAPPING BY CELL TYPE (starCAT):")
# Changed starCAT_plot to starCAT_label
tcr_by_type <- as.data.frame(table(CellType = sobj$starCAT_label, Has_TCR = !is.na(sobj$CTgene)))

tcr_summary <- tcr_by_type %>%
  pivot_wider(names_from = Has_TCR, values_from = Freq, names_prefix = "TCR_") %>%
  mutate(
    TCR_TRUE = replace_na(TCR_TRUE, 0),
    TCR_FALSE = replace_na(TCR_FALSE, 0),
    Total = TCR_TRUE + TCR_FALSE,
    Percent_Mapped = round((TCR_TRUE / Total) * 100, 1)
  ) %>%
  arrange(desc(Percent_Mapped))

print(tcr_summary)






sobj$Has_TCR <- ifelse(!is.na(sobj$CTgene), "Has VDJ", "No VDJ")

# Generate the UMAP
p_sanity <- DimPlot(sobj, group.by = "Has_TCR", reduction = "umap") +
  scale_color_manual(values = c("Has VDJ" = "#D55E00", "No VDJ" = "grey85")) +
  ggtitle("VDJ Integration Success", subtitle = "Cells with paired TCR Alpha/Beta chains") +
  theme_minimal()

# Save and print
ggsave("VDJ_SanityCheck_UMAP.png", plot = p_sanity, width = 8, height = 6, dpi = 300)
print(p_sanity)

