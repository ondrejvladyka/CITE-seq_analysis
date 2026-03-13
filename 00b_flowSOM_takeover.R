library(FlowSOM)
library(Seurat)
library(qs2)

message(">>> Loading existing GLM-PCA object...")
sobj <- qs_read("data/CD4_TargetCohorts_GLMPCA.qs2")

message("\n>>> HOSTILE TAKEOVER: Replacing Seurat clusters with FlowSOM Metaclusters...")
# 1. Extract the math we already calculated
pca_matrix <- Embeddings(sobj, "glmpca")[, 1:20]

# 2. Match the number of clusters
num_clusters <- length(unique(sobj$seurat_clusters))

# 3. Build FlowSOM (Takes a few seconds)
set.seed(42)
fSOM <- ReadInput(pca_matrix)
fSOM <- BuildSOM(fSOM, colsToUse = 1:20)
fSOM <- BuildMST(fSOM)
meta_clusters <- metaClustering_consensus(fSOM$map$codes, k = num_clusters)
cell_assignments <- meta_clusters[fSOM$map$mapping[, 1]]

# 4. Hijack the Seurat identity!
sobj$FlowSOM_clusters <- factor(cell_assignments) 
sobj$seurat_clusters <- factor(cell_assignments)  
Idents(sobj) <- sobj$seurat_clusters               

# 5. Overwrite the save file
qs_save(sobj, "data/CD4_TargetCohorts_GLMPCA.qs2")
message(">>> Takeover complete and saved! You just saved 45 minutes.")