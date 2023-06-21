#Converting h5ad file to Seurat object
#Keep metadata of interest
#Using Herring integrated datasets

library(Seurat)
library(ggplot2)
library(anndata)

#Use anndata to read in h5ad file
adata <- read_h5ad('herr_int_sc.h5ad')

#Turn counts into matrix and transpose
ad_t <- t(as.matrix(adata$X))

# Convert the count matrix to a Seurat object
int <- CreateSeuratObject(counts = ad_t, assay = "RNA")

#Add age metdata to Seurat object
int[["age"]] <- adata$obs$age

# Extract the "obsm" layer as a matrix of the calculated UMAP coordinates
obsm_m <- as.matrix(adata$obsm$X_umap)

#Add cell names to umap matrix
rownames(obsm_m) <- colnames(int)

#Add UMAP to Seurat reductions layer
int[["UMAP"]] <- CreateDimReducObject(embeddings = obsm_m, key = "UMAP")

#Plot UMAP
DimPlot(int, reduction = "UMAP", group.by = "age")