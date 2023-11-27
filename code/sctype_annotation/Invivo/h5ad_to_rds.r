# Convert h5ad file to rds
# Created by Fallon Ratner
# Example using Herring Day2 scRNAseq dataset

library(Seurat)
library(anndata)
library(reticulate)

#read in h5ad file
ad <- read_h5ad("2d_cleaned_count_matrices.h5ad")

# Transpose the count matrix
ad_t <- t(as.matrix(ad$X))

# Create Seurat object
d2 <- CreateSeuratObject(counts = ad_t, assay = "RNA")

#change rownames of meta data to cell id
rownames(d2@meta.data) <- colnames(d2)

#convert all factor vectors to character vectors
i <- sapply(d2@meta.data, is.factor)
d2@meta.data[i] <- lapply(d2@meta.data[i], as.character)

# save d2 as RDS file
saveRDS(d2, "day2.rds")