library(Seurat)
library(HGNChelper)
library(anndata)

#Read in the Herring data
ad <- read_h5ad("ga22_cleaned_count_matrices.h5ad")

# Transpose the count matrix
ad_t <- t(as.matrix(ad$X))

# Create Seurat object
ga22 <- CreateSeuratObject(counts = ad_t, assay = "RNA")

#change rownames of meta data to cell id
rownames(ga22@meta.data) <- colnames(ga22)

#convert all factor vectors to character vectors
i <- sapply(ga22@meta.data, is.factor)
ga22@meta.data[i] <- lapply(ga22@meta.data[i], as.character)


# normalize data
ga22[["percent.mt"]] <- PercentageFeatureSet(ga22, pattern = "^MT-")
ga22 <- NormalizeData(ga22, normalization.method = "LogNormalize", scale.factor = 10000)
ga22 <- FindVariableFeatures(ga22, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
ga22 <- ScaleData(ga22, features = rownames(ga22))
ga22 <- RunPCA(ga22, features = VariableFeatures(object = ga22))

# Check number of PC components (I chose 15 PCs for downstream analysis, based on Elbow plot)
ElbowPlot(ga22)

# cluster and visualize
ga22 <- FindNeighbors(ga22, dims = 1:15)
ga22 <- FindClusters(ga22, resolution = 0.8)
ga22 <- RunUMAP(ga22, dims = 1:15)
DimPlot(ga22, reduction = "umap", label = TRUE)


#Assign Cell types using scType

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

#Use my gs_list file
db_ = "gs_listv4.xlsx";
tissue = "Brain" 
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

#Assign cell type to the clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = ga22[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 


# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(ga22@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ga22@meta.data[ga22@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ga22@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

write.csv(sctype_scores, file = "ga22_cluster_sctype.csv", row.names = TRUE)
#Overlay cell type assignments on UMAP
ga22@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  ga22@meta.data$customclassif[ga22@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

dev.new()
DimPlot(ga22, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'customclassif')        

###############################################################################
#ga24
ad <- read_h5ad("ga24_cleaned_count_matrices.h5ad")

# Transpose the count matrix
ad_t <- t(as.matrix(ad$X))

# Create Seurat object
ga24 <- CreateSeuratObject(counts = ad_t, assay = "RNA")

#change rownames of meta data to cell id
rownames(ga24@meta.data) <- colnames(ga24)
#convert all factor vectors to character vectors
i <- sapply(ga24@meta.data, is.factor)
ga24@meta.data[i] <- lapply(ga24@meta.data[i], as.character)

# normalize data
ga24[["percent.mt"]] <- PercentageFeatureSet(ga24, pattern = "^MT-")
ga24 <- NormalizeData(ga24, normalization.method = "LogNormalize", scale.factor = 10000)
ga24 <- FindVariableFeatures(ga24, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
ga24 <- ScaleData(ga24, features = rownames(ga24))
ga24 <- RunPCA(ga24, features = VariableFeatures(object = ga24))

# Check number of PC components (I chose 15 PCs for downstream analysis, based on Elbow plot)
ElbowPlot(ga24)

# cluster and visualize
ga24 <- FindNeighbors(ga24, dims = 1:15)
ga24 <- FindClusters(ga24, resolution = 0.8)
ga24 <- RunUMAP(ga24, dims = 1:15)
DimPlot(ga24, reduction = "umap", label = TRUE)


#Assign cell type to the clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = ga24[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 


# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(ga24@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ga24@meta.data[ga24@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ga24@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

write.csv(sctype_scores, file = "ga24_cluster_sctype.csv", row.names = TRUE)
#Overlay cell type assignments on UMAP
ga24@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  ga24@meta.data$customclassif[ga24@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

dev.new()

DimPlot(ga24, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'customclassif')        
#############################################################################
#ga34
ad <- read_h5ad("ga34_cleaned_count_matrices.h5ad")

# Transpose the count matrix
ad_t <- t(as.matrix(ad$X))

# Create Seurat object
ga34 <- CreateSeuratObject(counts = ad_t, assay = "RNA")

#change rownames of meta data to cell id
rownames(ga34@meta.data) <- colnames(ga34)
#convert all factor vectors to character vectors
i <- sapply(ga34@meta.data, is.factor)
ga34@meta.data[i] <- lapply(ga34@meta.data[i], as.character)

# normalize data
ga34[["percent.mt"]] <- PercentageFeatureSet(ga34, pattern = "^MT-")
ga34 <- NormalizeData(ga34, normalization.method = "LogNormalize", scale.factor = 10000)
ga34 <- FindVariableFeatures(ga34, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
ga34 <- ScaleData(ga34, features = rownames(ga34))
ga34 <- RunPCA(ga34, features = VariableFeatures(object = ga34))

# Check number of PC components (I chose 15 PCs for downstream analysis, based on Elbow plot)
ElbowPlot(ga34)

# cluster and visualize
ga34 <- FindNeighbors(ga34, dims = 1:15)
ga34 <- FindClusters(ga34, resolution = 0.8)
ga34 <- RunUMAP(ga34, dims = 1:15)
DimPlot(ga34, reduction = "umap")

#Assign cell type to the clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = ga34[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(ga34@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ga34@meta.data[ga34@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ga34@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

write.csv(sctype_scores, file = "ga34_cluster_sctype.csv", row.names = TRUE)
#Overlay cell type assignments on UMAP
ga34@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  ga34@meta.data$customclassif[ga34@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

dev.new()
DimPlot(ga34, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'customclassif')        
