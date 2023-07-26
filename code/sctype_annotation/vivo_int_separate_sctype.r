#Created by Fallon Ratner on 24.07.23
library(Seurat)
library(dplyr)
library(HGNChelper)
library(ggplot2)
library(anndata)
library(reticulate)

###############################
setwd("C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/output/Integration/Vivo_Int_24_7")
#Convert h5ad to seurat obj
ad <- read_h5ad("pol_24_7.h5ad")
ad <- read_h5ad("liu_24_7.h5ad")
ad <- read_h5ad("herr22_24_7.h5ad")
ad <- read_h5ad("herr24_24_7.h5ad")
ad <- read_h5ad("herr34_24_7.h5ad")
ad <- read_h5ad("herr2d_24_7.h5ad")
ad <- read_h5ad("han11b1_24_7.h5ad")
ad <- read_h5ad("han11b2_24_7.h5ad")
ad <- read_h5ad("han12_24_7.h5ad")
ad <- read_h5ad("han13_24_7.h5ad")
ad <- read_h5ad("cou13_24_7.h5ad")
ad <- read_h5ad("cou17_24_7.h5ad")
ad <- read_h5ad("cou19_24_7.h5ad")
# Transpose the count matrix
ad_t <- t(as.matrix(ad$X))
# Convert the count matrix to a Seurat object
adata <- CreateSeuratObject(counts = ad_t, assay = "RNA")
#Add metdata to Seurat object
adata[["Source"]] <- ad$obs$Source
# Extract the "obsm" layer as a matrix
obsm_m <- as.matrix(ad$obsm$X_umap)
rownames(obsm_m) <- colnames(adata)
adata[["UMAP"]] <- CreateDimReducObject(embeddings = obsm_m, key = "UMAP")

#Filter
adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
adata <- subset(adata, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# Visualize QC metrics as a violin plot
VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#normalize
adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
# scale and run PCA
adata <- ScaleData(adata, features = rownames(adata))
adata <- RunPCA(adata, features = VariableFeatures(object = adata))


# cluster and visualize
adata <- FindNeighbors(adata, dims = 1:30)
adata <- FindClusters(adata, reduction.use = "UMAP", resolution = 0.8)
adata <- RunUMAP(adata, reduction.use = "UMAP", dims = 1:30)
pdf("umap_cou19_noannot24_7.pdf")
DimPlot(adata, reduction = "UMAP", label = FALSE)
dev.off()
#Assign Cell types using scType

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

#Use my gs_list file
db_ = "C:/Users/fallo/Documents/Internship_2023/tables/gs_listv4.xlsx";
tissue = "Brain" 
# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

#Assign cell type to the clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = adata[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(adata@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(adata@meta.data[adata@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(adata@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


#Overlay cell type assignments on UMAP
adata@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  adata@meta.data$customclassif[adata@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

dev.new()

pdf("sctype_umap_annot_cou19_24_7.pdf")
DimPlot(adata, reduction = "UMAP", label = FALSE, pt.size = 0.5, group.by = 'customclassif')        
dev.off()

pdf("sctype_umap_annot_eze10_20_7v2.pdf")
DimPlot(adata, reduction = "UMAP", label = TRUE, pt.size = 0.5, group.by = 'customclassif')        
dev.off()
#Repeat annotation on interneuron MGE and CGE clusters
#Subset the data
selected_cell_types <- c('MGE INs', 'CGE INs')

# Subset the Seurat object based on the selected cell types
gaba <- subset(adata, customclassif %in% selected_cell_types)
# cluster and visualize
gaba <- FindNeighbors(gaba, dims = 1:15)
gaba <- FindClusters(gaba, resolution = 0.8)
gaba <- RunUMAP(gaba, dims = 1:15)
pdf("umap_gaba_cou19_noannot24_7.pdf")
DimPlot(gaba, reduction = "UMAP", label = TRUE, group.by = 'customclassif')
DimPlot(gaba, reduction = "UMAP", label = FALSE)

dev.off()
#Use my gs_list file
db_2 = "C:/Users/fallo/Documents/Internship_2023/tables/gs_listv5.xlsx";
# prepare gene sets
gs_list2 = gene_sets_prepare(db_2, tissue)
#Assign cell type to the clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = gaba[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list2$gs_positive, gs2 = gs_list2$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(gaba@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(gaba@meta.data[gaba@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(gaba@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores2 = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores2$type[as.numeric(as.character(sctype_scores2$scores)) < sctype_scores2$ncells/4] = "Unknown"

#Overlay cell type assignments on UMAP
gaba@meta.data$customclassif = ""
for(j in unique(sctype_scores2$cluster)){
  cl_type = sctype_scores2[sctype_scores2$cluster==j,]; 
  gaba@meta.data$customclassif[gaba@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

dev.new()

pdf("sctype_umap_in_annot_cou19_24_7.pdf")
DimPlot(gaba, reduction = "UMAP", label = FALSE, pt.size = 0.5, group.by = 'customclassif')        
dev.off()


write.csv(sctype_scores, file = "cou19_sctype24_7.csv", row.names = TRUE)
write.csv(sctype_scores2, file = "cou19_sctype_INs_24_7.csv", row.names = TRUE)

