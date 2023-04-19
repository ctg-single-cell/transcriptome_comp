#Created by: Fallon Ratner on 12.04.23
#Project: scRNAseq Prototype Analysis with scType

#import libraries
library(Seurat)
library(data.table)
library(dplyr)
library(HGNChelper)
library(ggplot2)

# set file path for scRNAseqData
file_path <- "{Include path to your file}" # path to the scRNAseq data file

#Loading data if in matrix (mtx) format
load(file_path)
data_mat <- as.matrix(raw_counts_mat)

#Create a seurat oobject
myseurat <- CreateSeuratObject(counts = data_mat)

#Loading data if in txt format
mydata <- fread(file_path, header = TRUE, sep = "\t")
#make the genes the rownames
row.names(mydata) <- mydata$GENE

#create a seurat object
myseurat <- CreateSeuratObject(counts = mydata, project = "myproject", min.cells = 3, min.features = 200)

#Make sure that the rownames are genes before continuing
data <- myseurat[["RNA"]]@data
row.names(data)

#Use Seurat for Pre-processing (based on scType workflow)
# normalize data
myseurat[["percent.mt"]] <- PercentageFeatureSet(myseurat, pattern = "^MT-")
myseurat <- NormalizeData(myseurat, normalization.method = "LogNormalize", scale.factor = 10000)
myseurat <- FindVariableFeatures(myseurat, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
myseurat <- ScaleData(myseurat, features = rownames(myseurat))
myseurat <- RunPCA(myseurat, features = VariableFeatures(object = myseurat))

# Check number of PC components (we selected 10 PCs for downstream analysis, based on Elbow plot)
ElbowPlot(myseurat)

# cluster and visualize
myseurat <- FindNeighbors(myseurat, dims = 1:10)
myseurat <- FindClusters(myseurat, resolution = 0.8)
myseurat <- RunUMAP(myseurat, dims = 1:10)
DimPlot(myseurat, reduction = "umap")

#Assign Cell types using scType

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file: input cell marker file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

#Assign cell type to the clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = myseurat[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = NULL) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(myseurat@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(myseurat@meta.data[myseurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(myseurat@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

#save the annotations to csv for further processing via python
write.csv(sctype_scores, file = "{file name}.csv", row.names = TRUE)

#Overlay cell type assignments on UMAP
myseurat@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  myseurat@meta.data$customclassif[myseurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

dev.new()    

pdf("{plot name}.pdf")
DimPlot(myseurat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'customclassif')        

#can save the seurat object with new annotations
saveRDS(myseurat, file = "{object name}.rds")

