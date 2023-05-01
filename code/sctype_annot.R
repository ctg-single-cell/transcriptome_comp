#Created by: Fallon Ratner on 13.04.23
library(Seurat)
library(SeuratDisk)
library(data.table)
library(rhdf5)
library(dplyr)
library(HGNChelper)
library(ggplot2)
# set file path for scRNAseqData
file_path <- "C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/data/Velasco/GSE129519_expression_PGP1.3mon.txt/expression_PGP1.3mon.txt" # path to the scRNAseq data file
file_path <- "C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/data/Velasco/expression_PGP1.6mon.txt/expression_PGP1.6mon.txt" # path to the scRNAseq data file
file_path <- "C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/data/Polio/raw_counts_mat.rdata"


#Load the Polio data
load(file_path)
data_mat <- as.matrix(raw_counts_mat)
myseurat <- CreateSeuratObject(counts = data_mat)


#Load the Velasco data
mydata <- fread(file_path, header = TRUE, sep = "\t")
row.names(mydata) <- mydata$GENE

myseurat <- CreateSeuratObject(counts = mydata, project = "myproject", min.cells = 3, min.features = 200)

#check that the rownames are genes
data <- myseurat[["RNA"]]@data
row.names(data)

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
pdf("umap_annot_polio_noannot.pdf")
DimPlot(myseurat, reduction = "umap")
dev.off()
#Assign Cell types using scType

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file: input cell marker file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

#or use my own gene list
#Use my gs_list file
#db_2 = "gs_list_FR.xlsx";
#tissue = "Brain" 
# prepare gene sets
#gs_list2 = gene_sets_prepare(db_2, tissue)

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

write.csv(sctype_scores, file = "Polio_cluster_sctype.csv", row.names = TRUE)

#Overlay cell type assignments on UMAP
myseurat@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  myseurat@meta.data$customclassif[myseurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

dev.new()
DimPlot(myseurat, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')        

pdf("sctype_umap_annot_Polio.pdf")
DimPlot(myseurat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'customclassif')        

saveRDS(myseurat, file = "polio_sctype.rds")
