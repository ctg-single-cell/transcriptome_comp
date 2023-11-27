#Created by Fallon Ratner on 17.07.23
library(Seurat)
library(dplyr)
library(HGNChelper)
library(ggplot2)
library(data.table)

# Function to read and process Seurat object
read_and_process_seurat <- function(file_path, project_name, min_cells = 3, min_features = 200, dims = 1:15, resolution = 0.8) {
  # Detect file type based on extension
  file_extension <- tools::file_ext(file_path)
  
  if (file_extension == "rds") {
    # Read RDS file
    adata <- readRDS(file_path)
  } else if (file_extension == "txt") {
    # Read text file
    mydata <- read.table(file_path, header = TRUE, sep = " ")
    adata <- CreateSeuratObject(counts = mydata, project = project_name, min.cells = min_cells, min.features = min_features)
  } else if (file_extension == "filtered_gene_matrices") {
    # Read 10X data
    mydata <- Read10X(data.dir = file_path)
    adata <- CreateSeuratObject(counts = mydata, project = project_name, min.cells = min_cells, min.features = min_features)
  } else {
    stop("Unsupported file type.")
  }
  
  # Filter and preprocess
  adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
  adata <- subset(adata, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  
  # Cluster and visualize
  adata <- FindNeighbors(adata, dims = dims)
  adata <- FindClusters(adata, resolution = resolution)
  adata <- RunUMAP(adata, dims = dims)
  
  return(adata)
}

# Lookup table for datasets and clustering parameters
dataset_lookup_table <- data.frame(
  file_path = c("ga22.rds", "ga24.rds", "ga34.rds", "day2.rds", "liu.rds", "polio.rds", 
                "FetalBrain3.rawdge.txt", "FetalBrain4.rawdge.txt", "FetalBrain5.rawdge.txt", "FetalBrain6.rawdge.txt", 
                "HFA567_total.filtered_gene_matrices/", "HFA570_total.filtered_gene_matrices/", "HFA571_total.filtered_gene_matrices/"),
  project_name = c("ga22", "ga24", "ga34", "day2", "liu", "polio", "FetalBrain3", "FetalBrain4", "FetalBrain5", "FetalBrain6", "couturier567", "couturier570", "couturier571"),
  min_cells = rep(3, 13),
  min_features = rep(200, 13),
  dims = list(1:15, 1:15, 1:15, 1:10, 1:10, 1:15, 1:5, 1:5, 1:10, 1:5, 1:15, 1:15, 1:15),
  resolution = rep(c(0.8, 0.8, 0.8, 0.8, 0.5, 0.8, 0.5, 0.5, 0.5, 0.8, 0.5, 0.5, 0.5), 1)
)

# Loop through the lookup table and process each dataset
for (i in seq_len(nrow(dataset_lookup_table))) {
  file_path <- dataset_lookup_table$file_path[i]
  project_name <- dataset_lookup_table$project_name[i]
  min_cells <- dataset_lookup_table$min_cells[i]
  min_features <- dataset_lookup_table$min_features[i]
  dims <- dataset_lookup_table$dims[[i]]
  resolution <- dataset_lookup_table$resolution[i]
  
  adata <- read_and_process_seurat(file_path, project_name, min_cells, min_features, dims, resolution)
 
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

  #Overlay cell type assignments on UMAP
  adata@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    adata@meta.data$customclassif[adata@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }

  DimPlot(adata, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'customclassif')        

  #Repeat annotation on interneuron MGE and CGE clusters
  #Subset the data
  selected_cell_types <- c('MGE INs', 'CGE INs')
  # Subset the Seurat object based on the selected cell types
  gaba <- subset(adata, customclassif %in% selected_cell_types)
  # cluster and visualize
  gaba <- FindNeighbors(gaba, dims = dims)
  gaba <- FindClusters(gaba, resolution = 0.5)
  gaba <- RunUMAP(gaba, dims = 1:15)
 
  DimPlot(gaba, reduction = "umap", label = TRUE, group.by = 'customclassif')

  #scType annotation
  db_2 = "gs_listv5.xlsx";
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

  DimPlot(gaba, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'customclassif')        
 
  write.csv(sctype_scores, file = paste0(adata_)"sctype.csv", row.names = TRUE)
  write.csv(sctype_scores2, file = paste0(adata_)"sctype_INs.csv", row.names = TRUE)
}