library(Seurat)
library(HGNChelper)
library(anndata)

# List your datasets
datasets <- c("ga22","ga24","ga34")
# Loop through each dataset
for (dataset in datasets) {
  # Construct the filename for each dataset
  h5ad_file <- paste0(dataset, "_cleaned_count_matrices.h5ad")
   # Load h5ad file
  ad <- read_h5ad(h5ad_file)
  
  ad_t <- t(as.matrix(ad$X))
  #create seurat object
  adata <- CreateSeuratObject(counts = ad_t, assay = "RNA")
  #change rownames of meta data to cell id
  rownames(adata@meta.data) <- colnames(adata)
  #convert all factor vectors to character vectors
  i <- sapply(adata@meta.data, is.factor)
  adata@meta.data[i] <- lapply(adata@meta.data[i], as.character)
  #Filter and normalize
  adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")
  #normalize
  adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
  adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
  # scale and run PCA
  adata <- ScaleData(adata, features = rownames(adata))
  adata <- RunPCA(adata, features = VariableFeatures(object = adata))
  
  
  # cluster and visualize
  adata <- FindNeighbors(adata, dims = 1:15)
  adata <- FindClusters(adata, reduction.use = "UMAP", resolution = 0.8)
  adata <- RunUMAP(adata, reduction.use = "UMAP", dims = 1:15)
  
  #Assign Cell types using scType
  # load gene set preparation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  # load cell type annotation function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  #Use the curated gene list: gs_list file version4
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
  print(sctype_scores[,1:3])
  
  
  #Overlay cell type assignments on UMAP
  adata@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    adata@meta.data$customclassif[adata@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  DimPlot(adata, reduction = "UMAP", label = FALSE, pt.size = 0.5, group.by = 'customclassif') 
  # Define paths for outputs
  sctype_file <- paste0(dataset, "_sctype.csv")
  # Save output
  write.csv(sctype_scores, file = sctype_file, row.names = TRUE)
}
  
