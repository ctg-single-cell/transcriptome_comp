#!/usr/bin/env R
library(Seurat)
library(dplyr)
library(HGNChelper)
library(ggplot2)
library(anndata)
#######################################################################
# Set the working directory
#working_directory = '/home/rfallon/integration'
output_directory = '/home/rfallon/output'
#Convert h5ad to Seurat
adata <- read_h5ad('/home/rfallon/integration/vitro_int.h5ad')

#Turn into matrix and transpose
ad_t <- t(as.matrix(adata$X))

# Convert the count matrix to a Seurat object
int <- CreateSeuratObject(counts = ad_t, assay = "RNA")

#Add metdata to Seurat object
int[["Source"]] <- adata$obs$Source
# Extract the "obsm" layer as a matrix
obsm_m <- as.matrix(adata$obsm$X_umap)
rownames(obsm_m) <- colnames(int)
int[["UMAP"]] <- CreateDimReducObject(embeddings = obsm_m, key = "UMAP")
#DimPlot(int, reduction = "UMAP", group.by = "Source")


# Plot the UMAP
umap1 <- DimPlot(int, reduction = "UMAP", group.by = "Source")
# Save the plot as a PDF file
output_path <- "/home/rfallon/output/vitro_umap1.pdf"
#output_path <- "C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/output/Integration/vitro_umap1.pdf"
ggsave(output_path, plot = umap1)

# normalize data
int[["percent.mt"]] <- PercentageFeatureSet(int, pattern = "^MT-")
int <- NormalizeData(int, normalization.method = "LogNormalize", scale.factor = 10000)
int <- FindVariableFeatures(int, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
int <- ScaleData(int, features = rownames(int))

int <- RunPCA(int, features = VariableFeatures(object = int))

# cluster and visualize
int <- FindNeighbors(int, dims = 1:30)
int <- RunUMAP(int, reduction.use = "UMAP", dims = 1:30)
int <- FindClusters(int, reduction.use = "UMAP", resolution = 0.8, dims = 1:30)

umap2 <- DimPlot(int, reduction = "UMAP", label = TRUE)
# Save the plot as a PDF file
output_path <- "/home/rfallon/output/vitro_umap2.pdf"
ggsave(output_path, plot = umap2)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file: input cell marker file
db_ = "gs_listv4.xlsx";

tissue = "Brain" 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

#Assign cell type to the clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = int[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 


# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(int@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(int@meta.data[int@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(int@meta.data$seurat_clusters==cl)), 10)
}))


sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


#Overlay cell type assignments on UMAP
int@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  int@meta.data$customclassif[int@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}



umap3 <- DimPlot(int, reduction = "UMAP", label = FALSE, group.by = 'customclassif') 
# Save the plot as a PDF file
output_path <- "/home/rfallon/output/vitro_umap3.pdf"
ggsave(output_path, plot = umap3)

umap4 <- DimPlot(int, reduction = "UMAP", label = FALSE, group.by = 'customclassif', split.by = 'Source')        
# Save the plot as a PDF file
output_path <- "/home/rfallon/output/vitro_umap4.pdf"
ggsave(output_path, plot = umap4)

#Make a dataframe that also has the clusters and corresponding age
# Create an empty dataframe to store cell type, source, and ncells information
cell_type_source_df <- data.frame(cell_type = character(),
                               source = character(),
                               ncells = integer(),
                               stringsAsFactors = FALSE)

# Loop over unique clusters
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  cell_type <- as.character(cl_type$type[1])
  
  # Get unique source value for the cluster
  source <- unique(int@meta.data$Source[int@meta.data$seurat_clusters == j])
  
  # Get the number of cells in the cluster
  ncells <- sum(int@meta.data$seurat_clusters == j)
  
  # Append the cell type, age, and ncells information to the dataframe
  cell_type_source_df <- rbind(cell_type_source_df, data.frame(cell_type = cell_type, source = source, ncells = ncells))
}

#Save output 
# Create the file path for the CSV
csv_file_path <- paste0(output_directory, "vitro_int_sctype.csv")
write.csv(cell_type_source_df, file = csv_file_path, row.names = TRUE)