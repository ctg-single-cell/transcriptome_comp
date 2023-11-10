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
adata <- read_h5ad('/home/rfallon/vitro_int_09_08.h5ad')

#Turn into matrix and transpose
ad_t <- t(as.matrix(adata$X))
#removing this gene because causing error in creating seurat object
ad_t <- ad_t[rownames(ad_t) != "RP4-633O19--A.1", ]
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
#output_path <- "/home/rfallon/output/vitro_umap1.pdf"
#output_path <- "C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/output/Integration/vitro_umap1.pdf"
#ggsave(output_path, plot = umap1)

# normalize data
int[["percent.mt"]] <- PercentageFeatureSet(int, pattern = "^MT-")
int <- NormalizeData(int, normalization.method = "LogNormalize", scale.factor = 10000)
int <- FindVariableFeatures(int, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
int <- ScaleData(int, features = rownames(int))

int <- RunPCA(int, features = VariableFeatures(object = int))

# cluster and visualize
int <- FindNeighbors(int, dims = 1:40)
int <- RunUMAP(int, reduction.use = "UMAP", dims = 1:40)
int <- FindClusters(int, reduction.use = "UMAP", resolution = 0.8, dims = 1:40)

umap2 <- DimPlot(int, reduction = "UMAP", label = TRUE)
# Save the plot as a PDF file
#output_path <- "/home/rfallon/output/vitro_umap2.pdf"
#ggsave(output_path, plot = umap2)

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
output_path <- "/home/rfallon/output/vitro_umap3_1011.pdf"
ggsave(output_path, plot = umap3)

umap4 <- DimPlot(int, reduction = "UMAP", label = FALSE, group.by = 'customclassif', split.by = 'Source')        
# Save the plot as a PDF file
output_path <- "/home/rfallon/output/vitro_umap4_1011.pdf"
ggsave(output_path, plot = umap4)

#Make a dataframe that also has the clusters and corresponding age
cell_type_source_df <- data.frame(cell_type = character(),
                            source = character(),
                            ncells = integer(),
                            seurat_cluster = integer(),
                            stringsAsFactors = FALSE)

# Loop over unique clusters
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  cell_type <- as.character(cl_type$type[1])
  
  # Get the unique source values (samples) for the cluster
  sources <- unique(int@meta.data$Source[int@meta.data$seurat_clusters == j])
  
  for (source in sources) {
    # Get the number of cells in the cluster for the current sample
    ncells <- sum(int@meta.data$seurat_clusters == j & int@meta.data$Source == source)
    
    # Append the cell type, sample, ncells, and Seurat cluster number to the dataframe
    cell_type_source_df <- rbind(cell_type_source_df, data.frame(
      cell_type = cell_type,
      source = source,
      ncells = ncells,
      seurat_cluster = j  # Add the Seurat cluster number
    ))
  }
}
#Save output 
# Create the file path for the CSV
csv_file_path <- paste0(output_directory, "vitro_int_sctype1011.csv")
write.csv(cell_type_source_df, file = csv_file_path, row.names = TRUE)

##############################################################################################################
#Repeat annotation on interneuron MGE and CGE clusters
#Subset the data
selected_cell_types <- c('MGE INs', 'CGE INs')

# Subset the Seurat object based on the selected cell types
gaba <- subset(int, customclassif %in% selected_cell_types)
# cluster and visualize
gaba <- FindNeighbors(gaba, dims = 1:15)
gaba <- FindClusters(gaba, resolution = 0.5)
gaba <- RunUMAP(gaba, dims = 1:15)

gabamap <- DimPlot(gaba, reduction = "umap", label = TRUE, group.by = 'customclassif')
# Save the plot as a PDF file
output_path <- "/home/rfallon/output/vitro_gabamap1011.pdf"
ggsave(output_path, plot = gabamap)
gabamap2 <- DimPlot(gaba, reduction = "umap", label = TRUE)
# Save the plot as a PDF file
output_path <- "/home/rfallon/output/vitro_gabamap21011.pdf"
ggsave(output_path, plot = gabamap2)
#Use my gs_list file
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

gabamap3 <- DimPlot(gaba, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'customclassif')        
# Save the plot as a PDF file
output_path <- "/home/rfallon/output/vitro_gabamap3.pdf"
ggsave(output_path, plot = gabamap3)

#Make a dataframe that also has the clusters and corresponding age
# Create an empty dataframe to store cell type, source, and ncells information
cell_type_source_df2 <- data.frame(cell_type = character(),
                            source = character(),
                            ncells = integer(),
                            seurat_cluster = integer(),
                            stringsAsFactors = FALSE)

# Loop over unique clusters
for (j in unique(sctype_scores2$cluster)) {
  cl_type <- sctype_scores2[sctype_scores2$cluster == j, ]
  cell_type <- as.character(cl_type$type[1])
  
  # Get the unique source values (samples) for the cluster
  sources <- unique(gaba@meta.data$Source[gaba@meta.data$seurat_clusters == j])
  
  for (source in sources) {
    # Get the number of cells in the cluster for the current sample
    ncells <- sum(gaba@meta.data$seurat_clusters == j & gaba@meta.data$Source == source)
    
    # Append the cell type, sample, ncells, and Seurat cluster number to the dataframe
    cell_type_source_df2 <- rbind(cell_type_source_df2, data.frame(
      cell_type = cell_type,
      source = source,
      ncells = ncells,
      seurat_cluster = j  # Add the Seurat cluster number
    ))
  }
}

#Save output 
# Create the file path for the CSV
csv_file_path <- paste0(output_directory, "vitro_gaba_int_sctype1011.csv")
write.csv(cell_type_source_df2, file = csv_file_path, row.names = TRUE)

saveRDS(int, "int_1011.rds")
saveRDS(gaba, "gaba_int_1011.rds")