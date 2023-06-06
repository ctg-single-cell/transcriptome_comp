library(Seurat)
library(SeuratDisk)
library(data.table)
library(dplyr)
library(HGNChelper)
library(ggplot2)
library(anndata)
library(reticulate)
library(SeuratWrappers)

#read in the integrated file with the zellkonverter library
setwd('C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/data/Herring/')
library(zellkonverter)
adata <- readH5AD('herr_int_sc.h5ad')
adata_Seurat <- as.Seurat(adata, counts = "X", data = NULL) #this is very slow

#save as rds
saveRDS(adata_Seurat, file = "herr_int.rds")

# Plot the Scanorama umap 
pdf("scan_herr_int.pdf")
DimPlot(adata_Seurat, reduction = "Scanorama", group.by = "age")
DimPlot(adata_Seurat, reduction = "Scanorama", split.by = "age")


# normalize data
adata_Seurat[["percent.mt"]] <- PercentageFeatureSet(adata_Seurat, pattern = "^MT-")
adata_Seurat <- NormalizeData(adata_Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
adata_Seurat <- FindVariableFeatures(adata_Seurat, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
adata_Seurat <- ScaleData(adata_Seurat, features = rownames(adata_Seurat))

adata_Seurat <- RunPCA(adata_Seurat, features = VariableFeatures(object = adata_Seurat))

# Check number of PC components (I chose 15 PCs for downstream analysis, based on Elbow plot)
ElbowPlot(adata_Seurat)

# cluster and visualize
adata_Seurat <- FindNeighbors(adata_Seurat, dims = 1:15)
adata_Seurat <- RunUMAP(adata_Seurat, reduction.use = "Scanorama", dims = 1:15)
adata_Seurat <- FindClusters(adata_Seurat, reduction.use = "Scanorama", resolution = 0.8, dims = 1:15)

pdf("umap_herr_seurint_noannot.pdf")
DimPlot(adata_Seurat, reduction = "Scanorama", label = TRUE)
dev.off()

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
es.max = sctype_score(scRNAseqData = adata_Seurat[["originalexp"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(adata_Seurat@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(adata_Seurat@meta.data[adata_Seurat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(adata_Seurat@meta.data$seurat_clusters==cl)), 10)
}))


sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


#Overlay cell type assignments on UMAP
adata_Seurat@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  adata_Seurat@meta.data$customclassif[adata_Seurat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

dev.new()

pdf("sctype_umap_annot_herr_intv2.pdf")
DimPlot(adata_Seurat, reduction = "Scanorama", label = FALSE, group.by = 'customclassif', split.by = 'age')        
dev.off()


# Create an empty dataframe to store cell type, age, and ncells information
cell_type_age_df <- data.frame(cell_type = character(),
                               age = character(),
                               ncells = integer(),
                               stringsAsFactors = FALSE)

# Loop over unique clusters
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[sctype_scores$cluster == j, ]
  cell_type <- as.character(cl_type$type[1])
  
  # Get unique age value for the cluster
  age <- unique(adata_Seurat@meta.data$age[adata_Seurat@meta.data$seurat_clusters == j])
  
  # Get the number of cells in the cluster
  ncells <- sum(adata_Seurat@meta.data$seurat_clusters == j)
  
  # Append the cell type, age, and ncells information to the dataframe
  cell_type_age_df <- rbind(cell_type_age_df, data.frame(cell_type = cell_type, age = age, ncells = ncells))
}

write.csv(cell_type_age_df, "herr_sc_int_sctype.csv", row.names = TRUE)

herr <- read.csv("herr_sc_int_sctype.csv")

# Calculate the percentage for each subset
herr2 <- herr %>%
  group_by(age, cell_type) %>%
  summarize(total_cells = sum(ncells),
            .groups = "drop") %>%
  group_by(age) %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)


#the following code is to make sure all of the bars are the same width
df_all_combinations <- expand.grid(age = unique(herr2$age), cell_type = unique(herr2$cell_type))
# Merge the new data frame with the original data frame
combo_full <- merge(df_all_combinations, herr2, by = c("age", "cell_type"), all.x = TRUE)

combo_full <- replace(combo_full, is.na(combo_full), 0)
# Convert the 'group' column into an ordered factor
combo_full$age <- factor(combo_full$age, levels = c("ga22","ga24", "ga34"), ordered = TRUE)


library(stringr)
pdf("herr_int_sctype_2_6.pdf")
# Create the ggplot bar chart
bar5 <- ggplot(combo_full, aes(fill = age, y = percentage, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

bar5 + labs(title="scType Annotation on Scanorama Integration", subtitle="Herring et al., 2022", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values=c( 
                                                                     "ga22" = "#D3D3D3", 
                                                                     "ga24" = "#A9A9A9",
                                                                     "ga34"= "#696969"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

