library(Seurat)
library(dplyr)
library(HGNChelper)
library(ggplot2)
library(anndata)
library(zellkonverter)
library(stringr)
#######################################################################
#Going from Scanpy to R with zellkonverter
adata <- readH5AD('herr_int_sc.h5ad')
adata_Seurat <- as.Seurat(adata, counts = "X", data = NULL) #this is very slow

# Plot the UMAP
DimPlot(int, reduction = "X_umap", group.by = "age")
DimPlot(int, reduction = "X_umap", split.by = "age")

# normalize data
int[["percent.mt"]] <- PercentageFeatureSet(int, pattern = "^MT-")
int <- NormalizeData(int, normalization.method = "LogNormalize", scale.factor = 10000)
int <- FindVariableFeatures(int, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
int <- ScaleData(int, features = rownames(int))

int <- RunPCA(int, features = VariableFeatures(object = int))

# cluster and visualize
int <- FindNeighbors(int, dims = 1:30)
int <- RunUMAP(int, reduction.use = "X_umap", dims = 1:30)
int <- FindClusters(int, reduction.use = "X_umap", resolution = 0.8, dims = 1:30)

DimPlot(int, reduction = "X_umap", label = TRUE)


#Assign cell type to the clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = int[["originalexp"]]@scale.data, scaled = TRUE, 
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

dev.new()

DimPlot(int, reduction = "X_umap", label = FALSE, group.by = 'customclassif') 

DimPlot(int, reduction = "X_umap", label = FALSE, group.by = 'customclassif', split.by = 'age')        

#Make a dataframe that also has the clusters and corresponding age
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
  age <- unique(int@meta.data$age[int@meta.data$seurat_clusters == j])
  
  # Get the number of cells in the cluster
  ncells <- sum(int@meta.data$seurat_clusters == j)
  
  # Append the cell type, age, and ncells information to the dataframe
  cell_type_age_df <- rbind(cell_type_age_df, data.frame(cell_type = cell_type, age = age, ncells = ncells))
}

write.csv(cell_type_age_df, "herr_sc_int_sctype.csv", row.names = TRUE)

################################################################################################
# Calculate the percentage for each subset
herr <- cell_type_age_df %>%
  group_by(age, cell_type) %>%
  summarize(total_cells = sum(ncells),
            .groups = "drop") %>%
  group_by(age) %>%
  mutate(percentage = (total_cells / sum(total_cells)) * 100)


#the following code is to make sure all of the bars are the same width
df_all_combinations <- expand.grid(age = unique(herr$age), cell_type = unique(herr$cell_type))
# Merge the new data frame with the original data frame
combo_full <- merge(df_all_combinations, herr, by = c("age", "cell_type"), all.x = TRUE)

combo_full <- replace(combo_full, is.na(combo_full), 0)
# Convert the 'group' column into an ordered factor
combo_full$age <- factor(combo_full$age, levels = c("ga22","ga24", "ga34"), ordered = TRUE)

# Create the ggplot bar chart
bar <- ggplot(combo_full, aes(fill = age, y = percentage, x = cell_type)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.8)

bar + labs(title="scType Annotation on Scanorama Integration", subtitle="Herring et al., 2022", y = "Proportion (%)", fill ="scRNAseq Data") + 
  theme(axis.title.x = element_blank()) + scale_fill_manual(values=c( 
    "ga22" = "#D3D3D3", 
    "ga24" = "#A9A9A9",
    "ga34"= "#696969"))+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(size = 9)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
