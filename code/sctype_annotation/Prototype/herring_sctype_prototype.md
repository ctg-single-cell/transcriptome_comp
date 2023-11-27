scType Annotation Herring Prototype
================
Fallon Ratner
2023-11-24

# This script uses scType to annotate one of the Herring datasets

-   First the brain database gene markers from scType is used.
-   Then the custom gene markers are used with scType

# Setting up

-   Load libraries
-   Load Herring dataset and convert from h5ad to Seurat object
-   Perform pre-processing, normalization, and clustering
-   Assign Cell Types using scType

``` r
library(Seurat)
```

    ## The legacy packages maptools, rgdal, and rgeos, underpinning this package
    ## will retire shortly. Please refer to R-spatial evolution reports on
    ## https://r-spatial.org/r/2023/05/15/evolution4.html for details.
    ## This package is now running under evolution status 0

    ## Attaching SeuratObject

``` r
library(HGNChelper)
library(anndata)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

## Loading h5ad file and converting to seurat object

``` r
#set own wd
setwd("C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/data/Herring/")

# Load h5ad file
ad <- read_h5ad("ga22_cleaned_count_matrices.h5ad")

#Transpose the matrix
ad_t <- t(as.matrix(ad$X))

#create the seurat object
adata <- CreateSeuratObject(counts = ad_t, assay = "RNA")

#change rownames of meta data to cell id
rownames(adata@meta.data) <- colnames(adata)

#convert all factor vectors to character vectors
i <- sapply(adata@meta.data, is.factor)
adata@meta.data[i] <- lapply(adata@meta.data[i], as.character)
```

## Pre-processing, Normalizing, and Clustering

``` r
#Filter for mitochondrial genes
adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^MT-")

#normalize
adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
adata <- ScaleData(adata, features = rownames(adata))
```

    ## Centering and scaling data matrix

``` r
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
```

    ## PC_ 1 
    ## Positive:  RALYL, KCNQ5, KHDRBS2, MACROD2, LRRC4C, PRICKLE2, CNTN4, CALN1, RYR2, MEF2C 
    ##     NEBL, FAM189A1, LDB2, ZNF385B, PTPRK, SYBU, FGF12, MGAT4C, GRM7, PCDH11Y 
    ##     MARCH1, DOK5, KLHL1, FRMPD4, PCLO, GABRB1, ARHGAP20, SAMD5, CCBE1, TENM2 
    ## Negative:  QKI, ZBTB20, NPAS3, ERBB4, CHD7, NFIB, MSI2, DLX6-AS1, FBXL7, ZEB1 
    ##     RP11-436D23.1, ZNF536, SYNE2, VCAN, SOX6, SLC1A3, CDH4, GPR98, PDZRN3, TMTC2 
    ##     ADARB2, PAM, GLI3, ZFHX4, RP11-76I14.1, SOX2, MAGI1, SEMA5B, FGFR2, RFTN2 
    ## PC_ 2 
    ## Positive:  HS3ST4, TLE4, KAZN, GRIK3, OPCML, LRP1B, ZNF385D, RGS6, NEGR1, PLXDC2 
    ##     CLSTN2, COL12A1, MLIP, LRRC16A, ZFPM2, PDZD2, PCDH9, KHDRBS3, KIAA1217, SLIT3 
    ##     EPHA7, DLC1, TRPM3, CDH18, GAS7, DPP10, SORCS1, SULF1, PDE8B, FAM160A1 
    ## Negative:  CUX2, LINC01158, IQCJ-SCHIP1, GLIS3, PLXNA4, NDST3, FAM19A1, UNC5D, CACNA2D1, EPHA6 
    ##     EPHA3, CPNE8, GRM7, STK32B, MEIS2, ARL15, CNTN3, ZNF804A, CBLN2, L3MBTL4 
    ##     CCBE1, SDK1, VSTM2L, PRSS12, POU3F2, BHLHE22, LINC01102, R3HDM1, ZNF608, TMOD1 
    ## PC_ 3 
    ## Positive:  TNC, BCAN, CDH20, SLC1A3, ATP1A2, PTN, EEPD1, PARD3B, SPATA6, LRIG1 
    ##     SLCO1C1, SEZ6L, TMEM132D, ZEB1, PON2, RORA, NTNG1, GPC5, LTBP1, SEMA5A 
    ##     RP11-849I19.1, SLC24A3, PHLPP1, EGFR, LRRC4C, ADAMTS6, RP11-192P3.5, FBXL7, SALL3, HSPA1A 
    ## Negative:  PDE1A, HS3ST4, DSCAML1, DSCAM, ST18, CDH18, MLIP, FAM160A1, ASTN2, SLIT3 
    ##     SLC24A2, ZNF385D, HCRTR2, GRIK3, COL12A1, TLE4, KIAA1217, TMEM178A, TRPM3, HS6ST3 
    ##     NR4A2, EPHA7, NGEF, KIAA1456, CNTNAP4, BCL11B, CLSTN2, NFIB, ANKRD33B, PLXDC2 
    ## PC_ 4 
    ## Positive:  ERBB4, NRXN3, DLX6-AS1, ZNF536, PDZRN3, NPAS3, GAD1, NXPH1, THRB, ADARB2 
    ##     RP11-242P2.1, GAD2, ARX, VSTM2A, CNTNAP2, FAM65B, ATRNL1, PAM, NR3C2, LHFPL3 
    ##     SLC6A1, PLS3, PTPRT, NR2F2-AS1, ST8SIA5, TENM2, NPAS1, GRIK2, SORCS3, MAGI1 
    ## Negative:  MEIS2, DSCAM, HS6ST3, PTPRZ1, NLGN1, RP11-436D23.1, CDH4, KIAA1456, RP11-76I14.1, PDE1A 
    ##     SEMA3C, DAB1, ASTN2, SOX5, PALMD, LINC01158, RAI14, LMO3, SLC24A2, GPR98 
    ##     TLE4, POU3F2, ITPR2, CACNA2D1, HS3ST4, NRG1, FAM13A, MCTP1, DCC, ZFPM2 
    ## PC_ 5 
    ## Positive:  STK32B, ASIC2, CPNE8, FAM19A2, RP11-85M11.2, CACNA2D3, GALNTL6, EPHA6, CTC-340A15.2, KCNJ6 
    ##     MCTP1, DAB1, CNTN3, CUX2, LEPREL1, KIAA1377, SLIT2, ATP8A2, CNTNAP5, CTC-535M15.2 
    ##     HS6ST3, SPHKAP, PRKG1, NRXN3, NDST3, TRPM3, ZNF608, KIAA1217, GLIS3, KIF26B 
    ## Negative:  POU6F2, IL1RAPL2, ZNF804B, KLHL1, KCNH7, MYO10, PTPRK, ZNF385B, FOXP1, UNC5C 
    ##     ADAMTS3, LINGO2, INHBA, HTR7, FNDC1, NEUROD6, LHX2, PTPRT, PTPN13, AC007740.1 
    ##     DCC, GABRA4, FOXP2, KLHL14, SAMD5, KIAA1239, SYNDIG1, CNTNAP3B, RP11-307P5.1, CCDC85A

``` r
# cluster and visualize
adata <- FindNeighbors(adata, dims = 1:15)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
adata <- FindClusters(adata, reduction.use = "UMAP", resolution = 0.8)
```

    ## Suggested parameter: reduction instead of reduction.use

    ## Suggested parameter: reduction instead of reduction.use

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 10466
    ## Number of edges: 364317
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8973
    ## Number of communities: 24
    ## Elapsed time: 2 seconds

``` r
adata <- RunUMAP(adata, reduction.use = "UMAP", dims = 1:15)
```

    ## Suggested parameter: reduction instead of reduction.use

    ## 14:32:31 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 14:32:31 Read 10466 rows and found 15 numeric columns

    ## 14:32:31 Using Annoy for neighbor search, n_neighbors = 30

    ## 14:32:31 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 14:32:34 Writing NN index file to temp file C:\Users\fallo\AppData\Local\Temp\RtmpG6xiRU\file2d7871be30bc
    ## 14:32:34 Searching Annoy index using 1 thread, search_k = 3000
    ## 14:32:39 Annoy recall = 100%
    ## 14:32:39 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 14:32:40 Initializing from normalized Laplacian + noise (using irlba)
    ## 14:32:40 Commencing optimization for 200 epochs, with 427480 positive edges
    ## 14:32:53 Optimization finished

## Load scType Functions

``` r
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

## Assign Cell Types with scType Gene Markers

``` r
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain"

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = adata[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(adata@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(adata@meta.data[adata@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells =  sum(adata@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
```

    ## # A tibble: 24 × 3
    ## # Groups:   cluster [24]
    ##    cluster type                  scores
    ##    <fct>   <chr>                  <dbl>
    ##  1 0       Unknown                 206.
    ##  2 8       Neuroblasts             579.
    ##  3 22      Radial glial cells      355.
    ##  4 13      Radial glial cells     1419.
    ##  5 1       Glutamatergic neurons   373.
    ##  6 11      Dopaminergic neurons    217.
    ##  7 18      GABAergic neurons       461.
    ##  8 2       GABAergic neurons       495.
    ##  9 6       Mature neurons          318.
    ## 10 5       Unknown                 150.
    ## # ℹ 14 more rows

``` r
#Overlay cell type assignments on UMAP
adata@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  adata@meta.data$customclassif[adata@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(adata, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'customclassif') 
```

![](herring_sctype_prototype_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> -
There is a large unknown population when using the scType Brain Database

## Assign Cell Types with Custom Gene Markers

``` r
#Use the curated gene list: gs_list file version4
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
```

    ## # A tibble: 24 × 3
    ## # Groups:   cluster [24]
    ##    cluster type                    scores
    ##    <fct>   <chr>                    <dbl>
    ##  1 0       Upper Excitatory Layers   530.
    ##  2 8       CGE INs                  1224.
    ##  3 22      Immature Astrocytes      1162.
    ##  4 13      Glioblasts               1280.
    ##  5 1       Maturing Excitatory      1042.
    ##  6 11      Deep Excitatory Layers    534.
    ##  7 18      MGE INs                   480.
    ##  8 2       Upper Excitatory Layers   605.
    ##  9 6       Deep Excitatory Layers    732.
    ## 10 5       Immature Excitatory       516.
    ## # ℹ 14 more rows

``` r
#Overlay cell type assignments on UMAP
adata@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  adata@meta.data$customclassif[adata@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(adata, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'customclassif')
```

![](herring_sctype_prototype_files/figure-gfm/unnamed-chunk-6-1.png)<!-- --> -
The Unknown population is no longer identified when using the customized
gene markers.

## Annotate Interneuron Sub-type with the Custom Gene Markers

``` r
#Repeat annotation on interneuron MGE and CGE clusters
selected_cell_types <- c('MGE INs', 'CGE INs')
# Subset the Seurat object based on the selected cell types
gaba <- subset(adata, customclassif %in% selected_cell_types)
# cluster and visualize
gaba <- FindNeighbors(gaba, dims = 1:15)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
gaba <- FindClusters(gaba, resolution = 0.8)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1149
    ## Number of edges: 41537
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.7993
    ## Number of communities: 7
    ## Elapsed time: 0 seconds

``` r
gaba <- RunUMAP(gaba, dims = 1:15)
```

    ## 14:33:12 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 14:33:12 Read 1149 rows and found 15 numeric columns

    ## 14:33:12 Using Annoy for neighbor search, n_neighbors = 30

    ## 14:33:12 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 14:33:13 Writing NN index file to temp file C:\Users\fallo\AppData\Local\Temp\RtmpG6xiRU\file2d787b9c5a63
    ## 14:33:13 Searching Annoy index using 1 thread, search_k = 3000
    ## 14:33:13 Annoy recall = 100%
    ## 14:33:13 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 14:33:14 Initializing from normalized Laplacian + noise (using irlba)
    ## 14:33:14 Commencing optimization for 500 epochs, with 46390 positive edges
    ## 14:33:18 Optimization finished

``` r
DimPlot(gaba, reduction = "umap", label = TRUE, group.by = 'customclassif')
```

![](herring_sctype_prototype_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
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

DimPlot(gaba, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'customclassif')
```

![](herring_sctype_prototype_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->
