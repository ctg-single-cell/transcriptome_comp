# Annotate Cells with ScType
### Date: 09-06-2023   
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
### This README contains instructions on how to run scType on a local computer to annotate single cell/nuclei RNA-seq from in vivo brain tissue.
### _NOTE_: The annotation pipeline has been optimized for fetal brain tissue and may need further optimization for other tissue types.

# What is scType?
### scType is a tool to annotate cells based on unsupervised learning. It utilizes a combination of positive and negative gene markers to annotate cells and clusters based on a computed score. scType requires as input: a gene marker list, the scType functions, a normalized + scaled gene expression matrix, and Seurat computed clusters. After running the computation, scType will output a dataframe consisting of the cluster and the assigned cell label (which has the highest score). This can then be used for further processing and visualization.
### For more information about scType please refer to the orginal publication: 
### Ianevski, A., Giri, A.K. & Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations
### from single-cell transcriptomic data. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w
### Or visit the scType Github: https://github.com/IanevskiAleksandr/sc-type


# 0. Installation of software and packages
### scType pipeline requires the following
* [Python]: version 3.10.9
* Python packages (Scanpy)
* [R]: version 4.2.3
* R packages (Seurat; HGNChelper; anndata; dplyr; stringr)

# 1. Downloading the scRNAseq Data
### Download the Data
* Go to the link provided in the following text (Dataset Used) and download this file: 'Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad'
### Dataset Used:
 Herring et al., 2022
        snRNAseq of fetal Prefrontal Cortex (PFC)
        Number of Nuclei in ga22: 11,660 (10,466 cells kept)
        Total Genes ga22: 27,565
        Number of Nuclei in ga24: 10,809 (9,376 cells kept)
        Total Genes ga24: 27,737
        Number of Nuclei in ga34: 8,431 (6,738 cells kept)
        Total Genes ga34: 26,083
        Sequencing Technology: Illumina NextSeq 550
        Link:https://console.cloud.google.com/storage/browser/neuro-dev/Processed_data;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
        Accessed: 19-04-2023

# 2. Extracting the individual scRNAseq Datasets
### In the extract_herring.py script the downloaded file from the previous step, will be used to extract the individual datasets which have been cleaned by the original authors. The output will be a cleaned count matrix in the h5ad file format.

# 3. Annotate with scType
### In the sctype_annot_herring.R script, the files from the previous step will be used for annotation with the scType tool.
1. The h5ad file needs to be converted into a Seurat object which is done with the anndata and Seurat packages.
2. The data is pre-processed using a standard Seurat pipeline to normalize and scale the data which is required as input for Sctype.
3. The data is clustered using the standard Seurat pipeline which is also required as input for Sctype.
4. The scType functions and gene list are loaded into the environment.
    1. The gene list used is gs_listv4.xlsx, which consists of genes curated by scType and from the literature.
5. The sctype score is calculated for each cell and each cluster resulting in a dataframe with an assigned label for each cluster.
    1. The score df is saved as a text file for each format and can be used for further processing.
6. The cell labels are visualized in the Seurat UMAP. 

# 3. Visualize scType Ouput with Bar Chart
### In the plot_herr_sctype.R script, the files from the previous step are used to generate a bar chart.
1. Read in the files, combine into one dataframe, and calculate the proportion as a percentage.
2. Not all cell types are found in each dataset, so those values will be assigned to 0.
3. Visualize the percentage in a bar chart. 