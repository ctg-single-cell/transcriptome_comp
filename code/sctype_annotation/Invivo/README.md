# Annotate Cells with ScType
### Date: 09-06-2023   
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
### This README contains instructions on how to run scType on a local computer to annotate single cell/nuclei RNA-seq from in vivo brain tissue.

**NOTE**: The annotation pipeline has been optimized for fetal brain tissue and may need further optimization for other tissue types.

# 0. Installation of software and packages
### scType pipeline requires the following
* [Python]: version 3.10.9
* Python packages (Scanpy)
* [R]: version 4.2.3
* R packages (Seurat; HGNChelper; anndata; dplyr; stringr)

# 1. Downloading the scRNAseq Data
### Download the Data
* Go to the link provided in the following text (Dataset Used) and download this file: 'Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad'
### Datasets Used:
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
Polioudakis et al., 2019
    scRNAseq of cortical 17-18 GW
    Number of Cells: 40,000 (pre-qc)
    Link: http://solo.bmap.ucla.edu/shiny/webapp/
    Accessed: 28-3-2023
Han et al., 2020
    scRNAseq of 11-13GW
    Han 11GW rep1: 3,920 (pre-qc)
    Han 11GW rep2: 1,705 (pre-qc)
    Han 12GW: 5,096 (pre-qc)
    Han 13GW: 2,904 (pre-qc)
    Link:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134355
    Download:GSE134355_RAW.tar
    Accessed: 22-06-2023
Couturier et al., 2020
    scRNAseq at 13, 17, & 19 GW
    Couturier 13GW: 12,455 cells (pre-qc)
    Couturier 17GW: 7,434 cells (pre-qc)
    Couturier 19GW: 10,436 cells (pre-qc)
    Link:https://github.com/mbourgey/scRNA_GBM
    Download: GBM_cellranger_matrix.tar
    Accessed: 14-06-2023
Liu et al., 2023
    scRNAseq of cortical 17-19 GW
    Number of Cells: 12,455 (pre-qc)
    Link:https://www.synapse.org/#!Synapse:syn51201773/files/
    Download: hNSPC_raw_counts.h5ad
    Accessed: 19-07-2023

# 1. Convert h5ad files to rds with the h5ad_to_rds.r script. 
- H5ad files are compatible with the scanpy package in python. In order to make this single cell object compatible with the Seurat object in R, an rds object is generated using the anndata, reticulate, and Seurat packages.
- This code also transfers the metadata from the h5ad to rds object. 

# 2. Annotate Datasets Individually with scType
### In the sctype_annot_vivo.R script, the datasets will be used for annotation with the scType tool.
1. The individual datasets need to be converted into Seurat objects based on the file format provided (h5ad files are converted in the h5ad_to_rds script).
2. The data is pre-processed using a standard Seurat pipeline to normalize and scale the data which is required as input for sctype.
3. The data is clustered using the standard Seurat pipeline which is also required as input for Sctype.
    1. The Annotation_Guidelines.md which specifies the number of Principal Components (PCs) and Resolution used as this was optimized for each individual dataset. 
4. The scType functions and gene list are loaded into the environment.
    1. The gene list used for the first round of annotations is gs_listv4.xlsx, which consists of genes curated by scType and from the literature.
5. The sctype score is calculated for each cell and each cluster resulting in a dataframe with an assigned label for each cluster.
6. The cell labels are visualized in the Seurat UMAP. 
7. If the cells were labeled as MGE or CGE Interneurons, those cells will be subsetted into another Seurat object which undergoes a new round of clustering.
8. The interneuron subset will be annotated with the gs_listv5.xlsx which contains gene markers for interneuron sub-types. 
    1. Again scType will calculate the score and the ouput will be saved as a text file for further processing.
9. The interneuron cell labels are visualized in the Seurat UMAP.
10. The sctype annotations are saved as a csv file. 

