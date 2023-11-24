# Integrate Herring Datasets with Scanorama
- Date: 24-11-2023
- Author: Fallon Ratner (f.t.ratner@student.vu.nl)
- This README contains instructions on how to run Scanorama on a local computer to integrate Herring et al.,2022 single nuclei RNA-seq of fetal brain tissue.

**NOTE**: Integrating more datasets should be done on a cluster computer.

## What is Scanorama?
- Scanorama is a tool to integrate datasets and correct for batch effects. The method is based on using highly variable genes to compute mutual nearest neighbors among the datasets. Scanorama requires as input: a list of scanpy objects, a list of highly variables genes, and the Scanorama package. After running the computation, Scanorama will produce a corrected gene matrix which can be used to generate a UMAP. The output can then be used for further processing such as annotation with scType.
- For more information about Scanorama please refer to the orginal publication: 
    - Hie, B., Bryson, B. & Berger, B. Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. Nat Biotechnol 37, 685–691 (2019). https://doi.org/10.1038/s41587-019-0113-3
    - Or visit the Scanorama Github: https://github.com/brianhie/scanorama


## 0. Installation of software and packages
### scType pipeline requires the following
* [Python]: version 3.10.9
* Python packages (scanpy; scanorama; numpy)
* [R]: version 4.2.3
* R packages (Seurat; HGNChelper; anndata; dplyr; stringr)

## 1. Downloading the scRNAseq Data
### Download the Data
* Go to the link provided in the following text (Dataset Used) and download this file: 'Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad'
### Dataset Used:
* Herring et al., 2022
    - snRNAseq of fetal Prefrontal Cortex (PFC)
    - Number of Nuclei in ga22: 11,660 (10,466 cells kept)
    - Total Genes ga22: 27,565
    - Number of Nuclei in ga24: 10,809 (9,376 cells kept)
    - Total Genes ga24: 27,737
    - Number of Nuclei in ga34: 8,431 (6,738 cells kept)
    - Total Genes ga34: 26,083
    - Sequencing Technology: Illumina NextSeq 550
    - Link:https://console.cloud.google.com/storage/browser/neuro-dev/Processed_data;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
    - Accessed: 19-04-2023

## 2. Extracting the individual scRNAseq Datasets
- In the extract_herring.py script the downloaded file from the previous step will be used to extract the individual datasets which have been cleaned by the original authors. The output will be a cleaned count matrix in the h5ad file format.

## 3. Integrate with Scanorama
### In the scanorama_herring.ipynb script, the files from the previous step will be used for integration with Scanorama.
1. The h5ad files need to be combined and converted into a scanpy object using the scanpy package.
2. Perform Variance Stabilizing Transformation (VST) which stabilizes the variance of gene expression.
3. Next, highly variable genes are identified in the data.
4. Then a list of the individual datatsets is made with the variable genes calculated in the previous step.
5. Integration and batch correction using the scanorama package is performed.
6. The uncorrected and corrected data is visualized in a UMAP.
7. The output of this script is a scanpy object with the: corrected matrix, UMAP, and metadata. 



