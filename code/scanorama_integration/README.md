# Integrate Datasets with Scanorama
- Date: 09-06-2023
- Author: Fallon Ratner (f.t.ratner@student.vu.nl)
- This README contains instructions on how to run Scanorama on a local computer to integrate single cell/nuclei RNA-seq from in vivo brain tissue.
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
* R packages (Seurat; HGNChelper; zellkonverter; anndata; dplyr; stringr)

## 1. Downloading the scRNAseq Data
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

## 2. Extracting the individual scRNAseq Datasets
### In the extract_herring.py script the downloaded file from the previous step, will be used to extract the individual datasets which have been cleaned by the original authors. The output will be a cleaned count matrix in the h5ad file format.

## 3. Integrate with Scanorama
### In the scanorama_int_herring.py script, the files from the previous step will be used for integration with Scanorama.
1. The h5ad files needs to be combined and converted into a scanpy object using the scanpy package.
2. Perform Variance Stabilizing Transformation (VST) which stabilizes variance of gene expression.
3. Identify highly variable genes in the data.
4. Make a list of the individual datatsets with the variable genes calculated in the previous step.
5. Perform the integration and batch correction using the scanorama package.
6. Visualize the uncorrected and corrected data in a UMAP.
7. Save the scanpy objects as an h5ad file which will store the: corrected matrix, UMAP, and metadata. 

## 4. Annotate Integrated Data with scType and Visualize
### In the int_sctype_herr.R script, the files from the previous step are used to annotate cells and generate a bar chart.
1. Read in the h5ad file and convert to a Seurat object with the zellkonverter package.
2. Visualize the UMAP to check if it's the same as the one generated in Python.
3. The data is pre-processed using a standard Seurat pipeline to normalize and scale the data which is required as input for Sctype.
4. The data is clustered using the standard Seurat pipeline which is also required as input for Sctype.
5. The scType functions and gene list are loaded into the environment.
    1. The gene list used is gs_listv4.xlsx, which consists of genes curated by scType and from the literature.
6. The sctype score is calculated for each cell and each cluster resulting in a dataframe with an assigned label for each cluster.
7. The cell labels are visualized in the Scanroama generated UMAP.
8. The cluster annotation corresponding to the different ages is selected and stored in a new dataframe.
    1. This df is saved as a text file and can be used for further processing.
9. Calculate the proportion of cell types as a percentage for each cell in each age group.
10. Not all cell types are found in each dataset, so those values will be assigned to 0.
11. Visualize the percentage in a bar chart. 

## 5. Integrate More Datasets with Scanorama
### In the scripts: invivo_int2.py, invitro_int.py, vivo_vitro_int.py, and vivo_vitro_int2.py multiple datasets are integrated with Scanroma.
1. The different files need to be made into scanpy objects based on their file type. Then the type of dataset and first author is added to the metadata. 
2. Perform standard scanpy pre-processing to the individual datasets. 
4. Combine the datasets into one scanpy object and perform standard scanpy pre-processing.
4. Identify highly variable genes in the combined data.
5. Make a list of the individual datatsets with the variable genes calculated in the previous step.
6. Perform the integration and batch correction using the scanorama package. (This is done on snellius in the vivo_vitro_int2.py script)
7. Visualize the uncorrected and corrected data in a UMAP.
8. Save the scanpy objects as an h5ad file which will store the: corrected matrix, UMAP, and metadata.
