# Annotate Cells with ScType
### Date: 24-11-2023   
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
### This README contains instructions on how to run scType to annotate single cell/nuclei RNA-seq from in vivo brain tissue and brain organoids.
### _NOTE_: The annotation pipeline has been optimized for fetal brain tissue and brain organoids and may need further optimization for other tissue types.

# What is scType?
### scType is a tool to annotate cells based on unsupervised learning. It utilizes a combination of positive and negative gene markers to annotate cells and clusters based on a computed score. scType requires as input: a gene marker list, the scType functions, a normalized + scaled gene expression matrix, and Seurat computed clusters. After running the computation, scType will output a dataframe consisting of the cluster and the assigned cell label (which has the highest score). This can then be used for further processing and visualization.
### For more information about scType please refer to the orginal publication: 
    - Ianevski, A., Giri, A.K. & Aittokallio, T. Fully-automated and ultra-fast cell-type identification using specific marker combinations
from single-cell transcriptomic data. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w
    - Or visit the scType Github: https://github.com/IanevskiAleksandr/sc-type


# 0. Installation of software and packages
### scType pipeline requires the following
* [Python]: version 3.10.9
* Python packages (Scanpy)
* [R]: version 4.2.3
* R packages (Seurat; HGNChelper; anndata; dplyr; stringr)

# 1. scType requires gene markers in an excel file as input for the annotation scoring algorithm.
 - gs_listv4.xlsx contains cell types and gene markers from the scType databse as well as genes curated from the literature.This is used for the first round of annotations.
 - gs_listv5.xlsx contains interneuron subtype gene markers from the scType databse as well as genes curated from the literature.This is used for the second round of annotations. If MGE and CGE interneurons are classified in the first round they are then used to identify interneuron subtypes.

# Annotation Pipeline:
- In the Prototype folder scripts related to the scType annotation with and without the custom gene list of one of the Herring datasets is found. 
- In the Invivo folder, there are scripts related to the scType annotation of indidvidual fetal brain datasets.
- In the Invitro folder, there are scripts related to the scType annotation of the integrated brain organoid datasets.
