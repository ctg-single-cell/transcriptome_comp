# Assigning Cell Types to scRNAseq with Cell Typist
### Date: 09-06-2023   
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
### This README contains instructions on how to use Cell Typist to annotate cell types in scRNAseq.

# What is Cell Typist?
### Cell Typist is a a reference atlas of scRNAseq data from many different tissues. The fetal brain atlas consists of the following time points: 5-14 PCW. The input consists of: the developing brain atlas (129 cell types) and normalized count matrices (processed with scanpy). The ouput is the cell identity prediction and decision score computed by the Cell Typist algorithm. This can be used for further processing
### For more information about Cell Typist please refer to the orginal publication: 
### Xu C., Prete M., Webb S., Jardine L.,Stewart B., Hoo R., He P.,Teichmann S. Automatic cell type harmonization and integration across Human Cell Atlas datasets. bioRxiv 2023.05.01.538994; doi: https://doi.org/10.1101/2023.05.01.538994
### Or visit the Scanorama Github: https://github.com/Teichlab/celltypist
### The source for the developing brain model: Braun E., Danan-Gotthold M., Borm L., Vinsland E., Lee K., Lönnerberg P., Hu L., Li X., He X., Andrusivová Z., Lundeberg J., Arenas E., Barker R., Sundström E., Linnarsson . Comprehensive cell atlas of the first-trimester developing human brain. bioRxiv 2022.10.24.513487; doi: https://doi.org/10.1101/2022.10.24.513487

# 0. Installation of software and packages
### scType pipeline requires the following
* [Python]: version 3.10.9
* Python packages (scanpy; pandas; numpy; matplotlib; celltypist)
* [R]: version 4.2.3
* R packages (ggplot2; dplyr)

# 1. Downloading the scRNAseq Data
### Download the Data
### Download the Data
* Go to the links provided in the following text (Dataset Used) and download the specified files
### Dataset Used:
 Velasco et al., 2019
        scRNAseq of PGP1 (iPSC cell line) - derived 3 & 6 month brain organoids
        Number of Cells in 3 months: 17,774
        Number of Cells in 6 months: 21,213
        Sequencing Technology: Illumina NextSeq 500, 10x Genomics
        Link: https://singlecell.broadinstitute.org/single_cell/study/SCP282/reproducible-brain-organoids#study-download
        Download: 'expression_PGP1.3mon.txt.gz', 'expression_PGP1.6mon.txt.gz'
        Accessed: 28-03-2023
    
    Polioudakis et al., 2019
        scRNAseq of 17-18GW fetal neocortex
        Number of Cells: 40,000
        Sequencing Technology: Illumina HiSeq 2500, Drop-Seq
        Link:http://solo.bmap.ucla.edu/shiny/webapp/
        Download: Download the Data which will be az ipped folder contaning 'Polio_matrix.csv' 
        Accessed:28-03-2023
    
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
        Download this file: 'Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad'
        Accessed: 19-04-2023

# 2. Extracting the individual Herring Datasets
### In the extract_herring.py script the downloaded file from the previous step, will be used to extract the individual datasets which have been cleaned by the original authors. The output will be a cleaned count matrix in the h5ad file format.

# 3. Perform Cell Typist Annoation
### In the cell_typist_annot.py script, the referece model will be used to annotate the cells.
1. Select the developing brain model from cell typist to use as the reference.
2. Pre-process the data using a standard Scanpy pipeline.
3. Predict the cell types using the model and Cell Typist scoring algorithm.
4. Calculate the proportion of cells as a perecentage.
5. The output can be furtherd used for other processes like visualization. 

# 4. Visualize Cell Typist Annoation
### In the cell_typist_plot.R script, the output from the previous script will be used to visualize the annotations with a bar chart.
* 2 bar charts are generated: one with all cell types assigned, and one with the top 10 cell types