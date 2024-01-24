# Assigning Cell Types to scRNAseq using Author's provided cell labels
### Date: 22-11-2023   
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
### This README contains instructions on how the author's annotations from a few scRNAseq datasets was used. 

# 0. Installation of software and packages
### Extracting the author's annotations and visualizing the results requires the following:
* [Python]: version 3.10.9
* Python packages (scanpy; pandas; numpy; matplotlib)
* [R]: version 4.2.3
* R packages (ggplot2, dplyr, stringr, tidyr, VennDiagram, gridExtra, UpSetR)

# 1. Downloading the scRNAseq Data
### Download the Data
* Go to the links provided in the following text (Dataset Used) and download the specified files
### Dataset Used:
   Velasco et al., 2019
        scRNAseq of PGP1 (iPSC cell line) - derived 3 & 6 month brain organoids
        Number of Cells in 3 months: 17,774
        Number of Cells in 6 months: 21,213
        Sequencing Technology: Illumina NextSeq 500, 10x Genomics
        Link: https://singlecell.broadinstitute.org/single_cell/study/SCP282/reproducible-brain-organoids#study-download
        Download: 'expression_PGP1.3mon.txt.gz', 'expression_PGP1.6mon.txt.gz', 'meta_combined.txt'
        Accessed: 28-03-2023
    
    Polioudakis et al., 2019
        scRNAseq of 17-18GW fetal neocortex
        Number of Cells: 40,000
        Sequencing Technology: Illumina HiSeq 2500, Drop-Seq
        Link:http://solo.bmap.ucla.edu/shiny/webapp/
        Download: Download the Data which will be az ipped folder contaning 'Polio_matrix.csv' & 'cell_metadata.csv'
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

# 2. Read the files, assign the author's cell types, and calculate the proportion.
### In the metadata_annotation.py script the files will be loaded into python and the author's annotations are assigned to each cell. Next, the proportion of each cell type per dataset is calculated as a percentage.

# 3. Visualize the results.
### In the meta_plots.md file the output of the previous script is used to visualize the author's annotations.
