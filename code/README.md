# Project Title: transcriptome_comp
   
## This repo contains scripts relating to analyzing scRNAseq data.

### Creator: Fallon Ratner

### Date: 30-03-2023

## Datasets Used:
    
    Velasco et al., 2019
        scRNAseq of PGP1 (iPSC cell line) - derived 3 & 6 month brain organoids
        Number of Cells in 3 months: 17,774
        Number of Cells in 6 months: 21,213
        Sequencing Technology: Illumina NextSeq 500, 10x Genomics
        Link: https://singlecell.broadinstitute.org/single_cell/study/SCP282/reproducible-brain-organoids#study-download

    Polioudakis et al., 2019
        scRNAseq of 17-18GW fetal neocortex
        Number of Cells: 40,000
        Sequencing Technology: Illumina HiSeq 2500, Drop-Seq
        Link:http://solo.bmap.ucla.edu/shiny/webapp/

## Scripts:
    
    raw_annot_prototype28.03.py
        Python script to open and combine raw gene expression and metadata files from both datasets.
        The output processed files can be used in the following scripts.
    
    processed_annot_prototype29.03.23
        Python script to open and analyze processed data from both datasets.
        Data is further subsetted into dataframes
    
    bar_plot_30.03.23
        R script to plot cell types based on data processed in the previous Python scripts.
