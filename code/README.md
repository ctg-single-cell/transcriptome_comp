# Project Title: transcriptome_comp
## This repo contains scripts relating to analyzing scRNAseq data.
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
### Date: 09-06-2023

## Folders
### Scripts related to annotating cells with the scType tool are found in: code/sctye_annotation
### Scripts related to integrating cells with Scanorama are found in: code/scanorama_integration
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
    
    Herring et al., 2022
        snRNAseq of fetal Prefrontal Cortex (PFC)
        Number of Nuclei in ga22: 11,660 (10,466 cells kept)
        Total Genes ga22: 27,565
        Number of Nuclei in ga24: 10,809 (9,376 cells kept)
        Total Genes ga24: 27,737
        Number of Nuclei in ga34: 8,431 (6,738 cells kept)
        Total Genes ga34: 26,083
        Sequencing Technology: Illumina NextSeq 550
        Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168408
