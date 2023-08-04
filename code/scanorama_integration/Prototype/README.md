# Integrate Datasets with Scanorama
- Date: 04-08-2023
- Author: Fallon Ratner (f.t.ratner@student.vu.nl)
- This README contains instructions on how to run Scanorama  to integrate single cell/nuclei RNA-seq from in vivo brain tissue and in vitro organoids.
**NOTE**: Integrating these datasets was done on a cluster computer.

## What is Scanorama?
- Scanorama is a tool to integrate datasets and correct for batch effects. The method is based on using highly variable genes to compute mutual nearest neighbors among the datasets. Scanorama requires as input: a list of scanpy objects, a list of highly variables genes, and the Scanorama package. After running the computation, Scanorama will produce a corrected gene matrix which can be used to generate a UMAP. The output can then be used for further processing such as annotation with scType.
- For more information about Scanorama please refer to the orginal publication: 
- Hie, B., Bryson, B. & Berger, B. Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. Nat Biotechnol 37, 685–691 (2019). https://doi.org/10.1038/s41587-019-0113-3
- Or visit the Scanorama Github: https://github.com/brianhie/scanorama


## 0. Installation of software and packages
### Scanorama pipeline requires the following
* [Python]: version 3.10.9
* Python packages (scanpy; scanorama; numpy)


## 1. Downloading the scRNAseq Data
### Download the Data
* Go to the links provided in the following text (Dataset Used) and download the files: 
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
Velasco et al., 2019
    scRNAseq of PGP1 (iPSC cell line) - derived 3 & 6 month brain organoids
    Number of Cells in 3 months: 12,990 (post-qc)
    Number of Cells in 6 months: 17,150 (post-qc)
    Sequencing Technology: Illumina NextSeq 500, 10x Genomics
    Link: https://singlecell.broadinstitute.org/single_cell/study/SCP282/reproducible-brain-organoids#study-download
    Download: 'expression_PGP1.3mon.txt.gz', 'expression_PGP1.6mon.txt.gz'
    Accessed: 28-03-2023
Trujillo et al., 2019
    scRNAseq of iPSC derived cortical organoids - 1,3,6, and 10 months
    Number of Cells Trujillo 1 month: 2,639 (post-qc)
    Number of Cells Trujillo 3 months: 2,426 (post-qc)
    Number of Cells Trujillo 6 months: 4,629 (post-qc)
    Number of Cells Trujillo 10 month: 3,063 (post-qc)
    Link:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130238
    Download: GSE130238_RAW.tar
    Accessed: 23-06-2023
Giandomenico et al., 2019
    scRNAseq of ESC derived cerebral organoids - 75 days
    Number of Cells: 7,937 (post-qc)
    Link:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124174
    Download: GSE124174_org75_seuratdata.txt
    Accessed: 27-03-2023
Meng et al., 2022
    scRNAseq of iPSC derived forebrain organoids - 50 days, 2 cells lines: U1M, U2F
    Number of Cells U1M: 10,747 (post-qc)
    Number of Cells U2F: 4,814 (post-qc)
    Link:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186814
    Download: GSE186814_RAW.tar
    Accessed:27-06-2023
Xiang et al., 2017
    scRNAseq of iPSC/ESC derived organoids - 30, 72, and 79 days
    Number of Cells 30 days batch1: 10,023 (post-qc)
    Number of Cells 30 days batch2: 6,317 (post-qc)
    Number of Cells 72 days: 6,210 (post-qc)
    Number of Cells 79 days: 6,427 (post-qc)    
    Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97882
    Download: GSE98201_barcodes.tsv.gz
    Accessed: 23-06-2023
Author: Fair et al., 2020
    scRNAseq of iPSC-derived cerebral organoids - 93, 140 days
    Number of Cells 93 days: 329 (post-qc)
    Number of Cells 93 days: 644 (post-qc)
    Link:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157019
    Download: GSE157019_RAW.tar
    Accessed: 25-06-2023

## 2. Integrate In vivo Datasets or In vivo and In vitro Datasets with Scanorama
### In the invivo_int.py script (in vivo only) and vivo_vitro_int.py (in vivo and in vitro), the listed datasets are used for integration with Scanorama.
1. All of the files are made into scanpy objects based on the format provided by the Author's, some of the data needs to be wrangled to access the correct time points and conditions. The Source is added as an observed value to the metadata. 
2. Perform standard pre-processing on the individual datasets.
3. Combine all of the individual datasets into one scanpy object and repeat standard pre-processing steps.
4. Identify highly variable genes in the data.
5. Make a list of the individual datatsets with the variable genes calculated in the previous step.
6. Perform the integration and batch correction using the scanorama package.
7. Visualize the uncorrected and corrected data in a UMAP.
8. Save the scanpy objects as an h5ad file which will store the: corrected matrix, UMAP, and metadata. 


