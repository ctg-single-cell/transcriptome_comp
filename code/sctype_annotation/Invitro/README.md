# Annotate Cells with ScType
### Date: 09-06-2023   
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
### This README contains instructions on how to run scType on a local computer to annotate single cell/nuclei RNA-seq from in vivo brain tissue.
### _NOTE_: The annotation pipeline has been optimized for fetal brain tissue and may need further optimization for other tissue types.


# 0. Installation of software and packages
### scType pipeline requires the following
* [Python]: version 3.10.9
* Python packages (Scanpy)
* [R]: version 4.2.3
* R packages (Seurat; HGNChelper; anndata; dplyr; stringr)

# 1. Downloading the scRNAseq Data
### Download the Data
* Go to the link provided in the following text (Dataset Used) and download this file: 'Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad'
### Datasets Used: (cells reported as post-qc)
Velasco et al., 2019
    scRNAseq of PGP1 (iPSC cell line) - derived 3 & 6 month brain organoids
    Number of Cells in 3 months: 12,990
    Number of Cells in 6 months: 17,150
    Sequencing Technology: Illumina NextSeq 500, 10x Genomics
    Link: https://singlecell.broadinstitute.org/single_cell/study/SCP282/reproducible-brain-organoids#study-download
    Download: 'expression_PGP1.3mon.txt.gz', 'expression_PGP1.6mon.txt.gz'
    Accessed: 28-03-2023
Trujillo et al., 2019
    scRNAseq of iPSC derived cortical organoids - 1,3,6, and 10 months
    Number of Cells Trujillo 1 month: 2,639
    Number of Cells Trujillo 3 months: 2,426
    Number of Cells Trujillo 6 months: 4,629
    Number of Cells Trujillo 10 month: 3,063
    Link:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130238
    Download: GSE130238_RAW.tar
    Accessed: 23-06-2023
Giandomenico et al., 2019
    scRNAseq of ESC derived cerebral organoids - 75 days
    Number of Cells: 7,937
    Link:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124174
    Download: GSE124174_org75_seuratdata.txt
    Accessed: 27-03-2023
Madhavan et al., 2018
    scRNAseq of iPSC derived oligo-cortical spheroids - 12 weeks
    Number of Cells: 2,931
    Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110006
    Download: GSE110006_RAW.tar
    Accessed: 27-03-2021
Xiang et al., 2017
    scRNAseq of iPSC/ESC derived organoids - 30, 72, and 79 days
    Number of Cells 30 days batch1: 10,023
    Number of Cells 30 days batch2: 6,317
    Number of Cells 72 days: 6,210
    Number of Cells 79 days: 6,427   
    Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97882
    Download: GSE98201_barcodes.tsv.gz
    Accessed: 23-06-2023
Author: Fair et al., 2020
    scRNAseq of iPSC-derived cerebral organoids - 93, 140 days
    Number of Cells 93 days: 329
    Number of Cells 93 days: 644
    Link:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157019
    Download: GSE157019_RAW.tar
    Accessed: 25-06-2023
Author: Popova et al., 2021
    scRNAseq of iPSC-derived cortical organoids with or without primary microglia - 7 weeks
    Number of Cells batch1: 223
    Number of Cells batch2: 446
    Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180945
    Download: GSE180945_RAW.tar
    Accessed: 31-07-2023
Author: Bhaduri et al., 2020
    scRNAseq of ESC-derived forebrain organoid - 3,5,8,and 10 weeks (H1S=least directed, H1X=most directed)
    Number of cells H1S_3w: 4,383
    Number of cells H1X_3w: 5,545
    Number of cells H1S_5w: 2,217
    Number of cells H1X_5w: 4,648
    Number of cells H1S_8w: 6,670 
    Number of cells H1X_8w: 1,452
    Number of cells H1S_10w: 18,591
    Number of cells H1X_10w: 5,499
    Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132672
    Download: GSE132672_allorganoids_withnew_matrix.txt.gz
    Accessed: 31-07-2023

# 2. Annotate the Integaret Brain Organoid Datasets with scType
### In the vitro_int_sctype.R script, the Scanorama integrated data is annotated with the Sctype tool.
1. The h5ad file needs to be converted into a Seurat object which is done with the anndata and Seurat packages.
2. The data is pre-processed using a standard Seurat pipeline to normalize and scale the data which is required as input for Sctype.
3. The data is clustered using the standard Seurat pipeline which is also required as input for Sctype.
4. The scType functions and gene list are loaded into the environment.
    1. The gene list used is gs_listv4.xlsx, which consists of genes curated by scType and from the literature.
5. The sctype score is calculated for each cell and each cluster resulting in a dataframe with an assigned label for each cluster.
    1. The score df is saved as a text file for each format and can be used for further processing.
6. The cell labels are visualized in the Seurat UMAP.
7. The cells labeled as MGE or CGE Intenreurons are subsetted into another Seurat object which undergoes a new round of clustering.
8. The interneuron subset will be annotated with the gs_listv5.xlsx which contains gene markers for interneuron sub-types. 
    1. Again scType will calculate the score and the ouput will be saved as a text file for further processing.
9. The interneuron cell labels are visualized in the Seurat UMAP. 

# 3. Visualize scType Ouput with Bar Charts
### In the bar_plot_int.R script, the files from the previous step are used to generate bar charts.
1. Read in the output files and calculate the proportion as a percentage.
2. Not all cell types are found in each dataset, so those values will be assigned to 0.
3. Visualize the percentage of cell types for each dataset in a bar chart. 

 