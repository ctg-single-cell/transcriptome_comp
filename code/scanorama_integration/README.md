# Integrate Datasets with Scanorama
- Date: 04-08-2023
- Author: Fallon Ratner (f.t.ratner@student.vu.nl)
- This README contains instructions on how to run Scanorama to integrate single cell RNA-seq from brain organoids.
**NOTE**: Integrating these datasets was done on a cluster computer.

## What is Scanorama?
- Scanorama is a tool to integrate datasets and correct for batch effects. The method is based on using highly variable genes to compute mutual nearest neighbors among the datasets. Scanorama requires as input: a list of scanpy objects, a list of highly variables genes, and the Scanorama package. After running the computation, Scanorama will produce a corrected gene matrix which can be used to generate a UMAP. The output can then be used for further processing such as annotation with scType.
- For more information about Scanorama please refer to the orginal publication: 
- Hie, B., Bryson, B. & Berger, B. Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. Nat Biotechnol 37, 685–691 (2019). https://doi.org/10.1038/s41587-019-0113-3
- Or visit the Scanorama Github: https://github.com/brianhie/scanorama


## 0. Installation of software and packages
### Scanorama pipeline requires the following
* [Python]: version 3.10.9
* Python packages (scanpy; scanorama; numpy; pandas; scipy.io)


## 1. Downloading the scRNAseq Data
### Download the Data
* Go to the links provided in the following text (Dataset Used) and download the files: 
### Datasets Used: (Cells reported post-qc)
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


## 2. Integrate with Scanorama
### In the invitro_int.py script, the files from the previous step will be used for integration with Scanorama.
1. All of the files are made into scanpy objects based on the format provided by the Author's, some of the data needs to be wrangled to access the correct time points and conditions. The Source is added as an observed value to the metadata. 
2. Perform standard pre-processing on the individual datasets.
3. Combine all of the individual datasets into one scanpy object and repeat standard pre-processing steps.
4. Identify highly variable genes in the data.
5. Make a list of the individual datatsets with the variable genes calculated in the previous step.
6. Perform the integration and batch correction using the scanorama package.
7. Visualize the uncorrected and corrected data in a UMAP.
8. Save the scanpy objects as an h5ad file which will store the: corrected matrix, UMAP, and metadata. 
