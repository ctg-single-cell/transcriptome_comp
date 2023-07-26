Dataset Information:
#The authors use different filter methods so I will not have the same cell counts as reported for Han & Couturier
01. 
Author: Herring et al., 2022
Type: In Vivo  22, 24, 34 GA, and 2d
Number of Cells: 
22GA: 10,466 cells (15 PCs, 0.8 resolution)
24GA: 9,376 cells (15 PCs, 0.8 res)
34GA: 6,738 cells (15 PCs, 0.8 res)
2d: 9,369 cells (10 PCs, 0.8 res)
File: cleaned_count_matrices.h5ad
Raw or Log Transformed: default count matrix in RNA is raw, downsampled in layer ds_norm_cts
Meta Information: Yes ages included - need to add author name
Pre-Processing:
- Making a Scanpy object: sc.read({file})

02. 
Author: Polioudakis et al., 2019
Type: In Vivo 17-18 GW
Number of Cells: 33,986 post filtering is 33,185 (15 PCs, 0.5 resolution for unfiltered | 15 PCs, 0.8 res for filetered)
File: Polio_matrix.csv
Raw or Log Transformed: Raw counts
Meta Information: Available - need to add age and author info
Pre-processing:
- Making a Scanpy object: sc.read_csv("Polio_matrix.csv")
03. 
Author: Han et al., 2020 - Not using b/c of quality
Type: In Vivo 11-13 GW
brain 3 = 13 GW (M), brain4 = 11 GW (F), brain5 = 12 GW (F), brain 6 = 11 GW (F) 
4 -11GW(F): 2,013 Cells (5 PCs)
6 -11GW(F): 1,126 Cells (5 PCs)
12GW(F): 4,714 Cells (10 PCs) 
13GW: 678 Cells (5 PCs)
File: GSM4008678_Fetal-Brain3_dge.txt, GSM4008679_Fetal-Brain4_dge.txt, GSM4008680_Fetal-Brain5_dge.txt, GSM4008681_Fetal-Brain6_dge.txt
Raw or Log Transformed: Raw counts (dge = digital gene expression)
Meta Information: Not available - need to add age and author info
Link: https://figshare.com/articles/dataset/HCL_DGE_Data/7235471?file=23062979 --> dge_raw_data.tar.gz --> FetalBrain3.rawdge.txt
Pre-processing:
- Making a Scanpy object: sc.read_txt, transpose, need to filter for more than 500 UMIs to get same cells analyzed as authors
04. 
Author: Couturier et al., 2020
Type: In Vivo 13,17,19 GW
13 GW: 7434 --> 6,684 Cells (15 PCs)
17 GW: 10436 --> 8,546 Cells (15 PCs, 0.5 res)
19 GW: 3593 --> 2,689 Cells (15 PCs, 0.5 res)
File: GBM_cellranger_matrix.tar.gz --> HFA567_total.filtered_gene_matrices, HFA570_total.filtered_gene_matrices, HFA571_total.filtered_gene_matrices 
Raw or Log Transformed: Raw 
Meta Information: need to add age and author info
Pre-processing:
- Making a Scanpy object:
05. 
Author: Liu et al., 2023 - Not using b/c of quality
Type: In Vivo 22-23 GW
Number of Cells: 12,577 after filtering 544 cells (10 PCs)
File: hNSPC_raw_counts.h5ad, genes are in ENSEMBL Ids need to change to gene names
Raw or Log Transformed: Raw counts 
Meta Information: Not available - need to add age and author info
Pre-processing:
- Making a Scanpy object: 
06. 
Author: Zhong et al., 2018 - Not using b/c in TPM
Type: In Vivo 8-26 GW
Number of Cells: 2,309
File: GSE104276_all_pfc_2394_UMI_count_NOERCC.xls.gz --> this file does not work, want to try using the TPM file
Raw or Log Transformed: Raw counts (dont' work), try TPM
Meta Information: Not available - need to add age and author info
Pre-processing:
- Making a Scanpy object: First save xls as txt then use sc.read_text
07. 
Author: Fan et al., 2018 - Not using b/c in TPM
Type: In Vivo 22-23 GW
Number of Cells: 4,213
File: GSM2779721_BST_22WF_B1_gene_expression_TPM.txt (multiple files)
Raw or Log Transformed: TPM
Meta Information: Not available - need to add age and author info
Pre-processing:
- Making a Scanpy object: First save xls as txt then use sc.read_text
08. 
Author: Eze et al., 2021
Type: In Vivo 6-10 GW (reportes in CS - Carnegie Stages)
6 GW (CS12 & 13): 8210 --> 5640 Cells (post filter) (15 PCs)
6 GW (CS12): 2299 --> 2224 Cells (post filter) (10 PCs)
7 GW (CS14): 9100 --> 7286 Cells (post filter) (10 PCs)
8 GW (CS15): 9819 --> 9648 Cells (post filter)
9 GW (CS19 & 20): 8162 --> 6985 Cells (15 PCs)
10 GW (CS22): 9865 --> 5298 Cells (15 PCs)
File: counts_exprMatrix.tsv
Raw or Log Transformed: Counts
Meta Information: Available
Pre-processing: In extract_eze.py
08. 
Author: Velasco et al., 2019
Type: In Vitro 3 & 6 months
Number of Cells: 38,987
File: expression_PGP1.3mon.txt, expression_PGP1.6mon.txt
Raw or Log Transformed: Raw counts 
Meta Information: Yes but need to add age and author info

09. 
Author: Trujillo et al., 2019
Type: In Vitro 1, 3, 6, 10 months
Number of Cells: 15,990
File: GSE130238_RAW.tar
Raw or Log Transformed: raw counts
Meta Information: need to add age and author info

10. 
Author: Giandomenico et al., 2019
Type: In Vitro 75 days
Number of Cells: 13,280
File: GSE124174_org75_seuratdata.txt
Raw or Log Transformed: log transformed, maybe raw data in RNA assay - need to check
Meta Information: No - need to add age and author info

11. 
Author: Meng et al., 2022
Type: In Vitro 50 days
Number of Cells: 10,819
File: GSE186814_RAW --> barcodes, genes, and mtx --> u1 = U1M iPS line, u2 = U2F iPS line
Raw or Log Transformed: raw counts
Meta Information: No - need to add age and author info

12. 
Author: Xiang et al.,2017
Type: In Vitro 30, 72, 79 days
Number of Cells: 32,286
File: GSE98201_barcodes.tsv.gz --> barcodes, genes, and mtx
Raw or Log Transformed: Raw counts
Meta Information: No - need to add age and author info

13. 
Author: Fair et al., 2020
Type: In Vitro 93, 140 days on MEAs
Number of Cells: 8,678
File: GSE157019_RAW --> barcodes, genes, and mtx
Raw or Log Transformed: Raw counts
Meta Information: No - need to add age and author info