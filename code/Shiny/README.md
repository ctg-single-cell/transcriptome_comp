# Sharing Results in a Shiny App
### Date: 27-09-2023   
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
### This README contains instructions on how to share the results in a Shiny App.

# 0. Installation of software and packages
### scType pipeline requires the following
* [Python]: version 3.10.9
* Python packages (numpy; pandas; scanpy; os; itertools; scipy.stats)
* [R]: version 4.2.3
* R packages (Seurat; dplyr; stringr; ggplot2; grid; shiny; plotly; DT)

# 1. Compute Gene Signatures for each cell type in the fetal brain in the gene_signatures.py script.
* In this script, gene signatures for each cell type in each fetal dataset is calculated.
* First, the individual datasets are pre-processed and normalized.
* Next, the gene markers are determined for each cell type compared to the rest of the cell types in the dataset using the rank genes group function in scanpy. This is a Wilcoxon ranked sum test and it is only calculated if there are at least 10 cells in the given cell type.
* For the in vivo datasets, the results are averaged for datasets pertaining to specific developmental stages:
    * Early Development (11- 13 GW): Han 11GW (batc h1 &2), Han 12GW, Han 13GW, and Couturier 13GW
    * Mid Development (17 - 19 GW): Couturier 17GW, Couturier 19 GW, Liu 17-19 GW, and Polioudakis 17-18 GW
    * Late Development (22 GW - Day 2) = Herring 22GW,  Herring 24GW, Herring 34GW, Herring Day 2
* The gene signature is saved as csv for each developmental stage and the number of cells and genes per dataset is also saved as a csv.
* This process is repreated for the GABA sub-types found in both in vivo and in vitro. 

# 2. Compute Gene Signatures for each cell type in the organoid datasets in the vitro_gene_sig.py script.
* In this script, gene signatures for each cell type in each organoid dataset is calculated.
* First, the h5ad file is loaded with a csv file correpsonding to the Source observation.
* Next, the gene markers are determined for each cell type compared to the rest of the cell types in each unique dataset using the rank genes group function in scanpy. This is a Wilcoxon ranked sum test and it is only calculated if there are at least 10 cells in the given cell type.
* The gene signature is saved as csv for each dataset and the number of cells and genes per dataset is also saved as a csv.
* This process is repreated for the GABA sub-types found in the organoids. 

# 3. Compute Correlation between cell types in vivo aqnd in vitro in the pairwise_comparison.py script.
* A pairwise comaparison is conducted between each cell type from each developmental stage and each organoid dataset.
* First the csv output from the previous scripts are loaded.
* Then, the common genes between in vivo and in vitro are determined. This log fold change of the same genes are compared using the Spearman Correlation.
* The result of this statistical test is saved as a csv and can be found in the Shiny_App folder. 

# 4. Share the results with a Shiny App in the shiny_app.R script.
* In this script the comparison between the in vivo brain and in vitro brain organoids is visualized and shared.
* In the first tab, there is a guide on how to compare cell types between fetal brain and the brain organoids. There are also UMAPs to visualize the annotated cell types per dataset.
* In the side panel the cell type of interest and developmental stage of interest can be selected and the output is in the second tab.Here the reuslts of the pairwise comparison are presented in a table and as a heatmap.
* In the third tab more information about the datasets are presented. 