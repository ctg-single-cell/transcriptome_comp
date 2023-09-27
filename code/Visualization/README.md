# Project Title: Data Visualization
### Date: 11-08-2023
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
## This repo contains scripts relating to the visualization of annotated scRNAseq data.

# 0. Installation of software and packages
### 
* [R]: version 4.2.3
* R packages (ggplot; dplyr; stringr)
* [Python]: version 3.10.9
* Python packages (numpy; pandas; scanpy; matplotlib)

# 1. In the vitro_vivo_comparison.ipynb script, R is used to visualize the annotation output conducted with scType (see sctype_annotation)
* In the script the csv output from sctype conducted on the individual in vivo data and integrated in vitro data is utilized. 
* First, the in vivo csv files are loaded and the cell type proportion is calculated, then all datasets are combined into one dataframe.
* Next, the in vitro csv file is loaded and the cell type proportion is calculated.
* The in vitro and in vivo data is combined into one dataframe to be used for visualization.
* The annotations of each datasets are subsequently visualized in a bar chart, stacked area chart, and rose plot. 

# 2. In the dendro_clust.py script, python is used to generate a dendrogram and dendrogram matrix to investigate the simialrirties between the datsets.
* First, h5ad files of all in vivo datasets and all in vitro datasets are loaded into the program.
* Then, the objects are combined, pre-processed, and normalized.
* A PCA is conducted which is used to compute hierachical clustering.
* Lastly a dendrogram matrix is generated based on Pearson's correlation. 