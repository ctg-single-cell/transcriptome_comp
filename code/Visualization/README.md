# Project Title: Data Visualization
### Date: 11-08-2023
### Author: Fallon Ratner (f.t.ratner@student.vu.nl)
## This repo contains scripts relating to the visualization of annotated scRNAseq data.

# 0. Installation of software and packages
### 
* [R]: version 4.2.3
* R packages (ggplot; dplyr; stringr)

# 1. In the vitro_vivo_comparison.md script, R is used to visualize the annotation output conducted with scType (see sctype_annotation)
 ## In the script the csv output from sctype conducted on the individual in vivo data and integrated in vitro data is utilized. 
 ## First, the in vivo csv files are loaded and the cell type proportion is calculated, then all datasets are combined into one dataframe.
 ## Next, the in vitro csv file is loaded and the cell type proportion is calculated.
 ## The in vitro and in vivo data is combined into one dataframe to be used for visualization.
 ## The annotations of each datasets are subsequently visualized in a bar chart, stacked area chart, rose plot, scatter plot, and violin chart. Statistical analysis is performed using a Wilcoxon ranked sum test. 
