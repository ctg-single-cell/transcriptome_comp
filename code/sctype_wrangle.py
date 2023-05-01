"""
Created by: Fallon Ratner
Date: 01/05/23
Project: scRNAseq Prototype Analysis - using scType score output from Seurat
"""
#import libraries
import pandas as pd
import numpy as np

########
#Velasco & Polioudakis Datasets
# Read in the three count dataframes
cell_ct_3mon = pd.read_csv("Velasco_3mon_cluster_sctype.csv", index_col=0)
cell_ct_6mon = pd.read_csv("Velasco_6mon_cluster_sctype.csv", index_col=0)
cell_ct_Pol = pd.read_csv("Polio_cluster_sctype.csv", index_col=0)

# Add a new 'Group' column to each DataFrame
cell_ct_3mon['Group'] = '3months'
cell_ct_6mon['Group'] = '6months'
cell_ct_Pol['Group'] = '17-18 GW'

# Calculate the proportion of each cell type count for each dataframe
for df in [cell_ct_3mon, cell_ct_6mon, cell_ct_Pol]:
    total_count = df['ncells'].sum()
    proportions = df['ncells'] / total_count
    df['Proportion'] = proportions
    df['Percentage'] = (proportions * 100)

# Concatenate the three dataframes
combo = pd.concat([cell_ct_3mon, cell_ct_6mon, cell_ct_Pol])

# Save the new dataframe as a csv file
combo.to_csv('proportion_velasco_polio_sctype_cluster.csv')
######################################################
#Herring & Polioudakis Datasets
# Read in the three count dataframes
cell_ct_ga22 = pd.read_csv("ga22_sctype_v3.csv",index_col=0)
cell_ct_ga24 = pd.read_csv("ga24_sctype_v3.csv",index_col=0)
cell_ct_Pol = pd.read_csv("poliou_sctype_v3.csv",index_col=0)

# Add a new 'Group' column to each DataFrame
cell_ct_ga22['Group'] = '22 GA'
cell_ct_ga24['Group'] = '24 GA'
cell_ct_Pol['Group'] = '17-18 GW'

# Calculate the proportion of each cell type count for each dataframe
for df in [cell_ct_ga22, cell_ct_ga24, cell_ct_Pol]:
    total_count = df['ncells'].sum()
    proportions = df['ncells'] / total_count
    df['Proportion'] = proportions
    df['Percentage'] = (proportions * 100)

# Concatenate the dataframes
combo = pd.concat([cell_ct_ga22, cell_ct_ga24, cell_ct_Pol])
combo.head()

#save new df as csv
combo.to_csv('proportion_herring_polio_sctypev3.csv')

#################################################################################
#compare scType with herring gs_list to herring meta
cell_ct_ga22 = pd.read_csv("ga22_meta_qc_cluster_sctype.csv",index_col=0)

# Add a new 'Group' column to each DataFrame
cell_ct_ga22['Group'] = '22 GA scType'
# Calculate the proportion of each cell type count for each dataframe
for df in [cell_ct_ga22]:
    total_count = df['ncells'].sum()
    proportions = df['ncells'] / total_count
    df['Proportion'] = proportions
    df['Percentage'] = (proportions * 100)
total_count = cell_ct_ga22['ncells'].sum()

# read in metadata df
meta = pd.read_csv("proportion_herring_meta_annot.csv",index_col=0)

#subset the dataframe to ga22
meta = meta[meta['Group'] == '22 GA']

#rename the cell types to be the same in both dfs
# Immature PNs = PN_dev
# L5-6_THEMIS
# PNs L2-3 = L2-3_CUX2
# PNs L5/6 TLE4 = L5-6_TLE4
# Immature INs = MGE_dev & CGE_dev
# Mature Astrocytes & Immature Astrocytes = Astro
# PNs L4 = L4_RORB
# OPCs = OPC
# VIP = VIP
# ID2
# SST = SST
# Vas
# Oligo
# Micro
# Poor-Quality
# PV

cell_ct_ga22['type'] = cell_ct_ga22['type'].replace({'Immature PNs': 'PN_dev', 'PNs L2-3': 'L2-3_CUX2', 'PNs L5/6 TLE4': 'L5-6_TLE4', 'Immature INs': 'MGE_dev & CGE_dev','Mature Astrocytes': 'Astro', 'Immature Astrocytes': 'Astro', 'PNs L4': 'L4_RORB', 'OPCs': 'OPC'})

# Concatenate the dataframes
combo = pd.concat([meta, cell_ct_ga22])

#save new df as csv
combo.to_csv('proportion_herring_meta__sctype_annot.csv')