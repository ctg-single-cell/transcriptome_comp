"""
Created by: Fallon Ratner
Date: 27/03/23
Project: scRNAseq Prototype Analysis
"""
#import libraries
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

#In Vitro Dataset: Velasco et al., 2019 PGP 3mo organoids

# file path = 'C:/Users/fallo/OneDrive/Documents/Internship_2023/transcriptome_comp/data/Velasco/GSE129519_expression_PGP1.3mon.txt

# read the expression file into an AnnData object
adata = sc.read_csv('expression_PGP1.3mon.txt', delimiter='\t')

adata_transpose = adata.transpose()

##########################################################

# read the metadata file into an AnnData object
metadata = pd.read_csv('meta_combined_clean.txt', sep='\t')

#view the metadata
metadata.head()

# count the number of cells:
metadata.shape[0]

######################################################
# subset the metadata so that it contains only the cells for 3 months PGP1
# get a list of cells that are 3 months PGP1

expression_cell_id = adata_transpose.obs.index.tolist()
len(expression_cell_id)

# subset metadata
metadata_PGP1_3mon = metadata[metadata["NAME"].isin(expression_cell_id)]

metadata_PGP1_3mon.head()

metadata_PGP1_3mon.shape[0]

#Combine expression data and metadata into a scanpy object

adata_transpose.obs = metadata_PGP1_3mon

# access the expression matrix
adata_transpose.X

adata_transpose.X.shape

# access the gene names
adata_transpose.var

# access the metadata
adata_transpose.obs

#save new df as csv

metadata_PGP1_3mon.to_csv('metadata_PGP1_3mon.csv')
############################################################
#Visualize 3month Cell types with a bar chart

print(metadata_PGP1_3mon['CellType'].unique())


# Count the number of cell types and create a new DataFrame
counts = metadata_PGP1_3mon['CellType'].value_counts().reset_index()
counts.columns = ['CellType', 'Count']

# Display the new DataFrame
print(counts)

#quick bar chart
ax = counts.plot.bar(x='CellType', y='Count', rot=0)

##Matplotlib
# Set up the figure
fig, ax = plt.subplots(figsize=(8, 6))  # Set the size of the plot

# Create the bar chart
counts.plot.bar(x='CellType', y='Count', rot=0, ax=ax)

# Set the font size of the tick labels
labels = ['CPNs', 'Immature\nCPNs', 'CFuPNs', 'Immature\nPNs', 
          'IPCs/Immature\nPNs', 'Immature\nCFuPNs', 'oRG', 'RG', 'Cycling']
ax.set_xticklabels(labels, fontsize = 7)

# Add labels and title
ax.set_xlabel('Cell Type')
ax.set_ylabel('Count')
ax.set_title('Cell Type Counts Velasco 3mo Organoid')

# Display the chart
plt.show()
#################################################################
##################################################################
#In Vitro Dataset: Velasco et al., 2019 PGP 6mo organoids

# file path = C:\Users\fallo\OneDrive\Documents\Internship_2023\transcriptome_comp\data\Velasco\expression_PGP1.6mon.txt

# read the expression file into an AnnData object
adata = sc.read_csv('expression_PGP1.6mon.txt', delimiter='\t')

adata_transpose = adata.transpose()

##########################################################

# read the metadata file into an AnnData object
metadata = pd.read_csv('meta_combined_clean.txt', sep='\t')

#view the metadata
metadata.head()

# count the number of cells:
metadata.shape[0]

######################################################
# subset the metadata so that it contains only the cells for 6 months PGP1
# get a list of cells that are 6 months PGP1

expression_cell_id = adata_transpose.obs.index.tolist()
len(expression_cell_id)

# subset metadata
metadata_PGP1_6mon = metadata[metadata["NAME"].isin(expression_cell_id)]

metadata_PGP1_6mon.head()

metadata_PGP1_6mon.shape[0]

#Combine expression data and metadata into a scanpy object

adata_transpose.obs = metadata_PGP1_6mon

# access the expression matrix
adata_transpose.X

adata_transpose.X.shape

# access the gene names
adata_transpose.var

# access the metadata
adata_transpose.obs

#save new df as csv

metadata_PGP1_6mon.to_csv('metadata_PGP1_6mon.csv')
############################################################
#Visualize 6 motnhs Cell types with a bar chart

print(metadata_PGP1_6mon['CellType'].unique())


# Count the number of cell types and create a new DataFrame
counts2 = metadata_PGP1_6mon['CellType'].value_counts().reset_index()
counts2.columns = ['CellType', 'Count']

# Display the new DataFrame
print(counts2)

#quick bar chart
counts2.plot.bar(x='CellType', y='Count', rot=0)

##Matplotlib
# Set up the figure
fig, ax = plt.subplots(figsize=(8, 6))  # Set the size of the plot

# Create the bar chart
counts2.plot.bar(x='CellType', y='Count', rot=0, ax=ax)

# Set the font size of the tick labels
labels = ['oRG','Astroglia','Immature\nCPNs','Immature\nInterneurons',
          'Cycling', 'CPNs', 'Immature\nPNs', 'Ventral\nPrecursors', 
          'Unknown', 'RG', 'IPCs']
ax.set_xticklabels(labels, fontsize = 7)

# Add labels and title
ax.set_xlabel('Cell Type')
ax.set_ylabel('Count')
ax.set_title('Cell Type Counts Velasco 6mo Organoid')

# Display the chart
plt.show()

##############################################################
#Create a stacked bar chart for Velasco 3mo and 6mo

# Add a new 'Group' column to each DataFrame
counts['Group'] = '3months'
counts2['Group'] = '6months'

# Concatenate the two dataframes
combo = pd.concat([counts, counts2])
combo.head()

#save new df as csv

combo.to_csv('counts_3and6mon.csv')
# Pivot the dataframe to create a stacked bar chart
pivot = combo.pivot(index='CellType', columns='Group', values='Count')

# Create the stacked bar chart
ax = pivot.plot(kind='bar', stacked=True)

# Add labels and title
ax.set_xlabel('Cell Type')
ax.set_ylabel('Count')
ax.set_title('Cell Type Counts by Group Velasco et al., 2019')

# Display the chart
plt.show()
#################################################################
##################################################################

#In Vivo Dataset: Polioudakis et al., 2019 17-18 GW neocrotex

# file path = 'C:/Users/fallo/OneDrive/Documents/Internship_2023/transcriptome_comp/data/Polio/Polio_matrix.csv'


# Load expression data and metadata as pandas dataframes
expr_df = pd.read_csv("Polio_matrix.csv", index_col=0)
metadata_df = pd.read_csv("cell_metadata.csv", index_col=0)

expr_df_transposed = expr_df.T

# Find overlapping samples
common_samples = set(expr_df_transposed.index) & set(metadata_df.index)

# Subset expression data and metadata to match common samples
expr_df_transposed = expr_df_transposed.loc[common_samples]
metadata_df = metadata_df.loc[common_samples]

# Convert expression data and metadata to AnnData objects
adata_expr = sc.AnnData(expr_df_transposed)

#adata_metadata = sc.AnnData(metadata_df)

# Combine expression data and metadata into a single AnnData object
#adata = adata_expr.concatenate(adata_metadata, join="inner", batch_key="metadata")

###
expression_cell_id = adata_expr.obs.index.tolist()
len(expression_cell_id)


metadata.shape[0]

#Combine expression data and metadata into a scanpy object

adata_expr.obs = metadata_df

# access the expression matrix
adata_expr.X

adata_expr.X.shape

# access the gene names
adata_expr.var

# access the metadata
adata_expr.obs

#save new df as csv

metadata_df.to_csv('metadata_Polioudakis.csv')

############################################################
#Visualize 17-18GW Cell types with a bar chart

print(metadata_df['Cluster'].unique())


# Count the number of cell types and create a new DataFrame
count = metadata_df['Cluster'].value_counts().reset_index()
count.columns = ['Cluster', 'Count']

# Display the new DataFrame
print(count)

#Quick bar chart
count.plot.bar(x='Cluster', y='Count', rot=0)

##Matplotlib
# Set up the figure
fig, ax = plt.subplots(figsize=(8, 6))  # Set the size of the plot

# Create the bar chart
count.plot.bar(x='Cluster', y='Count', rot=0, ax=ax)

# Set the font size of the tick labels
labels = ['ExN','ExM','ExDp1','InMGE',
          'vRG', 'ExM-U', 'oRG', 'IP', 
          'PgS', 'InCGE', 'PgG2M','Per', 'End' ,'Mic', 'OPC', 'ExDp2']
ax.set_xticklabels(labels, fontsize = 7)

# Add labels and title
ax.set_xlabel('Cell Type')
ax.set_ylabel('Count')
ax.set_title('Cell Type Counts Poliodaskus 17-18 GW Neocortex')

# Display the chart
plt.show()

#######################################################################
#Combine Velasco Polioudakis into one df

# Add a new 'Group' column to Polio df
count['Group'] = '17-18 GW'

#Change cluster to cell type
count = count.rename(columns={'Cluster': 'CellType'})

# Concatenate the two dataframes
combo2 = pd.concat([count, counts, counts2])
combo2.head()

#save new df as csv
combo2.to_csv('counts_velasco_polio.csv')

# Pivot the dataframe 
pivot2 = combo2.pivot(index='CellType', columns='Group', values='Count')

# Create a stacked bar chart
ax = pivot2.plot(kind='bar', stacked=True)

# Add labels and title
ax.set_xlabel('Cell Type')
ax.set_ylabel('Count')
ax.set_title('Cell Type Counts by Group Velasco et al., 2019 & Polioudakis et al., 2019')

# Display the chart
plt.show()

#####################################################################
#Calculate the cell proportions as new column
# Calculate the sum of all cell counts
total_count = count['Count'].sum()

# Calculate the proportion of each cell type count
proportions = count['Count'] / total_count

# Add the proportions as a new column to the DataFrame
count['Proportion'] = proportions

# Convert the proportion column to a percentage
count['Percentage'] = (proportions * 100).apply(np.ceil).astype(int)

# Display the updated DataFrame
print(count)

# Calculate the sum of all cell counts
total_counts = counts['Count'].sum()

# Calculate the proportion of each cell type count
proportions2 = counts['Count'] / total_counts

# Add the proportions as a new column to the DataFrame
counts['Proportion'] = proportions2

# Convert the proportion column to a percentage
counts['Percentage'] = (proportions2 * 100).apply(np.ceil).astype(int)

# Display the updated DataFrame
print(counts)

# Calculate the sum of all cell counts
total_counts2 = counts2['Count'].sum()

# Calculate the proportion of each cell type count
proportions3 = counts2['Count'] / total_counts2

# Add the proportions as a new column to the DataFrame
counts2['Proportion'] = proportions3

# Convert the proportion column to a percentage
counts2['Percentage'] = (proportions3 * 100).apply(np.ceil).astype(int)

# Display the updated DataFrame
print(counts2)

# Concatenate the two dataframes
combo3 = pd.concat([count, counts, counts2])
combo3.head()

#save new df as csv
combo3.to_csv('proportion_velasco_polio.csv')

# Pivot the dataframe 
pivot3 = combo3.pivot(index='CellType', columns='Group', values='Percentage')

# Create a stacked bar chart
ax = pivot3.plot(kind='bar', stacked=True)

# Add labels and title
ax.set_xlabel('Cell Type')
ax.set_ylabel('Proportion (%)')
ax.set_title('Cell Type Percentage by Group Velasco et al., 2019 & Polioudakis et al., 2019')

##############################################################
#Plot only oRG (outer radial glia)
org = pd.DataFrame(pivot3.loc['oRG'])

# Create the stacked bar chart
ax = org.plot(kind='bar')

# Add labels and title
ax.set_xlabel('Cell Type')
ax.set_ylabel('Proportion (%)')
ax.set_title('Outer Radial Glia Percentage by Group Velasco et al., 2019 & Polioudakis et al., 2019')
