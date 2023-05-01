"""
Created by: Fallon Ratner
Date: 01/05/23
Project: scRNAseq Prototype Analysis
"""
#import libraries
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

#In Vitro Dataset: Velasco et al., 2019 PGP 3mo organoids

def process_scanpy_data(exp_file_3mon, exp_file_6mon, metadata_file, metadata_output_file):
    # read the expression data file for 3 months into a scanpy object
    adata_3mon = sc.read_csv(exp_file_3mon, delimiter='\t')
    # transpose the expression data
    adata_3mon_transpose = adata_3mon.transpose()

    # read the expression data file for 6 months into a scanpy object
    adata_6mon = sc.read_csv(exp_file_6mon, delimiter='\t')
    # transpose the expression data
    adata_6mon_transpose = adata_6mon.transpose()

    # read the metadata file into a pandas dataframe
    metadata = pd.read_csv(metadata_file, sep='\t')

    # subset the metadata to only include cells with 3 months or 6 months PGP1
    metadata_PGP1_3mon = metadata[metadata["NAME"].isin(adata_3mon_transpose.obs.index.tolist()]
    metadata_PGP1_6mon = metadata[metadata["NAME"].isin(adata_6mon_transpose.obs.index.tolist())]


    # assign the metadata to the scanpy objects
    adata_3mon_transpose.obs = metadata_PGP1_3mon
    adata_6mon_transpose.obs = metadata_PGP1_6mon

    # save the metadata as a CSV file
    metadata_PGP1_3mon.to_csv(metadata_output_file)
    metadata_PGP1_6mon.to_csv(metadata_output_file)
    
    # access the expression matrix for 3 months
    print("Expression matrix for 3 months:")
    print(adata_3mon_transpose.X)
    print(adata_3mon_transpose.X.shape)

    # access the expression matrix for 6 months
    print("Expression matrix for 6 months:")
    print(adata_6mon_transpose.X)
    print(adata_6mon_transpose.X.shape)

    # access the metadata
    print("Metadata 3months:")
    print(metadata_PGP1_3mon)
    print("Metadata 6months:")
    print(metadata_PGP1_6mon)
#Process the files
process_scanpy_data('expression_PGP1.3mon.txt', 'expression_PGP1.6mon.txt', 'meta_combined_clean.txt', 'metadata_PGP1_3mon_6mon.csv')


############################################################
#VisualizeCell types with a bar chart
def plot_cell_type_counts(metadata):
    # Get the unique cell types and their counts
    cell_type_counts = metadata['CellType'].value_counts().reset_index()
    cell_type_counts.columns = ['CellType', 'Count']

    # Create a bar chart of the cell type counts
    fig, ax = plt.subplots(figsize=(8, 6))
    cell_type_counts.plot.bar(x='CellType', y='Count', rot=0, ax=ax)

    # Set the font size of the tick labels
    ax.set_xticklabels(fontsize=7)

    # Add labels and title
    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Count')
    ax.set_title('Cell Type Counts')

    plt.show()
    
# Plot cell type counts for 3 months data
plot_cell_type_counts(metadata_PGP1_3mon)

# Plot cell type counts for 6 months data
plot_cell_type_counts(metadata_PGP1_6mon)

##################################################################

#In Vivo Dataset: Polioudakis et al., 2019 17-18 GW neocrotex


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

# Combine expression data and metadata into a single AnnData object

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
#Visualize Polioudakis Cell types with a bar chart
#Change cluster to cell type
metadata_df = metadata_df.rename(columns={'Cluster': 'CellType'})
# Count the number of cell types and create a new DataFrame
count = metadata_df['CellTyoe'].value_counts().reset_index()
count.columns = ['CellType', 'Count']# Add a new 'Group' column to Polio df
count['Group'] = '17-18 GW'

#Plot
plot_cell_type_counts(count)


#####################################################################
#Calculate proportions
def calculate_cell_type_percentages(metadata, group):
    # Count the number of each cell type
    count = metadata['CellType'].value_counts().reset_index()
    count.columns = ['CellType', 'Count']

    # Calculate the proportion of each cell type count
    total_count = count['Count'].sum()
    proportions = count['Count'] / total_count

    # Add the proportions as a new column to the DataFrame
    count['Proportion'] = proportions

    # Convert the proportion column to a percentage
    count['Percentage'] = (proportions * 100)

    # Add a column with the group identifier
    count['Group'] = group

    # Save the dataframe to a CSV file
    count.to_csv(f'proportion_{group}.csv', index=False)

    return count


# Call the function for each dataset
counts_PGP1_3mon = calculate_cell_type_percentages(metadata_PGP1_3mon, '3mon')
counts_PGP1_6mon = calculate_cell_type_percentages(metadata_PGP1_6mon, '6mon')
counts_Polio = calculate_cell_type_percentages(metadata_df, '17_18GW')

# Combine the dataframes into one
combined_counts = pd.concat([counts_PGP1_3mon, counts_PGP1_6mon, counts_Polio], ignore_index=True)

# Save the combined dataframe to a CSV file
combined_counts.to_csv('combined_proportions.csv', index=False)

############################################################################
#Select common variables
# Susbet IP (intermediate progenitors)
ip = combined_counts[combined_counts['CellType'].str.contains('IP')]
test =pd.DataFrame(ip, columns=('CellType', 'Group', 'Percentage'))

combo.to_csv('proportion_velasco_poliov2.csv')

#In Polioudakis combine Pg2M and PgS values & rename it Cycling
cy = combined_counts[combined_counts['CellType'].str.contains('Pg')]
#sum cycling variables
new_row = cy.sum(axis=0)
combo2 = combo.append(new_row, ignore_index=True)
combo2['CellType'] = combo2['CellType'].replace('PgSPgG2M', 'Cycling')

cy2 = combo2[combo2['CellType'].str.contains('Cycling')]

test2 =pd.DataFrame(cy2, columns=('CellType', 'Group', 'Percentage'))
test2['Group'] = test2['Group'].replace('17-18 GW17-18 GW', '17-18 GW')

#Save df as csv
test2.to_csv('proportion_velasco_polio_cycling.csv')

#####################################################################################################
# Herring et al., 2022
  #snRNAseq of fetal Prefrontal Cortex (PFC)
    #Number of Nuclei in ga22: 10,466 cells
    #Number of Nuclei in ga24: 9,376 cells 
    #Number of Nuclei in ga34: 6,738 cells

# read in adata
adata = sc.read("Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad")

#read in metadata
meta = pd.read_csv("Processed_data_RNA-all_BCs-meta-data.csv")

#merge dfs based on cell IDs
merged_df = pd.merge(adata, meta, on='Unnamed: 0')

#Look at names in batch
print(meta['batch'].unique())

#subset the dataframe to ga22,24, and 34
meta_22 = meta[meta['batch'] == 'RL2103_ga22_v3']
#save new df as csv
meta_22.to_csv('herring_meta_ga22.csv')

meta_24 = meta[meta['batch'] == 'RL2107_ga24_v3']
#save new df as csv
meta_24.to_csv('herring_meta_ga24.csv')

meta_34 = meta[meta['batch'] == 'RL2121_ga34_v3']
#save new df as csv
meta_34.to_csv('herring_meta_ga34.csv')

#subset dataframe to cell annotations and calculate the proportion
meta_dfs = [meta_22, meta_24, meta_34]

for meta_df in meta_dfs:
    meta_cell = meta_df.loc[:, ['Unnamed: 0', 'age', 'cell_type', 'major_clust', 'sub_clust']]
    # use the new df for calculations
    total_count = meta_cell['major_clust'].sum()
    proportions = meta_cell['major_clust'] / total_count
    meta_cell['Proportion'] = proportions
    meta_cell['Percentage'] = proportions * 100
    if meta_df is meta_22:
        meta_cell['Group'] = '22 GA'
    elif meta_df is meta_24:
        meta_cell['Group'] = '24 GA'
    elif meta_df is meta_34:
        meta_cell['Group'] = '34 GA'
     combo = pd.concat([combo, meta_cell])

combo.to_csv('proportion_herring_meta_annot.csv', index=False)




