"""
Created by: Fallon Ratner
Date: 01/05/23
Project: scRNAseq Prototype Analysis
"""
#import libraries
import pandas as pd
import scanpy as sc

def process_scanpy_data(exp_file_3mon, exp_file_6mon, metadata_file, metadata_output_file):
    adata_3mon = sc.read_csv(exp_file_3mon, delimiter='\t').transpose()
    adata_6mon = sc.read_csv(exp_file_6mon, delimiter='\t').transpose()

    metadata = pd.read_csv(metadata_file, sep='\t')

    metadata_PGP1_3mon = metadata[metadata["NAME"].isin(adata_3mon.obs.index.tolist())]
    metadata_PGP1_6mon = metadata[metadata["NAME"].isin(adata_6mon.obs.index.tolist())]

    adata_3mon.obs = metadata_PGP1_3mon
    adata_6mon.obs = metadata_PGP1_6mon

    metadata_PGP1_3mon.to_csv(metadata_output_file.format("3mon"), index=False)
    metadata_PGP1_6mon.to_csv(metadata_output_file.format("6mon"), index=False)

    print("Expression matrix for 3 months:")
    print(adata_3mon.X)
    print(adata_3mon.X.shape)

    print("Expression matrix for 6 months:")
    print(adata_6mon.X)
    print(adata_6mon.X.shape)

    print("Metadata 3 months:")
    print(metadata_PGP1_3mon)
    print("Metadata 6 months:")
    print(metadata_PGP1_6mon)

# Process the Velasco In Vitro Dataset
process_scanpy_data('expression_PGP1.3mon.txt', 'expression_PGP1.6mon.txt', 'meta_combined_clean.txt', 'metadata_PGP1_{}_6mon.csv')

# Process the Polioudakis In Vivo Dataset
expr_df = pd.read_csv("Polio_matrix.csv", index_col=0).T
metadata_df = pd.read_csv("cell_metadata.csv", index_col=0)

common_samples = set(expr_df.index) & set(metadata_df.index)
expr_df = expr_df.loc[common_samples]
metadata_df = metadata_df.loc[common_samples]

adata_expr = sc.AnnData(expr_df)
adata_expr.obs = metadata_df

metadata_df.to_csv('metadata_Polioudakis.csv', index=False)

# Calculate proportions
def calculate_cell_type_percentages(metadata, group):
    count = metadata['CellType'].value_counts().reset_index()
    count.columns = ['CellType', 'Count']
    total_count = count['Count'].sum()
    proportions = count['Count'] / total_count
    count['Proportion'] = proportions
    count['Percentage'] = proportions * 100
    count['Group'] = group
    count.to_csv(f'proportion_{group}.csv', index=False)
    return count

counts_PGP1_3mon = calculate_cell_type_percentages(metadata_PGP1_3mon, '3mon')
counts_PGP1_6mon = calculate_cell_type_percentages(metadata_PGP1_6mon, '6mon')
counts_Polio = calculate_cell_type_percentages(metadata_df, '17_18GW')

combined_counts = pd.concat([counts_PGP1_3mon, counts_PGP1_6mon, counts_Polio], ignore_index=True)
combined_counts.to_csv('proportion_meta_velasco_polio.csv', index=False)

# Process the Herring In vivo Datasets
adata = sc.read("Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad")
meta = pd.read_csv("Processed_data_RNA-all_BCs-meta-data.csv")
merged_df = pd.merge(adata, meta, on='Unnamed: 0')

meta_22 = meta[meta['batch'] == 'RL2103_ga22_v3']
meta_22.to_csv('herring_meta_ga22.csv')

meta_24 = meta[meta['batch'] == 'RL2107_ga24_v3']
meta_24.to_csv('herring_meta_ga24.csv')

meta_34 = meta[meta['batch'] == 'RL2121_ga34_v3']
meta_34.to_csv('herring_meta_ga34.csv')

meta_dfs = [meta_22, meta_24, meta_34]
combo = pd.DataFrame()

for meta_df in meta_dfs:
    meta_cell = meta_df.loc[:, ['Unnamed: 0', 'age', 'cell_type', 'major_clust', 'sub_clust']]
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




