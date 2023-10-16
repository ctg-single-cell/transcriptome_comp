#Created by Fallon Ratner on 20-09-23
import scanpy as sc
import pandas as pd
import numpy as np
# Read the h5ad file and source names
adata = sc.read("vitro_annot_0908.h5ad")
source_names = pd.read_csv("vitro_annot.csv", index_col=0)

# Ensure that the 'x' column is of type string
source_names['x'] = source_names['x'].astype(str)

# Replace the altered column with the original source names
adata.obs['Source'] = source_names['x'].values

# Unique datasets from the 'Source' column
unique_datasets = adata.obs['Source'].unique().tolist()

# List of cell types
cell_types = adata.obs['customclassif'].unique().tolist()

# Initialize an empty dataframe to store cell and gene counts
cell_gene_counts = pd.DataFrame(columns=['Cell_Type', 'Num_Cells', 'Num_Genes'])

for dataset in unique_datasets:
    
    # Subset the data for the current dataset
    dataset_subset = adata[adata.obs['Source'] == dataset].copy()  # Use copy() to ensure independent subsets
    
    for target_cell_type in cell_types:
        
        # Check if the subset contains the target cell type
        if target_cell_type in dataset_subset.obs['customclassif'].unique().tolist():
            
            
            # Create a temporary binary category for differential expression
            dataset_subset.obs['temp_binary_group'] = np.where(dataset_subset.obs['customclassif'] == target_cell_type, target_cell_type, 'other')
            
            target_cells_count = np.sum(dataset_subset.obs['temp_binary_group'] == target_cell_type)
            other_cells_count = np.sum(dataset_subset.obs['temp_binary_group'] == 'other')
            
            if target_cells_count > 10 and other_cells_count > 10:
                
                # Differential expression analysis
                sc.tl.rank_genes_groups(dataset_subset, groupby='temp_binary_group', groups=[target_cell_type], reference='other', method='wilcoxon')
                
                # Collecting results
                results_df = pd.DataFrame({
                    'genes': dataset_subset.uns['rank_genes_groups']['names'][target_cell_type],
                    'logfoldchanges': dataset_subset.uns['rank_genes_groups']['logfoldchanges'][target_cell_type],
                    'pvals': dataset_subset.uns['rank_genes_groups']['pvals'][target_cell_type],
                    'pvals_adj': dataset_subset.uns['rank_genes_groups']['pvals_adj'][target_cell_type]
                })
                
                # Capture number of cells and genes for this cell type
                selected_subset = adata[adata.obs['customclassif'] == target_cell_type]
                num_cells = selected_subset.shape[0]
                num_genes = selected_subset.shape[1]

                # Append to the dataframe using pd.concat
                new_row = pd.DataFrame({'Cell_Type': [target_cell_type], 'Num_Cells': [num_cells], 'Num_Genes': [num_genes]})
                cell_gene_counts = pd.concat([cell_gene_counts, new_row], ignore_index=True)

                # Save results
                output_filename = f"gene_signature_{dataset}_{target_cell_type}.csv"
                results_df.to_csv(output_filename, index=False)
            
            # Clean up the temporary binary group
            del dataset_subset.obs['temp_binary_group']

cell_gene_counts.to_csv("cell_gene_counts_organoids.csv") 
#########################################################################################################################################
# Read the h5ad file and source names
adata = sc.read("vitro_annot_gaba_0908.h5ad")
source_names = pd.read_csv("vitro_annot_gaba_.csv", index_col=0)

# Ensure that the 'x' column is of type string
source_names['x'] = source_names['x'].astype(str)

# Replace the altered column with the original source names
adata.obs['Source'] = source_names['x'].values

# Unique datasets from the 'Source' column
unique_datasets = adata.obs['Source'].unique().tolist()

# List of cell types
cell_types = adata.obs['customclassif'].unique().tolist()

# Initialize an empty dataframe to store cell and gene counts
cell_gene_counts = pd.DataFrame(columns=['Cell_Type', 'Num_Cells', 'Num_Genes'])

for dataset in unique_datasets:
    
    # Subset the data for the current dataset
    dataset_subset = adata[adata.obs['Source'] == dataset].copy()  # Use copy() to ensure independent subsets
    
    for target_cell_type in cell_types:
        
        # Check if the subset contains the target cell type
        if target_cell_type in dataset_subset.obs['customclassif'].unique().tolist():
            
            
            # Create a temporary binary category for differential expression
            dataset_subset.obs['temp_binary_group'] = np.where(dataset_subset.obs['customclassif'] == target_cell_type, target_cell_type, 'other')
            
            target_cells_count = np.sum(dataset_subset.obs['temp_binary_group'] == target_cell_type)
            other_cells_count = np.sum(dataset_subset.obs['temp_binary_group'] == 'other')
            
            if target_cells_count > 10 and other_cells_count > 10:
                
                # Differential expression analysis
                sc.tl.rank_genes_groups(dataset_subset, groupby='temp_binary_group', groups=[target_cell_type], reference='other', method='wilcoxon')
                
                # Collecting results
                results_df = pd.DataFrame({
                    'genes': dataset_subset.uns['rank_genes_groups']['names'][target_cell_type],
                    'logfoldchanges': dataset_subset.uns['rank_genes_groups']['logfoldchanges'][target_cell_type],
                    'pvals': dataset_subset.uns['rank_genes_groups']['pvals'][target_cell_type],
                    'pvals_adj': dataset_subset.uns['rank_genes_groups']['pvals_adj'][target_cell_type]
                })
                
                # Capture number of cells and genes for this cell type
                selected_subset = adata[adata.obs['customclassif'] == target_cell_type]
                num_cells = selected_subset.shape[0]
                num_genes = selected_subset.shape[1]

                # Append to the dataframe using pd.concat
                new_row = pd.DataFrame({'Cell_Type': [target_cell_type], 'Num_Cells': [num_cells], 'Num_Genes': [num_genes]})
                cell_gene_counts = pd.concat([cell_gene_counts, new_row], ignore_index=True)

                # Save results
                output_filename = f"gene_signature_{dataset}_{target_cell_type}.csv"
                results_df.to_csv(output_filename, index=False)
            
            # Clean up the temporary binary group
            del dataset_subset.obs['temp_binary_group']

cell_gene_counts.to_csv("cell_gene_counts_gaba_organoids.csv") 