#Created by Fallon Ratner on 20-09-23
#Gene Signature Correlation
#gene signature using GSEA for each devo stage and each organoid 

#Use marker genes to establish a gene signature
#I want to generate an average gene signature for each in vivo devo stage:
#Early Development (11-13GW) --> han11.1, han11.2, han12, han13, cou13
#Mid Development (17-19GW) --> cou17, cou19, polio, liu
#Late Development (22GW-Day2) --> herr22, herr24, herr34, herrday2
#and for each cell type
#I only include cell types if there are more than 10 
#I conduct a differential expression analysis of the selected cell types vs all other cell types in the dataset using the wilcoxon test
import numpy as np
import pandas as pd
import scanpy as sc



def process_datasets(cell_types, dataset_names, phase):
    """
    Process the datasets for the given cell types and phase.
    """
    cell_gene_counts = pd.DataFrame(columns=['Cell_Type', 'Num_Cells', 'Num_Genes'])
    for target_cell_type in cell_types:
        dataset_results = []
        for dataset_name in dataset_names:
            adata = sc.read(f"{dataset_name}_annot.h5ad")
             
             #If the current dataset is 'polio', rename the 'customclassif' column to 'annotation'
            if dataset_name == "polio" and "customclassif" in adata.obs.columns:
                adata.obs.rename(columns={"customclassif": "annotation"}, inplace=True)
            
            adata.X[np.isnan(adata.X)] = 0  # Set any NaN values to 0
            # Preprocessing steps
            sc.pp.filter_cells(adata, min_genes=200)
            sc.pp.filter_genes(adata, min_cells=3)
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            sc.pp.normalize_total(adata, target_sum=1e4)
            
            if target_cell_type in adata.obs['annotation'].unique():
                # Create a temporary binary category
                adata.obs['temp_binary_group'] = np.where(adata.obs['annotation'] == target_cell_type, target_cell_type, 'other')

                target_cells_count = np.sum(adata.obs['temp_binary_group'] == target_cell_type)
                other_cells_count = np.sum(adata.obs['temp_binary_group'] == 'other')

                if target_cells_count > 10 and other_cells_count > 10:
                    # Differential expression analysis
                    sc.tl.rank_genes_groups(adata, groupby='temp_binary_group', groups=[target_cell_type], reference='other', method='wilcoxon')

                    # Collecting results
                    results_df = pd.DataFrame({
                        'genes': adata.uns['rank_genes_groups']['names'][target_cell_type],
                        'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][target_cell_type],
                        'pvals': adata.uns['rank_genes_groups']['pvals'][target_cell_type],
                        'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][target_cell_type]
                    })
                    dataset_results.append(results_df)

                    # Capture number of cells and genes for this cell type
                    selected_subset = adata[adata.obs['annotation'] == target_cell_type]
                    num_cells = selected_subset.shape[0]
                    num_genes = selected_subset.shape[1]
                    new_row = pd.DataFrame({'Cell_Type': [target_cell_type], 'Num_Cells': [num_cells], 'Num_Genes': [num_genes]})
                    cell_gene_counts = pd.concat([cell_gene_counts, new_row], ignore_index=True)

                # Cleanup
                del adata.obs['temp_binary_group']

        # Aggregate results
        if dataset_results:
            combined = pd.concat(dataset_results)
            grouped = combined.groupby('genes').agg({
                'logfoldchanges': 'mean',
                'pvals': 'mean',
                'pvals_adj': 'mean'
            }).reset_index()

            # Save results
            output_filename = f"gene_signature_{target_cell_type}_{phase}.csv"
            grouped.to_csv(output_filename)

    # Save cell gene counts
    cell_gene_counts.to_csv(f"cell_gene_counts_{phase}.csv")

# Cell types
cell_types = ["outer radial glia", "ventral Radial Glia", "Intermediate Progenitors",
               "Neuroblasts", "Immature Excitatory", "Maturing Excitatory", "Migrating Excitatory",
               "Deep Excitatory Layers", "Upper Excitatory Layers", "Mixed",
               "Interneuron Precursors", "MGE INs", "CGE INs",
               "Glioblasts", "Immature Astrocytes", "Mature Astrocytes", "OPCs",
               "Oligodendrocytes", "Microglia", "Endothelial Cells", "Mural"]

# Dataset names
early_datasets = ["han11.1", "han11.2", "han12", "han13", "cou13"]
mid_datasets = ["cou17", "cou19", "liu", "polio"]
late_datasets = ["ga22", "ga24", "ga34", "day2"]

# Call the function for each phase
process_datasets(cell_types, early_datasets, 'early_devo')
process_datasets(cell_types, mid_datasets, 'mid_devo')
process_datasets(cell_types, late_datasets, 'late_devo')
            
###############################################################################
#For the Gaba Subtypes: only doing the cells found in both fetal brain and organoids
def process_gaba_datasets(cell_types, dataset_names, phase):
    """
    Process the datasets for the given cell types and phase.
    """
    cell_gene_counts = pd.DataFrame(columns=['Cell_Type', 'Num_Cells', 'Num_Genes'])
    for target_cell_type in cell_types:
        dataset_results = []
        for dataset_name in dataset_names:
            adata = sc.read(f"{dataset_name}_gaba_annot.h5ad")
            
            #If the current dataset is 'polio', rename the 'customclassif' column to 'annotation'
            if dataset_name == "polio" and "customclassif" in adata.obs.columns:
                adata.obs.rename(columns={"customclassif": "annotation"}, inplace=True)

            adata.X[np.isnan(adata.X)] = 0  # Set any NaN values to 0
            # Preprocessing steps
            sc.pp.filter_cells(adata, min_genes=200)
            sc.pp.filter_genes(adata, min_cells=3)
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            sc.pp.normalize_total(adata, target_sum=1e4)
            
            if target_cell_type in adata.obs['annotation'].unique():
                # Create a temporary binary category
                adata.obs['temp_binary_group'] = np.where(adata.obs['annotation'] == target_cell_type, target_cell_type, 'other')

                target_cells_count = np.sum(adata.obs['temp_binary_group'] == target_cell_type)
                other_cells_count = np.sum(adata.obs['temp_binary_group'] == 'other')

                if target_cells_count > 10 and other_cells_count > 10:
                    # Differential expression analysis
                    sc.tl.rank_genes_groups(adata, groupby='temp_binary_group', groups=[target_cell_type], reference='other', method='wilcoxon')

                    # Collecting results
                    results_df = pd.DataFrame({
                        'genes': adata.uns['rank_genes_groups']['names'][target_cell_type],
                        'logfoldchanges': adata.uns['rank_genes_groups']['logfoldchanges'][target_cell_type],
                        'pvals': adata.uns['rank_genes_groups']['pvals'][target_cell_type],
                        'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj'][target_cell_type]
                    })
                    dataset_results.append(results_df)

                    # Capture number of cells and genes for this cell type
                    selected_subset = adata[adata.obs['annotation'] == target_cell_type]
                    num_cells = selected_subset.shape[0]
                    num_genes = selected_subset.shape[1]
                    new_row = pd.DataFrame({'Cell_Type': [target_cell_type], 'Num_Cells': [num_cells], 'Num_Genes': [num_genes]})
                    cell_gene_counts = pd.concat([cell_gene_counts, new_row], ignore_index=True)

                # Cleanup
                del adata.obs['temp_binary_group']

        # Aggregate results
        if dataset_results:
            combined = pd.concat(dataset_results)
            grouped = combined.groupby('genes').agg({
                'logfoldchanges': 'mean',
                'pvals': 'mean',
                'pvals_adj': 'mean'
            }).reset_index()

            # Save results
            output_filename = f"gene_signature_{target_cell_type}_gaba_{phase}.csv"
            grouped.to_csv(output_filename)

    # Save cell gene counts
    cell_gene_counts.to_csv(f"cell_gene_counts_gaba_{phase}.csv")

# Cell types
cell_types = ["NDNF", "SST"]

# Dataset names
early_datasets = ["han11.2", "han12", "han13","cou13"]
mid_datasets = ["cou17", "cou19","polio"]
late_datasets = ["ga22", "ga24", "ga34", "day2"]


# Call the function for each phase
process_gaba_datasets(cell_types, early_datasets, 'early_devo')
process_gaba_datasets(cell_types, mid_datasets, 'mid_devo')
process_gaba_datasets(cell_types, late_datasets, 'late_devo')


    