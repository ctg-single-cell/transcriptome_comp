#Created by Fallon Ratner 
#Conduct pairwise comparison to correlate cell types in fetal brain vs organoids based on the logfc

import pandas as pd
import os
import itertools
from scipy.stats import spearmanr

datasets = ["Bhaduri_H1S_10w", "Bhaduri_H1S_10w", "Bhaduri_H1S_8w", "Bhaduri_H1X_8w",
            "Bhaduri_H1S_5w","Bhaduri_H1X_5w","Giandomenico_10w","Fair_20w","Madhavan_12w",
            "Popova_batch2_7w", "Trujillo_12w","Trujillo_24w","Trujillo_40w","Velasco_12w",
            "Velasco_24w","Xiang_10w","Xiang_11w","Xiangb1_4w","Xiangb2_4w"] 

cell_types = ["outer radial glia", "ventral Radial Glia", "Intermediate Progenitors",
               "Neuroblasts","Immature Excitatory", "Maturing Excitatory", "Upper Excitatory Layers", 
               "Interneuron Precursors", "MGE INs", "CGE INs","NDNF", "SST",
               "Immature Astrocytes",  "OPCs", "Endothelial Cells", "Mural"] 

developmental_stages = ["early_devo", "mid_devo", "late_devo"]

correlations = []

# Loop over developmental stages
for dev_stage in developmental_stages:

    # Loop over cell types
    for cell_type in cell_types:

       # Try to read the CSV
        try:
            dev_stage_df = pd.read_csv(f"gene_signature_{cell_type}_{dev_stage}.csv")

            # If the file is empty, skip this iteration
            if dev_stage_df.empty:
                raise ValueError(f"File for {cell_type} in {dev_stage} is empty.")

        except (FileNotFoundError, ValueError) as e:
            print(e)
            continue

        # Initialize list to store correlations for this cell type
        correlations_for_celltype = []


        
        # Compare the developmental stage with each dataset
        for dataset in datasets:
           
            # Try to read the CSV
            try:
                dataset_df = pd.read_csv(f"gene_signature_{dataset}_{cell_type}.csv")
                
                # If the file is empty, skip this iteration
                if dataset_df.empty:
                    raise ValueError(f"File for {cell_type} in {dataset} is empty.")

            except (FileNotFoundError, ValueError) as e:
                print(e)
                continue
         
            common_genes = set(dev_stage_df.genes).intersection(dataset_df.genes)
            
            subset_dev = dev_stage_df.set_index("genes").loc[common_genes].reset_index()
            subset_dataset = dataset_df.set_index("genes").loc[common_genes].reset_index()

            # Calculate correlation for logFC values
            correlation, _ = spearmanr(subset_dev['logfoldchanges'], subset_dataset['logfoldchanges'])

            correlations_for_celltype.append({
                "DevelopmentalStage": dev_stage,
                "Dataset": dataset,
                "Correlation": correlation
            })

        # Convert to DataFrame and save for the current cell type
        correlation_df = pd.DataFrame(correlations_for_celltype)
        output_filename = f"pairwise_correlations_{dev_stage}_vs_datasets_{cell_type}.csv"
        correlation_df.to_csv(output_filename, index=False)
