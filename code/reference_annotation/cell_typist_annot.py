# Cell Typist Annotation
# Created by Fallon Ratner


# Install celltypist
# pip install celltypist

import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models

# Show all available models that can be downloaded and used.
models.models_description()

# Download a specific model, for example, `Developing_Human_Brain.pkl`.
models.download_models(model='Developing_Human_Brain.pkl')

# Update the models by re-downloading the latest versions if they may be outdated.
models.download_models(model=['Developing_Human_Brain.pkl'], force_update=True)

# Show the local directory storing these models.
models.models_path  # 'C:\\Users\\fallo\\.celltypist\\data\\models'

# Get an overview of the models that are downloaded.
models.models_description(on_the_fly=True)

# Select the model from the list (default to 'Developing_Human_Brain.pkl' if not provided).
model = models.Model.load(model='Developing_Human_Brain.pkl')

# Display the model summary information.
model

# Examine cell types contained in the model.
model.cell_types

# Examine genes/features contained in the model.
model.features

# Function to calculate cell type percentages and save the results to a CSV file.
def calculate_cell_type_percentages(df, group):
    count = df['predicted_labels'].value_counts().reset_index()
    count.columns = ['predicted_labels', 'Count']
    total_count = count['Count'].sum()
    proportions = count['Count'] / total_count
    count['Proportion'] = proportions
    count['Percentage'] = proportions * 100
    count['Group'] = group
    count.to_csv(f'proportion_{group}.csv', index=False)
    return count

# Process datasets
def process_dataset(file_path, group_prefix):
    adata = sc.read(file_path).transpose()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    predictions = celltypist.annotate(adata, model='Developing_Human_Brain.pkl', transpose_input=True)

    predictions.to_table(folder='output', prefix=f'{group_prefix}_')

    annot = pd.read_csv(f"{group_prefix}_predicted_labels.csv")
    annot['Group'] = group_prefix

    # Call the function for each dataset
    calculate_cell_type_percentages(annot, group_prefix)

# Process datasets
process_dataset("Polio_matrix.csv", "17-18 GW")
process_dataset("C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/data/Herring/ga22_cleaned_count_matrices.h5ad", "22 GA")
process_dataset("ga24_cleaned_count_matrices.h5ad", "24 GA")
process_dataset("expression_PGP1.3mon.txt", "3 months")
process_dataset("expression_PGP1.6mon.txt", "6 months")

