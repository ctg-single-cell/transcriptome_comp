# Cell Typist Annotation
# Created by Fallon Ratner
# Date: 11/05/23

# First install cell typist
pip install celltypist

import numpy as np
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models

#Show all available models that can be downloaded and used.
models.models_description()
#Download a specific model, for example, `Developing_Human_Brain.pkl`.
models.download_models(model = 'Developing_Human_Brain.pkl')
#Update the models by re-downloading the latest versions if you think they may be outdated.
models.download_models(model = ['Developing_Human_Brain.pkl'], force_update = True)
#Show the local directory storing these models.
models.models_path #'C:\\Users\\fallo\\.celltypist\\data\\models'
#Get an overview of the models that are downloaded in `1.2.`.
#By default (`on_the_fly = False`), all possible models (even those that are not downloaded) are shown.
models.models_description(on_the_fly = True)
#Select the model from the above list. If the `model` argument is not provided, will default to `Developing_Human_Brain.pkl`.
model = models.Model.load(model = 'Developing_Human_Brain.pkl')
#The model summary information.
model
#Examine cell types contained in the model.
model.cell_types
#Examine genes/features contained in the model.
model.features

###########################################################################################################
#First Dataset Polioudakis 
input_file = sc.read_csv("Polio_matrix.csv")
adata = input_file.transpose()
#need normalized expression for prediction
#pre-processing with Scanpy
#Filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
#Quality Control
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
#Filtering and slicing
adata = adata[adata.obs.pct_counts_mt < 5, :]
#Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
#Logarithmize
sc.pp.log1p(adata)
adata.raw = adata
#Predict the identity of each input cell.
predictions = celltypist.annotate(adata, model = 'Developing_Human_Brain.pkl')

#In case your input file is a gene-by-cell mtx file.
predictions = celltypist.annotate(adata, model = 'Developing_Human_Brain.pkl', transpose_input = True)

#Summary information for the prediction result.
predictions
#Examine the prediction results.
predictions.predicted_labels
#Examine the matrix representing the decision score of each cell belonging to a given cell type.
predictions.decision_matrix
#Examine the matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
predictions.probability_matrix

#Query cell will get the label of 'Unassigned' if it fails to pass the probability cutoff in each cell type.
#Query cell will get multiple label outputs (concatenated by '|') if more than one cell type passes the probability cutoff.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl', mode = 'prob match', p_thres = 0.5)

#Export the three results to csv tables.
predictions.to_table(folder = 'C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/output', prefix = 'Polio_')

annot = pd.read_csv("Polio_predicted_labels.csv")

annot['Group'] = '17-18 GW'

#Calculate proportions
def calculate_cell_type_percentages(df, group):
    # Count the number of each cell type
    count = df['predicted_labels'].value_counts().reset_index()
    count.columns = ['predicted_labels', 'Count']

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
counts_Polio = calculate_cell_type_percentages(annot, '17-18 GW')
######################################################################################
input_file =sc.read('C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/data/Herring/ga22_cleaned_count_matrices.h5ad')

#need normalized expression for prediction
#pre-processing with Scanpy
#Filtering
sc.pp.filter_cells(input_file, min_genes=200)
sc.pp.filter_genes(input_file, min_cells=3)
#Quality Control
input_file.var['mt'] = input_file.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(input_file, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
input_file = input_file[input_file.obs.n_genes_by_counts < 2500, :]
#Filtering and slicing
input_file = input_file[input_file.obs.pct_counts_mt < 5, :]
#Normalize
sc.pp.normalize_total(input_file, target_sum=1e4)
#Logarithmize
sc.pp.log1p(input_file)
input_file.raw = input_file
#Predict the identity of each input cell.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl')

#In case your input file is a gene-by-cell table.
#predictions = celltypist.annotate(adata, model = 'Developing_Human_Brain.pkl', transpose_input = True)
#In case your input file is a gene-by-cell mtx file.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl', transpose_input = True)

#Summary information for the prediction result.
predictions
#Examine the prediction results.
predictions.predicted_labels
#Examine the matrix representing the decision score of each cell belonging to a given cell type.
predictions.decision_matrix
#Examine the matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
predictions.probability_matrix

#Query cell will get the label of 'Unassigned' if it fails to pass the probability cutoff in each cell type.
#Query cell will get multiple label outputs (concatenated by '|') if more than one cell type passes the probability cutoff.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl', mode = 'prob match', p_thres = 0.5)

#Export the three results to csv tables.
predictions.to_table(folder = 'output', prefix = 'GA22_')

annot = pd.read_csv("GA22_predicted_labels.csv")

annot['Group'] = '22 GA'

#Calculate proportions

# Call the function for each dataset
counts_ga22 = calculate_cell_type_percentages(annot, '22 GA')
######################################################################################
input_file =sc.read('ga24_cleaned_count_matrices.h5ad')

#need normalized expression for prediction
#pre-processing with Scanpy
#Filtering
sc.pp.filter_cells(input_file, min_genes=200)
sc.pp.filter_genes(input_file, min_cells=3)
#Quality Control
input_file.var['mt'] = input_file.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(input_file, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
input_file = input_file[input_file.obs.n_genes_by_counts < 2500, :]
#Filtering and slicing
input_file = input_file[input_file.obs.pct_counts_mt < 5, :]
#Normalize
sc.pp.normalize_total(input_file, target_sum=1e4)
#Logarithmize
sc.pp.log1p(input_file)
input_file.raw = input_file
#Predict the identity of each input cell.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl')

#In case your input file is a gene-by-cell table.
#predictions = celltypist.annotate(adata, model = 'Developing_Human_Brain.pkl', transpose_input = True)
#In case your input file is a gene-by-cell mtx file.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl', transpose_input = True)

#Summary information for the prediction result.
predictions
#Examine the prediction results.
predictions.predicted_labels
#Examine the matrix representing the decision score of each cell belonging to a given cell type.
predictions.decision_matrix
#Examine the matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
predictions.probability_matrix

#Query cell will get the label of 'Unassigned' if it fails to pass the probability cutoff in each cell type.
#Query cell will get multiple label outputs (concatenated by '|') if more than one cell type passes the probability cutoff.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl', mode = 'prob match', p_thres = 0.5)

#Export the three results to csv tables.
predictions.to_table(folder = 'output', prefix = 'GA24_')

annot = pd.read_csv("GA24_predicted_labels.csv")

annot['Group'] = '24 GA'

#Calculate proportions
# Call the function for each dataset
counts_ga24 = calculate_cell_type_percentages(annot, '24_GA')    
######################################################################################
input_file =sc.read('ga34_cleaned_count_matrices.h5ad')

#need normalized expression for prediction
#pre-processing with Scanpy
#Filtering
sc.pp.filter_cells(input_file, min_genes=200)
sc.pp.filter_genes(input_file, min_cells=3)
#Quality Control
input_file.var['mt'] = input_file.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(input_file, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
input_file = input_file[input_file.obs.n_genes_by_counts < 2500, :]
#Filtering and slicing
input_file = input_file[input_file.obs.pct_counts_mt < 5, :]
#Normalize
sc.pp.normalize_total(input_file, target_sum=1e4)
#Logarithmize
sc.pp.log1p(input_file)
input_file.raw = input_file
#Predict the identity of each input cell.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl')

#In case your input file is a gene-by-cell mtx file.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl', transpose_input = True)

#Summary information for the prediction result.
predictions
#Examine the prediction results.
predictions.predicted_labels
#Examine the matrix representing the decision score of each cell belonging to a given cell type.
predictions.decision_matrix
#Examine the matrix representing the probability each cell belongs to a given cell type (transformed from decision matrix by the sigmoid function).
predictions.probability_matrix

#Query cell will get the label of 'Unassigned' if it fails to pass the probability cutoff in each cell type.
#Query cell will get multiple label outputs (concatenated by '|') if more than one cell type passes the probability cutoff.
predictions = celltypist.annotate(input_file, model = 'Developing_Human_Brain.pkl', mode = 'prob match', p_thres = 0.5)

#Export the three results to csv tables.
predictions.to_table(folder = 'output', prefix = 'GA34_')

annot = pd.read_csv("GA34_predicted_labels.csv")

annot['Group'] = '34 GA'

#Calculate proportions
# Call the function for each dataset
counts_ga34 = calculate_cell_type_percentages(annot, '34_GA')    

####################################################################################\
input_file = sc.read_csv('expression_PGP1.3mon.txt', delimiter='\t')
# transpose the expression data
adata = input_file.transpose()

#need normalized expression for prediction
#pre-processing with Scanpy
#Filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
#Quality Control
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
#Filtering and slicing
adata = adata[adata.obs.pct_counts_mt < 5, :]
#Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
#Logarithmize
sc.pp.log1p(adata)
adata.raw = adata
#Predict the identity of each input cell.
predictions = celltypist.annotate(adata, model = 'Developing_Human_Brain.pkl')

#In case your input file is a gene-by-cell mtx file.
predictions = celltypist.annotate(adata, model = 'Developing_Human_Brain.pkl', transpose_input = True)
#Export the three results to csv tables.
predictions.to_table(folder = 'output', prefix = 'Velasco_3mon_')

annot = pd.read_csv("Velasco_3mon_predicted_labels.csv")

annot['Group'] = '3 months'
# Calculate the proportions
counts_3mon = calculate_cell_type_percentages(annot, '3 months')
    
####################################################################################
input_file = sc.read_csv('expression_PGP1.6mon.txt', delimiter='\t')
# transpose the expression data
adata = input_file.transpose()

#need normalized expression for prediction
#pre-processing with Scanpy
#Filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
#Quality Control
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
#Filtering and slicing
adata = adata[adata.obs.pct_counts_mt < 5, :]
#Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
#Logarithmize
sc.pp.log1p(adata)
adata.raw = adata
#Predict the identity of each input cell.
predictions = celltypist.annotate(adata, model = 'Developing_Human_Brain.pkl')

#In case your input file is a gene-by-cell mtx file.
predictions = celltypist.annotate(adata, model = 'Developing_Human_Brain.pkl', transpose_input = True)
#Export the three results to csv tables.
predictions.to_table(folder = 'output', prefix = 'Velasco_6mon_')

annot = pd.read_csv("Velasco_6mon_predicted_labels.csv")

annot['Group'] = '6 months'

# Calculate the proportions
counts_6mon = calculate_cell_type_percentages(annot, '6 months')
