"""
# Scanorama Data Integration Prototype In vivo Datasets

#Created by: Fallon Ratner 14/06/23
"""
import pandas as pd
import numpy as np
import scanpy as sc
import anndata
import scanorama
from PIL import Image
working_directory = '/home/rfallon/data/'
# Set output directory
#output_dir = '/home/rfallon/output'
#Herring Data: 22, 24, 34 GA
herr22 = sc.read(working_directory +"ga22_cleaned_count_matrices.h5ad")
herr24 = sc.read(working_directory +"ga24_cleaned_count_matrices.h5ad")
herr34 = sc.read(working_directory +"ga34_cleaned_count_matrices.h5ad")

# Assign the age & author as a new column in the metadata as Source
herr22.obs['Source'] = 'Herring_GA22'
herr24.obs['Source'] = 'Herring_GA24'
herr34.obs['Source'] = 'Herring_GA34'

#Polioudakis: 17-18 GW 
pol = sc.read_csv(working_directory +"Polio_matrix.csv")
pol = pol.transpose()

# Assign the age & author as a new column in the metadata as Source
pol.obs['Source'] = 'Polioudakis_GW17-18'

#Han Data: 11, 12, 13 GW, authors removed cells with less then 500 UMI counts
han11b1 = sc.read_text(working_directory +'FetalBrain4.rawdge.txt')

#han11b1 = sc.read_text('GSM4008679_Fetal-Brain4_dge.txt')
han11b1 = han11b1.transpose()

##
han11b2 = sc.read_text(working_directory +'FetalBrain6.rawdge.txt')
#han11b2 = sc.read_text('GSM4008681_Fetal-Brain6_dge.txt')
han11b2 = han11b2.transpose()

###
han12 = sc.read_text(working_directory +'FetalBrain5.rawdge.txt')
#han12 = sc.read_text('GSM4008680_Fetal-Brain5_dge.txt')
han12 = han12.transpose()

##
han13 = sc.read_text(working_directory +'FetalBrain3.rawdge.txt')
#han13 = sc.read_text('GSM4008678_Fetal-Brain3_dge.txt')
han13 = han13.transpose()


# Assign the age & author as a new column in the metadata as Source
han11b1.obs['Source'] = 'Han_GW11.1'
han11b2.obs['Source'] = 'Han_GW11.2'
han12.obs['Source'] = 'Han_GW12'
han13.obs['Source'] = 'Han_GW13'



#Couturier Data: 13, 17, 19 GW, authors removed cells with less then a 1000 cells
cou13 = sc.read_10x_mtx(working_directory +'HFA567_total.filtered_gene_matrices', var_names='gene_symbols', cache=True)
cou17 = sc.read_10x_mtx(working_directory +'HFA570_total.filtered_gene_matrices', var_names='gene_symbols', cache=True)
cou19 = sc.read_10x_mtx(working_directory +'HFA571_total.filtered_gene_matrices', var_names='gene_symbols', cache=True)

# Assign the age & author as a new column in the metadata as Source
cou13.obs['Source'] = 'Couturier_GW13'
cou17.obs['Source'] = 'Couturier_GW17'
cou19.obs['Source'] = 'Couturier_GW19'

#Fan Data: 22-23 GW
liu = sc.read(working_directory +'hNSPC_raw_counts.h5ad')


# Assign the age & author as a new column in the metadata as Source
liu.obs['Source'] = 'Liu_GW17-19'



#Zhong Data: 8-26 GW
#zho = sc.read_text('Zhong_matrix.txt', delimiter='\t')
#zho = zho.transpose()

#Assign the age & author as a new column in the metadata as Source
#import re

# Extract the string before the underscore from the cell names
#age_info = []  # List to store age information

#for cell_name in zho.obs_names:
#    match = re.search(r'(.+?)_', cell_name)
#    if match:
#        age_info.append(match.group(1) + "_Zhong")
#    else:
#        age_info.append('Unknown')

# Create a new column in the metadata and assign the extracted age information to it
#zho.obs['Source'] = age_info

#need to remove some values from var
# Define the regular expression pattern
#pattern = r'\d{1,2}-[A-Za-z]{3}'

# Create a mask to identify variable names that match the pattern
#mask = np.array([re.match(pattern, name) is not None for name in zho.var_names])

# Remove the variable names that match the pattern from the Scanpy object
#zho = zho[:, ~mask]



dataset_list = [herr22, herr24, herr34, pol, cou13, cou17, cou19, han11b1, han11b2, han12, han13, liu]  


for adata in dataset_list:
    # Preprocessing steps for each adata object
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata, copy=False)
      
####################################################################################

# Convert all data to numpy.float32
for adata in dataset_list:
    adata.X = adata.X.astype(np.float32)
    
# for each element  make var names unique
for adata in dataset_list:
    adata.var_names_make_unique()

# Make observation names unique in each `adata` object
for adata in dataset_list:
    adata.obs_names_make_unique()
    
test = sc.concat(dataset_list, axis = 0, join = 'outer')

# Perform preprocessing steps
sc.pp.filter_cells(test, min_genes=200)
sc.pp.filter_genes(test, min_cells=3)
test.var['mt'] = test.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(test, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
test = test[test.obs.n_genes_by_counts < 2500, :]
test = test[test.obs.pct_counts_mt < 5, :]
sc.pp.normalize_total(test, target_sum=1e4)
sc.pp.log1p(test, copy=False)

# Identify highly variable genes using batch information
sc.pp.highly_variable_genes(test, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='Source')

# Detect Variable Genes
var_genes_all = test.var.highly_variable
print("Highly variable genes: %d" % sum(var_genes_all))

print("Highly variable genes intersection: %d" % sum(test.var.highly_variable_intersection))
print("Number of batches where gene is variable:")
print(test.var.highly_variable_nbatches.value_counts())

var_genes_batch = test.var.highly_variable_nbatches > 0

# Calculate the variable genes (sum) for all batches
print("Any batch var genes: %d" % sum(var_genes_batch))
print("All data var genes: %d" % sum(var_genes_all))

# Identify variable genes that are the same between batches 
print("Overlap: %d" % sum(var_genes_batch & var_genes_all))
#This should be the number of batches
print("Variable genes in all batches: %d" % sum(test.var.highly_variable_nbatches == 12))
print("Overlap batch intersection and all: %d" % sum(var_genes_all & test.var.highly_variable_intersection))

var_select = test.var.highly_variable_nbatches > 2
var_genes = var_select.index[var_select]
len(var_genes)

#Data Integration
#First we need to create individual AnnData objects from each of the datasets.
# split per batch into new objects.
batches = test.obs['Source'].cat.categories.tolist()
alldata = {}
for batch in batches:
    alldata[batch] = test[test.obs['Source'] == batch,]

alldata
#subset the individual dataset to the variable genes we defined at the beginning
alldata2 = dict()
for ds in alldata.keys():
    print(ds)
    alldata2[ds] = alldata[ds][:,var_genes]

#convert to list of AnnData objects
adatas = list(alldata2.values())

# Run Integration and batch correction.
corrected = scanorama.correct_scanpy(adatas, return_dimred=True)

#scanorama adds the corrected matrix to adata.obsm in each of the datasets in adatas.
corrected[0].obsm['X_scanorama'].shape

# Get all the integrated matrices.
scanorama_int = [ad.obsm['X_scanorama'] for ad in corrected]
# make into one matrix.
all_s = np.concatenate(scanorama_int)
print(all_s.shape)

# add to the AnnData object, create a new object first
adata_sc = test.copy()
adata_sc.obsm["Scanorama"] = all_s

#perform umap
sc.pp.neighbors(adata_sc, n_pcs =30, use_rep = "Scanorama")
sc.tl.umap(adata_sc)
#map umap
sc.pl.umap(adata_sc, color="Source", title="Scanorama (integrated & batch corrected) umap")
sc.pl.umap(test, color="Source", title="Uncorrected")
#save adata file
adata_sc.write_h5ad('vivo2_int.h5ad')



