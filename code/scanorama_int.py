"""
# Scanorama Data Integration Tutorial

#Created by: Fallon Ratner 24/05/23
"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as plt
import anndata
#Read in Herring files
data_path = "C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/data/Herring/"
ga22 = sc.read("ga22_cleaned_count_matrices.h5ad")
ga24 = sc.read("ga24_cleaned_count_matrices.h5ad")
ga34 = sc.read("ga34_cleaned_count_matrices.h5ad")

## Combine the AnnData objects into one
combo_adata = sc.AnnData.concatenate(ga22, ga24, ga34)

#Dont perform variable gene selection on scaled data matrix
adata2 = combo_adata.raw
adata2.raw = None
# check that the matrix looks like noramlized counts
print(adata2.X[1:10,1:10])
#Detect Variable Genes
var_genes_all = combo_adata.var.highly_variable

print("Highly variable genes: %d"%sum(var_genes_all))
#Variance Stabilizing Transformation - VST
sc.pp.log1p(combo_adata, copy=False)
sc.pp.highly_variable_genes(combo_adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'batch')

print("Highly variable genes intersection: %d"%sum(combo_adata.var.highly_variable_intersection))

print("Number of batches where gene is variable:")
print(combo_adata.var.highly_variable_nbatches.value_counts())

var_genes_batch = combo_adata.var.highly_variable_nbatches > 0

#Calculate the the variable genes (sum) for all batches
print("Any batch var genes: %d"%sum(var_genes_batch))
print("All data var genes: %d"%sum(var_genes_all))
#Which variable genes are the same between batches
print("Overlap: %d"%sum(var_genes_batch & var_genes_all))
print("Variable genes in all batches: %d"%sum(combo_adata.var.highly_variable_nbatches == 6))
print("Overlap batch instersection and all: %d"%sum(var_genes_all & combo_adata.var.highly_variable_intersection))


var_select = combo_adata.var.highly_variable_nbatches > 2
var_genes = var_select.index[var_select]
len(var_genes)

#Data Integration
#First we need to create individual AnnData objects from each of the datasets.
# split per batch into new objects.
batches = combo_adata.obs['batch'].cat.categories.tolist()
alldata = {}
for batch in batches:
    alldata[batch] = combo_adata[combo_adata.obs['batch'] == batch,]

alldata


import scanorama


#subset the individual dataset to the variable genes we defined at the beginning
alldata2 = dict()
for ds in alldata.keys():
    print(ds)
    alldata2[ds] = alldata[ds][:,var_genes]

#convert to list of AnnData objects
adatas = list(alldata2.values())

# run scanorama.integrate
#scanorama.integrate_scanpy(adatas, dimred = 50)

# Run Integration and batch correction.
corrected = scanorama.correct_scanpy(adatas, return_dimred=True)


#scanorama adds the corrected matrix to adata.obsm in each of the datasets in adatas.

corrected[0].obsm['X_scanorama'].shape

# Get all the integrated matrices.
#scanorama_int = [corrected.obsm['X_scanorama'] for corrected in adatas]
scanorama_int = [ad.obsm['X_scanorama'] for ad in corrected]
# make into one matrix.
all_s = np.concatenate(scanorama_int)
print(all_s.shape)

# add to the AnnData object, create a new object first
adata_sc = combo_adata.copy()
adata_sc.obsm["Scanorama"] = all_s

#umap
sc.pp.neighbors(adata_sc, n_pcs =30, use_rep = "Scanorama")
sc.tl.umap(adata_sc)

sc.pl.umap(adata_sc, color="age", title="Scanorama (integrated & batch corrected) umap")
sc.pl.umap(combo_adata, color="age", title="Uncorrected")

#save adata file
adata_sc.write_h5ad('C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/data/Herring/herr_int_sc.h5ad')

adata = sc.read("herr_int_sc.h5ad")

