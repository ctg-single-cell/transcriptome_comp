"""
# Scanorama Data Integration Prototype In vitro Datasets

#Created by: Fallon Ratner 26/06/23
"""
import pandas as pd
import numpy as np
import scanpy as sc
import scanorama
working_directory = '/home/rfallon/data/'
#Read in combined vitro and vivo scnapy object
test = sc.read(working_directory +'vitro_vivo.h5ad')

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
print("Variable genes in all batches: %d" % sum(test.var.highly_variable_nbatches == 21))
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


#save adata file
adata_sc.write_h5ad('vitro_vivo_int.h5ad')