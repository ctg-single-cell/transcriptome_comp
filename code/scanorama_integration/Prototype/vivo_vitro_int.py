"""
# Scanorama Data Integration Prototype In vitro Datasets

#Created by: Fallon Ratner 25/06/23
"""
import pandas as pd
import numpy as np
import scanpy as sc
import scanorama
#Herring Data: 22, 24, 34 GA
herr22 = sc.read("ga22_cleaned_count_matrices.h5ad")
herr24 = sc.read("ga24_cleaned_count_matrices.h5ad")
herr34 = sc.read("ga34_cleaned_count_matrices.h5ad")

# Assign the age & author as a new column in the metadata as Source
herr22.obs['Source'] = 'Herring_GA22'
herr24.obs['Source'] = 'Herring_GA24'
herr34.obs['Source'] = 'Herring_GA34'
#################################################################################
#Polioudakis: 17-18 GW 
pol = sc.read_csv("Polio_matrix.csv")
pol = pol.transpose()

# Assign the age & author as a new column in the metadata as Source
pol.obs['Source'] = 'Polioudakis_GW17-18'
#################################################################################

#Han Data: 11, 12, 13 GW, authors removed cells with less then 500 UMI counts
han11b1 = sc.read_text('FetalBrain4.rawdge.txt')

han11b1 = han11b1.transpose()

##
han11b2 = sc.read_text('FetalBrain6.rawdge.txt')

han11b2 = han11b2.transpose()

###
han12 = sc.read_text('FetalBrain5.rawdge.txt')

han12 = han12.transpose()

##
han13 = sc.read_text('FetalBrain3.rawdge.txt')

han13 = han13.transpose()

# Assign the age & author as a new column in the metadata as Source
han11b1.obs['Source'] = 'Han_GW11.1'
han11b2.obs['Source'] = 'Han_GW11.2'
han12.obs['Source'] = 'Han_GW12'
han13.obs['Source'] = 'Han_GW13'
#################################################################################
#Couturier Data: 13, 17, 19 GW, authors removed cells with less then a 1000 cells
cou13 = sc.read_10x_mtx('HFA567_total.filtered_gene_matrices', var_names='gene_symbols', cache=True)
cou17 = sc.read_10x_mtx('HFA570_total.filtered_gene_matrices', var_names='gene_symbols', cache=True)
cou19 = sc.read_10x_mtx('HFA571_total.filtered_gene_matrices', var_names='gene_symbols', cache=True)

# Assign the age & author as a new column in the metadata as Source
cou13.obs['Source'] = 'Couturier_GW13'
cou17.obs['Source'] = 'Couturier_GW17'
cou19.obs['Source'] = 'Couturier_GW19'
#################################################################################
#Fan Data: 22-23 GW
liu = sc.read('hNSPC_raw_counts.h5ad')

# Assign the age & author as a new column in the metadata as Source
liu.obs['Source'] = 'Liu_GW17-19'
#################################################################################
#Velasco 3 and 6 months
vel3 = sc.read_csv('expression_PGP1.3mon.txt', delimiter='\t')
vel3 = vel3.transpose()

vel6 = sc.read_csv('expression_PGP1.6mon.txt', delimiter='\t')
vel6 = vel6.transpose()

# Assign the age & author as a new column in the metadata as Source
vel3.obs['Source'] = 'Velasco_3mo'
vel6.obs['Source'] = 'Velasco_6mo'
#################################################################################
#Trujillo 1, 3, 6, 10 months
truj1 = sc.read_10x_mtx('Trujillo/1_month', var_names='gene_symbols', cache=True)
truj3 = sc.read_10x_mtx('Trujillo/3_months', var_names='gene_symbols', cache=True)
truj6 = sc.read_10x_mtx('Trujillo/6_months', var_names='gene_symbols', cache=True)
truj10 = sc.read_10x_mtx('Trujillo/10_months', var_names='gene_symbols', cache=True)

# Assign the age & author as a new column in the metadata as Source
truj1.obs['Source'] = 'Trujillo_1mo'
truj3.obs['Source'] = 'Trujillo_3mo'
truj6.obs['Source'] = 'Trujillo_6mo'
truj10.obs['Source'] = 'Trujillo_10mo'
#################################################################################
#Giandomenico 75 days  #log transformed, maybe raw data in RNA assay - need to check
gian = sc.read_csv('GSE124174_org75_seuratdata.txt', delimiter=' ')
gian = gian.transpose()

# Assign the age & author as a new column in the metadata as Source
gian.obs['Source'] = 'Giandomenico_2.5mo'
#################################################################################
#Fair: 90 & 140 days (3mo & 4.5mo)
fair90 = sc.read_10x_mtx('Fair/90_days', var_names='gene_symbols', cache=True)
fair140 = sc.read_10x_mtx('140_days', var_names='gene_symbols', cache=True)

# Assign the age & author as a new column in the metadata as Source
fair90.obs['Source'] = 'Fair_3mo'
fair140.obs['Source'] = 'Fair_4.5mo'
#################################################################################

#Meng 50 days
mengu1 = sc.read_10x_mtx('Meng/U1M', var_names='gene_symbols', cache=True)
mengu2 = sc.read_10x_mtx('Meng/U2F', var_names='gene_symbols', cache=True)

# Assign the age & author as a new column in the metadata as Source
mengu1.obs['Source'] = 'Meng_U1M_1.5mo'
mengu2.obs['Source'] = 'Meng_U2F_1.5mo'
#################################################################################
dataset_list = [herr22, herr24, herr34, pol, cou13, cou17, cou19, han11b1, han11b2, han12, han13, liu, vel3, vel6, truj1, truj3, truj6, truj10, gian, mengu1, mengu2]  

for adata in dataset_list:
    # Print the number of cells and genes before filtering
    print(f"Before filtering: {adata.shape[0]} cells and {adata.shape[1]} genes")
    
    # Preprocessing steps for each adata object
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Print the number of cells and genes after filtering
    print(f"After cell and gene filtering: {adata.shape[0]} cells and {adata.shape[1]} genes")

    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    
    # Print the number of cells and genes after quality control filtering
    print(f"After quality control filtering: {adata.shape[0]} cells and {adata.shape[1]} genes")

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

# Print the number of cells and genes before filtering
print(f"Before filtering: {test.shape[0]} cells and {test.shape[1]} genes")
# Perform preprocessing steps
sc.pp.filter_cells(test, min_genes=200)
sc.pp.filter_genes(test, min_cells=3)
# Print the number of cells and genes after filtering
print(f"After cell and gene filtering: {test.shape[0]} cells and {test.shape[1]} genes")
test.var['mt'] = test.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(test, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
test = test[test.obs.n_genes_by_counts < 2500, :]
test = test[test.obs.pct_counts_mt < 5, :]
# Print the number of cells and genes after quality control filtering
print(f"After quality control filtering: {test.shape[0]} cells and {test.shape[1]} genes")
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
adata_sc.write_h5ad('C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/vitro_vivo_int.h5ad')