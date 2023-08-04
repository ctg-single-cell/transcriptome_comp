"""
# Scanorama Data Integration Prototype In vitro Datasets

#Created by: Fallon Ratner 14/06/23
"""
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.io import mmread
import scanorama

#Velasco 3 and 6 months
vel3 = sc.read_csv('expression_PGP1.3mon.txt', delimiter='\t')
vel3 = vel3.transpose()

vel6 = sc.read_csv('expression_PGP1.6mon.txt', delimiter='\t')
vel6 = vel6.transpose()

# Assign the age & author as a new column in the metadata as Source
vel3.obs['Source'] = 'Velasco_12w'
vel6.obs['Source'] = 'Velasco_24w'
###########################################################################################
#Trujillo 1, 3, 6, 10 months
truj1 = sc.read_10x_mtx('Trujillo/1_month', var_names='gene_symbols', cache=True)
truj3 = sc.read_10x_mtx('Trujillo/3_months', var_names='gene_symbols', cache=True)
truj6 = sc.read_10x_mtx('Trujillo/6_months', var_names='gene_symbols', cache=True)
truj10 = sc.read_10x_mtx('Trujillo/10_months', var_names='gene_symbols', cache=True)

# Assign the age & author as a new column in the metadata as Source
truj1.obs['Source'] = 'Trujillo_4w'
truj3.obs['Source'] = 'Trujillo_12w'
truj6.obs['Source'] = 'Trujillo_24w'
truj10.obs['Source'] = 'Trujillo_40w'
###########################################################################################
#Giandomenico 75 days  #log transformed, maybe raw data in RNA assay - need to check
gian = sc.read_csv('GSE124174_org75_seuratdata.txt', delimiter=' ')
gian = gian.transpose()

# Assign the age & author as a new column in the metadata as Source (10.7weeks rounded down)
gian.obs['Source'] = 'Giandomenico_10w'
###########################################################################################
#Fair: 90 & 140 days (3mo & 4.5mo)
fair90 = sc.read_10x_mtx('Fair/90_days', var_names='gene_symbols', cache=True)
fair140 = sc.read_10x_mtx('Fair/140_days', var_names='gene_symbols', cache=True)

# Assign the age & author as a new column in the metadata as Source
fair90.obs['Source'] = 'Fair_12w'
fair140.obs['Source'] = 'Fair_20w'
###########################################################################################
#Xiang 30, 72, 79 days
xiang = sc.read_10x_mtx('Xiang_data', var_names='gene_symbols', cache=True)

#Data label	Data Name	GEO ID
#1	early hCO rep1	GSM2589130    day30
#2	early hMGEO rep1	GSM2589129
#3	early hMGEO rep2	GSM2589133
#4	early hCO rep2	GSM2589134    day30 
#5	late hMGEO rep1	GSM2589131
#6	late hCO rep1	GSM2589132   day72 
#7	late hMGEO rep2	GSM2589135
#8	late hCO rep2	GSM2589136   day79

#use hCOs only: 1,4,6,8
# Create a mask for cells 
mask = xiang.obs_names.str.endswith('-1')
mask2 = xiang.obs_names.str.endswith('-4')
mask3 = xiang.obs_names.str.endswith('-6')
mask4 = xiang.obs_names.str.endswith('-8')
# Subset the AnnData object based on this mask
xiang_day30_b1 = xiang[mask]
xiang_day30_b2 = xiang[mask2]
xiang_day72 = xiang[mask3]
xiang_day79 = xiang[mask4]

xiang_day30_b1.obs['Source'] = 'Xiangb1_4w'
xiang_day30_b2.obs['Source'] = 'Xiangb2_4w'
xiang_day72.obs['Source'] = 'Xiang_10w'
xiang_day79.obs['Source'] = 'Xiang_11w'
###########################################################################################
#Madhavan data:12 weeks, oligo-cortical
mad = sc.read_10x_mtx('Madhavan/GSM3243667_H7noid', var_names='gene_symbols', cache=True)
# Assign the age & author as a new column in the metadata as Source
mad.obs['Source'] = 'Madhavan_12w'
###########################################################################################
#Popova data: 7 weeks, primary microglia & organoids
# load the batch1 data into pandas DataFrames
# load the data into pandas DataFrames
barcodes = pd.read_csv('barcodes.tsv', header=None, sep='\t')
genes = pd.read_csv('features.tsv', header=None, sep='\t')

# Read the sparse matrix data
# Note: use skiprows=3 to skip the header and the dimensions line
sparse_matrix = mmread('matrix.tsv')

# Adjust the format of sparse_matrix
sparse_matrix = sparse_matrix.T

# Convert to AnnData
popb1 = sc.AnnData(X=sparse_matrix,
                   obs=pd.DataFrame(index=barcodes[0]),
                   var=pd.DataFrame(index=genes[0]))
popb1.X = popb1.X.toarray()
# load the batch2 data into pandas DataFrames
barcodes = pd.read_csv('barcodes.tsv', header=None, sep='\t')
genes = pd.read_csv('features.tsv', header=None, sep='\t')

# Read the sparse matrix data
# Note: use skiprows=3 to skip the header and the dimensions line
sparse_matrix = mmread('matrix.tsv')

# Adjust the format of sparse_matrix
sparse_matrix = sparse_matrix.T

# Convert to AnnData
popb2 = sc.AnnData(X=sparse_matrix,
                   obs=pd.DataFrame(index=barcodes[0]),
                   var=pd.DataFrame(index=genes[0]))
popb2.X = popb2.X.toarray()
# Assign the age & author as a new column in the metadata as Source
popb1.obs['Source'] = 'Popova_batch1_7w'
popb2.obs['Source'] = 'Popova_batch2_7w'
###########################################################################################
#Bhaduri data
bhad = sc.read_csv('GSE132672_allorganoids_withnew_matrix.txt', delimiter='\t')
bhad = bhad.transpose()
# Extract week information and create a new column
bhad.obs['week'] = bhad.obs_names.str.extract(r'(Week\d+)', expand=False)
#split the object based on week time points
bhad_w3 = bhad[bhad.obs['week'] == 'Week3']
bhad_w5 = bhad[bhad.obs['week'] == 'Week5']
bhad_w8 = bhad[bhad.obs['week'] == 'Week8']
bhad_w10 = bhad[bhad.obs['week'] == 'Week10']
bhad_w15 = bhad[bhad.obs['week'] == 'Week15']
bhad_w24 = bhad[bhad.obs['week'] == 'Week24']

#Select for protocol type
bhad_w3.obs['protocol'] = bhad_w3.obs_names.str.extract(r'(H1S|H1X)', expand=False)

# Now you can split the AnnData object based on this new 'type' column
bhad_w3_H1S = bhad_w3[bhad_w3.obs['protocol'] == 'H1S']
bhad_w3_H1X = bhad_w3[bhad_w3.obs['protocol'] == 'H1X']

bhad_w5.obs['protocol'] = bhad_w5.obs_names.str.extract(r'(H1S|H1X)', expand=False)

# Now you can split the AnnData object based on this new 'type' column
bhad_w5_H1S = bhad_w5[bhad_w5.obs['protocol'] == 'H1S']
bhad_w5_H1X = bhad_w5[bhad_w5.obs['protocol'] == 'H1X']

bhad_w8.obs['protocol'] = bhad_w8.obs_names.str.extract(r'(H1S|H1X)', expand=False)

# Now you can split the AnnData object based on this new 'type' column
bhad_w8_H1S = bhad_w8[bhad_w8.obs['protocol'] == 'H1S']
bhad_w8_H1X = bhad_w8[bhad_w8.obs['protocol'] == 'H1X']

bhad_w10.obs['protocol'] = bhad_w10.obs_names.str.extract(r'(H1S|H1X)', expand=False)

# Now you can split the AnnData object based on this new 'type' column
bhad_w10_H1S = bhad_w10[bhad_w10.obs['protocol'] == 'H1S']
bhad_w10_H1X = bhad_w10[bhad_w10.obs['protocol'] == 'H1X']

bhad_w3_H1S.obs['Source'] = 'Bhaduri_H1S_3w'
bhad_w3_H1X.obs['Source'] = 'Bhaduri_H1X_3w'
bhad_w5_H1S.obs['Source'] = 'Bhaduri_H1S_5w'
bhad_w5_H1X.obs['Source'] = 'Bhaduri_H1X_5w'
bhad_w8_H1S.obs['Source'] = 'Bhaduri_H1S_8w'
bhad_w8_H1X.obs['Source'] = 'Bhaduri_H1X_8w'
bhad_w10_H1S.obs['Source'] = 'Bhaduri_H1S_10w'
bhad_w10_H1X.obs['Source'] = 'Bhaduri_H1X_10w'
###########################################################################################

dataset_list = [vel3, vel6, truj1, truj3, truj6, truj10, gian,
                mad, fair90, fair140, popb1, popb2, bhad_w3_H1S, bhad_w3_H1X, bhad_w5_H1S, bhad_w5_H1X,
                bhad_w8_H1S, bhad_w8_H1X, bhad_w10_H1S, bhad_w10_H1X,
                xiang_day30_b1, xiang_day30_b2, xiang_day72, xiang_day79]  


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
print("Variable genes in all batches: %d" % sum(test.var.highly_variable_nbatches == 9))
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
adata_sc.write_h5ad('C:/Users/fallo/Documents/Internship_2023/transcriptome_comp/vitro_int.h5ad')

