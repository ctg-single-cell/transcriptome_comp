#Created by: Fallon Ratner 15-09-23
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

#Read in files
vivo = sc.read("vivo_noint_1108.h5ad")
vitro = sc.read("vitro_int_09_08.h5ad")

dataset_list = [vivo, vitro] 
combined = sc.concat(dataset_list, axis=0, join='outer')
combined.X[np.isnan(combined.X)] = 0
 
# Perform preprocessing steps
sc.pp.filter_cells(combined, min_genes=200)
sc.pp.filter_genes(combined, min_cells=3)
combined.var['mt'] = combined.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(combined, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.normalize_total(combined, target_sum=1e4)
sc.tl.pca(combined)

# compute hierarchical clustering using PCs 
sc.tl.dendrogram(combined, 'Source')
plt.figure(figsize=(10, 8))
ax = sc.pl.dendrogram(combined, 'Source')
plt.savefig("dendrogram_vivvit_1509.pdf")

#compute Pearson's correlation
ax = sc.pl.correlation_matrix(combined, 'Source', figsize=(10,10))
plt.savefig("matrix_vivvit_1509.pdf")