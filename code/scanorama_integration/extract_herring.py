"""
Created by: Fallon Ratner
Date: 20/04/23
Project: snRNAseq Annotation Prototype Analysis
"""

import scanpy as sc


adata = sc.read("Processed_data_RNA-all_full-counts-and-downsampled-CPM.h5ad")

adata


# Subset cells from ga22
adata_batch1 = adata[adata.obs["batch"] == "RL2103_ga22_v3", :].copy()
adata_batch1.obs['CellID'] = adata_batch1.obs.index
adata_batch1.obs["batch"] = "RL2103_ga22_v3"
adata_batch1.write('ga22_cleaned_count_matrices.h5ad')

# Subset cells from ga24
adata_batch2 = adata[adata.obs["batch"] == "RL2107_ga24_v3", :].copy()
adata_batch2.obs['CellID'] = adata_batch2.obs.index
adata_batch2.obs["batch"] = "RL2107_ga24_v3"
adata_batch2.write('ga24_cleaned_count_matrices.h5ad')

# Subset cells from ga34
adata_batch3 = adata[adata.obs["batch"] == "RL2121_ga34_v3", :].copy()
adata_batch3.obs['CellID'] = adata_batch3.obs.index
adata_batch3.obs["batch"] = "RL2121_ga34_v3"
adata_batch3.write('ga34_cleaned_count_matrices.h5ad')

