
#Created by Fallon Ratner
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

#read in integrated file
combo = sc.read("vivo_int.h5ad")    
#Save the datasets as a separate scanpy object
sources = combo.obs['Source'].unique().tolist()
adata_dict = {source: combo[combo.obs['Source'] == source].copy() for source in sources}
herr22 = adata_dict['Herring_GA22']
herr22.write_h5ad('herr22_24_7.h5ad')
herr24 = adata_dict['Herring_GA24']
herr24.write_h5ad('herr24_24_7.h5ad')
herr34 = adata_dict['Herring_GA34']
herr34.write_h5ad('herr34_24_7.h5ad')
herr2d = adata_dict['Herring_Day2']
herr2d.write_h5ad('herr2d_24_7.h5ad')
pol = adata_dict['Polioudakis_GW17-18']
pol.write_h5ad('pol_24_7.h5ad')
han11b1 = adata_dict['Han_GW11.1']
han11b1.write_h5ad('han11b1_24_7.h5ad')
han11b2 = adata_dict['Han_GW11.2']
han11b2.write_h5ad('han11b2_24_7.h5ad')
han12 = adata_dict['Han_GW12']
han12.write_h5ad('han12_24_7.h5ad')
han13 = adata_dict['Han_GW13']
han13.write_h5ad('han13_24_7.h5ad')
cou13 = adata_dict['Couturier_GW13']
cou13.write_h5ad('cou13_24_7.h5ad')
cou17 = adata_dict['Couturier_GW17']
cou17.write_h5ad('cou17_24_7.h5ad')
cou19 = adata_dict['Couturier_GW19']
cou19.write_h5ad('cou19_24_7.h5ad')
liu = adata_dict['Liu_GW17-19']
liu.write_h5ad('liu_24_7.h5ad')

###############################################################################
#read in integrated file
vitro = sc.read('vitro_int.h5ad')
#Save the datasets as a separate scanpy object
sources = vitro.obs['Source'].unique().tolist()
adata_dict = {source: vitro[vitro.obs['Source'] == source].copy() for source in sources}
truj1mo = adata_dict['Trujillo_1mo']
truj1mo.write_h5ad('truj1mo_26_7.h5ad')
truj3mo = adata_dict['Trujillo_3mo']
truj3mo.write_h5ad('truj3mo_26_7.h5ad')
truj6mo = adata_dict['Trujillo_6mo']
truj6mo.write_h5ad('truj6mo_26_7.h5ad')
truj10mo = adata_dict['Trujillo_10mo']
truj10mo.write_h5ad('truj10mo_26_7.h5ad')
gian = adata_dict['Giandomenico_2.5mo']
gian.write_h5ad('gian_26_7.h5ad')
mengu1 = adata_dict['Meng_U1M_1.5mo']
mengu1.write_h5ad('mengu1m_26_7.h5ad')
mengu2 = adata_dict['Meng_U2F_1.5mo']
mengu2.write_h5ad('mengu2f_26_7.h5ad')
vel3mo = adata_dict['Velasco_3mo']
vel3mo.write_h5ad('vel3mo_26_7.h5ad')
vel6mo = adata_dict['Velasco_6mo']
vel6mo.write_h5ad('vel6mo_26_7.h5ad')




