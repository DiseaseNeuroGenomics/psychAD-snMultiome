import os
import numpy as np
import scanpy as sc
import seaborn as sns
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib.colors as mtpltcl
import glob
import re
import scanpy.external as sce
import pandas as pd
import sys
import subprocess
import pegasus as pg
sys.path.append("/sc/arion/projects/roussp01a/pengfei/tools/scanpy_myextra/")
import myutsc
import h5py
import muon as mu
os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files')
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc._settings.ScanpyConfig.n_jobs=16
sc.settings.set_figure_params(dpi=80, facecolor='white')

Ma2022_markers = pd.read_csv('Ma_2022_DLPFC_markers_clean.csv')
Ma2022_markers['Pericyte']=Ma2022_markers['PC']
Ma2022_markers.drop(columns='PC',inplace=True)
ress=[0.1,0.2,0.5,1]
nnb=30
prefix='RNA_2nd'
cluster_names=['2leiden'+str(res).replace('.', '_') for res in ress]

N_hvg=3000
npc=10

for c in 'EN','IN','Vascular','Immune','Astro','OPC','Oligo':
    adata=sc.read_h5ad(f'RNA_2nd_subset_{c}_raw.h5ad')
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
        
    z=adata.obs.newID.value_counts()    
    z=sc.pp.highly_variable_genes(adata[adata.obs.manual_batch.eq('other')*adata.obs.newID.isin(z.loc[z>10].index.values),], min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key='newID',inplace=False)
    adata.var['filtered']=adata.var['autosome']*adata.var['protein_coding']*(~adata.var['hb'])*(~adata.var['ribo'])*(~adata.var['mt'])
    adata.var.highly_variable=adata.var.index.isin(z.loc[z.index.isin(adata.var.loc[adata.var.filtered].index)*(z.means>0.0125)*(z.means<3)].sort_values('dispersions_norm').index[-N_hvg:])
    adata_raw=adata.copy()
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps =50)
    sc.pl.pca_variance_ratio(adata,50)
    adata.obsm['X_pca']=adata.obsm['X_pca'][:,:npc]
    sce.pp.harmony_integrate(adata, 'newID',max_iter_harmony = 40)
    sc.pp.neighbors(adata, use_rep = "X_pca_harmony",n_neighbors=nnb,metric='l2')
    sc.tl.umap(adata,maxiter=50)
    for i in range(len(ress)):
        res=ress[i]
        cluster_name=cluster_names[i]
        sc.tl.leiden(adata, resolution=res,key_added=cluster_name)
        if adata.obs[cluster_name].value_counts().shape[0]==1:
            continue
        mycolor=sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, as_cmap=True)
        sns.set_style("white")
        zdf=adata.obs.groupby(cluster_name).agg(dict(zip(Ma2022_markers.columns.tolist(),['mean']*len(Ma2022_markers.columns.tolist()))))
        sns.clustermap(zdf,  figsize=(15,5),cmap=mycolor,edgecolor="none",linewidths=0)


    adata_raw=adata_raw[adata.obs.index ,]
    adata_raw.obs=adata.obs
    adata_raw.uns=adata.uns
    adata_raw.obsm=adata.obsm
    adata_raw.obsp=adata.obsp
    adata_raw.var['highly_variable']=adata_raw.var.index.isin(adata.var.index)
    adata=adata_raw
    for cells in  [Ma2022_markers.keys()[i:i+3] for i in range(0, len(Ma2022_markers.keys()), 3)]:
        sc.pl.umap(adata,color=list(cells))
    for i in range(len(ress)):
        if adata.obs[cluster_name].value_counts().shape[0]==1:
            continue
        res=ress[i]
        cluster_name=cluster_names[i]
        pg.de_analysis(adata, cluster=cluster_name,de_key='de_'+cluster_name)
        marker_dict = pg.markers(adata,de_key='de_'+cluster_name)
        #for x in adata.obs[cluster_name].unique():
        #    pg.volcano(adata, cluster_id = x, dpi=200,de_key='de_'+cluster_name)
        pg.write_results_to_excel(marker_dict, prefix+'_'+c+'_'+cluster_name+".de.xlsx")
    adata.write( f'RNA_2nd_subset_{c}.h5ad',compression='gzip' 
