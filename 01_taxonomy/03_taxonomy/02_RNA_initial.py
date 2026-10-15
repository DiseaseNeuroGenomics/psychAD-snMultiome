#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 14:10:53 2026

@author: pengfeidong
"""

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
from gtfparse import read_gtf
os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files')
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc._settings.ScanpyConfig.n_jobs=1
sc.settings.set_figure_params(dpi=80, facecolor='white')
N_hvg=5000
ngb_min=500
ngb_max=10000
ntotal_max=100000
pct_mt=1
n_neighbor=100
N_PC=30
prefix=f'RNA_all_initial_cutoff_ngene_{ngb_min}_{ngb_max}_ntotal_{ntotal_max}_pct_mito{pct_mt}_n_neighbors{n_neighbor}_Nhvg{N_hvg}_NPC{N_PC}'
########################################
adata=sc.read_h5ad('RNA_all_raw.h5ad')
adata.var['percent_cell']=adata.var.n_cells_by_counts/adata.X.shape[0]
adata.obs['log10_total_counts']=np.log10(adata.obs['total_counts'])
adata.obs['scDblFinder_atac']=adata.obs['scDblFinder_atac'].values.reorder_categories(['singlet','doublet'])

print('N cell raw: ' +str(adata.obs.shape[0])+'\n')
adata = adata[adata.obs.total_counts< ntotal_max, :]
print(f'N cell total < {ntotal_max}: ' +str(adata.obs.shape[0])+'\n')
adata = adata[adata.obs.n_genes_by_counts< ngb_max, :]
print(f'N cell total < {ngb_max}: ' +str(adata.obs.shape[0])+'\n')

adata = adata[adata.obs.n_genes_by_counts> ngb_min, :]
print(f'N cell N_gene_counts > {ngb_min}: ' +str(adata.obs.shape[0])+'\n')

adata = adata[adata.obs.pct_counts_mt < pct_mt, :]
print(f'N cell pct_counts_mt <{pct_mt}: ' +str(adata.obs.shape[0])+'\n')
print('N cell: ' +str(adata.obs.shape[0])+'\n')
adata=adata[adata.obs.scDblFinder_gex.eq(1)*adata.obs.scDblFinder_atac.eq('singlet'),:]
print('doublet: ' +str(adata.obs.shape[0])+'\n')


# highly variable genes
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.var['filtered']=adata.var['autosome']*adata.var['protein_coding']*(~adata.var['hb'])*(~adata.var['ribo'])*(~adata.var['mt'])
z=sc.pp.highly_variable_genes(adata[adata.obs.manual_batch.eq('other'),], min_mean=0.0125, max_mean=3, min_disp=0.5,batch_key='newID',inplace=False)

sns.lmplot(x='means',y='dispersions',hue='highly_variable', data=z, fit_reg=False)
sns.lmplot(x='means',y='dispersions_norm',hue='highly_variable', data=z, fit_reg=False)

adata.var['filtered']=adata.var['autosome']*adata.var['protein_coding']*(~adata.var['hb'])*(~adata.var['ribo'])*(~adata.var['mt'])
adata.var.highly_variable=adata.var.index.isin(z.loc[z.index.isin(adata.var.loc[adata.var.filtered].index)*(z.means>0.0125)*(z.means<3)].sort_values('dispersions_norm').index[-N_hvg:])

#PCA KNN UMAP
cell_cycle_genes = [x.strip() for x in open('/sc/arion/projects/roussp01a/pengfei/ref/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

adata_raw=adata.copy()
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

sc.set_figure_params(scanpy=True, dpi=150, dpi_save=300)
sc.tl.pca(adata, svd_solver='arpack', n_comps =50)
sc.pl.pca_variance_ratio(adata, log=True,n_pcs =50)
adata.obsm['X_pca']=adata.obsm['X_pca'][:,:N_PC]
sc.pp.neighbors(adata,n_neighbors=n_neighbor,metric='l2')
sc.tl.umap(adata,maxiter=50)
# harmony 
sce.pp.harmony_integrate(adata, 'newID',max_iter_harmony = 40)
sc.pp.neighbors(adata, use_rep = "X_pca_harmony",n_neighbors=n_neighbor,metric='l2')
sc.tl.umap(adata,maxiter=50)
sc.tl.leiden(adata, resolution=0.5,key_added='leiden0_5')
sc.tl.leiden(adata, resolution=0.2,key_added='leiden0_2')
sc.tl.leiden(adata, resolution=0.1,key_added='leiden0_1')
# marker genes
sc.pl.umap(adata, color=['leiden0_1','leiden0_2'])
sc.pl.umap(adata, color=['leiden0_5'])
brain_cell_markers=myutsc.brain_cell_markers
marker_genes_in_data=dict()

adata_raw=adata_raw[adata.obs.index ,]
adata_raw.obs=adata.obs
adata_raw.uns=adata.uns
adata_raw.obsm=adata.obsm
adata_raw.obsp=adata.obsp
adata_raw.var['highly_variable']=adata_raw.var.index.isin(adata.var.index)
adata=adata_raw

for ct, markers in brain_cell_markers.items():
    markers_found = list()
    for marker in markers:
        if marker in adata.var.index:
            markers_found.append(marker)
    if len(markers_found)>0:
        marker_genes_in_data[ct] = markers_found

for ct in marker_genes_in_data.keys():
    print(ct)  # print cell subtype name
    sc.pl.umap(adata,color=marker_genes_in_data[ct],
        vmin=0,vmax="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
        sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
        frameon=False,
        cmap="Reds",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    )
    print("\n\n\n") 
for cutoff in ['0_1','0_2']:
    print('leiden'+cutoff)
    sc.pl.dotplot(adata,
    groupby="leiden"+cutoff,
    var_names=marker_genes_in_data,
    standard_scale="var")
#module score
Ma2022_markers = pd.read_csv('Ma_2022_DLPFC_markers_clean.csv')
Ma2022_markers['Pericyte']=Ma2022_markers['PC']
Ma2022_markers.drop(columns='PC',inplace=True)
for x in Ma2022_markers.keys():
    sc.tl.score_genes(adata,Ma2022_markers[x][Ma2022_markers[x].isin(adata.var.index)][:100],score_name=x,ctrl_size=400)
    sc.pl.umap(adata, color=[x,'leiden0_1','leiden0_2'])
anno_map=dict()
anno_map['0']='Oligo'
anno_map['1']='EN'
anno_map['2']='Astro'
anno_map['3']='Immune'
anno_map['4']='IN'
anno_map['5']='OPC'
anno_map['6']='IN'
anno_map['7']='Vascular'
anno_map['8']='EN'
anno_map['9']='Vascular'
anno_map['10']='Immune'
adata.obs['cell_type']=adata.obs.leiden0_1.map(anno_map)
adata.write(f'{prefix}_qc_01.h5ad',compression='gzip')
