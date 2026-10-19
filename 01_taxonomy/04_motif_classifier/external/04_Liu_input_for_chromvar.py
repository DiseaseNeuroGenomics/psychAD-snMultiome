#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 18:27:02 2025

@author: pengfeidong
"""
import muon as mu
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import statistics as stat
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
import os
import sys
from scipy import sparse
from joblib import Parallel, delayed
from pathlib import Path
import h5py
import pyranges as pr
def get_pseduobulk(adata,colname,mean=False):
    z=pd.get_dummies(adata.obs[colname])    
    mat=pd.get_dummies(z).to_numpy()
    if mean:
        mat=mat/mat.sum(axis=0)
    mat=sparse.csr_matrix(mat.T)
    psb=mat.dot(adata.X)
    psb=pd.DataFrame(psb.T.todense(),columns=z.columns,index=adata.var.index)   
    return psb
def process_one_celltype(celltype):

    cell_inds = np.where(meta["celltype"] == celltype)[0]
    target_n = int(min(meta.iloc[cell_inds]["RiP"].sum() / threshold*1.5, 1000))

    pbulk_centers_local = []
    barcode_groups_local = []

    rng = np.random.default_rng()

    while len(pbulk_centers_local) < target_n:

        # random center
        new_center = rng.choice(cell_inds)

        # nearest neighbors
        distances, indices = nbrs.kneighbors([embed[new_center, :]])
        nbr_inds = indices[0]

        # cumulative depth
        cumulative_depth = np.cumsum(cell_depth[nbr_inds])
        over = np.where(cumulative_depth > threshold)[0]
        if len(over) == 0:
            continue
        n_members = over[0]
        pbulk_members = nbr_inds[:n_members]

        # purity
        celltype_labels = meta["celltype"].values[pbulk_members]
        uniq, counts = np.unique(celltype_labels, return_counts=True)
        purity = counts.max() / counts.sum()

        if purity < 0.9:
            continue

        # success: store one pseudobulk
        pbulk_centers_local.append(new_center)
        group_id = f"{celltype}_pbulk_{len(pbulk_centers_local)}"

        barcode_groups_local.append(
            pd.DataFrame({
                "barcode": cell_barcodes[pbulk_members],
                "group": group_id
            })
        )

    return pbulk_centers_local, barcode_groups_local
def write_10X_h5(adata, file):
    if '.h5' not in file:
        file = f'{file}.h5'
    if Path(file).exists():
        raise FileExistsError(f"There already is a file `{file}`.")
    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    def str_max(x):
        return max([len(i) for i in x])
    if os.path.exists(file):
        os.remove(file)
    w = h5py.File(file, 'w')
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset("feature_type", data=np.repeat('Gene Expression',adata.var.shape[0]).astype('S15'))
    ftrs.create_dataset("genome", data=np.repeat('GRCh38',adata.var.shape[0]).astype('S6'))
    ftrs.create_dataset("id", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1], dtype=f'<i{int_max(adata.X.shape)}'))
cell=sys.argv[1]
map1={'EN':'Exc','IN':'Inh'}
cell2=map1[cell]
os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSHM_ATAC')
output_dir=f'/sc/arion/projects/CommonMind/roussp01a/snmulti/motif/files/class/{cell}'

adata=sc.read_h5ad(f'MultiRegion_cortical_ATAC_{cell2}.h5ad')
if cell=='EN':
    adata=adata[~adata.obs['ATAC.Subtype'].isin(['Exc CA1 pyramidal cells','Exc NXPH1 RNF220','Exc CA2, CA3 pyramidal cells']),]
if cell=='Micro_PVM':
    adata=adata[~adata.obs['ATAC.Subtype'].isin(['T']),]

adata.X.data[adata.X.data>1]=1
adata.var['cvg']=np.array(adata.X.sum(axis=0)).flatten()/adata.obs.shape[0]


z=pd.read_csv(f'/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/{cell}/peak_h5ad.csv',index_col=0)

z.columns=['Chromosome','Start','End']+z.columns[3:].tolist()
ref=pr.PyRanges(z)
peak=pr.PyRanges(adata.var)
peak=peak.overlap(ref)
sel=(adata.var.interval.isin(peak.interval.values)|(adata.var.cvg>0.015))*(adata.var.cvg>0.001)
adata=adata[:,sel]

adata.obs['celltype']=adata.obs['ATAC.Subtype']
adata.obs['RiP']=np.array(adata.X.sum(axis=1)).flatten()
z=adata.obs.celltype.value_counts()
adata=adata[adata.obs.celltype.isin(z.loc[z>100].index),]

embedding=pd.DataFrame(adata.obsm['X_umap'],
                       index=adata.obs_names,
                       columns=['UMAP1','UMAP2'])
meta=adata.obs
sh=adata.obs_names

cell_depth=meta['RiP'].values
threshold=meta.RiP.mean().round()*20

pbulk_centers = []
barcode_groups = []
cell_barcodes = sh
embed=embedding.values
# Fit a KNN on our cells

nbrs = NearestNeighbors(n_neighbors=embed.shape[0], algorithm='ball_tree').fit(embed)
meta['celltype']=meta['ATAC.Subtype']
# Go through each cell type
results = Parallel(n_jobs=16, backend="loky")(
    delayed(process_one_celltype)(ct)
    for ct in np.unique(meta["celltype"])
)
barcode_groups = pd.concat([df for _, dfs in results for df in dfs], ignore_index=True)
barcode_groups.to_csv(f'{output_dir}/RUSHMR_barcode_groups.csv.gz',compression='gzip')
embedding.to_csv(f'{output_dir}/RUSHMR_raw_embed.csv.gz',compression='gzip')


z=pd.get_dummies(barcode_groups.group)
mat=z.values
mat=mat/mat.sum(axis=0)
embed=pd.DataFrame(mat.T.dot(embedding.loc[barcode_groups.barcode.values]),
                   columns=['UMAP1','UMAP2'])
embed.index=z.columns.values
cov=meta.loc[barcode_groups.barcode.values,
             ['TSSEnrichment','NucleosomeRatio','BlacklistRatio','nFrags']]
cov=pd.DataFrame(mat.T.dot(cov),columns=['TSSEnrichment','NucleosomeRatio',
                                         'BlacklistRatio','nFrags'])
cov.index=z.columns.values



ad1=adata[barcode_groups.barcode.values,]
ad1.obs['pb']=barcode_groups['group'].values
pb=get_pseduobulk(ad1,'pb')
obs=pd.concat([embed,cov.loc[embed.index]],axis=1).loc[pb.columns]

adata=ad.AnnData(sparse.csr_matrix(pb.values.T),
              obs=obs,var=adata.var)

def write_input4R(adata,output_prefix):
    if os.path.exists(f"{output_prefix}.h5"):
        os.remove(f"{output_prefix}.h5")
    write_10X_h5(adata, f"{output_prefix}.h5")
    adata.obs.to_csv(f'{output_prefix}.obs.csv.gz', compression='gzip')
    adata.var.to_csv(f'{output_prefix}.var.csv.gz', compression='gzip')

write_input4R(adata, f'{output_dir}/RUSHMR_metacell_20')    

