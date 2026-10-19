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

        if purity < pu_cutoff:
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
os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/')
output_dir=f'/sc/arion/projects/CommonMind/roussp01a/snmulti/motif/files/class/{cell}'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
pu_cutoff=0.95
if cell=='EN':
    mdata=mu.read(f'{cell}_mudata_2.h5mu')
elif cell=='Micro_PVM':
    mdata=sc.read_h5ad('RNA_3rd_subset_Micro_PVM.h5ad')
    pu_cutoff=0.9
elif cell=='Astro_Epend':
    mdata=sc.read_h5ad('RNA_2nd_subset_Astro_Epend_HVG5000.h5ad')
    pu_cutoff=0.9
elif cell in ['Oligo','OPC']:
    mdata=sc.read_h5ad(f'RNA_2nd_subset_{cell}.h5ad')
    pu_cutoff=0.9
else:
    mdata=mu.read(f'{cell}_mudata_2_update.h5mu')

embedding=pd.DataFrame(mdata.obsm['X_umap'],
                       index=mdata.obs_names,
                       columns=['UMAP1','UMAP2'])

adata=sc.read_h5ad(f'ArchR/{cell}/merged.h5ad')
adata.X.data[adata.X.data>1]=1

anno=pd.read_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno_05092025_fixed3.csv',index_col=0)
anno=anno.loc[anno.inATAC*anno.inRNA]
anno=anno.loc[~anno.subclass.eq('Epend')]
if cell =='Vascular':
    anno=anno.loc[anno['class'].isin(['Mural','Endo'])]
else:
    anno=anno.loc[anno['class'].eq(cell)]

meta=pd.read_csv(f'ArchR/{cell}/cellmetadata.csv',index_col=0)
sh=np.intersect1d(embedding.index,meta.index)
sh=np.intersect1d(sh,anno.index)

embedding=embedding.loc[sh]
meta=meta.loc[sh]
adata=adata[sh,]
f=np.array(adata.X.sum(axis=1)).flatten()
meta['RiP']=f
cell_depth=meta['RiP'].values
threshold=meta.RiP.mean()*20

pbulk_centers = []
barcode_groups = []
cell_barcodes = sh
embed=embedding.values
# Fit a KNN on our cells

nbrs = NearestNeighbors(n_neighbors=embed.shape[0], algorithm='ball_tree').fit(embed)
anno=anno.loc[sh]
if cell in ['Micro_PVM','Astro_Epend','Oligo','OPC']:
    meta['celltype']=anno['subtype']
else:
    meta['celltype']=anno['subclass']
# Go through each cell type
results = Parallel(n_jobs=8, backend="loky")(
    delayed(process_one_celltype)(ct)
    for ct in np.unique(meta["celltype"])
)
barcode_groups = pd.concat([df for _, dfs in results for df in dfs], ignore_index=True)
barcode_groups.to_csv(f'{output_dir}/barcode_groups.csv.gz',compression='gzip')
embedding.to_csv(f'{output_dir}/raw_embed.csv.gz',compression='gzip')


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

write_input4R(adata, f'{output_dir}/metacell_20')
