import os
import sys
import pandas as pd
import numpy as np
from scipy import sparse
import scanpy as sc
import anndata as ad
import sys
import lightgbm as lgb
import pickle

def normalize_adata(adata):
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
def read_genes():
    genes=pd.read_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/metadata/cellranger_gene.csv')
    genes.loc[genes['gene_name'].duplicated(),'gene_name']=genes.loc[genes['gene_name'].duplicated(),'gene_name']+'-1'
    genes['mt'] = genes.gene_name.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    genes['ribo'] = genes.gene_name.str.startswith('RPL') | genes.gene_name.str.startswith('RPS')
    genes["hb"] = genes.gene_name.str.contains(("^HB[^(P)]"))
    genes['protein_coding'] = genes.gene_type.eq('protein_coding')
    genes['autosomeX']=genes.seqnames.isin(set([f'chr{x+1}' for x in range(22)]+['chrX']))
    genes['filtered']=genes.protein_coding & (~genes.mt) & (~genes.ribo) & (genes.autosomeX)
    genes=genes.loc[genes.filtered]
    return genes
os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate3') 
genes=read_genes()

input_fold='/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSHC_multiome/processed'
Obs=pd.read_csv(f'{input_fold}/rosmap_diversity_snmultiome_celltypes_metadata_20230410.csv')
Obs['barcodekey']=Obs.batch+'_'+Obs.barcode
Obs.set_index('barcodekey',inplace=True)


map1=dict(zip(['Astro','Vascular','Epend','IN','EN','Micro_PVM','Oligo','OPC','Adaptive'],
              ['Astro','Endo','Epend','Gaba','Glut','Micro','Olig','Opc','OtherImmune']))

map2=dict(zip(["caudate nucleus","dorsolateral prefrontal cortex","superior temporal gyrus"],
              ['CN','PFC','STG']))

#adata.obs['Region']=adata.obs.region.map(map1)
#adata.obs['class']=adata.obs.broadcellclass.map(map1)

dfs=[]
for cell in ['Astro','Micro_PVM','OPC','EN','IN','Oligo']:    
    cell2=map1[cell]
    adata=sc.read_10x_h5(f'{input_fold}/counts_{cell2}_filtered.h5')
    obs=Obs.loc[Obs.index.isin(adata.obs_names)]
    obs['SubID']=Obs.loc[obs.index,'individualID']
    obs['class']=cell
    obs['subclass']=obs['subcellclass']
    obs['Region']=obs.region.map(map2)
    obs=obs.loc[~obs.Region.eq('CN')]

    avg=pd.read_csv(f'{cell}/gene_avg2.csv',index_col=0)
    train_genes = avg.index.to_numpy()
    train_means = avg.iloc[:, 0].to_numpy(dtype=np.float32)

    shared_genes = adata.var_names[adata.var_names.isin(train_genes)]
    adata = adata[obs.index, shared_genes].copy()
    adata.obs=obs.loc[adata.obs_names]
    
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    shared_mask = np.isin(train_genes, adata.var_names)
    shared_idx_train = np.where(shared_mask)[0]
    shared_genes_ordered = train_genes[shared_mask]
    adata_idx = pd.Index(adata.var_names).get_indexer(shared_genes_ordered)


    class_models = []
    reg_models = []
    for x in range(20):
        with open(f'{cell}/lgb_classify_redo/model_fold{x}.pkl', 'rb') as f:
            class_models.append(pickle.load(f))
        with open(f'{cell}/lgb_regression_redo/model_fold{x}.pkl', 'rb') as f:
            reg_models.append(pickle.load(f))

    n_cells = adata.n_obs
    n_folds = 20
    n_classes = 3
    
    Probs = np.full((n_folds, n_cells, n_classes), np.nan, dtype=np.float32)
    Preds = np.full((n_cells, n_folds), np.nan, dtype=np.float32)
    
    chunk_size = 10000
    
    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        n = end - start

        X_chunk = np.tile(train_means, (n, 1)).astype(np.float32)
        X_chunk[:, shared_idx_train] = adata.X[start:end, :][:, adata_idx].toarray().astype(np.float32)

    # ensemble prediction
        for x in range(n_folds):
            Probs[x, start:end, :] = class_models[x].predict(X_chunk)
            Preds[start:end, x] = reg_models[x].predict(X_chunk)


    mean_probs = np.nanmean(Probs, axis=0)
    y_pred_classes = np.argmax(mean_probs, axis=1)
    adata.obs['Dx_cat_cell']=np.array(['control','early','late'])[y_pred_classes]
    adata.obs['Dx_cell']=np.nanmean(Preds,axis=1)
    dfs.append(adata.obs)
df=pd.concat(dfs)
df=df[['broadcellclass', 'subcellclass', 'cogdx_rush', 'braak_rush', 'cerad_rush', 'sex',
 'SubID', 'class', 'subclass', 'Region', 'Dx_cat_cell', 'Dx_cell']]
df.to_csv('RUSHC_Dx_cell.csv.gz',compression='gzip')

