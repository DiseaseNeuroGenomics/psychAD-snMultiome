import os
import scanpy as sc
import pandas as pd
import numpy as np
import lightgbm as lgb
import pickle
import argparse
parser = argparse.ArgumentParser(description='frozen cell Dx_cell projection')
parser.add_argument('--input', help='Input h5ad file', required=True)
parser.add_argument('--model_dir', help='frozen model dir', required=True)
parser.add_argument('--output_path', help='output file path', required=True)
parser.add_argument('--obs_output', help='qc obs file path', required=True)
args = parser.parse_args()
raw_path=args.input
model_dir=args.model_dir
output_path=args.output_path
obs_output=args.obs_output

raw=sc.read_h5ad(raw_path)
raw.obs['class']=raw.obs['Class'].astype('str')
raw.obs.loc[raw.obs['Class'].eq('Non-neuronal and Non-neural'),'class']=raw.obs.loc[raw.obs['Class'].eq('Non-neuronal and Non-neural'),'Subclass'].astype('str')
name_map={'Neuronal: Glutamatergic':'EN','Neuronal: GABAergic':'IN',
'Oligodendrocyte':"Oligo",'Astrocyte':'Astro',
'Microglia-PVM':'Micro_PVM','OPC':'OPC','VLMC':'Mural',
'Endothelial':'Endo'}

raw.obs['class']=raw.obs['class'].map(name_map)
raw.obs['class'].value_counts()
mask=raw.var.gene_ids.str.contains('ENSG')
raw.var.loc[mask,'gene_ids']=raw.var.loc[mask,'gene_ids'].str.replace(' ENSG.*','-1',regex=True)
raw.var.set_index('gene_ids',inplace=True) 

dfs=[]
for cell in 'EN','IN','Astro','Micro_PVM','Oligo','OPC':
    avg=pd.read_csv(f'{model_dir}/{cell}/ref_avg_expression.csv',index_col=0) # the average expression of the training data
    train_genes = avg.index.to_numpy()
    train_means = avg.iloc[:, 0].to_numpy(dtype=np.float32)

    shared_genes = raw.var_names[raw.var_names.isin(train_genes)]    
    adata = raw[raw.obs['class'].eq(cell),shared_genes].copy()
    adata.X = adata.layers["UMIs"]    
    del adata.layers
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    ###
    shared_mask = np.isin(train_genes, adata.var_names)
    shared_idx_train = np.where(shared_mask)[0]
    shared_genes_ordered = train_genes[shared_mask]
    adata_idx = pd.Index(adata.var_names).get_indexer(shared_genes_ordered)
    class_models = []
    reg_models = []
    n_cells = adata.n_obs
    n_folds = 20
    n_classes = 3
    
    for x in range(n_folds):
        with open(f'{model_dir}/{cell}/lgb_classify/model_fold{x}.pkl', 'rb') as f:
            class_models.append(pickle.load(f))
        with open(f'{model_dir}/{cell}/lgb_regression/model_fold{x}.pkl', 'rb') as f:
            reg_models.append(pickle.load(f))

    Probs = np.full((n_folds, n_cells, n_classes), np.nan, dtype=np.float32)
    Preds = np.full((n_cells, n_folds), np.nan, dtype=np.float32)
    
    chunk_size = 5000
    
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
df.to_csv(f'{output_path}',compression='gzip')
raw.X=raw.layers['UMIs']
raw.var['mt'] = raw.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
raw.var['ribo'] = raw.var_names.str.startswith('RPL') | raw.var_names.str.startswith('RPS')
raw.var["hb"] = raw.var_names.str.contains(("^HB[^(P)]"))
sc.pp.calculate_qc_metrics(raw, qc_vars=['mt','ribo','hb'], percent_top=None,inplace=True)
raw.obs.to_csv(obs_output,compression='gzip')