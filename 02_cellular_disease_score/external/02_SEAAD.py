#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:03:56 2026

@author: pengfeidong
"""

import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import lightgbm as lgb
from sklearn.metrics import mean_squared_error
from sklearn.metrics import accuracy_score
import pickle

def filter_genes(adata,includeX=True):
    adata.var['protein_coding']=adata.var.gene_type.eq('protein_coding')
    adata.var['mt']=adata.var.index.str.startswith('MT-')
    adata.var['ribo']=adata.var.index.str.startswith('RPL') | adata.var.index.str.startswith('RPS')
    adata.var['AS']= adata.var.index.str.endswith('-AS')
    adata.var['autosome']=adata.var.chr.isin(set([f'chr{x+1}' for x in range(22)]))
    if 'PC' in adata.var.columns:
        adata.var.drop(columns=['PC'],inplace=True)
    if includeX:
        adata.var['filtered']=adata.var.protein_coding & (~adata.var.mt) & (~adata.var.ribo) & (adata.var.autosome | adata.var.chr.eq('chrX'))
    else:
        adata.var['filtered']=adata.var.protein_coding & (~adata.var.mt) & (~adata.var.ribo) & adata.var.autosome

def process_adata(adata,includeX=True,min_cell_percent=0.001):
    min_cell=int(adata.obs.shape[0]*min_cell_percent)
    sc.pp.filter_genes(adata, min_cells=min_cell)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)


os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files')
raw=sc.read_h5ad('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/SEA-AD/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.h5ad')
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
    avg=pd.read_csv(f'ML_cellstate3/{cell}/gene_avg.csv',index_col=0) # the average expression of the training data
    adata = raw[raw.obs['class'].eq(cell),raw.var_names.isin(avg.index)].copy()
    adata.X = adata.layers["UMIs"]    
    del adata.layers
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    ###
    X_train=adata[:,avg.index].X
    Probs = np.full([20,adata.obs.shape[0],3], np.nan)
    Preds = np.full([adata.obs.shape[0],20], np.nan)
    for x in range(20):
        with open(f'ML_cellstate3/{cell}/lgb_classify_redo/model_fold{x}.pkl', 'rb') as f:
            model=pickle.load(f)
        
        y_proba = model.predict(X_train)
        Probs[x,:,:]=y_proba
        with open(f'ML_cellstate3/{cell}/lgb_regression_redo/model_fold{x}.pkl', 'rb') as f:
            model=pickle.load(f)
        
        Preds[:,x]=model.predict(X_train)
    y_pred_classes = np.argmax(np.nanmean(Probs,axis=0),axis=1)
    adata.obs['Dx_cat_cell']=np.array(['control','early','late'])[y_pred_classes]
    adata.obs['Dx_cell']=np.nanmean(Preds,axis=1)
    dfs.append(adata.obs)
df=pd.concat(dfs)
df.to_csv('ML_cellstate3/SEA_PFC_Dx_cell.csv.gz',compression='gzip')


