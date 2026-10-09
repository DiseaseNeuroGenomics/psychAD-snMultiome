#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:57:21 2026

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
################
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

def subset_by_cat(stratify,cat2,cellfrac,target=5000, T=10):
    np.random.seed(42)
    N = stratify.shape[0]
    cat_vals = np.unique(stratify)
    n_cat=cat_vals.shape[0]
    celltarget=np.ceil(cellfrac/cellfrac.sum()*target/n_cat).astype(int)
    SEL=[[] for i in range(T)]
    for c in cat_vals:
        for xx in celltarget.index.values:
            ind=np.arange(N)[(stratify==c)*(cat2==xx)]
            p=np.ones(ind.shape,dtype=np.float64)
            p=p/np.sum(p)
            for x in range(T):
                sel=np.random.choice(ind.shape[0],celltarget[xx],replace=False,p=p)
                p[sel]*=0.0001
                p=p/np.sum(p)
                SEL[x].append(ind[sel])
    SEL=[np.sort(np.concatenate(x)) for x in SEL]
    return SEL
def find_neg(pos,L):
    f=np.ones(L)
    f[pos]=0
    return np.where(f)[0]

def lgb_regress(X_train,X_valid,y_train,y_valid):
    train_data = lgb.Dataset(X_train, label=y_train)
    valid_data = lgb.Dataset(X_valid, label=y_valid, reference=train_data)
    params = {
    "num_threads":36,
    "objective": "huber", #"regression",
    "alpha": 0.9,  # Controls the threshold for outlier residuals
    "metric": "mae", #"rmse"
    "boosting_type": "gbdt",  # Gradient Boosting Decision Tree
    "num_leaves": 31,  # Number of leaves in each tree
    "learning_rate": 0.05,  # Learning rate
    "feature_fraction": 0.8,  # Fraction of features to use per iteration
    "bagging_fraction": 0.8,  # Fraction of data to use per iteration
    "bagging_freq": 5,  # Perform bagging every 5 iterations
    "verbose": -1}
    model = lgb.train(
    params,
    train_data,
    num_boost_round=1000,  # Maximum number of boosting iterations
    valid_sets=[train_data, valid_data],
    callbacks=[lgb.early_stopping(stopping_rounds=50), lgb.log_evaluation(50)])
    return model
if __name__ == '__main__':
    cell=sys.argv[1]
    level='class'
    os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files')
    os.makedirs(f'ML_cellstate3/{cell}/lgb_regression',exist_ok =True)
    anno=pd.read_csv('./cell_anno_05092025_fixed3.csv',index_col=0)
    anno=anno.loc[anno.inRNA]
    anno=anno.loc[~anno.subtype.eq('ependA')]
    
    anno['celltype']=anno['subtype']
    #sel=anno.level1_5.isin(['Astro_Epend','Oligo','OPC','Micro_PVM','IN_MGE'])
    #anno.loc[sel,'celltype']=anno.loc[sel,'subtype']
    
    anno=anno.loc[~anno.BrainRegion.isna()]
    anno=anno.loc[~anno.Dx_donor.isna()]

    anno=anno.loc[anno[level].eq(cell)]
    input_path=f'class_subset_{cell}_raw.h5ad'
    adata=sc.read_h5ad(input_path)
    sh=np.intersect1d(adata.obs_names,anno.index)
    adata=adata[sh,]
    for col in anno.columns:
        adata.obs[col]=anno.loc[sh,col]
    
    adata.obs['DxPC1_clust3']=adata.obs['Dx_cat_donor'].map({'control':1,'early':2,'late':3})
    adata.obs['DxPC1']=adata.obs['Dx_donor']
    
    filter_genes(adata,includeX=True)
    adata=adata[:,adata.var.filtered]
    process_adata(adata)
    
    M=adata.X.mean(axis=0).A1 #colmeans
    df_refavg=pd.DataFrame({'mean':M,
                            'std':np.sqrt(adata.X.power(2).mean(axis=0).A1 - M**2)})
    df_refavg.index=adata.var.index.values
    df_refavg.to_csv(f'ML_cellstate3/{cell}/ref_avg_expression.csv')

    
    adata.obs['gp']=adata.obs.celltype+'_'+adata.obs.BrainRegion
    
    f=pd.crosstab(adata.obs.gp,
                  adata.obs.Dx_cat_donor)
    cellfrac=f.min(axis=1)
    cellfrac=cellfrac.loc[cellfrac>0]
    n_class=f.shape[1]
    
    target=round(cellfrac.sum()*0.9*n_class)
    target=np.min((target,10000*n_class))
    if target < 2500:
        sys.exit(f'Stop low number of cells {cell}')
    fold=20
    #SEL=subset_by_cat(adata.obs.DxPC1_clust3.values,target=target, T=fold)
    SEL=subset_by_cat(adata.obs.Dx_cat_donor.values,
                      adata.obs.gp.values,cellfrac,target=target, T=fold)
    z=[adata.obs_names[x].values for x in SEL]
    with open(f'ML_cellstate3/{cell}/sel_cell.pkl', 'wb') as f:
        pickle.dump(z, f)
    Importance2=np.zeros(shape=[adata.var.shape[0],fold])
    Preds = np.full([adata.obs.shape[0],fold], np.nan)
    RMSE=np.zeros(fold)
    for x in range(fold):
        sel=SEL[x]
        validate=find_neg(sel,adata.obs.shape[0])
        X_train=adata.X[sel,:]
        X_valid=adata.X[validate,:]        
        
        y_train=adata.obs.DxPC1.values[sel]
        y_valid=adata.obs.DxPC1.values[validate]
        # reserve 10% cells for early stop
        rng = np.random.default_rng()
        insplit = rng.uniform(0,1,sel.shape[0])<0.9
        X_train0=X_train[insplit,:]
        X_train1=X_train[~insplit,:]
        y_train0=y_train[insplit]
        y_train1=y_train[~insplit]        
        model=lgb_regress(X_train0,X_train1,y_train0,y_train1)
        y_pred = model.predict(X_valid)
        Preds[validate,x]=y_pred
        Importance2[:,x] = model.feature_importance(importance_type="gain")
        RMSE[x] = np.sqrt(mean_squared_error(y_valid, y_pred))
        with open(f'ML_cellstate3/{cell}/lgb_regression/model_fold{x}.pkl', 'wb') as f:
            pickle.dump(model, f)    
    
    np.savetxt(f'ML_cellstate3/{cell}/lgb_regression/RMSE',RMSE)
    adata.obs['DxPC1_predict']=np.nanmean(Preds,axis=1)
    f=np.average(Importance2, weights=1/RMSE,axis=1)
    Importance2=pd.DataFrame({'importance':f,'gene':adata.var.index.values})
    Importance2.set_index('gene',inplace=True)
    Importance2.to_csv(f'ML_cellstate3/{cell}/lgb_regression/Importance.csv')
    zdf=adata.obs[['DxPC1','DxPC1_predict']]
    zdf.to_csv(f'ML_cellstate3/{cell}/lgb_predict.csv')

