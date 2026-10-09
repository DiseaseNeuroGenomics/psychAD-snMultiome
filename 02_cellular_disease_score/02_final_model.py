#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 18:00:07 2026

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

def lgb_classify(X_train,X_valid,y_train,y_valid):
    train_data = lgb.Dataset(X_train, label=y_train)
    valid_data = lgb.Dataset(X_valid, label=y_valid, reference=train_data)
    params = {
    "num_threads":36,
    "objective": "multiclass",
    "metric": "multi_logloss",  # Multi-class log loss
    "num_class": 3,   # Number of classes
    "boosting_type": "gbdt",
    "num_leaves": 31,
    "learning_rate": 0.05,
    "feature_fraction": 0.8}
    model = lgb.train(
    params,
    train_data,
    num_boost_round=1000,
    valid_sets=[train_data, valid_data],  # Validation datasets
    valid_names=["train", "valid"],      # Names for datasets
    callbacks=[lgb.early_stopping(stopping_rounds=50), lgb.log_evaluation(50)])
    return model

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
def de_scale(qu, ref):
    qu_mean = np.nanmean(qu)
    qu_std = np.nanstd(qu)
    ref_mean = np.nanmean(ref)
    ref_std = np.nanstd(ref)
    z = (qu - qu_mean) / qu_std * ref_std + ref_mean
    return z

if __name__ == '__main__':
    cell=sys.argv[1]
    level='class'
    os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files')
    os.makedirs(f'ML_cellstate3/{cell}/lgb_classify_final',exist_ok =True)
    os.makedirs(f'ML_cellstate3/{cell}/lgb_regression_final',exist_ok =True)
    
    anno=pd.read_csv('./cell_anno_05092025_fixed3.csv',index_col=0)
    anno=anno.loc[anno.inRNA]
    anno=anno.loc[~anno.subtype.eq('ependA')]
    anno['celltype']=anno['subtype']
    
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
    #########
    zdf=pd.read_csv(f'ML_cellstate3/{cell}/lgb_predict.csv',index_col=0)
    zdf=zdf.loc[~zdf.DxPC1_predict.isna()]
    zdf['diff']=zdf['DxPC1']-de_scale(zdf.DxPC1_predict, zdf.DxPC1)
    cutoff=np.min((zdf.DxPC1.std(),zdf['diff'].std()))
    zdf['consistent']=zdf['diff'].abs()<cutoff
    adata=adata[zdf.index,]
    for x in ['DxPC1_predict','diff','consistent']:
        adata.obs[x]=zdf[x].values
        
    filter_genes(adata,includeX=True)
    adata=adata[:,adata.var.filtered]
    process_adata(adata)

    zdf=pd.DataFrame({'mean':np.array(adata.X.mean(axis=0)).flatten()})
    zdf.index=adata.var.index.values
    zdf.to_csv(f'ML_cellstate3/{cell}/gene_avg.csv')

    adata.obs['gp']=adata.obs.celltype+'_'+adata.obs.BrainRegion
    
    f=pd.crosstab(adata.obs.loc[adata.obs.consistent].gp,
                  adata.obs.loc[adata.obs.consistent].Dx_cat_donor)
    
    cellfrac=f.min(axis=1)
    cellfrac=cellfrac.loc[cellfrac>0]
    n_class=f.shape[1]
    
    target=round(cellfrac.sum()*0.9*n_class)
    target=np.min((target,10000*n_class))
    if target < 2500:
        sys.exit(f'Stop low number of cells {cell}')
    fold=20
    SEL=subset_by_cat(adata.obs.loc[adata.obs.consistent].Dx_cat_donor.values,
                      adata.obs.loc[adata.obs.consistent].gp.values,
                      cellfrac,target=target, T=fold)
    SEL=[np.where(adata.obs.consistent)[0][x] for x in SEL]
    z=[adata.obs_names[x].values for x in SEL]
    with open(f'ML_cellstate3/{cell}/sel_cell_final.pkl', 'wb') as f:
        pickle.dump(z, f)
    
    n_class2=adata.obs.DxPC1_clust3.unique().shape[0]
    Importance1=np.zeros(shape=[adata.var.shape[0],fold])
    Importance2=np.zeros(shape=[adata.var.shape[0],fold])
    Probs = np.full([fold,adata.obs.shape[0],n_class2], np.nan)
    Preds = np.full([adata.obs.shape[0],fold], np.nan)
    ACC=np.zeros(fold)
    ACC2=np.zeros(fold)
    RMSE=np.zeros(fold)
    RMSE2=np.zeros(fold)
    for x in range(fold):
        sel=SEL[x]
        rest=find_neg(sel,adata.obs.shape[0])
        validate=rest[adata.obs.consistent[rest]]
       
        X_train=adata.X[sel,:]
        X_valid=adata.X[validate,:]
        X_rest=adata.X[rest,:]
        y_train=adata.obs.DxPC1_clust3.values[sel]-1
        y_valid=adata.obs.DxPC1_clust3.values[validate]-1
        y_rest=adata.obs.DxPC1_clust3.values[rest]-1
       # train classify model
        rng = np.random.default_rng()
        insplit = rng.uniform(0,1,sel.shape[0])<0.9
        X_train0=X_train[insplit,:]
        X_train1=X_train[~insplit,:]
        y_train0=y_train[insplit]
        y_train1=y_train[~insplit]
        model=lgb_classify(X_train0,X_train1,y_train0,y_train1)

        y_pred_proba = model.predict(X_valid)
        y_pred_classes = np.argmax(y_pred_proba, axis=1)
        y_rest_proba = model.predict(X_rest)
        y_rest_classes = np.argmax(y_rest_proba, axis=1)
        
        Probs[x,rest,:]=y_rest_proba
        Importance1[:,x] = model.feature_importance(importance_type="gain")
        ACC[x]= accuracy_score(y_valid, y_pred_classes)
        ACC2[x]= accuracy_score(y_rest, y_rest_classes)
        with open(f'ML_cellstate3/{cell}/lgb_classify_final/model_fold{x}.pkl', 'wb') as f:
           pickle.dump(model, f)
       
        # regression model
        y_train=adata.obs.DxPC1.values[sel]
        y_valid=adata.obs.DxPC1.values[validate]
        y_rest0=adata.obs.DxPC1.values[rest]
        
        y_train0=y_train[insplit]
        y_train1=y_train[~insplit]
        model=lgb_regress(X_train0,X_train1,y_train0,y_train1)
        
        y_pred = model.predict(X_valid)
        y_rest = model.predict(X_rest)
        Preds[rest,x]=y_rest
        Importance2[:,x] = model.feature_importance(importance_type="gain")
        RMSE[x] = np.sqrt(mean_squared_error(y_valid, y_pred))
        RMSE2[x] = np.sqrt(mean_squared_error(y_rest0, y_rest))
        with open(f'ML_cellstate3/{cell}/lgb_regression_final/model_fold{x}.pkl', 'wb') as f:
            pickle.dump(model, f)

    np.savetxt(f'ML_cellstate3/{cell}/lgb_classify_final/ACC_consistent',ACC)
    np.savetxt(f'ML_cellstate3/{cell}/lgb_regression_final/RMSE_consistent',RMSE)
    np.savetxt(f'ML_cellstate3/{cell}/lgb_classify_final/ACC_all',ACC2)
    np.savetxt(f'ML_cellstate3/{cell}/lgb_regression_final/RMSE_all',RMSE2)
    adata.obs['DxPC1_clust3_predict2']=np.argmax(np.nanmean(Probs,axis=0),axis=1)
    f=np.average(Importance1, weights=ACC,axis=1)
    Importance1=pd.DataFrame({'importance':f,'gene':adata.var.index.values})
    Importance1.set_index('gene',inplace=True)
    Importance1.to_csv(f'ML_cellstate3/{cell}/lgb_classify_final/Importance.csv')

    adata.obs['DxPC1_predict2']=np.nanmean(Preds,axis=1)
    f=np.average(Importance2, weights=1/RMSE,axis=1)
    Importance2=pd.DataFrame({'importance':f,'gene':adata.var.index.values})
    Importance2.set_index('gene',inplace=True)
    Importance2.to_csv(f'ML_cellstate3/{cell}/lgb_regression_final/Importance.csv')

    zdf=adata.obs[['DxPC1_clust3','DxPC1_clust3_predict2',
                 'DxPC1','DxPC1_predict','DxPC1_predict2']]
    zdf.to_csv(f'ML_cellstate3/{cell}/lgb_predict_final.csv')
