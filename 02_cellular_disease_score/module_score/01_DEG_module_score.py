#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 19:28:03 2025

@author: pengfeidong
"""

import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

def read_anno(path,level,cell):
    anno=pd.read_csv(path,index_col=0)
    anno=anno.loc[anno.inRNA]
    anno=anno.loc[anno[level].eq(cell)]
    return anno

def filter_genes(adata,includeX=True,min_cell_percent=0.001):
    min_cell=int(adata.obs.shape[0]*min_cell_percent)
    sc.pp.filter_genes(adata, min_cells=min_cell)
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
    adata=adata[:,adata.var.filtered]
    return adata

def normalize_adata(adata):
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

def add_gene_score(adata,deg,prefix,topN=50):
    deg=deg.loc[deg.ID.isin(adata.var_names)]
    deg = deg.assign(abst=deg['t'].abs()).sort_values('abst', ascending=False)
    up=deg.loc[deg.logFC>0].ID.values[:topN]
    down=deg.loc[deg.logFC<0].ID.values[:topN]
    sc.tl.score_genes(adata,up,score_name=f'{prefix[0]}',ctrl_size=2*topN)
    sc.tl.score_genes(adata,down,score_name=f'{prefix[1]}',ctrl_size=2*topN)


def read_lgb_anno(fh,adata):
    lgb=pd.read_csv(fh,index_col=0)
    lgb=lgb.loc[~lgb.DxPC1_predict2.isna()]
    lgb=lgb.loc[~lgb.DxPC1_clust3_predict2.isna()]
    adata=adata[lgb.index,]
    for x in lgb.columns:
        adata.obs[x]=lgb[x]
    adata.obs['DxPC1_clust3_predict2']=adata.obs['DxPC1_clust3_predict2'].astype('int')
    adata.obs['DxPC1_clust3']=adata.obs['DxPC1_clust3'].astype('int')
    adata.obs['cat']=adata.obs['DxPC1_clust3'].astype('str')+'_'+adata.obs['DxPC1_clust3_predict2'].astype('str')
    return adata


if __name__ == '__main__':
    cell=sys.argv[1] # 
    level='class'
    root='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files'
    os.chdir(root)

    anno=read_anno('./cell_anno_05092025_fixed3.csv',level,cell)
    anno=anno.loc[~anno.subtype.eq('ependA')]
    DxPC1=anno.loc[(~anno.SubID.duplicated())*(~anno.Dx_donor.isna())].Dx_donor.values
    input_path=f'{level}_subset_{cell}_raw.h5ad'
    adata=sc.read_h5ad(input_path)
    sh=np.intersect1d(anno.index,adata.obs_names)
    adata=adata[sh,]
    for xx in anno.columns:
        adata.obs[xx]=anno.loc[sh,xx]
    
    adata=read_lgb_anno(f'ML_cellstate3/{cell}/lgb_predict_redo.csv',adata)
    adata=filter_genes(adata,includeX=True)
    normalize_adata(adata)
    # the DEG file is avaliable at the supplementary files of https://doi.org/10.1101/2024.10.31.24316513
    DEG=pd.read_csv('/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/analysis/results/MSSM/MSSM_AD__controls/SubID/class/topTable.tsv.gz',index_col=0)
    
    if(cell=='Micro_PVM'):
        deg=DEG.loc[DEG.assay.eq('Immune')]
    else:
        deg=DEG.loc[DEG.assay.eq(cell)]
    prefix=['PsychAD_up','PsychAD_down']
    add_gene_score(adata,deg,prefix,topN=50)

    sel_col=['class', 'subclass', 'subtype','level1_5', 'DxPC1', 'DxPC1_clust3',
             'BrainRegion', 'SubID', 'Age',
             'DxPC1_clust3_predict2', 'DxPC1_predict','DxPC1_predict2',
             'PsychAD_up','PsychAD_down',
             'mito_genes','mito_ribo', 'ribo_genes','n_genes_by_counts',
             'total_counts','pct_counts_mt']
    adata.obs[sel_col].to_csv(f'ML_cellstate3/{cell}/DxPC1_predict2_DEG_module_score.csv.gz',
                              compression='gzip')
    #####
