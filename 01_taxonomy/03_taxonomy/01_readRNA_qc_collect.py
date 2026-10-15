# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 16:12:18 2023

@author: dongp
"""

import pandas as pd
import os
import scanpy as sc
import anndata as ad
import numpy as np
import pegasus as pg
import rpy2.robjects as robjects 
from rpy2.robjects import pandas2ri
from gtfparse import read_gtf
pandas2ri.activate()


def read_sample(sample):
    os.chdir('/'.join([inputDir,sample,'outs']))
    if os.path.isfile('soupX.h5'):
        ad1=sc.read_10x_h5('soupX.h5')
        soupX=True
    else:
        ad1=sc.read_10x_h5('filtered_feature_bc_matrix.h5')
        ad1.var_names_make_unique()
        soupX=False
    robjects.r['load']('barcodeinfo.Rdata')
    meta=robjects.r('barcodeinfo')
    meta=meta.drop(['gex_barcode','atac_barcode','is_cell','excluded_reason'],axis=1)
    
    meta['soupX']=soupX
    allInfo0=allInfo.loc[sample]
    for x in allInfo0.keys():
        meta[x]=allInfo0[x]
    ad1.obs=meta.loc[ad1.obs.index]
    ad1.obs_names=allInfo0['order']+'#'+ad1.obs_names
    return ad1
    
inputDir='/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/files'
outDir='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files'
allInfo_path='/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/files/metaInfo_11102023.csv'
robjects.r['load'](allInfo_path)

allInfo=pd.read_csv(allInfo_path,index_col=0)
allInfo_extra=pd.read_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/files/allInfo_extra.csv')
allInfo=allInfo.loc[allInfo['Estimated_number_of_cells'].between(200,20000,inclusive='left')]
samples=np.r_[allInfo['newID'].to_numpy(),np.unique(allInfo_extra['newID'])]
adatas=[read_sample(x) for x in samples]
adata = ad.concat(adatas, index_unique=None)


genes=read_gtf('/hpc/packages/minerva-centos7/cellranger-arc/2.0.0/refdata-cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz')
genes = genes[genes["feature"] == "gene"]
genes.loc[genes['gene_name'].duplicated(),'gene_name']=genes.loc[genes['gene_name'].duplicated(),'gene_name']+'-1'
genes=genes.set_index('gene_name')
adata.var['gene_type']=genes.loc[adata.var.index,'gene_type']
#mt ribo PC
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
adata.var['ribo'] = adata.var_names.str.startswith('RPL') | adata.var_names.str.startswith('RPS')
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
adata.var['protein_coding'] = adata.var.gene_type.eq('protein_coding')
#
mitocarta = pd.read_csv('/sc/arion/projects/CommonMind/leed62/ref/MitoCarta/Human.MitoCarta3.0.csv')
adata.var['mitocarta'] = adata.var_names.isin(mitocarta.Symbol)

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)
###################################
#pegasus
for signature in ["cell_cycle_human", "gender_human", "mitochondrial_genes_human", "ribosomal_genes_human","apoptosis_human"]:
    pg.calc_signature_score(adata,signature)
adata=adata.to_anndata()
#############################################
if not os.path.isdir(outDir):
    os.makedirs(outDir)  
os.chdir(outDir)
adata.write_h5ad('RNA_all_raw.h5ad',compression='gzip') 
