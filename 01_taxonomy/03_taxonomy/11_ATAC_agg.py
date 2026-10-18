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
import pegasus as pg
import sys
import subprocess
import h5py
import anndata as ad
from scipy import sparse
import muon as mu
from muon import atac as ac
sys.path.append("/sc/arion/projects/roussp01a/pengfei/tools/scanpy_myextra/")
import myutsc
os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/')
sc.settings.verbosity = 2             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc._settings.ScanpyConfig.n_jobs=8
sc.settings.set_figure_params(dpi=80, facecolor='white')

nnb=30
npc=10
atac_path='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save2/fine_cluster2/'
fhs=glob.glob(f'{atac_path}/*/meta_PCA.h5')
cells=[os.path.basename(os.path.dirname(x)) for x in fhs]
cells.sort(reverse=True)
cells

ARI=[]
CR=[]
dfs=[]
for i in range(len(cells)):
    cell=cells[i]

    RNA=sc.read_h5ad(f'RNA_3nd_subset_{cell}.h5ad').obs
    fh=h5py.File(atac_path+cell+'/meta_PCA.h5')
    atac_PCA=fh['PCA'][:][:(npc-1)]
    atac_HPCA=fh['Harmony_PCA'][:(npc-1)]
    atac_meta=pd.read_csv(atac_path+cell+'/meta.csv',index_col=0)
    adata=ad.AnnData(obs=atac_meta,
                  X=sparse.csr_matrix((atac_meta.shape[0], 1)))
    adata.obsm['X_pca']=atac_PCA.T
    adata.obsm['X_pca_harmony']=atac_HPCA.T
    rna=RNA.loc[RNA.index.isin(adata.obs.index)]
    rna['mask']=np.random.uniform(size=rna.shape[0])>0.2 # reserve for validation
    adata.obs['ref']=np.nan
    adata.obs.loc[rna.index,'ref']=rna.subtype
    adata.obs['validation']=adata.obs.index.isin(rna.loc[~rna['mask']].index)
    adata.obs['subtype']=adata.obs['ref']
    adata.obs.loc[adata.obs.validation,'subtype']=np.nan
    ad_query=adata[adata.obs.subtype.isna(),]
    ad_ref=adata[~adata.obs.subtype.isna(),]
    sc.pp.neighbors(ad_ref, use_rep = "X_pca_harmony",n_neighbors=nnb)
    pred=myutsc.injest_from_same_cohort(ad_ref,ad_query,use_rep='X_pca_harmony',n_pcs=npc,
                            epsilon=0.1,obs='subtype')
    confusion_matrix=pd.crosstab(pred.obs.loc[pred.obs.validation].ref,pred.obs.loc[pred.obs.validation].pred_subtype)
    from sklearn.metrics import adjusted_rand_score as ari
    ariv=ari(pred.obs.loc[pred.obs.validation].ref,
            pred.obs.loc[pred.obs.validation].pred_subtype)
    ind=np.intersect1d(confusion_matrix.index.values,confusion_matrix.columns.values)
    correct_ratio=np.sum(confusion_matrix.loc[ind,ind].values.diagonal())/confusion_matrix.sum().sum()
    ARI.append(ariv)
    CR.append(correct_ratio)
    dfa=pred.obs[['ref','pred_subtype']]
    dfb=ad_ref.obs[['ref','subtype']]
    dfa=dfa.rename(columns={'pred_subtype':'subtype'})
    dfa.loc[~dfa.ref.isna(),'subtype']=dfa.loc[~dfa.ref.isna(),'ref']
    dfa
    df=pd.concat((dfa,dfb))
    dfs.append(df)
ATAC_anno=pd.concat(dfs)   
