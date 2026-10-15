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
sys.path.append("/sc/arion/projects/roussp01a/pengfei/tools/scanpy_myextra/")
import myutsc
import h5py
import muon as mu
import anndata as ad
from scipy.sparse import csr_matrix
def compare_cluster(obs_df,refc,clusterc):
    import seaborn as sns
    from sklearn.metrics import adjusted_rand_score as ari
    import matplotlib.pyplot as plt
    ariv=ari(obs_df.loc[~(pd.isnull(obs_df[refc])|pd.isnull(obs_df[clusterc]))][refc],
            obs_df.loc[~(pd.isnull(obs_df[refc])|pd.isnull(obs_df[clusterc]))][clusterc])
    confusion_matrix=pd.crosstab(obs_df.loc[~(pd.isnull(obs_df[refc])|pd.isnull(obs_df[clusterc]))][clusterc],
                                obs_df.loc[~(pd.isnull(obs_df[refc])|pd.isnull(obs_df[clusterc]))][refc]).values
    purity=(np.sum(np.max(confusion_matrix, axis=1)) / np.sum(confusion_matrix))
    
    ncluster=obs_df.loc[~pd.isnull(obs_df[refc])].groupby(clusterc).size().reset_index(name="ncluster").set_index(clusterc)
    df_frac=obs_df.groupby([clusterc,refc]).size().reset_index(name="n").set_index(clusterc).join(ncluster).reset_index()
    df_frac = df_frac.assign(frac = lambda x: x.n / x.ncluster)
    ### plot heatmap
    mycolor=sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, as_cmap=True)
    df_wide = df_frac.set_index(clusterc).pivot(columns=refc, values="frac").T
    sns.set_style("white")
    annot_labels = np.empty_like(df_wide, dtype=np.str_)
    annot_mask = df_wide > 0.1
    labels=df_wide.to_numpy()[annot_mask]
    labels=np.array(["{:0.0f}".format(x*100) for x in labels])
    annot_labels=annot_labels.astype(labels.dtype)
    annot_labels[annot_mask] =  labels
    plt.figure(figsize = (5+0.12*df_wide.shape[1],4+0.12*df_wide.shape[0]))
    sns.heatmap(df_wide, annot=annot_labels, fmt='',cmap=mycolor,edgecolor="none",linewidths=0,
                yticklabels=True,xticklabels=True).set(title=f'ari{ariv:.3f} purity{purity:.3f}')
    
    plt.show()
    return df_wide 


os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files')
sc.settings.verbosity = 2             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc._settings.ScanpyConfig.n_jobs=1

RNA_prefix="RNA_all_initial_cutoff_ngene_500_10000_ntotal_100000_pct_mito1_n_neighbors100_Nhvg5000_NPC30"
RNA=sc.read_h5ad(f'{RNA_prefix}_qc_01.h5ad')
## ATAC
fh=h5py.File('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save2/initial_meta_PCA.h5')
atac_PCA=fh['PCA'][:]
atac_HPCA=fh['Harmony_PCA'][:]
atac_meta=pd.read_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save2/initial_meta.csv',index_col=0)
adata=ad.AnnData(obs=atac_meta,var=RNA.var.iloc[:1],# fake ATAC-seq h5ad var and X
                  X=csr_matrix((atac_meta.shape[0], 1)))
adata.obsm['X_pca']=atac_PCA.T
adata.obsm['X_pca_harmony']=atac_HPCA.T

sc.pp.neighbors(adata, use_rep = "X_pca_harmony",n_neighbors=100)
sc.tl.umap(adata,maxiter=50)
sc.tl.leiden(adata, resolution=0.1,key_added='leiden0_1')
sc.tl.leiden(adata, resolution=0.2,key_added='leiden0_2')
adata.obs.newID=adata.obs.newID.astype('category')

####
sel_obs=np.intersect1d(RNA.obs.index,adata.obs.index)
atac=adata.obs.loc[sel_obs]
rna=RNA.obs.loc[sel_obs]
atac['rna_level1']=rna.cell_type
df_wide=compare_cluster(atac,'rna_level1','leiden0_1')
anno_map=dict()
anno_map['0']='Oligo'
anno_map['1']='EN'
anno_map['2']='Astro'
anno_map['3']='Immune'
anno_map['4']='IN'
anno_map['5']='IN'
anno_map['6']='OPC'
anno_map['7']='Vascular'
anno_map['8']='Immune'
adata.obs['cell_type']=adata.obs.leiden0_1.map(anno_map)
adata.write('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save2/qc_00.h5ad')

### WNN 
ad1=RNA[sel_obs,].copy()
ad2=adata[sel_obs,].copy()
mdata = mu.MuData({'rna': ad1, 'atac': ad2})
sc.pp.neighbors(mdata['rna'], use_rep = "X_pca_harmony",n_neighbors=50)
sc.pp.neighbors(mdata['atac'], use_rep = "X_pca_harmony",n_neighbors=50)
mu.pp.neighbors(mdata, key_added='wnn')
mu.tl.umap(mdata, neighbors_key='wnn', random_state=10)
mdata.obsm["X_wnn_umap"] = mdata.obsm["X_umap"]
sc.tl.leiden(mdata, resolution=0.1,neighbors_key='wnn',key_added='wnnleiden0_1')
sc.tl.leiden(mdata, resolution=0.2,neighbors_key='wnn',key_added='wnnleiden0_2')
df_wide=compare_cluster(mdata.obs,'rna:rna_level1','level1')
anno_map=dict()
anno_map['0']='Oligo'
anno_map['1']='EN'
anno_map['2']='Astro'
anno_map['3']='Immune'
anno_map['4']='IN'
anno_map['5']='OPC'
anno_map['6']='IN'
anno_map['7']='EN'
anno_map['8']='Vascular'
anno_map['9']='Vascular'
anno_map['10']='Immune'
adata.obs['cell_type']=adata.obs.leiden0_1.map(anno_map)
mdata.write("mudata_initial.h5mu")
