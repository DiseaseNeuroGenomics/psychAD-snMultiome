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
from sklearn.metrics import adjusted_rand_score as ari
sys.path.append("/sc/arion/projects/roussp01a/pengfei/tools/scanpy_myextra/")
import myutsc
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

os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/')
sc.settings.verbosity = 2             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc._settings.ScanpyConfig.n_jobs=32
sc.settings.set_figure_params(dpi=80, facecolor='white')

nnb=30
ress=[0.1,0.3]
Ma2022_markers = pd.read_csv('Ma_2022_DLPFC_markers_clean.csv')
Ma2022_markers['Pericyte']=Ma2022_markers['PC']
Ma2022_markers.drop(columns='PC',inplace=True)
atac_folder='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save2/'

   
ress=[0.1,0.2,0.5,1]
cluster_names=['2leiden'+str(res).replace('.', '_') for res in ress]
npc=10 
ARI=[]
CR=[]
for cell in 'Vascular','EN','IN','Immune','OPC','Astro','Oligo':
    rna=sc.read_h5ad(f'RNA_2nd_subset_{cell}.h5ad')
    fh=h5py.File(atac_folder+'/fine_cluster/'+cell+'/meta_PCA.h5')
    atac_PCA=fh['PCA'][:][:(npc-1)]
    atac_HPCA=fh['Harmony_PCA'][:(npc-1)]
    atac_meta=pd.read_csv(atac_folder+'/fine_cluster/'+cell+'/meta.csv',index_col=0)
    adata=ad.AnnData(obs=atac_meta,var=rna.var.iloc[:1],
                  X=sparse.csr_matrix((atac_meta.shape[0], 1)))
    adata.obsm['X_pca']=atac_PCA.T
    adata.obsm['X_pca_harmony']=atac_HPCA.T

    #adata=adata[adata.obs[f'LGlevel1_{cell}']>0.4,]
    sc.pp.neighbors(adata, use_rep = "X_pca_harmony",n_neighbors=30)
    sc.tl.umap(adata,maxiter=50)
    for i in range(len(ress)):
        res=ress[i]
        cluster_name=cluster_names[i]
        sc.tl.leiden(adata, resolution=res,key_added=cluster_name)

    adata.obs.newID=adata.obs.newID.astype('category')
    adata.obs['log10_nFrags']=np.log10(adata.obs['nFrags'])
    ################
    sel_obs=np.intersect1d(rna.obs.index,adata.obs.index)
    ad1=rna[sel_obs,].copy()
    ad2=adata[sel_obs,].copy()
    mdata = mu.MuData({'rna': ad1, 'atac': ad2})
    for i in range(len(ress)):
        res=ress[i]
        cluster_name=cluster_names[i]
        sns.set_style("white")
        z=compare_cluster(mdata.obs,f'rna:{cluster_name}',f'atac:{cluster_name}')
        z=compare_cluster(mdata.obs,f'atac:{cluster_name}',f'rna:{cluster_name}')
        mu.pl.embedding(mdata, basis="atac:X_umap", color=[f'rna:{cluster_name}',f'atac:{cluster_name}'])
        mu.pl.embedding(mdata, basis="rna:X_umap", color=[f'rna:{cluster_name}',f'atac:{cluster_name}'])
    print('WNN')
    sc.pp.neighbors(mdata['rna'], use_rep = "X_pca_harmony",n_neighbors=nnb)
    sc.pp.neighbors(mdata['atac'], use_rep = "X_pca_harmony",n_neighbors=nnb)
    mu.pp.neighbors(mdata, key_added='wnn')
    mu.tl.umap(mdata, neighbors_key='wnn', random_state=10)
    mdata.obsm["X_wnn_umap"] = mdata.obsm["X_umap"]
    #mu.pl.embedding(mdata, basis="X_wnn_umap", color=["rna:predicted_Ma", "rna:predicted_NPSADsubtype"])
    mu.pl.embedding(mdata, basis="X_wnn_umap", color=['atac:2leiden0_1','rna:2leiden0_1'])
    mu.pl.embedding(mdata, basis="X_wnn_umap", color=['atac:2leiden0_2','rna:2leiden0_2'])
    for i in range(len(ress)):
        res=ress[i]
        cluster_name=cluster_names[i]
        sc.tl.leiden(mdata, resolution=res,neighbors_key='wnn',key_added='wnn'+cluster_name)
        sns.set_style("white")

        rna.obs[f'wnn{cluster_name}']=np.nan
        rna.obs.loc[mdata.obs.index,f'wnn{cluster_name}']=mdata.obs[f'wnn{cluster_name}']
        
        query=rna[rna.obs[f'wnn{cluster_name}'].isna(),]
        ref=rna[~rna.obs[f'wnn{cluster_name}'].isna(),]
        sc.pp.neighbors(ref, use_rep = "X_pca_harmony",n_neighbors=nnb)
        pred=myutsc.injest_from_same_cohort(ref,query,use_rep='X_pca_harmony',n_pcs=npc,
                            epsilon=0.1,obs='wnn{cluster_name}')
        rna.obs.loc[query.obs_names,f'wnn{cluster_name}']=pred.obs[f'pred_wnn{cluster_name}']

        adata.obs[f'wnn{cluster_name}']=np.nan
        adata.obs.loc[mdata.obs.index,f'wnn{cluster_name}']=mdata.obs[f'wnn{cluster_name}']
        query=adata[adata.obs[f'wnn{cluster_name}'].isna(),]
        ref=adata[~adata.obs[f'wnn{cluster_name}'].isna(),]
        sc.pp.neighbors(ref, use_rep = "X_pca_harmony",n_neighbors=nnb)
        pred=myutsc.injest_from_same_cohort(ref,query,use_rep='X_pca_harmony',n_pcs=npc,
                            epsilon=0.1,obs=f'wnn{cluster_name}')
        adata.obs.loc[query.obs_names,f'wnn{cluster_name}']=pred.obs[f'pred_wnn{cluster_name}']
    
    adata.write(atac_folder+'/'+cell+'_qc01.h5ad')
    rna.write(f'RNA_2nd_subset_{cell}_wnn.h5ad')
    mdata.write(cell+'_mudata_2.h5mu')    
    # select 20% of cells for validation
    i=1
    res=ress[i]
    cluster_name=cluster_names[i]
    
    adata.obs['inRNA']=adata.obs.index.isin(rna.obs.index)
    adata.obs['ref']=adata.obs.inRNA*(np.random.uniform(size=adata.obs.shape[0])>0.2)
    adata.obs['validation']=adata.obs.inRNA*(~adata.obs.ref)
    adata.obs[f'wnn{cluster_name}']=np.nan
    adata.obs.loc[adata.obs.ref,f'wnn{cluster_name}']=mdata.obs.loc[adata.obs.loc[adata.obs.ref].obs_names,f'wnn{cluster_name}']
    
    query=adata[adata.obs[f'wnn{cluster_name}'].isna(),]
    ref=adata[~adata.obs[f'wnn{cluster_name}'].isna(),]
    sc.pp.neighbors(ref, use_rep = "X_pca_harmony",n_neighbors=nnb)
    pred=myutsc.injest_from_same_cohort(ref,query,use_rep='X_pca_harmony',n_pcs=npc,
                            epsilon=0.1,obs=f'wnn{cluster_name}')
    confusion_matrix=pd.crosstab(mdata.obs.loc[pred[pred.obs.validation,].obs_names,f'wnn{cluster_name}'],
             pred.obs.loc[pred[pred.obs.validation,].obs_names,f'pred_wnn{cluster_name}'])
    ind=np.intersect1d(confusion_matrix.index.values,confusion_matrix.columns.values)
    correct_ratio=np.sum(confusion_matrix.loc[ind,ind].values.diagonal())/confusion_matrix.sum().sum()
    ariv=ari(mdata.obs.loc[pred[pred.obs.validation,].obs_names,f'wnn{cluster_name}'],
             pred.obs.loc[pred[pred.obs.validation,].obs_names,f'pred_wnn{cluster_name}'])

    ARI.append(ariv)
    CR.append(correct_ratio)   

met={'ARI':ARI,'correct_ratio':CR,'class':['Vascular','EN','IN','Immune','OPC','Astro','Oligo']}
met=pd.DataFrame(met)
met.to_csv('class_ATAC_ARI_correct_ratio.csv')

