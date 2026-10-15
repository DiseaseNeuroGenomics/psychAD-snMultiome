import os
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib.colors as mtpltcl
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm
def sel_celltype_rm_criteria(obs1,mod,col):
    obs1=obs1.sort_values(col)
    obs1['rank']=np.arange(obs1.shape[0])//200
    df1=obs1[['rank',col,'same_withTO']].groupby('rank').agg('mean')
    plt.plot(df1.same_withTO)
    plt.title(f'{celltype} {mod} {col}')
    st=np.rint(df1.shape[0]*0.2).astype('int')
    ed=np.rint(df1.shape[0]*0.8).astype('int')
    consistency_cutoff=df1[st:ed].same_withTO.mean()-3*df1[st:ed].same_withTO.std()

    lowess = sm.nonparametric.lowess
    z = lowess(df1.same_withTO, df1.index.values,frac=0.05)
    plt.plot(z[:,0],z[:,1])
    plt.show()
    sel=np.where(z[:,1]>consistency_cutoff)[0][0]
    if sel>0:
        rm0={"celltype":celltype,"col":col,"cutoff":df1.loc[sel,col]}
        return rm0
    return None

os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files')
RNA_prefix="RNA_all_initial_cutoff_ngene_500_10000_ntotal_100000_pct_mito1_n_neighbors100_Nhvg5000_NPC30"
RNA=sc.read_h5ad(f'{RNA_prefix}_qc_01.h5ad')
adata=sc.read_h5ad('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save2/qc_00.h5ad')
sel_obs=np.intersect1d(RNA.obs.index,adata.obs.index)
#
RNA.obs['complexity']=np.nan
RNA.obs['mlog10_total_counts']=-1*RNA.obs['log10_total_counts']
adata.obs['mBlacklistRatio']=-1*adata.obs.BlacklistRatio
adata.obs['mlog10_nFrags']=-1*adata.obs['log10_nFrags']
sel_obs=np.intersect1d(RNA.obs.index,adata.obs.index)
atac=adata.obs.loc[sel_obs]
rna=RNA.obs.loc[sel_obs]
rna['same_withTO']=atac['same_withTO']=rna.cell_type==atac.cell_type


RNA_rm=[]
ATAC_rm=[]
for celltype in atac.level1.unique():
    atac1=atac.loc[atac.level1==celltype,]
    rna1=rna.loc[rna.level1==celltype,]
    ## add rna complexity'
    sel=rna1.same_withTO
    reg = LinearRegression().fit(rna1.loc[sel].percent_dup.values.reshape(-1, 1),rna1.loc[sel].log10_total_counts.values)
    sel=RNA.obs.level1==celltype
    RNA.obs.loc[sel,'complexity']=RNA.obs.loc[sel].log10_total_counts.values-reg.predict(RNA.obs.loc[sel].percent_dup.values.reshape(-1, 1))
    rna1['complexity']=RNA.obs.loc[rna1.index].complexity.values
    ###
    rm=sel_celltype_rm_criteria(rna1,'RNA','complexity')
    if rm is not None:
        RNA_rm.append(rm)
        rna1=rna1.loc[rna1[rm['col']]>rm['cutoff']]
    rm=sel_celltype_rm_criteria(rna1,'RNA','mlog10_total_counts')
    if rm is not None:
        RNA_rm.append(rm)

    rm=sel_celltype_rm_criteria(atac1,'ATAC','FRiP') 
    if rm is not None:
        ATAC_rm.append(rm)
        atac1=atac1.loc[atac1[rm['col']]>rm['cutoff']]
    rm=sel_celltype_rm_criteria(atac1,'ATAC','TSSEnrichment')
    if rm is not None:
        ATAC_rm.append(rm)
        atac1=atac1.loc[atac1[rm['col']]>rm['cutoff']]
    rm=sel_celltype_rm_criteria(atac1,'ATAC','log10_nFrags')
    if rm is not None:
        ATAC_rm.append(rm)
        atac1=atac1.loc[atac1[rm['col']]>rm['cutoff']]
    rm=sel_celltype_rm_criteria(atac1,'ATAC','mBlacklistRatio') 
    if rm is not None:
        ATAC_rm.append(rm)
        atac1=atac1.loc[atac1[rm['col']]>rm['cutoff']]
    rm=sel_celltype_rm_criteria(atac1,'ATAC','mlog10_nFrags')
    if rm is not None:
        ATAC_rm.append(rm)

RNA_rm=pd.DataFrame(RNA_rm)
ATAC_rm=pd.DataFrame(ATAC_rm)
RNA_rm.to_csv('RNA_remove_cutoff.csv',index=False)
ATAC_rm.to_csv('ATAC_remove_cutoff.csv',index=False)
############# adaptive filter and remove inconsistent cells 
RNA=sc.read_h5ad(f'{RNA_prefix}_qc_01.h5ad')
adata=sc.read_h5ad('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save2/qc_00.h5ad')
raw=sc.read_h5ad('RNA_all_raw.h5ad')
RNA_exclude=[]
ATAC_exclude=[]
for x in range(RNA_rm.shape[0]):
    cell=RNA_rm.loc[x,'celltype']
    col=RNA_rm.loc[x,'col']
    cutoff=RNA_rm.loc[x,'cutoff']
    RNA_exclude.append(RNA.obs.loc[RNA.obs.cell_type.eq(cell)*(RNA.obs[col]<cutoff)].index.values)
RNA_exclude=np.unique(np.concatenate(RNA_exclude))
for x in range(ATAC_rm.shape[0]):
    cell=ATAC_rm.loc[x,'celltype']
    col=ATAC_rm.loc[x,'col']
    cutoff=ATAC_rm.loc[x,'cutoff']
    ATAC_exclude.append(adata.obs.loc[adata.obs.cell_type.eq(cell)*(adata.obs[col]<cutoff)].index.values)
ATAC_exclude=np.unique(np.concatenate(ATAC_exclude))
adata=adata[~adata.obs_names.isin(ATAC_exclude),]
RNA=RNA[~RNA.obs_names.isin(RNA_exclude),]
sel_obs=np.intersect1d(RNA.obs.index,adata.obs.index)
RNA.obs['ATAC_celltype']=np.nan
RNA.obs.loc[sel_obs,'ATAC_celltype']=adata.obs[sel_obs,'cell_type']
mask = (RNA.obs['ATAC_celltype'].notna() &
    RNA.obs['ATAC_celltype'].ne(RNA.obs['cell_type']))

rm = RNA.obs.loc[mask]
RNA=RNA[~RNA.obs_names.isin(rm),]
adata=adata[~adata.obs_names.isin(rm),]
adata.write('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save2/qc_01.h5ad')
for cell in RNA.obs.cell_type.unique():
    rna=RNA[RNA.obs.cell_type.eq(cell),]
    ad1=raw[rna.obs_names,]
    ad1.obs=rna.obs
    ad1.write(f'RNA_2nd_subset_{cell}_raw.h5ad')

