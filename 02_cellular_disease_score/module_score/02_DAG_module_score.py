mport scanpy as sc
import pandas as pd
import os

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
    
DEG=pd.read_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/class_genescore/all/DxPC1.csv',index_col=0)
anno=pd.read_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno_05092025_fixed3.csv',index_col=0)
anno=anno.loc[anno.inATAC]
os.chdir("/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save3")
adata=sc.read_h5ad('mergedgene.h5ad')
normalize_adata(adata)
adata=adata[adata.obs_names.isin(anno.index),]
for x in 'class','subclass','subtype':
    adata.obs[x]=anno.loc[adata.obs_names,x]

dfs=[]
for cell in ['Astro','EN','IN','Micro_PVM','Oligo','OPC','Mural']:
    ad1=adata[adata.obs['class'].eq(cell),]
    prefix=['DAG_up','DAG_down']
    deg=DEG.loc[DEG.assay.eq(cell)]
        
    add_gene_score(ad1,deg,prefix,topN=100)
    dfs.append(ad1.obs[prefix+'class'])
df=pd.concat(dfs,axis=0)
df.to_csv('gene_DAG_score.csv.gz',compression='gzip')
