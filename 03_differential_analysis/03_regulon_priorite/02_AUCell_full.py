# conda activate /sc/arion/projects/roussp01a/pengfei/tools/anaconda3/envs/scenicplus
import os
import scanpy as sc
import numpy as np
import pandas as pd
from pyscenic.aucell import GeneSignature, derive_auc_threshold
from pyscenic.aucell import aucell4r
os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/GRN/files/filter')
def notsocreate_rankings(ex_mtx: pd.DataFrame,ind, seed=None) -> pd.DataFrame:
    DTYPE = "uint32"
    return (
        ex_mtx[ind]
        .rank(axis=1, ascending=False, method="first", na_option="bottom")
        .astype(DTYPE)- 1
    )
def aucell2(
    exp_mtx,
    signatures,
    ind,
    auc_threshold: float = 0.05,
    noweights: bool = False,
    normalize: bool = False,
    seed=None,
    num_workers: int = 1):
    return aucell4r(
        notsocreate_rankings(exp_mtx,ind,seed),
        signatures,
        auc_threshold,
        noweights,
        normalize,
        num_workers,
    )
#_corr_p_filtered_grn.csv
for cell in 'Astro','Oligo','OPC','Micro_PVM','EN','IN':
    for prefix in ['p']:
        grn=pd.read_csv(f'{cell}_corr_{prefix}_filtered_grn.csv',index_col=0)
        with open(f'{cell}_{prefix}_tmp.gmt','w') as f:
            for xx in grn.TF.unique():
                l=f'{xx}\t{xx}_{prefix}_p\t'+'\t'.join(grn.loc[grn.TF.eq(xx)].target.values)+'\n'
                f.write(l)
        signatures = GeneSignature.from_gmt(f'{cell}_{prefix}_tmp.gmt', field_separator='\t', gene_separator='\t')
        input_path=f'/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/RNA_subset_{cell}_raw.h5ad'
        raw=sc.read_h5ad(input_path)
        rd=raw.var.index.values.copy()
        np.random.seed(42)
        np.random.shuffle(rd)
        
        obs1=raw.obs.sample(n=10000,replace=False).index.values
        mat1=pd.DataFrame(raw[obs1,].X.todense(),index=obs1,columns=raw.var.index)
        percentiles = derive_auc_threshold(mat1)
        step=int(1e4)
        aucs_mtx=[]
        for st in range(0,raw.obs.shape[0],step):
            ed=min(st+step,raw.obs.shape[0])
            obs1=raw.obs_names[st:ed]
            mat1=pd.DataFrame(raw[obs1,].X.todense(),index=obs1,columns=raw.var.index)
            aucs_mtx1 = aucell2(mat1, signatures,rd, auc_threshold=percentiles[0.01], num_workers=16)
            aucs_mtx.append(aucs_mtx1)
        aucs_mtx=pd.concat(aucs_mtx)
        aucs_mtx.to_csv(f'{cell}_{prefix}_AUCcell.csv.gz',compression='gzip')
        del aucs_mtx
        del raw

