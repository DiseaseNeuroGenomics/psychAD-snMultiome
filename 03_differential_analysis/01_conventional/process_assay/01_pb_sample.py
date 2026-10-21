
import numpy as np
import scanpy as sc
import pandas as pd
import os
from scipy import sparse
import argparse
parser = argparse.ArgumentParser(description='h5ad to pseudobulk')

parser.add_argument('--input', help='Input h5ad file, X is raw counts', required=True)
parser.add_argument('--cellanno', help='cell type anno csv file including class subclass subtype', required=True)
parser.add_argument('--outDir', help='Output dir', required=True)
parser.add_argument('--outPrefix', help='output prefix', required=True)
parser.add_argument('--sampleCol', help='sampleID column', default='order')
parser.add_argument('--ATAC', action='store_true', default=False) 
args = parser.parse_args()
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
##
def get_combine_col(obs,cellanno,sample,newcol):
    obs[newcol]=obs[cellanno].astype('str')+'#'+obs[sample].astype('str')

def get_pseduobulk(adata,colname,mean=False,X=True):
    z=pd.get_dummies(adata.obs[colname])
    columns=z.columns
    mat=pd.get_dummies(z).to_numpy()
    ncell=mat.sum(axis=0)
    del z 
    if mean:
        mat=mat/mat.sum(axis=0)
    mat=sparse.csr_matrix(mat.T)
    if X:
        psb=mat.dot(adata.X)
    else:
        psb=mat.dot(adata.counts)
    psb=pd.DataFrame(psb.T.todense(),columns=columns,index=adata.var.index)   
    return psb,ncell
def get_by_cell(pseduobulk,ncells,newid):
    z=pd.DataFrame([a.replace(newid+'_','').split('#') for a in pseduobulk.columns])
    z.columns=['Cell','Sample']    
    pbs={}
    ncell={}
    for x in z.Cell.unique():
        pb1=pseduobulk[pseduobulk.columns[z.Cell.eq(x)]]
        pb1.columns=z.Sample.loc[z.Cell.eq(x)]
        pbs[x]=pb1
        ncell[x]=pd.DataFrame({'sample':z.Sample.loc[z.Cell.eq(x)],
                               'ncell':ncells[z.Cell.eq(x)]})
        ncell[x].set_index('sample',inplace=True)
    return pbs,ncell
####
raw_path=args.input
wnn_path=args.cellanno
outputdir=args.outDir
output_prefix=args.outPrefix
adata=sc.read_h5ad(raw_path)
wnn_cluster=pd.read_csv(wnn_path,index_col=0)
if args.ATAC:
    wnn_cluser=wnn_cluster.loc[wnn_cluster.inATAC]
else:
    wnn_cluster=wnn_cluster.loc[wnn_cluster.inRNA] 
###########################
sel=np.intersect1d(adata.obs.index,wnn_cluster.index)
adata=adata[sel,]

adata.obs['subclass']=wnn_cluster.loc[sel,'subclass']
adata.obs['subtype']=wnn_cluster.loc[sel,'subtype']
adata.obs['class']=wnn_cluster.loc[sel,'class']
adata.obs['sample']=adata.obs[args.sampleCol].astype('str') #order is sample ID
################

for level in 'class','subclass','subtype':
    get_combine_col(adata.obs,f'{level}','sample',f'{level}person')
    pb_subclass,ncells=get_pseduobulk(adata,[f'{level}person'])    
    pbs,ncells=get_by_cell(pb_subclass,ncells, f'{level}person')

    if os.path.exists(f'{outputdir}/{output_prefix}_{level}_pb_persample.h5'):
        os.remove(f'{outputdir}/{output_prefix}_{level}_pb_persample.h5')
    with pd.HDFStore(f'{outputdir}/{output_prefix}_{level}_pb_persample.h5') as store:
    # Write DataFrames to the file
        for x in ncells.keys():
            store.put(x+'/count_matrix',pbs[x])
            store.put(x+'/ncell',ncells[x])

    logging.info(f"{output_prefix} pseudobulk at {outputdir}/{output_prefix}_{level}_pb_persample.h5")

