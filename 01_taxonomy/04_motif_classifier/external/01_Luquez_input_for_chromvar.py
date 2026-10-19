import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import os
import sys
from joblib import Parallel, delayed
from pathlib import Path
import h5py
def get_pseduobulk(adata,colname,mean=False):
    z=pd.get_dummies(adata.obs[colname])    
    mat=pd.get_dummies(z).to_numpy()
    if mean:
        mat=mat/mat.sum(axis=0)
    mat=sparse.csr_matrix(mat.T)
    psb=mat.dot(adata.X)
    psb=pd.DataFrame(psb.T.todense(),columns=z.columns,index=adata.var.index)   
    return psb

def process_one_celltype_simple(celltype):
    """
    Creates meta-cells by sequentially aggregating cells of the same celltype 
    and donor_id until a target read depth (threshold) is reached. No shuffling.
    """
    # Filter metadata for the target cell type
    ct_meta = meta[meta["celltype"] == celltype]
    
    barcode_groups_local = []
    pbulk_count = 0
    
    # Get unique donors in the order they appear
    unique_donors = ct_meta["individualID"].unique()
    rng = np.random.default_rng()
    rng.shuffle(unique_donors)
    for donor in unique_donors:
        # Filter down to the specific donor
        donor_meta = ct_meta[ct_meta["individualID"] == donor]
        
        # Extract barcodes in their existing order
        barcodes = donor_meta.index.values
        
        current_pbulk_members = []
        current_depth = 0
        
        for bc in barcodes:
            current_pbulk_members.append(bc)
            current_depth += donor_meta.loc[bc, "RiP"]
            
            # Store the meta-cell once threshold is reached
            if current_depth >= threshold:
                pbulk_count += 1
                group_id = f"{celltype}_pbulk_{pbulk_count}"
                
                barcode_groups_local.append(
                    pd.DataFrame({
                        "barcode": current_pbulk_members,
                        "group": group_id
                    })
                )
                
                # Reset accumulators for the next meta-cell
                current_pbulk_members = []
                current_depth = 0
                
    # Combine all valid meta-cells for this celltype
    if len(barcode_groups_local) > 0:
        return pd.concat(barcode_groups_local, ignore_index=True)
    else:
        return pd.DataFrame(columns=["barcode", "group"])
    
def write_10X_h5(adata, file):
    if '.h5' not in file:
        file = f'{file}.h5'
    if Path(file).exists():
        raise FileExistsError(f"There already is a file `{file}`.")
    def int_max(x):
        return int(max(np.floor(len(str(int(max(x)))) / 4), 1) * 4)
    def str_max(x):
        return max([len(i) for i in x])
    if os.path.exists(file):
        os.remove(file)
    w = h5py.File(file, 'w')
    grp = w.create_group("matrix")
    grp.create_dataset("barcodes", data=np.array(adata.obs_names, dtype=f'|S{str_max(adata.obs_names)}'))
    grp.create_dataset("data", data=np.array(adata.X.data, dtype=f'<i{int_max(adata.X.data)}'))
    ftrs = grp.create_group("features")
    # this group will lack the following keys:
    # '_all_tag_keys', 'feature_type', 'genome', 'id', 'name', 'pattern', 'read', 'sequence'
    ftrs.create_dataset("feature_type", data=np.repeat('Gene Expression',adata.var.shape[0]).astype('S15'))
    ftrs.create_dataset("genome", data=np.repeat('GRCh38',adata.var.shape[0]).astype('S6'))
    ftrs.create_dataset("id", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    ftrs.create_dataset("name", data=np.array(adata.var.index, dtype=f'|S{str_max(adata.var.index)}'))
    grp.create_dataset("indices", data=np.array(adata.X.indices, dtype=f'<i{int_max(adata.X.indices)}'))
    grp.create_dataset("indptr", data=np.array(adata.X.indptr, dtype=f'<i{int_max(adata.X.indptr)}'))
    grp.create_dataset("shape", data=np.array(list(adata.X.shape)[::-1], dtype=f'<i{int_max(adata.X.shape)}'))
cell=sys.argv[1]

os.chdir('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSHC_multiome/processed/ArchR/')
output_dir=f'/sc/arion/projects/CommonMind/roussp01a/snmulti/motif/files/class/{cell}'

adata=sc.read_h5ad(f'{cell}/merged.h5ad')
anno=pd.read_csv('save1/cell_metadata_full.csv.gz',index_col=0)
#anno.index=anno.index.str.replace('#','_')
anno=anno.loc[anno.TSSEnrichment>3]
anno=anno.loc[~anno.subcellclass.str.startswith('Gaba_AC')]
anno=anno.loc[~anno.region.eq('caudate nucleus')]
adata=adata[adata.obs_names.isin(anno.index),]
adata.obs=anno.loc[adata.obs_names]
adata.X.data[adata.X.data>1]=1
z = adata.var.index.str.extract(r"([^:]+):(\d+)-(\d+)")
#[["Chromosome", "Start", "End"]]
adata.var["Chromosome"]=z[0].values
adata.var["Start"] = z[1].astype(int).values
adata.var["End"] = z[2].astype(int).values
adata.var['cvg']=np.array(adata.X.sum(axis=0)).flatten()/adata.obs.shape[0]


adata.obs['celltype']=adata.obs['subcellclass']
adata.obs['RiP']=np.array(adata.X.sum(axis=1)).flatten()
z=adata.obs.celltype.value_counts()
adata=adata[adata.obs.celltype.isin(z.loc[z>100].index)*(adata.obs.RiP>100),]


meta=adata.obs
sh=adata.obs_names

cell_depth=meta['RiP'].values
threshold=meta.RiP.mean().round()*20

pbulk_centers = []
barcode_groups = []
cell_barcodes = sh
#embed=embedding.values

# Go through each cell type
results = Parallel(n_jobs=16, backend="loky")(
    delayed(process_one_celltype_simple)(ct)
    for ct in np.unique(meta["celltype"])
)
barcode_groups = pd.concat(results, ignore_index=True)
barcode_groups.to_csv(f'{output_dir}/RUSHC2_barcode_groups.csv.gz',compression='gzip')



z=pd.get_dummies(barcode_groups.group)
mat=z.values
mat=mat/mat.sum(axis=0)

cov=meta.loc[barcode_groups.barcode.values,
             ['TSSEnrichment','NucleosomeRatio','BlacklistRatio','nFrags']]
cov=pd.DataFrame(mat.T.dot(cov),columns=['TSSEnrichment','NucleosomeRatio',
                                         'BlacklistRatio','nFrags'])
cov.index=z.columns.values


ad1=adata[barcode_groups.barcode.values,]
ad1.obs['pb']=barcode_groups['group'].values
pb=get_pseduobulk(ad1,'pb')
obs=cov.loc[pb.columns]

adata=ad.AnnData(sparse.csr_matrix(pb.values.T),
              obs=obs,var=adata.var)

def write_input4R(adata,output_prefix):
    if os.path.exists(f"{output_prefix}.h5"):
        os.remove(f"{output_prefix}.h5")
    write_10X_h5(adata, f"{output_prefix}.h5")
    adata.obs.to_csv(f'{output_prefix}.obs.csv.gz', compression='gzip')
    adata.var.to_csv(f'{output_prefix}.var.csv.gz', compression='gzip')

write_input4R(adata, f'{output_dir}/RUSHC_metacell_20')    

