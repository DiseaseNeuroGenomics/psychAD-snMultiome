# taxonomy
### Contents
```text
├──01_readRNA_qc_collect.py        #collect QC information and combine all h5ad files
├──02_RNA_initial.py               #major class level RNA-seq cell clustering
├──03_ATAC_initial.R               #major class level ATAC-seq LSI
├──04_ATA_WNN_initial.py           #class level ATAC-seq taxonomy and WNN
├──05_adaptive_qc.py               #major cell class level cell type specific filtering
├──06_RNA_2nd.py RNA-seq           #cell clustering within major cell class
├──07_ATAC_2nd.R ATAC-seq          #LSI within major cell class
├──08_ATAC_WNN_2nd.py              #within major cell class ATAC-seq taxonomy and WNN
├──09_RNA_3rd.py                   #cell clustering within major cell subclass
├──10_ATAC_3rd.R                   #LSI within major cell subclass
├──11_ATAC_agg.py                  #ATAC-seq single modality cells annotation 
├──12_ATAC_callPeak_major_class.R  #call peak for each major cell class & export h5ad         
├──13_consensus_peakset.R          #consensus peaks across cell types export peak/gene-score h5ad
└──helper                          #custom helper functions in R and python
```


