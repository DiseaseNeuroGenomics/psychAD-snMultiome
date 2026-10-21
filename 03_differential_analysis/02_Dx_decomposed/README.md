# Dx-cell decomposed differential analysis for different molecular features
### Contents
```text
├──01_GEX.R                                #Dx-cell decomposed differential gene expression
├──02_GeneScore.R                          #gene activity score
├──03_OCR.R                                #ATAC-seq peaks, run each major cell type at a time
├──04_motif.R                              #chromVar motif across samples, run each major cell type at a time
├──process_assay/.                         #dreamlet aggregateToPseudoBulk rewrap (can use processAssays instead)
│   ├──01_pb_sampleDx.py                   #get psedobulk from h5ad 
│   └──02_aggregateToPseudoBulk_rewrap.R   #prepare PseudoBulk for voom processing.               
```
Run process_assay/ to get pseudo bulk before differential analysis. Can also be done with dreamlet aggregateToPseudoBulk function, details https://diseaseneurogenomics.github.io/dreamlet/articles/h5ad_on_disk.html. 
