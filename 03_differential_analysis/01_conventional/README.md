# Conventional differential analysis with dreamlet
### Contents
```text
├──01_GEX.R                                #differential gene expression at class level
├──02_GeneScore.R                          #gene activity score
├──03_GEX_stack.R                          #Stack assays(assay means subclass) from pseudobulk to perform analysis across cell types  
├──04_GeneScore_stack                      #Stack assays(subclass) for gene activity score
├──05_peak_crosscelltype.R                 #input for chromvar stack assays analysis, peak across major cell types
├──06_chromVar_stack.R                     #Stack assays(subclass) for chromvar motif activity
├──07_stack_varPart_across_assay.R         #compare varianceParition result between gene expression and chromatin accessbility
├──process_assay/                          #dreamlet aggregateToPseudoBulk rewrap (can use processAssays instead)
    ├──01_pb_sampleDx.py                   #get psedobulk from h5ad 
    └──02_aggregateToPseudoBulk_rewrap.R   #prepare PseudoBulk for voom processing.    
```
Run process_assay/ to get pseudo bulk before differential analysis. Can also be done with dreamlet aggregateToPseudoBulk
 function, details https://diseaseneurogenomics.github.io/dreamlet/articles/h5ad_on_disk.html.

