.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(dreamlet)
library(SingleCellExperiment)
library(tidyverse)
library(zenith)
library(qvalue)
library(cowplot)
library(splines)
# ATAC-seq peak cross cell type
# python process_assay/01_pb_sampleDx.py  --input all.h5ad --cellanno cellanno.csv --outDir OUTDIR --outPrefix ATAC_{cell} --ATAC --binarize
# Rscript process_assay/02_aggregateToPseudoBulk_rewrap.R all.h5ad cellanno.csv  path_metaInfo OUTDIR ATAC_{cell} ATAC

read_big_csv <- function(fh,row1asrowname=T){
  library(data.table)
  dat=fread(fh)
  dat=as.data.frame(dat)
  if(row1asrowname){
    rownames(dat)=dat[,1]
    dat[,1]=NULL
  }
  return(dat)
}
cdmk_folder <- function(folder){
  if(!dir.exists(folder)){
    dir.create(folder,recursive = T)  
  }
  setwd(folder)
}
#########################################################
cell='all'
level='subclass'
BR='all'
OUTDIR='/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/'
prefix=paste0('ATAC_',cell)
sampleCol='order'
NCPU=8 # less CPU to save memory
bpparam = MulticoreParam(NCPU)
############################################################################

setwd(file.path(OUTDIR,level))
load(paste0('pb4dreamlet_',prefix,'_',level,'.Rdata'))
#voom
form=~(1|BrainRegion)+(1|SubID)+scale(TSSEnrichment)+scale(nFrags)+(1|Sex)+scale(Age)+scale(Dx_cell)+scale(Dx_donor)
res.proc <- processAssays(pb1, form,min.cells = 20, BPPARAM = MulticoreParam(NCPU)) # peak is sparse
z=reshape2::melt(do.call(rbind,pb1@int_colData$n_cells))
rownames(z)=with(z,paste(Var1,Var2,sep='_'))
proc@metadata$ncell=z[paste0(proc@metadata$sampleDx,'_',unlist(proc@metadata[,level])),]$value
cdmk_folder(prefix)
save(proc,file='proc.Rdata')
