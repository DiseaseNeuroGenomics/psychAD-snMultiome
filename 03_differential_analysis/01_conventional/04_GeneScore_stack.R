.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(dreamlet)
library(SingleCellExperiment)
library(tidyverse)
library(zenith)
library(qvalue)
library(cowplot)
library(splines)

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
level='subclass'
BR='all'
OUTDIR='/sc/arion/projects/CommonMind/roussp01a/snmulti/files/DE/files/'
prefix='GS'
cellanno_path='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno_05092025_fixed3.csv'
sampleCol='order'
NCPU=8
bpparam = MulticoreParam(NCPU)
setwd(file.path(OUTDIR,level))
##########################stacked assay #######################################
load(paste0('pb4dreamlet_',prefix,'_',level,'.Rdata'))
pb.stack <- stackAssays(pb1)
form=~(1|BrainRegion)+(1|SubID)+(1|Sex)+scale(Age)+scale(Dx_donor)+scale(PMI..min.)+scale(TSSEnrichment)+scale(log(nFrags))+stackedAssay
bpparam = MulticoreParam(8) # use less cpu
res.proc <- processAssays(pb.stack, form, min.count = 10,min.prop=0.01, BPPARAM = bpparam)
form=~(1|Dx_cat_donor)+(1|Sex)+scale(Age)+scale(PMI..min.)+(1|SubID)+(1|BrainRegion)+(1|class)+ (1|class:stackedAssay)+scale(TSSEnrichment)+scale(nFrags)
varPart=fitVarPart(res.proc,form, BPPARAM = bpparam)
save(varPart,file='stackedassay_varpart.Rdata')


