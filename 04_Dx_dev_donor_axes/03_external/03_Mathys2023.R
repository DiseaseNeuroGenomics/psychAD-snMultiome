.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(crumblr)
library(ggrepel)
library(rhdf5)
library(variancePartition)

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
setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/donor/new')

Info1=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSHM_ATAC/rush_dataset_865_cross-sectional_12-06-2024.csv')
Info2=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSHM_ATAC/ROSMAP_clinical.csv')
lgb=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate2/RUSHM_Dx_cell.csv.gz')
z=as.character(unique(lgb$projid))
metaInfo=data.frame(Info1[z,],Info2[z,setdiff(colnames(Info2),colnames(Info1))])
rownames(metaInfo)=metaInfo$individualID
Obs=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSHM_ATAC/RNAprocessed/obs_qc.csv.gz')
Dx_dev <- list()
Dx_dev_adj <- list()
for(cell in c('Astro','EN','Oligo','OPC','IN','Micro_PVM')){
  lgb1=lgb[lgb$class==cell,]
  agg=as.data.frame(lgb1%>%group_by(SubID)%>%summarise(Dx_cell=mean(Dx_cell)))
  rownames(agg)=agg$SubID
  agg$SubID=NULL
  
  covs=c('n_genes_by_counts','pct_counts_mt','pct_counts_ribo')
  covmetaInfo=as.data.frame(data.frame(Obs[rownames(lgb1),c(covs)],SubID=lgb1$SubID)%>%group_by(SubID)%>%summarise(n_genes_by_counts=mean(n_genes_by_counts),
                                                                                                pct_counts_mt=mean(pct_counts_mt),
                                                                                                pct_counts_ribo=mean(pct_counts_ribo)))
  rownames(covmetaInfo)=covmetaInfo$SubID
  form1=~scale(age_death)+msex+scale(gpath)+scale(pmi)
  zmetaInfo=data.frame(covmetaInfo[rownames(agg),],metaInfo[rownames(agg),])
  Dx_dev[[cell]]=t(residuals(dream(t(agg),form1,zmetaInfo)))
  form=~scale(age_death)+msex+scale(gpath)+scale(pmi)+scale(n_genes_by_counts)+scale(pct_counts_mt)+scale(pct_counts_ribo)
  Dx_dev_adj[[cell]]=t(residuals(dream(t(agg),form,zmetaInfo)))
}
f=unique(unlist(lapply(Dx_dev,rownames)))
Dx_dev=do.call(cbind,lapply(Dx_dev,function(x)data.frame(x)[f,1]))
rownames(Dx_dev)=f
Dx_dev_adj=do.call(cbind,lapply(Dx_dev_adj,function(x)data.frame(x)[f,1]))
rownames(Dx_dev_adj)=f
write.csv(Dx_dev,'external/Mathys23_Dx_dev.csv')
write.csv(Dx_dev_adj,'external/Mathys23_Dx_dev_adj.csv')

zInfo=metaInfo[rownames(Dx_dev_adj),]
form=~scale(age_death)+msex+scale(gpath)+scale(cogng_demog_slope)+scale(pmi)+(1|apoe_genotype)
zfit1=eBayes(dream(t(scale(Dx_dev_adj)),form,zInfo[rownames(Dx_dev_adj),]))
z1=topTable(zfit1,coef='scale(cogng_demog_slope)')
form=~scale(age_death)+msex+scale(gpath)+scale(cogng_path_slope)+scale(pmi)+(1|apoe_genotype)
zfit1=eBayes(dream(t(scale(Dx_dev_adj)),form,zInfo[rownames(Dx_dev_adj),]))
z2=topTable(zfit1,coef='scale(cogng_path_slope)')
form=~scale(age_death)+msex+scale(gpath)+scale(cts_mmse30_lv)+scale(pmi)+(1|apoe_genotype)
zfit1=eBayes(dream(t(scale(Dx_dev_adj)),form,zInfo[rownames(Dx_dev_adj),]))
z3=topTable(zfit1,coef='scale(cts_mmse30_lv)')
write.csv(z1,file='external/Mathys23_cogng_demog_slope_Dx_dev.csv')
write.csv(z3,file='external/Mathys23_cts_mmse30_lv_Dx_dev.csv')
