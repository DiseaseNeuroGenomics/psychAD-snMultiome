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
read_obsvar <- function(path,sub='obs'){
  dat = h5read(path,sub)
  l=sapply(dat,length)
  dat1=do.call(data.frame,lapply(dat[l==2],function(x)c(NA,x$categories)[x$codes+2]))
  dat2=do.call(data.frame,dat[l!=2])
  dat=data.frame(dat1,dat2)
  #rownames(dat)=dat$X_index
  #rownames(dat)=dat[,grepl('key$',colnames(dat))]
  dat
}
setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/donor/new')
lgb=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate3/RUSHC_Dx_cell.csv.gz')
metaInfo=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSH_metadata.csv')

metaInfo$Sex=metaInfo$sex
metaInfo$BrainRegion=metaInfo$Region
Obs=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSHC_multiome/processed/obs.csv.gz')
#rownames(Obs)=paste0(Obs$batch,'_',Obs$barcode)
Obs=Obs[Obs$region!='caudate nucleus' & Obs$broadcellclass!='Epend',]
Obs$sample=gsub('_ST','_STG',gsub('_MF','_PFC',Obs$specimenID))
Obs$subclass=Obs$subcellclass
Obs[grepl('Gaba_AC_',Obs$subclass),]$subclass=NA
######
lgb$sample=paste0(lgb$SubID,'_',lgb$Region)
cobjs=lapply(ncells,crumblr)
###############Dx_dev######
Dx_dev <- list()
Dx_dev_adj <- list()
for(cell in c('Astro','EN','Oligo','OPC','IN','Micro_PVM')){
  lgb1=lgb[lgb$class==cell,]
  agg=as.data.frame(lgb1%>%group_by(sample)%>%summarise(Dx_cell=mean(Dx_cell)))
  rownames(agg)=agg$sample
  agg$sample=NULL
  
  covs=c('n_genes_by_counts','mt_percent','pct_counts_ribo')
  covmetaInfo=as.data.frame(Obs[rownames(lgb1),c(covs,'sample')]%>%group_by(sample)%>%summarise(n_genes_by_counts=mean(n_genes_by_counts),
                                                                                            pct_counts_mt=mean(mt_percent),
                                                                                            pct_counts_ribo=mean(pct_counts_ribo)))
  rownames(covmetaInfo)=covmetaInfo$sample
  form1=~scale(Age)+Sex+scale(CERAD)+scale(BRAAK_AD)+scale(PMI)
  zmetaInfo=data.frame(covmetaInfo[rownames(agg),],metaInfo[rownames(agg),])
  Dx_dev[[cell]]=t(residuals(dream(t(agg),form1,zmetaInfo)))
  form=~scale(Age)+Sex+scale(CERAD)+scale(BRAAK_AD)+scale(PMI)+scale(n_genes_by_counts)+scale(pct_counts_mt)+scale(pct_counts_ribo)
  Dx_dev_adj[[cell]]=t(residuals(dream(t(agg),form,zmetaInfo)))
}
f=unique(unlist(lapply(Dx_dev,rownames)))
Dx_dev=do.call(cbind,lapply(Dx_dev,function(x)data.frame(x)[f,1]))
rownames(Dx_dev)=f
Dx_dev_adj=do.call(cbind,lapply(Dx_dev_adj,function(x)data.frame(x)[f,1]))
rownames(Dx_dev_adj)=f
write.csv(Dx_dev,'external/Luquez_Dx_dev.csv')
write.csv(Dx_dev_adj,'external/Luquez_Dx_dev_adj.csv')
Dx_dev=Dx_dev_adj
#########projection###################################
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSHC_multiome/processed/DE/class_proc.Rdata')
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/RUSHC_multiome/processed/DE/scale(BRAAK_AD)_scale(CERAD).Rdata')
GEX_Dobj=fit
LGexs=lapply(c('Astro','Oligo','OPC','Micro_PVM','EN','IN'),function(cell){
  ###
  E=res.proc[[cell]]$E
  fit1=GEX_Dobj[[cell]]
  rowsh=names(fit1$Amean)
  colsh=rownames(fit1$design)
  cov=colnames(fit1$coefficients)
  R=E[rowsh,colsh]-fit1$coefficients[rowsh,cov]%*%t(fit1$design[,cov])  
  
  rownames(R)=paste0(cell,":",rownames(R))
  R
})

f=table(unlist(lapply(LGexs,colnames)))
f=names(f)[f==length(LGexs)]
LGex=do.call(rbind,lapply(LGexs,function(x)x[,f]))
colnames(LGex)=gsub('_ST','_STG',gsub('_MF','_PFC',colnames(LGex)))
load('MOFA_top_25_10_agg.Rdata')
W=agg$consensus_weights
sh=intersect(rownames(LGex),rownames(W))
W1=W[sh,]
z=do.call(rbind,strsplit(rownames(W1),'\\:'))
Wgli=W1[!z[,1]%in%c("EN",'IN'),]
Wgli2=W1[!z[,1]%in%c("EN",'IN','Micro_PVM'),]
Fp=scale(t(LGex[rownames(W1),])%*%W1)
Fgli=data.frame(scale(t(LGex[rownames(Wgli),])%*%Wgli))
Fgli$Fmeta=scale(Fgli[,2])-scale(Fgli[,3])
Fgli2=data.frame(scale(t(LGex[rownames(Wgli2),])%*%Wgli2))
Fgli2$Fmeta=scale(Fgli2[,2])-scale(Fgli2[,3])

write.csv(Fp,'external/Luquez_Factors.csv')
write.csv(Fgli,'external/Luquez_Factors_glia.csv')
write.csv(Fgli2,'external/Luquez_Factors_glia_nomicro.csv')
F_sel=paste0('Factor',1:3)

#######correlation ######
pfc=rownames(Fp)[metaInfo[rownames(Fp),]$BrainRegion=='PFC']
pfc=Fp[pfc,]
rownames(pfc)=metaInfo[rownames(pfc),]$SubID
br='PFC'
zInfo=metaInfo[metaInfo$BrainRegion==br,]
zInfo$Dx_dev=rowMeans(scale(data.frame(Dx_dev)[rownames(zInfo),]),na.rm=T)
zInfo=zInfo[!duplicated(zInfo$SubID),]
rownames(zInfo)=zInfo$SubID
zInfo=data.frame(zInfo,data.frame(pfc)[rownames(zInfo),])
xx0=zInfo[,c('Dx_dev','BrainRegion',colnames(Fp))]
zdf=do.call(rbind,lapply(c(paste0('Factor',1:3)),
                         function(x)summary(lm(as.formula(sprintf('Dx_dev~scale(CERAD)+scale(BRAAK_AD)+Age+Sex+PMI+scale(%s)',x)),zInfo))$coefficient[sprintf('scale(%s)',x),]))
rownames(zdf)=paste0('Factor',1:3)
colnames(zdf)=c('FC','SE','t','P')
zdf=data.frame(zdf)
zdf1=zdf
br='STG'
zInfo=metaInfo[metaInfo$BrainRegion==br,]
zInfo$Dx_dev=rowMeans(scale(data.frame(Dx_dev)[rownames(zInfo),]),na.rm=T)
zInfo=zInfo[!duplicated(zInfo$SubID),]
rownames(zInfo)=zInfo$SubID
zInfo=data.frame(zInfo,data.frame(pfc)[rownames(zInfo),])
xx=zInfo[,c('Dx_dev','BrainRegion',colnames(Fp))]
zdf=do.call(rbind,lapply(c(paste0('Factor',1:3)),
                         function(x)summary(lm(as.formula(sprintf('Dx_dev~scale(CERAD)+scale(BRAAK_AD)+Age+Sex+PMI+scale(%s)',x)),zInfo))$coefficient[sprintf('scale(%s)',x),]))
rownames(zdf)=paste0('Factor',1:3)
colnames(zdf)=c('FC','SE','t','P')
zdf=data.frame(zdf)
zdf$BrainRegion='STG'
zdf1$BrainRegion='PFC'
zdf$coef=rownames(zdf)
zdf1$coef=rownames(zdf1)
zdf=rbind(zdf1,zdf)
zdf$q=p.adjust(zdf$P,method='BH')
write.csv(zdf,file='external/Luquez_brainregion_factor_Dx_dev_asso.csv')
write.csv(rbind(xx0,xx),file='external/Luquez_brainregion_F123_Dx_dev.csv')
