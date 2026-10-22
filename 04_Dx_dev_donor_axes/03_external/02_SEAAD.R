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
Obs=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/class_SEAAD/SeaAD_obs.csv.gz')
lgb1=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate3/SEA_MTG_Dx_cell.csv.gz')
lgb2=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate3/SEA_PFC_Dx_cell.csv.gz')

lgb1$sample=paste0(lgb1[,'Donor ID'],'_','MTG')
lgb2$sample=paste0(lgb2[,'Donor ID'],'_','PFC')
lgb=rbind(lgb1[,c('sample','Dx_cell','Dx_cat_cell','class','Subclass')],
          lgb2[,c('sample','Dx_cell','Dx_cat_cell','class','Subclass')])
sh=intersect(rownames(lgb),rownames(Obs))
lgb=lgb[sh,]
#Obs=Obs[sh,]
Info=unique(lgb1[,c('Donor ID','Continuous Pseudo-progression Score',
                    'Sex','Age at Death','PMI','Cognitive Status','APOE Genotype','Severely Affected Donor','CERAD score','Braak')])
colnames(Info)=c('SubID','CPS','Sex','Age','PMI','Cogn','apoe','SeverDonor','CERAD','BRAAK_AD')
Info$BRAAK_AD=as.numeric(factor(Info$BRAAK_AD,levels=c("Braak 0","Braak II","Braak III","Braak IV","Braak V","Braak VI"),
                               ordered=T))
Info$CERAD=as.numeric(factor(Info[,'CERAD'],levels=c("Absent","Sparse","Moderate","Frequent"),
                            ordered=T))
Info$PMI=as.numeric(Info$PMI)
a=b=Info
a$BrainRegion='MTG'
b$BrainRegion='PFC'
rownames(a)=a$sample=paste0(a$SubID,'_','MTG')
rownames(b)=b$sample=paste0(b$SubID,'_','PFC')
metaInfo=rbind(a,b)
####
genes=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/metadata/cellranger_gene.csv')
genes[duplicated(genes$gene_name),'gene_name']=paste0(genes[duplicated(genes$gene_name),'gene_name'],'-1')
genes$mt = grepl('^MT-',genes$gene_name)
genes$ribo = grepl('^RPL',genes$gene_name) |grepl('^RPS',genes$gene_name) # annotate the group of mitochondrial genes as 'mt'
genes$hb=grepl('^HB',genes$gene_name)
genes$protein_coding = genes$gene_type=='protein_coding'
genes$filtered=genes$protein_coding & (!genes$ribo)& (!genes$mt) & genes$seqnames%in%paste0('chr',c(1:22,'X'))
#genes=genes[genes$filtered,]
rownames(genes)=genes$gene_name
##
#
Neu=lgb[lgb$class%in%c('EN'),]
ncells=table(Neu$sample,Neu$Subclass)
ncells=apply(ncells,2,c)
df=data.frame(vul=ncells[,c('L2/3 IT')],
              nvul=rowSums(ncells)-ncells[,'L2/3 IT'])
cobj=crumblr(df)
vul1=cobj
Neu=lgb[lgb$class%in%c('IN'),]
ncells=table(Neu$sample,Neu$Subclass)
ncells=apply(ncells,2,c)
df=data.frame(vul=ncells[,c('Sst')],
              nvul=rowSums(ncells)-ncells[,'Sst'])
cobj=crumblr(df)
vul2=cobj

Neu=lgb[lgb$class%in%c('EN','IN'),]
ncells=table(Neu$sample,Neu$Subclass)
ncells=apply(ncells,2,c)
df=data.frame(vul=rowSums(ncells[,c('L2/3 IT','Sst')]),
              nvul=rowSums(ncells)-rowSums(ncells[,c('L2/3 IT','Sst')]))
cobj=crumblr(df)
vul=cobj
sh=intersect(colnames(vul1),colnames(vul2))
Vul=new("EList")
Vul$E <- rbind(vul$E['vul',sh],
               vul1$E['vul',sh],
               vul2$E['vul',sh])
Vul$weights <- rbind(vul$weights['vul',sh],
                     vul1$weights['vul',sh],
                     vul2$weights['vul',sh])
rownames(Vul$E)=rownames(Vul$weights)=c('Neu','EN','IN')
#######subtype#######
Obs$celltype=Obs$Supertype
Obs[Obs$class%in%c('EN'),]$celltype=Obs[Obs$class%in%c('EN'),]$Subclass
ncells=lapply(split(Obs,Obs$class),function(x)table(x$Sample,x$celltype))
ncells=lapply(ncells,function(x)apply(x,2,c))
z1=apply(table(Obs$Sample,Obs$class),2,c)
z2=apply(table(Obs$Sample,Obs$Subclass),2,c)
z=apply(table(Obs[Obs$class%in%c('IN',"EN"),]$Sample,Obs[Obs$class%in%c('IN',"EN"),]$Supertype),2,c)
ncells$class=z1
ncells$subclass=z2
ncells$neuron=z3
cobjs=lapply(ncells,crumblr)
####Dx_dev
Dx_dev <- list()
Dx_dev_adj <- list()
for(cell in c('Astro','EN','Oligo','OPC','IN','Micro_PVM')){
  lgb1=lgb[lgb$class==cell,]
  agg=as.data.frame(lgb1%>%group_by(sample)%>%summarise(Dx_cell=mean(Dx_cell)))
  rownames(agg)=agg$sample
  agg$sample=NULL

  covs=c('n_genes_by_counts','pct_counts_mt','pct_counts_ribo')
  covInfo=as.data.frame(Obs[rownames(lgb1),c(covs,'Sample')]%>%group_by(Sample)%>%summarise(n_genes_by_counts=mean(n_genes_by_counts),
                                                                                            pct_counts_mt=mean(pct_counts_mt),
                                                                                            pct_counts_ribo=mean(pct_counts_ribo)))
  rownames(covInfo)=covInfo$Sample
  form1=~scale(Age)+Sex+scale(CERAD)+scale(BRAAK_AD)+scale(PMI)
  zInfo=data.frame(covInfo[rownames(agg),],metaInfo[rownames(agg),])
  Dx_dev[[cell]]=t(residuals(dream(t(agg),form1,zInfo)))
  form=~scale(Age)+Sex+scale(CERAD)+scale(BRAAK_AD)+scale(PMI)+scale(n_genes_by_counts)+scale(pct_counts_mt)+scale(pct_counts_ribo)
  Dx_dev_adj[[cell]]=t(residuals(dream(t(agg),form,zInfo)))
}
f=unique(unlist(lapply(Dx_dev,rownames)))
Dx_dev=do.call(cbind,lapply(Dx_dev,function(x)data.frame(x)[f,1]))
rownames(Dx_dev)=f
Dx_dev_adj=do.call(cbind,lapply(Dx_dev_adj,function(x)data.frame(x)[f,1]))
rownames(Dx_dev_adj)=f
write.csv(Dx_dev,'external/SEAAD_Dx_dev.csv')
write.csv(Dx_dev_adj,'external/SEAAD_Dx_dev_adj.csv')
Dx_dev=Dx_dev_adj # use adjusted Dx_dev
#######projection#######################
F_sel=paste0('Factor',1:3)
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/class_SEAAD/all.proc.Rdata')
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/class_SEAAD/all/scale(BRAAK_AD)_scale(CERAD).Rdata')
GEX_Dobj=fit
SGexs=lapply(c('Astro','Oligo','OPC','Micro_PVM','EN','IN'),function(cell){
  ###
  E=proc[[cell]]$E
  fit1=GEX_Dobj[[cell]]
  rowsh=names(fit1$Amean)
  colsh=rownames(fit1$design)
  cov=colnames(fit1$coefficients)
  R=E[rowsh,colsh]-fit1$coefficients[rowsh,cov]%*%t(fit1$design[,cov])  
  
  rownames(R)=paste0(cell,":",rownames(R))
  R
})
load('MOFA_top_25_10_agg.Rdata')
W=agg$consensus_weights
f=table(unlist(lapply(SGexs,colnames)))
f=names(f)[f==length(SGexs)]
SGex=do.call(rbind,lapply(SGexs,function(x)x[,f]))
sh=intersect(rownames(SGex),rownames(W))
W1=W[sh,]
z=do.call(rbind,strsplit(rownames(W1),'\\:'))
Wgli=W1[!z[,1]%in%c("EN",'IN'),]
Wgli2=W1[!z[,1]%in%c("EN",'IN','Micro_PVM'),]
Fp=scale(t(SGex[rownames(W1),])%*%W1)
Fgli=data.frame(scale(t(SGex[rownames(Wgli),])%*%Wgli))
Fgli$Fmeta=scale(Fgli[,2])-scale(Fgli[,3])
Fgli2=data.frame(scale(t(SGex[rownames(Wgli2),])%*%Wgli2))
Fgli2$Fmeta=scale(Fgli2[,2])-scale(Fgli2[,3])
write.csv(Fp,'external/SEAAD_Factors.csv')
write.csv(Fgli,'external/SEAAD_Factors_glia.csv')
write.csv(Fgli2,'external/SEAAD_Factors_glia_nomicro.csv')
#######################correlation with cell composition #####
##
pfc=rownames(Fp)[metaInfo[rownames(Fp),]$BrainRegion=='PFC']
pfc=Fp[pfc,]
rownames(pfc)=metaInfo[rownames(pfc),]$SubID
br='PFC'
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
zdf$q=p.adjust(zdf$P,method='BH')
zdf1=zdf
xx0=xx
br='MTG'
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
zdf$BrainRegion='MTG'
zdf1$BrainRegion='PFC'
zdf1$coef=rownames(zdf1)
zdf$coef=rownames(zdf)
zdf=rbind(zdf1,zdf)
zdf$q=p.adjust(zdf$P,method='BH')
write.csv(zdf,file='external/SEAAD_brainregion_factor_Dx_dev_asso.csv')
write.csv(rbind(xx0,xx),file='external/SEAAD_brainregion_F123_Dx_dev.csv')
zInfo=data.frame(Dx_dev=rowMeans(scale(Dx_dev),na.rm=T)[rownames(Fp)],
                 metaInfo[rownames(Fp),c('Age','DxPC1','CPS')])

###########################
sh=intersect(intersect(colnames(cobjs$neuron),rownames(Fgli)), # use MTG only
             rownames(metaInfo[metaInfo$BrainRegion=='MTG',]))
cobj=cobjs$neuron[,sh]

zInfo=data.frame(metaInfo[sh,],data.frame(Fgli)[sh,])
zInfo$Dx_dev=rowMeans(scale(Dx_dev[,c('Astro','Oligo','OPC','Micro_PVM')]),na.rm=T)[rownames(zInfo)]
form1=~scale(Age)+Sex+scale(CPS)+scale(PMI)
zfit1=eBayes(dream(cobj,form1,zInfo))
z1=topTable(zfit1,'scale(CPS)',number=Inf)
write.csv(z1,file='external/SEAAD_CPS_subtype_crumblr.csv')

for(xx in c('Factor1',"Factor2",'Factor3','Dx_dev')){
  form=as.formula(sprintf('~scale(Age)+Sex+scale(CPS)+scale(PMI)+scale(%s)',xx))
  zfit=eBayes(dream(cobj,form,zInfo))
  z2=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  z2$se=z2$logFC/z2$t
  write.csv(z2,file=sprintf('external/SEAAD_%s_subtype_crumblr.csv',xx))
}

############################vulnerbale neurons######
agg1=lgb%>%group_by(sample,class)%>%summarise(Dx_cell=mean(Dx_cell))
agg1=reshape2::acast(agg1,sample~class,value.var='Dx_cell')

sh=intersect(intersect(colnames(Vul),rownames(Fgli)),
             rownames(metaInfo[metaInfo$BrainRegion=='MTG',]))
cobj=Vul[,sh]
zInfo=data.frame(metaInfo[sh,],data.frame(Fgli)[sh,])
zInfo$Dx_dev=rowMeans(scale(Dx_dev[,c('Astro','Oligo','OPC','Micro_PVM')]),na.rm=T)[rownames(zInfo)]
zInfo$Dx_cell=rowMeans(scale(agg1[,c('Astro','Oligo','OPC','Micro_PVM')]),na.rm=T)[rownames(zInfo)]

form0=~scale(Age)+Sex+scale(Dx_cell)+scale(PMI)
zfit=eBayes(dream(cobj,form0,zInfo))
z1=topTable(zfit,'scale(Dx_cell)',number=Inf)
write.csv(z1,file='external/SEAAD_Dx_cell_vulneuron_crumblr.csv')

form0=~scale(Age)+Sex+scale(CPS)+scale(PMI)
zfit=eBayes(dream(cobj,form0,zInfo))
z1=topTable(zfit,'scale(CPS)',number=Inf)
write.csv(z1,file='external/SEAAD_CPS_vulneuron_crumblr.csv')
form1=~scale(Age)+Sex+scale(CPS)+scale(PMI)+scale(Dx_dev)
zfit=eBayes(dream(cobj,form1,zInfo))
z2=topTable(zfit,'scale(Dx_dev)',number=Inf)
z2$se=z2$logFC/z2$t
write.csv(z2,file=sprintf('external/SEAAD_Dx_dev_vulneuron_crumblr.csv',xx))
###########DAM###########

sh=intersect(intersect(colnames(cobjs$Micro_PVM),rownames(Fgli2)),
             rownames(metaInfo[metaInfo$BrainRegion=='MTG',]))
cobj=cobjs$Micro_PVM[,sh]
zInfo=data.frame(metaInfo[sh,],data.frame(Fgli)[sh,])
zInfo$Dx_dev=rowMeans(scale(Dx_dev[,c('Astro','Oligo','OPC')]),na.rm=T)[rownames(zInfo)]
zInfo$Dx_cell=rowMeans(scale(agg1[,c('Astro','Oligo','OPC')]),na.rm=T)[rownames(zInfo)]

form0=~scale(Age)+Sex+scale(Dx_cell)+scale(PMI)
zfit=eBayes(dream(cobj,form0,zInfo))
z1=topTable(zfit,'scale(Dx_cell)',number=Inf)
write.csv(z1,file='external/SEAAD_Dx_cell_DAM_crumblr.csv')

form0=~scale(Age)+Sex+scale(CPS)+scale(PMI)
zfit=eBayes(dream(cobj,form0,zInfo))
z1=topTable(zfit,'scale(CPS)',number=Inf)
write.csv(z1,file='external/SEAAD_Dx_donor_DAM_crumblr.csv')

form0=~scale(Age)+Sex+scale(CPS)+scale(PMI)+scale(Dx_dev)
zfit=eBayes(dream(cobj,form0,zInfo))
z1=topTable(zfit,'scale(Dx_dev)',number=Inf)
write.csv(z1,file='external/SEAAD_Dx_dev_DAM_crumblr.csv')

