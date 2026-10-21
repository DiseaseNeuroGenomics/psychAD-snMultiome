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
NPSADcontrasts=local({load('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/metadata/contrasts_for_dreamlet.Rdata');environment()})

## read information 
genes=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/metadata/cellranger_gene.csv')
genes[duplicated(genes$gene_name),'gene_name']=paste0(genes[duplicated(genes$gene_name),'gene_name'],'-1')
genes$mt = grepl('^MT-',genes$gene_name)
genes$ribo = grepl('^RPL',genes$gene_name) |grepl('^RPS',genes$gene_name) # annotate the group of mitochondrial genes as 'mt'
genes$hb=grepl('^HB',genes$gene_name)
genes$protein_coding = genes$gene_type=='protein_coding'
genes$filtered=genes$protein_coding & (!genes$ribo)& (!genes$mt) & genes$seqnames%in%paste0('chr',c(1:22,'X'))
#genes=genes[genes$filtered,]
rownames(genes)=genes$gene_name

meta=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/metadata/psychAD_meta_PC1_pathonly.csv')
MDx_dev=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate3/Dx_dev.csv')
MInfo=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/files/allInfo_11102023.csv')
rownames(MInfo)=MInfo$order
Mlgb=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ML_cellstate3/lgb_anno.csv')
lgb=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate3/NPSAD_Dx_cell3.csv.gz')
lgb=lgb[lgb$subclass!='Adaptive',]
Obs_raw=read_obsvar('/sc/arion/projects/psychAD/NPS-AD/freeze2_proc/240124_PsychAD_freeze3_FULL_clean.h5ad')
rownames(Obs_raw)=Obs_raw$barcodekey
Obs_raw$barcodekey=NULL
Nick=do.call(rbind,lapply(c('Astro','Immune','Oligo','OPC','EN','IN'),function(cell){
  read_obsvar(sprintf('/sc/arion/projects/psychAD/NPS-AD/freeze2_proc/PsychADxD/7_trajectory/processed_data/model_output/redo/%s_barcode.h5ad',cell))
})) ##### PsychAD deep learning cell score results
rownames(Nick)=Nick$barcodekey
Nick$barcodekey=NULL
#lgb=lgb[lgb$SubID%in%Nick$SubID,]
#
metaInfo=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/metadata/clinical_metadata.csv')
metaInfo$group=c('RUSH','MSSM nonoverlap','MSSM overlap')[(metaInfo$Brain_bank=='MSSM')+(rownames(metaInfo)%in%MInfo$SubID)+1]
f=c(meta2$DxPC1,meta1$DxPC1)
names(f)=c(meta2$SubID,meta1$SubID)
metaInfo$Dx_donor=metaInfo$DxPC1=f[rownames(metaInfo)]
Obs=Obs_raw[Obs_raw$SubID%in%rownames(metaInfo),]
metaInfo$Cognitive_ResilienceUse=metaInfo$Cognitive_Resilience
Dx_cutoff=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/metadata/psychAD_meta_cutoff.csv',row.names=1)

metaInfo$Cognitive_ResilienceUse[metaInfo$DxPC1<Dx_cutoff$x[1]]=NA

########## cell composition change #################
#########neuron loss #####################################
Neu=Obs[Obs$class%in%c('EN'),]
ncells=table(Neu$SubID,Neu$subclass)
ncells=apply(ncells,2,c)
df=data.frame(vul=ncells[,c('EN_L2_3_IT')],
              nvul=rowSums(ncells)-ncells[,'EN_L2_3_IT'])
cobj=crumblr(df)
vul1=cobj
######
Neu=Obs[Obs$class%in%c('IN'),]
ncells=table(Neu$SubID,Neu$subclass)
ncells=apply(ncells,2,c)
df=data.frame(vul=ncells[,c('IN_SST')],
              nvul=rowSums(ncells)-ncells[,'IN_SST'])
cobj=crumblr(df)
vul2=cobj
###########
Neu=Obs[Obs$class%in%c('EN','IN'),]
ncells=table(Neu$SubID,Neu$subclass)
ncells=apply(ncells,2,c)
df=data.frame(vul=rowSums(ncells[,c('EN_L2_3_IT','IN_SST')]),
              nvul=rowSums(ncells)-rowSums(ncells[,c('EN_L2_3_IT','IN_SST')]))
cobj=crumblr(df)
vul=cobj
sh=intersect(colnames(vul1),colnames(vul2))
vul_df=zdf=data.frame(EN=(vul1$E['vul',sh]),IN=(vul2$E['vul',sh]),
                      Neu=(vul$E['vul',sh]))
cor(zdf,use='complete.obs',method='spearman')
Vul=new("EList")
Vul$E <- rbind(vul$E['vul',sh],
               vul1$E['vul',sh],
               vul2$E['vul',sh])
Vul$weights <- rbind(vul$weights['vul',sh],
                     vul1$weights['vul',sh],
                     vul2$weights['vul',sh])
rownames(Vul$E)=rownames(Vul$weights)=c('Neu','EN','IN')
####################subtypes#############
Obs$celltype=Obs$subtype
Obs[Obs$class%in%c('EN','IN','Mural'),]$celltype=Obs[Obs$class%in%c('EN','IN','Mural'),]$subclass
ncells=lapply(split(Obs,Obs$class),function(x)table(x$SubID,x$celltype))
ncells=lapply(ncells,function(x)apply(x,2,c))
z1=apply(table(Obs$SubID,Obs$class),2,c)
z2=apply(table(Obs$SubID,Obs$subclass),2,c)
z3=apply(table(Obs[Obs$class%in%c('EN','IN'),]$SubID,Obs[Obs$class%in%c('EN','IN'),]$subclass),2,c)
ncells$class=z1
ncells$subclass=z2
ncells$neu=z3
#### PsychAD Microglia annotation from https://www.nature.com/articles/s41588-026-02716-6
Micro_anno=read_obsvar('/sc/arion/projects/CommonMind/roussp01a/snmulti/external/230905_PsychAD_freeze3_MicroPVM_ADAM_taxonomy_final.h5ad')
rownames(Micro_anno)=Micro_anno$barcodekey
Micro_anno$barcodekey=NULL
z3=apply(table(Micro_anno$SubID,Micro_anno$subclass),2,c)
z4=apply(table(Micro_anno$SubID,Micro_anno$subtype),2,c)
ncells$Micro_subclass=z3
ncells$Micro_subtype=z4
cobjs=lapply(ncells,crumblr)

##################analysis##############################
#################Dx_dev#################
Dx_devs <- list()
Dx_dev_adjs <- list()
Vp=data.frame()
for(cell in unique(lgb$class)){
  lgb1=lgb[lgb$class==cell,]
  agg1=as.data.frame(lgb1%>%group_by(SubID)%>%summarise(Dx_cell=mean(Dx_cell)))
  rownames(agg1)=agg1$order;agg1$order=NULL
  
  nick=Nick[Nick$class==cell,]
  agg2=as.data.frame(nick%>%group_by(SubID)%>%summarise(pred_BRAAK_AD=mean(pred_BRAAK_AD),
                                                        pred_Dementia=mean(pred_Dementia),
                                                        BRAAK_AD=mean(BRAAK_AD),
                                                        CERAD=mean(CERAD)))
  rownames(agg1)=agg1$SubID
  rownames(agg2)=agg2$SubID
  
  zdf=inner_join(agg1,agg2,by='SubID')
  agg=zdf[,c('SubID','pred_BRAAK_AD','pred_Dementia','Dx_cell')]
  rownames(agg)=agg$SubID
  agg$SubID=NULL
  covs=c('n_genes','mito_genes','ribo_genes','mito_ribo')
  
  covInfo=as.data.frame(Obs[rownames(lgb1),c(covs,'SubID')]%>%group_by(SubID)%>%summarise(n_genes=mean(n_genes),
                                                                                          mito_genes=mean(mito_genes),
                                                                                          ribo_genes=mean(ribo_genes),
                                                                                          mito_ribo=mean(mito_ribo)))
  rownames(covInfo)=covInfo$SubID
  zInfo=data.frame(covInfo[rownames(agg),],metaInfo[rownames(agg),])
  form1=~scale(Age)+Sex+scale(CERAD)+scale(BRAAK_AD)+scale(PMI)+Brain_bank
  form=~scale(Age)+Sex+scale(CERAD)+scale(BRAAK_AD)+scale(PMI)+Brain_bank+scale(n_genes)+scale(mito_genes)+scale(ribo_genes)+scale(mito_ribo)
  vform=~scale(Age)+(1|Sex)+(1|ApoE_gt)+scale(prs_scaled_AD_Bellenguez)+scale(CERAD)+scale(BRAAK_AD)+scale(PMI)+(1|Brain_bank)+scale(n_genes)+scale(mito_genes)+scale(ribo_genes)+scale(mito_ribo)
  vp=fitExtractVarPartModel(t(agg),vform,zInfo)
  vp=data.frame(vp,ID=rownames(vp),assay=cell)
  Vp=rbind(Vp,vp)
  Dx_devs[[cell]]=t(residuals(dream(t(agg),form1,zInfo)))
  Dx_dev_adjs[[cell]]=t(residuals(dream(t(agg),form,zInfo)))
}
f=unique(unlist(lapply(Dx_devs,rownames)))
Dx_dev=lapply(colnames(Dx_devs[[1]]),function(col){data.frame(do.call(cbind,lapply(Dx_devs,function(x)data.frame(x)[f,col])),row.names = f)})
names(Dx_dev)=colnames(Dx_devs[[1]])
Dx_dev_adj=lapply(colnames(Dx_dev_adjs[[1]]),function(col){data.frame(do.call(cbind,lapply(Dx_dev_adjs,function(x)data.frame(x)[f,col])),row.names = f)})
names(Dx_dev_adj)=colnames(Dx_dev_adjs[[1]])
write.csv(Dx_dev$Dx_cell,'external/NPSAD_Dx_dev.csv')
write.csv(Dx_dev$pred_BRAAK_AD,'external/NPSAD_Dx_dev_NickBraak.csv')
write.csv(Dx_dev$pred_Dementia,'external/NPSAD_Dx_dev_NickDementia.csv')
write.csv(Dx_dev_adj$Dx_cell,'external/NPSAD_Dx_dev_adj.csv')
write.csv(Dx_dev_adj$pred_BRAAK_AD,'external/NPSAD_Dx_dev_NickBraak_adj.csv')
write.csv(Dx_dev_adj$pred_Dementia,'external/NPSAD_Dx_dev_NickDementia_adj.csv')
Dx_dev=Dx_dev_adj
f=table(unlist(lapply(Dx_dev,function(x)rownames(x))))
f=names(f)[f==3]
zdf=do.call(cbind,lapply(Dx_dev,function(x)rowMeans(scale(x),na.rm=T)[f]))
zdf=data.frame(zdf,group=metaInfo[rownames(zdf),]$group)
zdf1=reshape2::melt(zdf,id.vars=c('Dx_cell','group'))
zdf=zdf%>%group_by(group)%>%summarise(rho1=cor(Dx_cell,pred_BRAAK_AD,use='complete.obs',method='spearman'),rho2=cor(Dx_cell,pred_Dementia,use='complete.obs',method='spearman'))
zdf=reshape2::melt(zdf,id.vars='group')
print(zdf)
zdf=zdf[zdf$group!='RUSH',]
pdf('external/Dx_dev_vs_NPSADmodel.pdf',width=5.5,height = 5.5)
print(ggplot(zdf1[zdf1$group!='RUSH',])+
        geom_point(aes(x=Dx_cell,y=value,color=group))+
        geom_smooth(aes(x=Dx_cell,y=value,color=group),method='lm',se=F)+facet_wrap(~variable)+
        theme_bw()+
        scale_color_manual(values=c("#E41A1C","#377EB8"))+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
print(ggplot(zdf)+geom_bar(aes(y=variable,x=value,fill=group),width=0.7,
                           stat='identity',position='dodge')+
        theme_bw()+
        scale_fill_manual(values=c("#E41A1C","#377EB8"))+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
dev.off()
zdf=data.frame(M=rowMeans(scale(MDx_dev),na.rm=T),MInfo[rownames(MDx_dev),c('SubID','BrainRegion2')])
zdf=zdf[!duplicated(zdf[,c('SubID','BrainRegion2')]),]
zdf=data.frame(reshape2::acast(zdf,SubID~BrainRegion2,value.var='M'))
zdf$NPSAD=rowMeans(Dx_dev$Dx_cell)[rownames(zdf)]
zdf$pred_BRAAK=rowMeans(Dx_dev$pred_BRAAK_AD)[rownames(zdf)]
zdf$pred_Dementia=rowMeans(Dx_dev$pred_Dementia)[rownames(zdf)]
zdf=zdf[order(rowMeans(zdf[,1:3],na.rm=T)),]

Mord=rownames(zdf)
pdf('external/Dx_dev_vsNPSAD.pdf',width=6,height = 6)
cmat=cor(zdf,use='complete.obs',method='spearman')
corrplot::corrplot.mixed(cmat[c('pred_Dementia','pred_BRAAK','NPSAD','PHG','STG','PFC'),c('pred_Dementia','pred_BRAAK','NPSAD','PHG','STG','PFC')],
                         upper.col = rev(corrplot::COL2()),upper='square')
corrplot::corrplot(cmat[c('pred_Dementia','pred_BRAAK','NPSAD','PHG','STG','PFC'),c('pred_Dementia','pred_BRAAK','NPSAD','PHG','STG','PFC')],
                   col = rev(corrplot::COL2()),type = 'upper',method  = 'square',diag = FALSE)
corrplot::corrplot(cmat[c('NPSAD','PHG','STG','PFC'),c('NPSAD','PHG','STG','PFC')],
                   col = rev(corrplot::COL2()),type = 'upper',method  = 'square',diag = FALSE)
zdf1=pivot_longer(data.frame(scale(zdf)),cols=c('pred_Dementia','pred_BRAAK','NPSAD','PHG','STG'))
zdf1$Var2=factor(zdf1$name,levels = c('pred_Dementia','pred_BRAAK','NPSAD','PHG','STG'))

zdf1=pivot_longer(data.frame(scale(zdf)),cols=c('pred_Dementia','pred_BRAAK','NPSAD','PHG','STG'))

zdf1$Var2=factor(zdf1$name,levels = c('pred_Dementia','pred_BRAAK','NPSAD','PHG','STG'))
print(ggplot(zdf1)+
        geom_point(aes(x=PFC,y=value,color=Var2,color=Var2))+
        scale_color_manual(values=c("#E41A1C","#253494","#1B9E77","#7570B3","#D95F02"))+
        theme_bw()+facet_wrap(~grepl('pred',name))+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
print(ggplot(zdf1)+
        geom_point(aes(x=PFC,y=value,color=Var2,color=Var2),alpha=0.5)+
        geom_smooth(aes(x=PFC,y=value,color=Var2,color=Var2),method='lm',se=F)+
        scale_color_manual(values=c("#E41A1C","#253494","#1B9E77","#7570B3","#D95F02"))+
        theme_bw()+facet_wrap(~grepl('pred',name),scales='free')+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),legend.position = 'bottom',
              panel.grid.minor = element_blank(),aspect.ratio = 1))
print(ggplot(zdf1[!grepl('pred',zdf1$name),])+
        geom_point(aes(x=PFC,y=value,color=Var2,color=Var2),alpha=0.5)+
        geom_smooth(aes(x=PFC,y=value,color=Var2,color=Var2),method='lm',se=F)+
        scale_color_manual(values=c("#1B9E77","#7570B3","#D95F02"))+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
dev.off()
################projection. ##################################
F_sel=paste0('Factor',1:3)
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/MSSM/class/DxPC1.Rdata')
NcovInfo=do.call(rbind,lapply(c('Astro','Oligo','OPC','Immune','EN','IN'),function(cell){z=fit[[cell]]$data;
data.frame(sample=rownames(z),assay=cell,scale(z[,c('n_genes','mito_genes','ribo_genes')]))}))
NcovInfos=lapply(list('EN','IN',c("EN",'IN'),c("Astro",'Oligo','OPC','Immune'),c("Astro",'Oligo','OPC','Immune',"EN",'IN')),function(x){
  z=NcovInfo[NcovInfo$assay%in%x,]
  z=data.frame(z%>%pivot_wider(names_from=assay,values_from=c(n_genes,ribo_genes,mito_genes)))
  rownames(z)=z[,1]
  z[,1]=NULL
  z[rowSums(is.na(z))==0,]
})
names(NcovInfos)=c('EN','IN','Neu','Gli','All')


NGexs=lapply(c('Astro','Oligo','OPC','Immune','EN','IN'),function(cell){
  R=residuals(fit[[cell]])
  rownames(R)=paste0(gsub('Immune','Micro_PVM',cell),":",rownames(R))
  R
})
f=table(unlist(lapply(NGexs,colnames)))
f=names(f)[f==length(NGexs)]
NGex=do.call(rbind,lapply(NGexs,function(x)x[,f]))
load('MOFA_top_25_10_agg.Rdata')
####
W=agg$consensus_weights
Fs=agg$consensus_scores
sh=intersect(rownames(NGex),rownames(W))
NFp=scale(t(NGex[sh,])%*%W[sh,])
sh=intersect(rownames(NFp),rownames(Fs))
z=do.call(rbind,strsplit(rownames(W),'\\:'))
Wgli=W[!z[,1]%in%c("EN",'IN'),]
sh=intersect(rownames(NGex),rownames(Wgli))
NFgli=scale(t(NGex[sh,])%*%W[sh,])

Wgli2=W[!z[,1]%in%c("EN",'IN','Micro_PVM'),]
sh=intersect(rownames(NGex),rownames(Wgli2))
NFgli2=scale(t(NGex[sh,])%*%W[sh,])


write.csv(NFp,'external/NPSAD_Factors.csv')
write.csv(NFgli,'external/NPSAD_Factors_glia.csv')
write.csv(NFgli2,'external/NPSAD_Factors_glia_nomicro.csv')
######################resilience#####################
z=rowMeans(scale(Dx_dev$Dx_cell),na.rm=T)
zInfo=data.frame(Dx_cell=rowMeans(scale(agg1),na.rm=T)[names(z)],Dx_dev=z,metaInfo[names(z),])
sel1=intersect(rownames(agg1),rownames(zInfo[zInfo$group=='MSSM overlap',]))
sel2=intersect(rownames(agg1),rownames(zInfo[zInfo$group=='MSSM nonoverlap',]))
form=~scale(Age)+Sex+scale(DxPC1)+scale(PMI)+scale(Cognitive_ResilienceUse)+(1|ApoE_gt)
zfit1=eBayes(dream(t(scale(Dx_dev$Dx_cell)),form,zInfo[rownames(Dx_dev$Dx_cell),])) # all samples
z1=topTable(zfit1,coef='scale(Cognitive_ResilienceUse)') 
zfit2=eBayes(dream(t(scale(Dx_dev$Dx_cell[sel1,])),form,zInfo[sel1,])) ### overlap
z2=topTable(zfit2,coef='scale(Cognitive_ResilienceUse)') 
zfit3=eBayes(dream(t(scale(Dx_dev$Dx_cell[sel2,])),form,zInfo[sel2,])) #### not overlap
z3=topTable(zfit3,coef='scale(Cognitive_ResilienceUse)')
zdf=data.frame(z1,cell=rownames(z1),se=z1$logFC/z1$t)
zdf$sig=c('','+')[(zdf$adj.P.Val<0.05)+1]
zdf$cell=factor(zdf$cell,levels=c('Astro','Oligo','OPC','Immune','EN','IN'))

print(ggplot(zdf)+
        geom_linerange(aes(x=cell,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
        geom_point(aes(x=cell,y=logFC,size=-log10(adj.P.Val)))+
        geom_text(aes(x=cell,y=logFC,label=sig),color='white')+
        #scale_color_manual(values=colormap_class)+
        scale_size(limits=c(1,max(-log10(zdf$adj.P.Val))))+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))
zdf1=zdf;zdf1$group='all';zdf1$sig=NULL
zdf=data.frame(rbind(z2,z3),cell=c(rownames(z2),rownames(z3)),
               se=c(z2$logFC/z2$t,z3$logFC/z3$t),group=rep(c('overlap','no'),each=nrow(z2)))
write.csv(zdf,file='external/NPSAD_Dx_dev_resilience.csv')
##########
zInfo$Dx_dev=rowMeans(scale(Dx_dev$Dx_cell),na.rm=T)[rownames(zInfo)]
meta2=read.csv('/sc/arion/projects/psychAD/NPS-AD/freeze2_rc/metadata/clinical_metadata_full.csv',row.names=1)
zInfo$resilience=c('other','resilience')[(rownames(zInfo)%in%rownames(NPSADcontrasts$metadata[NPSADcontrasts$metadata$c28x=='AD_resilient',]))+1]#
zdf=zInfo[which(zInfo$DxPC1>Dx_cutoff$x[1] & zInfo$Brain_bank=='MSSM'),] #
pdf('NPSAD_resilience.pdf')
print(ggplot(zdf)+geom_smooth(aes(x=DxPC1,y=Dx_cell,color=resilience))+
        geom_vline(xintercept = Dx_cutoff$x,linetype=2)+
        scale_color_manual(values=rev(c("#E41A1C","#377EB8")))+
        theme_bw()+facet_wrap(~group)+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),legend.position = 'bottom',
              panel.grid.minor = element_blank(),aspect.ratio = 1))
print(ggplot(zdf)+geom_smooth(aes(x=DxPC1,y=Dx_cell,color=resilience))+
        geom_vline(xintercept = Dx_cutoff$x,linetype=2)+
        scale_color_manual(values=rev(c("#E41A1C","#377EB8")))+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
dev.off()

################subclass#####################
sh=intersect(colnames(cobjs$subclass),rownames(NFgli))
sh2=setdiff(sh,MInfo$SubID)
cobj=cobjs$subclass[,sh]
NFg=data.frame(NFgli)
zInfo=data.frame(metaInfo[sh,],data.frame(NFg)[sh,])
zInfo$Dx_dev=rowMeans(scale(Dx_dev$Dx_cell[,c('Astro','Oligo','OPC','Immune')]),na.rm=T)[rownames(zInfo)]

form1=~scale(Age)+Sex+scale(DxPC1)+scale(PMI)
zfit1=eBayes(dream(cobj,form1,zInfo))
z1=topTable(zfit1,'scale(DxPC1)',number=Inf)
zfit1=eBayes(dream(cobj[,sh2],form1,zInfo[sh2,]))
z1b=topTable(zfit1,'scale(DxPC1)',number=Inf)
zdf=data.frame(rbind(z1,z1b),cell=c(rownames(z1),rownames(z1b)),
               se=c(z1$logFC/z1$t,z1b$logFC/z1b$t),group=rep(c('all','no'),each=nrow(z1)))
write.csv(zdf,file='external/NPSAD_Dx_donor_subclass_crumblr.csv')

form2=~scale(Age)+Sex+scale(DxPC1)+scale(PMI)+scale(Dx_dev)
zfit2=eBayes(dream(cobj,form2,zInfo))
z2=topTable(zfit2,'scale(Dx_dev)',number=Inf)
zfit2=eBayes(dream(cobj[,sh2],form2,zInfo[sh2,]))
z2b=topTable(zfit2,'scale(Dx_dev)',number=Inf)
zdf=data.frame(rbind(z2,z2b),cell=c(rownames(z2),rownames(z2b)),
               se=c(z2$logFC/z2$t,z2b$logFC/z2b$t),group=rep(c('all','no'),each=nrow(z2)))
write.csv(zdf,file='external/NPSAD_Dx_dev_subclass_crumblr.csv')

z1$se=z1$logFC/z1$t
z2$se=z2$logFC/z2$t

sh2=setdiff(sh,MInfo$SubID)
form2=~scale(Age)+Sex+scale(DxPC1)+scale(PMI)+scale(Dx_dev)
zfit2=eBayes(dream(cobj[,sh2],form2,zInfo[sh2,]))

z3=topTable(zfit2,'scale(Dx_dev)',number=Inf)
z3$se=z3$logFC/z3$t

zdf=data.frame(rbind(z2,z3),cell=c(rownames(z2),rownames(z3)),
               se=c(z2$logFC/z2$t,z3$logFC/z3$t),group=rep(c('all','no'),each=nrow(z2)))
write.csv(zdf,file='external/NPSAD_Dx_dev_subclass_crumblr.csv')
###

for(xx in c('Factor1',"Factor2",'Factor3')){
  form=as.formula(sprintf('~scale(Age)+Sex+scale(DxPC1)+scale(PMI)+scale(%s)',xx))
  
  zfit=eBayes(dream(cobj,form,zInfo))
  z2=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  z2$se=z2$logFC/z2$t
  zfit2=eBayes(dream(cobj[,sh2],form,zInfo[sh2,]))
  
  z3=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  z3$se=z3$logFC/z3$t
  zdf=data.frame(rbind(z2,z3),cell=c(rownames(z2),rownames(z3)),
                 se=c(z2$logFC/z2$t,z3$logFC/z3$t),group=rep(c('all','no'),each=nrow(z2)))
  write.csv(zdf,file=sprintf('external/NPSAD_%s_subclass_crumblr.csv',xx))
}
dev.off()
############################vulnerbale neurons######
sh=intersect(colnames(Vul),rownames(NFgli))
cobj=Vul[,sh]
zInfo=data.frame(metaInfo[sh,],data.frame(NFg)[sh,])
zInfo$Dx_dev=rowMeans(scale(Dx_dev$Dx_cell[,c('Astro','Oligo','OPC','Immune')]),na.rm=T)[rownames(zInfo)]
zInfo$Dx_cell=rowMeans(scale(agg1[,c('Astro','Oligo','OPC','Immune')]),na.rm=T)[rownames(zInfo)]
sh2=setdiff(sh,MInfo$SubID)
form0=~scale(Age)+Sex+scale(Dx_cell)+scale(PMI)
zfit=eBayes(dream(cobj,form0,zInfo))
z0=topTable(zfit,'scale(Dx_cell)',number=Inf)
zfit=eBayes(dream(cobj[,sh2],form0,zInfo[sh2,]))
z1=topTable(zfit,'scale(Dx_cell)',number=Inf)
zdf=data.frame(rbind(z0,z1),cell=c(rownames(z0),rownames(z1)),
               se=c(z0$logFC/z0$t,z1$logFC/z1$t),group=rep(c('all','no'),each=nrow(z1)))
write.csv(zdf,file='external/NPSAD_Dx_cell_vulneuron_crumblr.csv')
form0=~scale(Age)+Sex+scale(DxPC1)+scale(PMI)
zfit=eBayes(dream(cobj,form0,zInfo))
z1=topTable(zfit,'scale(DxPC1)',number=Inf)
zfit=eBayes(dream(cobj[,sh2],form0,zInfo[sh2,]))
z1b=topTable(zfit,'scale(DxPC1)',number=Inf)
zdf=data.frame(rbind(z1,z1b),cell=c(rownames(z1),rownames(z1b)),
               se=c(z1$logFC/z1$t,z1b$logFC/z1b$t),group=rep(c('all','no'),each=nrow(z1)))
write.csv(zdf,file='external/NPSAD_Dx_donor_vulneuron_crumblr.csv')
z1$se=z1$logFC/z1$t
zdf1=data.frame(cell=rownames(z1),logFC=z1$logFC,se=z1$se,FDR=z1$adj.P.Val,
                group='Dx_donor',use='all')
df=zdf1
pdf('external/NPSAD_Vul_neuron_crumblr.pdf',width=5,height = 5)
for(xx in c('Dx_dev','Factor1','Factor2','Factor3')){
  form1=as.formula(sprintf('~scale(Age)+Sex+scale(DxPC1)+scale(PMI)+scale(%s)',xx))  
  zfit=eBayes(dream(cobj,form1,zInfo))
  z2=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  #write.csv(z2,file=sprintf('%s_vulneuron_crumblr.csv',xx))
  z2$se=z2$logFC/z2$t
  zfit=eBayes(dream(cobj[,sh2],form1,zInfo[sh2,]))
  z3=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  #write.csv(z3,file=sprintf('%s_vulneuron_crumblr.csv',xx))
  z3$se=z3$logFC/z3$t
  
  zdf=data.frame(rbind(z2,z3),cell=c(rownames(z2),rownames(z3)),
                 se=c(z2$logFC/z2$t,z3$logFC/z3$t),group=rep(c('all','no'),each=nrow(z2)))
  write.csv(zdf,file=sprintf('external/NPSAD_%s_vulneuron_crumblr.csv',xx))
}
dev.off()
#############DAM######################################
NFg2=data.frame(NFgli2)

sh=intersect(colnames(cobjs$Micro_subclass),rownames(NFgli2))
cobj=cobjs$Micro_subclass[,sh]
zInfo=data.frame(metaInfo[sh,],data.frame(NFg2)[sh,])
zInfo$Dx_dev=rowMeans(scale(Dx_dev$Dx_cell[,c('Astro','Oligo','OPC')]),na.rm=T)[rownames(zInfo)]
zInfo$Dx_cell=rowMeans(scale(agg1[,c('Astro','Oligo','OPC')]),na.rm=T)[rownames(zInfo)]
sh2=setdiff(sh,MInfo$SubID)

form0=~scale(Age)+Sex+scale(Dx_cell)+scale(PMI)
zfit=eBayes(dream(cobj,form0,zInfo))
z0=topTable(zfit,'scale(Dx_cell)',number=Inf)
zfit=eBayes(dream(cobj[,sh2],form0,zInfo[sh2,]))
z0b=topTable(zfit,'scale(Dx_cell)',number=Inf)
zdf=data.frame(rbind(z0,z0b),cell=c(rownames(z0),rownames(z0b)),
               se=c(z0$logFC/z0$t,z0b$logFC/z0b$t),group=rep(c('all','no'),each=nrow(z0)))
write.csv(zdf,file='external/NPSAD_Dx_cell_DAM_crumblr.csv')

form0=~scale(Age)+Sex+scale(DxPC1)+scale(PMI)
zfit=eBayes(dream(cobj,form0,zInfo))
z1=topTable(zfit,'scale(DxPC1)',number=Inf)
zfit=eBayes(dream(cobj[,sh2],form0,zInfo[sh2,]))
z1b=topTable(zfit,'scale(DxPC1)',number=Inf)
zdf=data.frame(rbind(z1,z1b),cell=c(rownames(z1),rownames(z1b)),
               se=c(z1$logFC/z1$t,z1b$logFC/z1b$t),group=rep(c('all','no'),each=nrow(z1)))
write.csv(zdf,file='external/NPSAD_Dx_donor_DAM_crumblr.csv')
z1$se=z1$logFC/z1$t
zdf1=data.frame(cell=rownames(z1),logFC=z1$logFC,se=z1$se,FDR=z1$adj.P.Val,
                group='Dx_donor',use='all')

for(xx in c('Dx_dev')){
  form1=as.formula(sprintf('~scale(Age)+Sex+scale(DxPC1)+scale(PMI)+scale(%s)',xx))  
  zfit=eBayes(dream(cobj,form1,zInfo))
  z2=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  #write.csv(z2,file=sprintf('%s_vulneuron_crumblr.csv',xx))
  z2$se=z2$logFC/z2$t

  zfit=eBayes(dream(cobj[,sh2],form1,zInfo[sh2,]))
  z3=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  #write.csv(z3,file=sprintf('%s_vulneuron_crumblr.csv',xx))
  z3$se=z3$logFC/z3$t
  
  zdf=data.frame(rbind(z2,z3),cell=c(rownames(z2),rownames(z3)),
                 se=c(z2$logFC/z2$t,z3$logFC/z3$t),group=rep(c('all','no'),each=nrow(z2)))
  write.csv(zdf,file=sprintf('external/NPSAD_%s_DAM_crumblr.csv',xx))
}

