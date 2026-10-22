library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(crumblr)
library(ggrepel)
library(readxl)
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
args=commandArgs(T)
lgb_path=args[1]
obs_path=args[2]
sccoda_res_path=args[3]
lgb=read_big_csv(lgb_path)
Obs=read_big_csv(obs_path)

Info=unique(lgb[,c('Donor ID','Continuous Pseudo-progression Score',
                    'Sex','Age at Death','PMI','Cognitive Status','APOE Genotype','Severely Affected Donor','CERAD score','Braak')])
colnames(Info)=c('SubID','CPS','Sex','Age','PMI','Cogn','apoe','SeverDonor','CERAD','BRAAK_AD')
Info$BRAAK_AD=as.numeric(factor(Info$BRAAK_AD,levels=c("Braak 0","Braak II","Braak III","Braak IV","Braak V","Braak VI"),
                                ordered=T))
Info$CERAD=as.numeric(factor(Info[,'CERAD'],levels=c("Absent","Sparse","Moderate","Frequent"),
                             ordered=T))
Info$PMI=as.numeric(Info$PMI)
rownames(Info)=Info$SubID
metaInfo=Info
##############
ncell=apply(table(Obs[Obs$class%in%c('IN',"EN"),]$`Donor ID`,Obs[Obs$class%in%c('IN',"EN"),]$Supertype),2,c)
cobj=crumblr(ncell)
##########
Dx_dev <- list()
Obs$SubID=Obs$`Donor ID`
lgb$SubID=lgb$`Donor ID`
for(cell in unique(lgb$class)){
  lgb1=lgb[lgb$class==cell,]
  agg=as.data.frame(lgb1%>%group_by(SubID)%>%summarise(Dx_cell=mean(Dx_cell)))
  rownames(agg)=agg$SubID
  agg$SubID=NULL
  
  covs=c('n_genes_by_counts','pct_counts_mt','pct_counts_ribo')
  covInfo=as.data.frame(Obs[rownames(lgb1),c(covs,'SubID')]%>%group_by(SubID)%>%summarise(n_genes_by_counts=mean(n_genes_by_counts),
                                                                                            pct_counts_mt=mean(pct_counts_mt),
                                                                                            pct_counts_ribo=mean(pct_counts_ribo)))
  rownames(covInfo)=covInfo$SubID
  zInfo=data.frame(covInfo[rownames(agg),],metaInfo[rownames(agg),])
  form=~scale(Age)+Sex+scale(CPS)+scale(PMI)+scale(n_genes_by_counts)+scale(pct_counts_mt)+scale(pct_counts_ribo)
  Dx_dev[[cell]]=t(residuals(dream(t(agg),form,zInfo)))
}
f=unique(unlist(lapply(Dx_dev,rownames)))
Dx_dev=do.call(cbind,lapply(Dx_dev,function(x)data.frame(x)[f,1]))
rownames(Dx_dev)=f
write.csv(Dx_dev,file='SEAAD-MTG-Dx_dev.csv')
sh=intersect(rownames(metaInfo),rownames(Dx_dev))
zInfo=metaInfo[sh,]
zInfo$Dx_dev=rowMeans(scale(Dx_dev[,c('Astro','Oligo','OPC','Micro_PVM')]),na.rm=T)[rownames(zInfo)]
####
form1=~scale(Age)+Sex+scale(CPS)+scale(PMI)
zfit1=eBayes(dream(cobj[,sh],form1,zInfo))
z1=topTable(zfit1,'scale(CPS)',number=Inf)
form2=~scale(Age)+Sex+scale(CPS)+scale(PMI)+scale(Dx_dev)
zfit2=eBayes(dream(cobj[,sh],form2,zInfo))
z2=topTable(zfit2,'scale(Dx_dev)',number=Inf)
write.csv(z1,file='SEAAD_CPS_subtype_crumblr.csv')
write.csv(z2,file='SEAAD_Dx_dev_subtype_crumblr.csv')
z1$se=z1$logFC/z1$t
z2$se=z2$logFC/z2$t
dat=read_excel(sccoda_res_path,'MTG_RNAseq Neu CPS')# (https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/Supplementary+Information/Supplementary+Table+4.xlsx)

dat <- dat %>%
  filter(Covariate == "Continuous_Pseudo-progression_Score") %>%
  transmute(
    reference = `Reference Cell Type`,
    cell = `Cell Type`,
    beta = as.numeric(`Final Parameter`),
    pip = as.numeric(`Inclusion probability`)
  )

# Check that all reference/cell combinations are present
stopifnot(!anyDuplicated(dat[c("reference", "cell")]))

# Average across reference models
res <- dat %>%
  group_by(cell) %>%
  summarise(
    mean_effect = mean(beta, na.rm = TRUE),
    sem = sd(beta, na.rm = TRUE) /
      sqrt(sum(!is.na(beta))),
    mean_pip = mean(pip, na.rm = TRUE),
    n_references = n_distinct(reference),
    .groups = "drop"
  ) %>%
  mutate(
    credible = mean_pip > 0.8
  )
write.csv(res, "SEAAD_scCODA_CPS_summary.csv",
          row.names = FALSE)
res=data.frame(res)
rownames(res)=res$cell
zdf=data.frame(cell=rownames(z1),
               Dx_donor=z1$logFC,donor_se=z1$se,
               Dx_dev=z2[rownames(z1),]$logFC,
               dev_se=z2[rownames(z1),]$se,
               ref=res[rownames(z1),]$mean_effect,
               ref_se=res[rownames(z1),]$sem,
               FDR=z2[rownames(z1),]$adj.P.Val)

zdf$subclass=gsub('_[0-9]+$','',zdf$cell)
zdf$cell2='other'
zdf$cell2[zdf$subclass%in%c('Sst','L2/3 IT')]=zdf$subclass[zdf$subclass%in%c('Sst','L2/3 IT')]
zdf$cell2=factor(zdf$cell2,c('Sst','L2/3 IT','other'))
pdf('SEAAD_MTG_neuron_cell_composition_change_Dx_dev_vs_CPS.pdf')
print(ggplot(zdf)+
        geom_smooth(aes(x=Dx_donor,y=Dx_dev),method='lm',linewidth = 0.7)+
        geom_point(aes(x=Dx_donor,y=Dx_dev,size=-log10(FDR),color=cell2))+
        geom_linerange(aes(x=Dx_donor,ymin=Dx_dev-dev_se,ymax=Dx_dev+dev_se,color=cell2),linewidth = 0.35)+
        geom_linerange(aes(y=Dx_dev,xmin=Dx_donor-donor_se,xmax=Dx_donor+donor_se,color=cell2),linewidth = 0.35)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        geom_vline(xintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('logFC (CPS)')+ylab('logFC (Dx_dev)')+
        ggtitle(paste0('rho ',round(cor(zdf$Dx_donor,zdf$Dx_dev,use='complete.obs',method='spearman'),2)))+
        scale_color_manual(values=c("#E41A1C","#377EB8",'#66666680'))+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
print(ggplot(zdf)+
        geom_smooth(aes(x=Dx_donor,y=Dx_dev),method='lm',linewidth = 0.7)+
        geom_point(aes(x=Dx_donor,y=Dx_dev,size=-log10(FDR),color=cell2))+
        geom_linerange(aes(x=Dx_donor,ymin=Dx_dev-dev_se,ymax=Dx_dev+dev_se,color=cell2),linewidth = 0.35)+
        geom_linerange(aes(y=Dx_dev,xmin=Dx_donor-donor_se,xmax=Dx_donor+donor_se,color=cell2),linewidth = 0.35)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        geom_vline(xintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('logFC (CPS)')+ylab('logFC (Dx_dev)')+
        ggtitle(paste0('rho ',round(cor(zdf$Dx_donor,zdf$Dx_dev,use='complete.obs',method='spearman'),2)))+
        geom_text_repel(aes(x=Dx_donor,y=Dx_dev,label=cell),
                        color = "black",seed = 123, box.padding = 0.45,
                        point.padding = 0.25,force = 2,force_pull = 0.1,
                        max.overlaps = Inf,
                        max.time = 5, min.segment.length = 0,
                        segment.color = "grey35",linewidth = 0.3)+
        scale_color_manual(values=c("#E41A1C","#377EB8",'#66666680'))+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
print(ggplot(zdf)+
        geom_smooth(aes(x=ref,y=Dx_donor),method='lm',linewidth = 0.7)+
        geom_point(aes(x=ref,y=Dx_donor,color=cell2))+
        geom_linerange(aes(x=ref,ymin=Dx_donor-donor_se,ymax=Dx_donor+donor_se,color=cell2),linewidth = 0.35)+
        geom_linerange(aes(y=Dx_donor,xmin=ref-ref_se,xmax=ref+ref_se,color=cell2),linewidth = 0.35)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        geom_vline(xintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('logFC (CPS) scCODA')+ylab('logFC (CPS) crumblr')+
        ggtitle(paste0('rho ',round(cor(zdf$Dx_donor,zdf$ref,use='complete.obs',method='spearman'),2)))+
        scale_color_manual(values=c("#E41A1C","#377EB8",'#66666680'))+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
print(ggplot(zdf)+
        geom_smooth(aes(x=ref,y=Dx_donor),method='lm',linewidth = 0.7)+
        geom_point(aes(x=ref,y=Dx_donor,color=cell2))+
        geom_linerange(aes(x=ref,ymin=Dx_donor-donor_se,ymax=Dx_donor+donor_se,color=cell2),linewidth = 0.35)+
        geom_linerange(aes(y=Dx_donor,xmin=ref-ref_se,xmax=ref+ref_se,color=cell2),linewidth = 0.35)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        geom_vline(xintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('logFC (CPS) scCODA')+ylab('logFC (CPS) crumblr')+
        ggtitle(paste0('rho ',round(cor(zdf$Dx_donor,zdf$ref,use='complete.obs',method='spearman'),2)))+
        geom_text_repel(aes(x=ref,y=Dx_donor,label=cell),
                        color = "black",seed = 123, box.padding = 0.45,
                        point.padding = 0.25,force = 2,force_pull = 0.1,
                        max.overlaps = Inf,
                        max.time = 5, min.segment.length = 0,
                        segment.color = "grey35",linewidth = 0.3)+
        scale_color_manual(values=c("#E41A1C","#377EB8",'#66666680'))+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
dev.off()