.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
options(stringsAsFactors = F)
library(data.table)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)
de_scale <- function(qu,ref){
  z=(qu-mean(qu,na.rm=T))/sd(qu,na.rm=T)*sd(ref,na.rm=T)+mean(ref,na.rm=T)
  return(z)
}
read_big_csv <- function(fh,row1asrowname=T){
  dat=fread(fh)
  dat=as.data.frame(dat)
  if(row1asrowname){
    rownames(dat)=dat[,1]
  }
  return(dat)
}
DAG=read_big_csv('gene_DAG_score.csv.gz')

setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ML_cellstate3')
bugn=colorRampPalette(RColorBrewer::brewer.pal(9,'BuGn'))
cell_anno=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno_05092025_fixed3.csv',row.names=1)
LGB=data.frame()
df1=data.frame()
df2=data.frame()
for(cell in c('Astro','EN','IN','Micro_PVM','Oligo','OPC')){
  
  lgb=read.csv(paste0(cell,'/lgb_predict_redo.csv'),row.names=1)
  sel=which(is.na(lgb$DxPC1_predict2))
  if(length(sel)>0){
    lgb[sel,'DxPC1_predict2']=lgb[sel,]$DxPC1
  }
  sel=which(is.na(lgb$DxPC1_clust3_predict2))
  if(length(sel)>0){
    lgb[sel,'DxPC1_clust3_predict2']=lgb[sel,]$DxPC1_clust3-1
  }
  #lgb=lgb[!is.na(lgb$DxPC1_predict2),]
  #lgb=lgb[!is.na(lgb$DxPC1_clust3_predict2),]
  lgb=data.frame(lgb,cell_anno[rownames(lgb),])
  
  lgb$Dx_cell=de_scale(lgb$DxPC1_predict2,lgb$DxPC1)
  lgb$Dx_donor=lgb$DxPC1
  lgb$Dx_cell_raw=lgb$DxPC1_predict2
  lgb$Dx_cat_cell=c('control','early','late')[lgb$DxPC1_clust3_predict2+1]
  lgb$Dx_cat_donor=c('control','early','late')[lgb$DxPC1_clust3]
  lgb$delta=lgb$Dx_cell-lgb$Dx_donor
  lgb$abs_delta=abs(lgb$delta)
  lgb$abs_delta_rankbin100=as.numeric(cut(lgb$abs_delta,
                                          c(-Inf,quantile(lgb$abs_delta,seq(0,1,0.01))[-1])))/100
  lgb$consistent=lgb$abs_delta < quantile(lgb$abs_delta,0.7)
  lgb$cell=cell
  LGB=rbind(LGB,lgb)
  pdf(paste0(cell,'_summary.pdf'))
  z=melt(lgb[,c('Dx_cell','Dx_donor','Dx_cat_cell','Dx_cat_donor')],
         value.name=c('Dx'),variable.name=c('Dx_type'))
  z=melt(z,id.var=c('Dx','Dx_type'),
         value.name=c('Dx_cat'),variable.name=c('Dx_cat_type'))
  print(ggplot(z)+geom_boxplot(aes(x=Dx_type,y=Dx,fill=Dx_cat))+
          theme_bw()+ylab('Dx')+facet_wrap(~Dx_cat_type)+
          theme(axis.text.x = element_text(colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 1))
  if(cell %in% c('IN','EN','Mural')){
    lgb$celltype=lgb$subclass
  }else{
    lgb$celltype=lgb$subtype
  }
  print(ggplot(lgb)+geom_boxplot(aes(x=celltype,y=Dx_cell,fill=Dx_cat_donor),position='dodge')+
          facet_grid(.~BrainRegion)+theme_bw()+
          theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 1))
  print(ggplot(lgb)+geom_boxplot(aes(x=celltype,y=Dx_cell,fill=BrainRegion),position='dodge')+
          facet_grid(.~Dx_cat_donor)+theme_bw()+
          theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 1))
  
  
  cell_deg_score=read_big_csv(paste0('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ML_cellstate3/',cell,'/DxPC1_predict2_DEG_module_score.csv.gz'))
  cell_deg_score$PsychAD_DE=scale(cell_deg_score$PsychAD_up)-scale(cell_deg_score$PsychAD_down)
  sh=intersect(rownames(lgb),rownames(cell_deg_score))
  zdf=data.frame(lgb[sh,c('Dx_cell','Dx_donor','abs_delta_rankbin100','consistent','BrainRegion','cell')],
                 cell_deg_score[sh,c('PsychAD_DE'),drop=F])
  zdf=as.data.frame(tidyr::pivot_longer(zdf,cols=c('Dx_cell','Dx_donor'),
                                        names_to = "Dx_type",values_to = "Dx"))
  zdf1=result <- zdf %>%
    group_by(consistent,Dx_type,cell) %>%
    summarise(Dx_cor = cor(Dx, PsychAD_DE, use = "complete.obs",method='spearman'), .groups = "drop")
  df1=rbind(df1,zdf1)
  zdf2=result <- zdf %>%
    group_by(abs_delta_rankbin100,Dx_type) %>%
    summarise(Dx_cor = cor(Dx, PsychAD_DE, use = "complete.obs",method='spearman'), .groups = "drop")
  
  print(ggplot(zdf2,aes(x=abs_delta_rankbin100,y=Dx_cor,color=Dx_type))+
          geom_point()+geom_smooth()+geom_vline(xintercept=70/100,linewidth=0.5,linetype = "longdash")+
          xlab('Delta bin (Dx_cell-Dx_donor)')+ylab('Cor(Dx score (PsychAD),Dx)')+
          ggtitle(cell)+theme_bw()+
          theme(axis.text.x = element_text(colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 1))
  zdf2<- zdf %>%
    group_by(abs_delta_rankbin100,Dx_type,BrainRegion) %>%
    summarise(Dx_cor = cor(Dx, PsychAD_DE, use = "complete.obs",method='spearman'), .groups = "drop")
  
  print(ggplot(zdf2,aes(x=abs_delta_rankbin100,y=Dx_cor,color=Dx_type))+
          geom_point()+geom_smooth()+geom_vline(xintercept=70/100,linewidth=0.5,linetype = "longdash")+
          xlab('Delta bin (Dx_cell-Dx_donor)')+ylab('Cor(Dx score (PsychAD),Dx)')+
          ggtitle(cell)+theme_bw()+facet_wrap(~BrainRegion)+
          theme(axis.text.x = element_text(colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                legend.position = 'bottom',
                panel.grid.minor = element_blank(),aspect.ratio = 1))
  zdf2<- zdf %>%
    group_by(Dx_type,consistent) %>%
    summarise(Dx_cor = cor(Dx, PsychAD_DE, use = "complete.obs",method='spearman'), 
              .groups = "drop")
  xs=quantile(zdf$Dx,c(0.01,0.99))
  ys=max(abs(quantile(zdf$PsychAD_DE,c(0.01,0.99),na.rm=T)))*c(-1,1)
  a=min(xs)*0.8
  b=max(ys)*0.9
  print(ggplot(zdf,aes(x=Dx,y=PsychAD_DE))+
          geom_density_2d_filled(contour_var = "ndensity") + 
          geom_smooth(method='lm')+xlim(xs[1],xs[2])+ylim(ys[1],ys[2])+
          geom_text(data=zdf2,x=a,y=b,aes(label=paste0('SCC :',round(Dx_cor,3))))+
          facet_grid(consistent~Dx_type)+
          ggtitle(cell)+theme_bw()+scale_fill_manual(values=bugn(10))+
          theme(axis.text.x = element_text(colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 1))

  cell_dac_score=DAG[DAG$cell==cell,]  
  cell_dac_score$DAG=scale(cell_dac_score$DAG_up)-scale(cell_dac_score$DAG_down)
  sh=intersect(rownames(lgb),rownames(cell_dac_score))
  zdf=data.frame(lgb[sh,c('Dx_cell','Dx_donor','abs_delta_rankbin100','consistent','BrainRegion','cell')],
                 cell_dac_score[sh,c('DAG'),drop=F])
  zdf=as.data.frame(tidyr::pivot_longer(zdf,cols=c('Dx_cell','Dx_donor'),
                                        names_to = "Dx_type",values_to = "Dx"))
  
  zdf2<- zdf %>%
    group_by(consistent,Dx_type,cell) %>%
    summarise(Dx_cor = cor(Dx, DAG, use = "complete.obs",method='spearman'), 
              .groups = "drop")
  df2=rbind(df2,zdf2)
  zdf2<- zdf %>%
    group_by(abs_delta_rankbin100,Dx_type) %>%
    summarise(Dx_cor = cor(Dx, DAG, use = "complete.obs",method='spearman'), 
              .groups = "drop")
  
  print(ggplot(zdf2,aes(x=abs_delta_rankbin100,y=Dx_cor,color=Dx_type))+
          geom_point()+geom_smooth()+geom_vline(xintercept=70/100,linewidth=0.5,linetype = "longdash")+
          xlab('Delta bin (Dx_cell-Dx_donor)')+ylab('Cor(Dx score (DAG),Dx)')+
          ggtitle(cell)+theme_bw()+
          theme(axis.text.x = element_text(colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 1))
  zdf2<- zdf %>%
    group_by(abs_delta_rankbin100,Dx_type,BrainRegion) %>%
    summarise(Dx_cor = cor(Dx, DAG, use = "complete.obs",method='spearman'), .groups = "drop")
  
  print(ggplot(zdf2,aes(x=abs_delta_rankbin100,y=Dx_cor,color=Dx_type))+
          geom_point()+geom_smooth()+geom_vline(xintercept=70/100,linewidth=0.5,linetype = "longdash")+
          xlab('Delta bin (Dx_cell-Dx_donor)')+ylab('Cor(Dx score (DAG),Dx)')+
          ggtitle(cell)+theme_bw()+facet_wrap(~BrainRegion)+
          theme(axis.text.x = element_text(colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                legend.position = 'bottom',
                panel.grid.minor = element_blank(),aspect.ratio = 1))
  zdf2<- zdf %>%
    group_by(Dx_type,consistent) %>%
    summarise(Dx_cor = cor(Dx, DAG, use = "complete.obs",method='spearman'), 
              .groups = "drop")
  ys=max(abs(quantile(zdf$DAG,c(0.01,0.99),na.rm=T)))*c(-1,1)
  a=min(xs)*0.8
  b=max(ys)*0.9
  print(ggplot(zdf,aes(x=Dx,y=DAG))+
          xlim(xs[1],xs[2])+ylim(ys[1],ys[2])+
          geom_density_2d_filled(contour_var = "ndensity") + 
          geom_smooth(method='lm')+
          geom_text(data=zdf2,x=a,y=b,aes(label=paste0('SCC :',round(Dx_cor,3))))+
          facet_grid(consistent~Dx_type)+
          ggtitle(cell)+theme_bw()+scale_fill_manual(values=bugn(10))+
          theme(axis.text.x = element_text(colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 1))
  dev.off()
  
}
write.csv(LGB,file='lgb_anno.csv')
#######################################################################################
pdf('cor_DEG_module_score_Dx.pdf',width=5,height=5)
df1=as.data.frame(df1);df1$cell=factor(df1$cell)
print(ggplot(df1,aes(x=cell,y=Dx_cor,color=Dx_type))+
        geom_point(shape=19)+geom_line(aes(group=Dx_type),linewidth=0.5)+
        facet_grid(consistent~.)+
        theme_bw()+ylab('spearman cor(DEG module score,Dx)')+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 0.5))
df2=as.data.frame(df2);df2$cell=factor(df2$cell)
print(ggplot(df2,aes(x=cell,y=Dx_cor,color=Dx_type))+
        geom_point(shape=19)+geom_line(aes(group=Dx_type),linewidth=0.5)+
        facet_grid(consistent~.)+
        theme_bw()+ylab('spearman cor(DAG module score,Dx)')+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 0.5))
dev.off()
##################



