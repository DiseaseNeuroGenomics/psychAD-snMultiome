.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(ggplot2)
library(dplyr)
library(tidyr)
library(crumblr)
library(ggrepel)
library(dplyr)
setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/donor/new')
z1=read.csv('Dx_donor_subclass_crumblr.csv',row.names=1)
z2=read.csv('Dx_dev_subclass_crumblr.csv',row.names=1)
z1$se=z1$logFC/z1$t
z2$se=z2$logFC/z2$t
df1=data.frame(cell=rownames(z1),
               Dx_donor=z1$logFC,donor_se=z1$se,
               Dx_dev=z2[rownames(z1),]$logFC,
               dev_se=z2[rownames(z1),]$se,
              donor_fdr=z1$adj.P.Val,
              dev_fdr=z2[rownames(z1),]$adj.P.Val,group='current')
z1=read.csv('external/NPSAD_Dx_donor_subclass_crumblr.csv',row.names=1)
z2=read.csv('external/NPSAD_Dx_dev_subclass_crumblr.csv',row.names=1)
z2=z2[z2$group=='no',];rownames(z2)=z2$cell;z2$cell=NULL
z1=z1[z1$group=='no',];rownames(z1)=z1$cell;z1$cell=NULL
z1$se=z1$logFC/z1$t
z2$se=z2$logFC/z2$t
df2=data.frame(cell=rownames(z1),
               Dx_donor=z1$logFC,donor_se=z1$se,
               Dx_dev=z2[rownames(z1),]$logFC,
               dev_se=z2[rownames(z1),]$se,
               donor_fdr=z1$adj.P.Val,
               dev_fdr=z2[rownames(z1),]$adj.P.Val,group='psychAD')
z1=read.csv('external/SEAAD_CPS_subtype_crumblr.csv',row.names=1)
z2=read.csv('external/SEAAD_Dx_dev_subtype_crumblr.csv',row.names=1)
z1$se=z1$logFC/z1$t
z2$se=z2$logFC/z2$t
df3=data.frame(cell=rownames(z1),
               Dx_donor=z1$logFC,donor_se=z1$se,
               Dx_dev=z2[rownames(z1),]$logFC,
               dev_se=z2[rownames(z1),]$se,
               donor_fdr=z1$adj.P.Val,
               dev_fdr=z2[rownames(z1),]$adj.P.Val,group='SEAAD')
df3$cell2=gsub('_[0-9]+$','',df3$cell)
df3$cell3='other'
df3$cell3[df3$cell2=='Sst']='IN_SST'
df3$cell3[df3$cell2=='L2/3 IT']='EN_L2_3_IT'
df=rbind(df1,df2)
df$cell2=df$cell
df$cell3='other'
df$cell3[df$cell2%in%c('EN_L2_3_IT','IN_SST')]=df$cell2[df$cell2%in%c('EN_L2_3_IT','IN_SST')]
df=rbind(df,df3)
pdf('subclass_dx_donor_vsdx_dev.pdf',width=6,height = 6)
print(ggplot(df)+
        geom_point(data=df[df$cell3=='other',],aes(x=Dx_donor,y=Dx_dev,color=cell3))+ #size=-log10(donor_fdr)
        geom_point(data=df[df$cell3!='other',],aes(x=Dx_donor,y=Dx_dev,color=cell3))+ #size=-log10(donor_fdr)
        geom_linerange(aes(x=Dx_donor,ymin=Dx_dev-dev_se,ymax=Dx_dev+dev_se,color=cell3),linewidth = 0.35)+
        geom_linerange(aes(y=Dx_dev,xmin=Dx_donor-donor_se,xmax=Dx_donor+donor_se,color=cell3),linewidth = 0.35)+
        geom_smooth(aes(x=Dx_donor,y=Dx_dev),method='lm',linewidth = 0.7,se=F)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        geom_vline(xintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('logFC (Dx_donor)')+ylab('logFC (Dx_dev)')+facet_wrap(~group)+
        #scale_color_manual(values=colormap_class)+
        scale_color_manual(values=c("#E41A1C","#377EB8",'#66666680'))+
       
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
dev.off()
a=c('Dx_cell','Dx_donor','Dx_dev')
b=c('Dx_cell','Dx_donor','Dx_dev')
zdf=data.frame()
for(i in 1:length(a)){
  xx=a[i]
  yy=b[i]
  z1=read.csv(sprintf('%s_vulneuron_crumblr.csv',yy),row.names=1)  
  z1=z1[rownames(z1)!='Neu',]
  z1$se=NULL
  zdf=rbind(zdf,data.frame(z1,cell=rownames(z1),coef=xx,group='current'))
  if(xx=='Dx_dev'){
    z2=read.csv(sprintf('%s_DAM_crumblr.csv',yy),row.names=1)
  }else{
    z2=read.csv(sprintf('%s_Micro_PVM_crumblr.csv',yy),row.names=1)
  }
  z2$se=NULL
  zdf=rbind(zdf,data.frame(z2['DAM',],cell='Micro_PVM',coef=xx,group='current'))
}
zdf$z.std=NULL
for(i in 1:length(a)){
  xx=a[i]
  yy=b[i]
  z1=read.csv(sprintf('external/NPSAD_%s_vulneuron_crumblr.csv',yy),row.names=1)  
  z1=z1[rownames(z1)!='Neu' & z1$cell!='Neu' &z1$group=='no',]
  z1$se=NULL;z1$group=NULL;z1$se.1=NULL
  zdf=rbind(zdf,data.frame(z1,coef=xx,group='psychAD'))
  z2=read.csv(sprintf('external/NPSAD_%s_DAM_crumblr.csv',yy),row.names=1)
  z2=z2[z2$group=='no' & z2$cell=='ADAM',]
  z2$se=NULL;z2$group=NULL;z2$se.1=NULL;z2$cell=NULL
  zdf=rbind(zdf,data.frame(z2,cell='Micro_PVM',coef=xx,group='psychAD'))
}
for(i in 1:length(a)){
  xx=a[i]
  yy=b[i]
  z1=read.csv(sprintf('external/SEAAD_%s_vulneuron_crumblr.csv',yy),row.names=1)  
  z1=z1[rownames(z1)!='Neu' ,]
  z1$z.std=NULL
  z1$se=NULL;z1$group=NULL;z1$se.1=NULL
  zdf=rbind(zdf,data.frame(z1,cell=rownames(z1),coef=xx,group='SEAAD'))
  z2=read.csv(sprintf('external/SEAAD_%s_DAM_crumblr.csv',yy),row.names=1)
  z2$se=NULL;z2$group=NULL;z2$se.1=NULL;z2$cell=NULL
  z2$z.std=NULL
  zdf=rbind(zdf,data.frame(z2['Micro-PVM_3-SEAAD',],cell='Micro_PVM',coef=xx,group='SEAAD'))
}
zdf$se=zdf$logFC/zdf$t
zdf$FDR=zdf$adj.P.Val
zdf$sig=c('','+')[(zdf$FDR<0.05)+1]
zdf$coef=factor(zdf$coef,levels=c('Dx_donor','Dx_dev','Dx_cell'),ordered=T)
pdf('vulnerable_cell_across_cohort.pdf',width = 4,height = 4)
print(ggplot(zdf)+
        geom_linerange(aes(x=group,ymin=logFC-se,ymax=logFC+se,color=coef),
                       position=position_dodge(0.7),linewidth = 0.35)+
        geom_point(aes(x=group,y=logFC,color=coef,size=-log10(FDR)),
                   position=position_dodge(0.7))+
        geom_text(aes(x=group,y=logFC,group=coef,label=sig),color='white',
                  position=position_dodge(0.7))+
        scale_color_manual(values=c("#FDE725FF","#21908CFF","#472D7BFF"))+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        theme_bw()+facet_grid(.~cell,scale='free_x')+coord_flip()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 2))
print(ggplot(zdf)+
        geom_linerange(aes(x=group,ymin=logFC-se,ymax=logFC+se,color=coef),
                       position=position_dodge(0.7),linewidth = 0.35)+
        geom_point(aes(x=group,y=logFC,color=coef,size=-log10(FDR)),
                   position=position_dodge(0.7))+
        geom_text(aes(x=group,y=logFC,group=coef,label=sig),color='white',
                  position=position_dodge(0.7))+
        scale_color_manual(values=c("#FDE725FF","#21908CFF","#472D7BFF"))+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        scale_size(range=c(2,4))+
        theme_bw()+facet_grid(cell~.,scale='free_y')+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 0.4))
dev.off()
z1=read.csv('Dx_devCognitive_Resilience.csv',row.names=1)
z2=read.csv('external/NPSAD_Dx_dev_resilience.csv',row.names = 1)
z2=z2[z2$group=='no',]
z1$cell=gsub('_Epend','',rownames(z1))
z1$se=z1$logFC/z1$t
z2$B=NULL;z2$group=NULL
z2$cell=gsub('Immune','Micro_PVM',z2$cell)
z1$group='current'
z2$group='psychAD'
zdf=rbind(z1,z2)
zdf$cell=factor(zdf$cell,levels=c('Astro','Oligo','OPC','Micro_PVM','EN','IN'))
zdf$FDR=zdf$adj.P.Val
zdf$sig=c('','+')[(zdf$FDR<0.05)+1]
pdf('resilience_summary.pdf',width=4,height = 4)
print(ggplot(zdf)+
        geom_linerange(aes(x=cell,ymin=logFC-se,ymax=logFC+se,color=group),
                       position=position_dodge(0.6),linewidth = 0.35)+
        geom_point(aes(x=cell,y=logFC,color=group,size=-log10(FDR)),
                   position=position_dodge(0.6))+
        geom_text(aes(x=cell,y=logFC,group=group,label=sig),color='white',
                  position=position_dodge(0.6))+
        scale_color_manual(values=c("#E41A1C","#377EB8"))+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        scale_size(range=c(2,4))+theme_bw()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 0.5))
dev.off()

z1=read.csv('external/Mathys23_cts_mmse30_lv_Dx_dev.csv',row.names=1)
z2=read.csv('external/Mathys23_cogng_demog_slope_Dx_dev.csv',row.names = 1)
z1$cell=gsub('_Epend','',rownames(z1))
z1$se=z1$logFC/z1$t

z2$cell=gsub('_Epend','',rownames(z2))
z2$se=z2$logFC/z2$t
z1$group='MMSE'
z2$group='decline slope'

zdf=rbind(z1,z2)
zdf$cell=factor(zdf$cell,levels=c('Astro','Oligo','OPC','Micro_PVM','EN','IN'))
zdf$FDR=zdf$adj.P.Val
zdf$sig=c('','+')[(zdf$FDR<0.05)+1]
zdf$group=factor(zdf$group,levels=c('MMSE','decline slope'))
pdf('Mathys_resilience_summary.pdf',width=4,height = 4)
print(ggplot(zdf)+
        geom_linerange(aes(x=cell,ymin=logFC-se,ymax=logFC+se,color=group),
                       position=position_dodge(0.6),linewidth = 0.35)+
        geom_point(aes(x=cell,y=logFC,color=group,size=-log10(FDR)),
                   position=position_dodge(0.6))+
        geom_text(aes(x=cell,y=logFC,group=group,label=sig),color='white',
                  position=position_dodge(0.6))+
        scale_color_manual(values=c("#E41A1C","#377EB8"))+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        scale_size(range=c(2,4))+theme_bw()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 0.5))
dev.off()
