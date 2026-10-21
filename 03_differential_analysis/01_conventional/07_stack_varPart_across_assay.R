.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(variancePartition)
library(dreamlet)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(scattermore)
library(patchwork) 
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/subclass/stackedassay_varpart.Rdata')
zdf1=data.frame(varPart2[,c("gene","BrainRegion","class","class.stackedAssay","Dx_cat_donor","Sex","SubID","scale.Age.","Residuals")],
                Tech=rowSums(data.frame(varPart2[,c("scale.PMI..min..","scale.mito_genes.","scale.n_genes_by_counts.","scale.mito_ribo.")])))
zdf1=reshape2::melt(zdf1,id.vars='gene')
zdf1$assay='GEX'
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/subclass_binarize_chromvar/stackedassay_varpart.Rdata')
zdf2=data.frame(varPart2[,c("gene","BrainRegion","class","class.stackedAssay","Dx_cat_donor","Sex","SubID","scale.Age.","Residuals")],
                Tech=rowSums(data.frame(varPart2[,c("scale.PMI..min..","scale.TSSEnrichment.","scale.nFrags..")])))
zdf2=reshape2::melt(zdf2,id.vars='gene')
zdf2$assay='ATAC(Motif)'

load('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/subclass_genescoreinpeak/stackedassay_varpart.Rdata')
zdf3=data.frame(varPart2[,c("gene","BrainRegion","class","class.stackedAssay","Dx_cat_donor","Sex","SubID","scale.Age.","Residuals")],
                Tech=rowSums(data.frame(varPart2[,c("scale.PMI..min..","scale.TSSEnrichment.","scale.nFrags..")])))
zdf3=reshape2::melt(zdf3,id.vars='gene')
zdf3$assay='ATAC(GS)'
ord=c("class","class.stackedAssay",
      "SubID","Sex",
      "BrainRegion",
      "Dx_cat_donor", "scale.Age.","Tech","Residuals")
zdf=rbind(zdf1,zdf2,zdf3)
zdf$assay=factor(zdf$assay,levels=c('GEX','ATAC(GS)','ATAC(Motif)'),ordered = T)
zdf$variable=factor(zdf$variable,levels=ord,
                    ordered = T)
setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/')
pdf('stackedAssay_varPart.pdf')
zdf$value=zdf$value*100
print(ggplot(zdf, aes(x = variable, y = value)) +
        geom_violin(aes(fill = assay), 
                    position = position_dodge(0.8), scale = 'width') +
        geom_boxplot(aes(color = assay, outlier.shape = assay),
                     fill = "white",
                     outlier.stroke = 0.2, 
                     width = 0.2, 
                     position = position_dodge(0.8)) +
        geom_boxplot(data=zdf,aes(x = variable, y = value,outlier.color = assay), 
                     color = "black", 
                     fill = NA,          
                     outlier.shape = NA, 
                     width = 0.2, 
                     position = position_dodge(0.8)) +
        scale_fill_manual(values = c("#440154FF","#21908CFF","#FDE725FF")) +
        scale_color_manual(values = c("#440154FF","#21908CFF","#FDE725FF")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio = 0.618))
print(ggplot(zdf, aes(x = variable, y = value)) +
        geom_violin(aes(fill = assay), 
                    position = position_dodge(0.8), scale = 'width') +
        geom_boxplot(aes(color = assay, outlier.shape = assay),
                     fill = "white",
                     outlier.stroke = 0.2, 
                     width = 0.2, 
                     position = position_dodge(0.8)) +
        geom_boxplot(data=zdf,aes(x = variable, y = value,outlier.color = assay), 
                     color = "black", 
                     fill = NA,          
                     outlier.shape = NA, 
                     width = 0.2, 
                     position = position_dodge(0.8)) +
        scale_fill_manual(values = c("#440154FF","#21908CFF","#FDE725FF")) +
        scale_color_manual(values = c("#440154FF","#21908CFF","#FDE725FF")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              aspect.ratio = 0.5))
dev.off()

pdf('ATAC_vs_GEX_varPart.pdf',width=8)
plot_list <- list()
for(x in ord){
  a = zdf1[zdf1$variable==x,]
  b = zdf3[zdf3$variable==x,]
  rownames(a) = a$gene
  rownames(b) = b$gene
  sh = intersect(rownames(a),rownames(b))
  
  zdf = data.frame(gene=sh, GEX=a[sh,'value'], ATAC=b[sh,'value'])
  
  g = union(head(zdf[order(-1*zdf$GEX),]$gene,3),
            head(zdf[order(-1*zdf$ATAC),]$gene,3))
  
  c1 = cor(zdf$GEX, zdf$ATAC)
  c2 = cor(zdf$GEX, zdf$ATAC, method='spearman')

  p <- ggplot(zdf) +
    geom_scattermore(aes(x=GEX,y=ATAC),pointsize = 2,alpha = 0.3,color='grey40') +
    geom_smooth(aes(x=GEX,y=ATAC),method='lm',alpha=0.8) +
    geom_point(data=zdf[zdf$gene%in%g,], aes(x=GEX,y=ATAC),size=0.5,color='darkred') +
    geom_text_repel(data=zdf[zdf$gene%in%g,], aes(x=GEX,y=ATAC,label=gene)) +
    theme_bw() + xlab('')+ylab('')+
    #ggtitle(paste0(x, "\nr: ", round(c1, 3), " \u03c1",': ', round(c2, 3))) +
    labs(
      title = x,
      subtitle = bquote("r:" ~ .(round(c1, 3)) ~ ~ rho * ":" ~ .(round(c2, 3)))) +
    #ggtitle(paste0(x,'\nr: ',round(c1,3),expression(rho),': ', round(c2,3))) +
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_rect(fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1)

  plot_list[[as.character(x)]] <- p
}

wrap_plots(plot_list,ncol=4)
dev.off()

zdf=do.call(rbind,
            lapply(setdiff(unique(zdf1$variable),'Residuals'),function(x){
              a=zdf1[zdf1$variable==x,]
              b=zdf3[zdf3$variable==x,]
              rownames(a)=a$gene
              rownames(b)=b$gene
              sh=intersect(rownames(a),rownames(b))
              zdf=data.frame(gene=sh,GEX=a[sh,'value'],
                             ATAC=b[sh,'value'])
              data.frame(Var=x,pcc=cor(zdf$GEX,zdf$ATAC),
                         scc=cor(zdf$GEX,zdf$ATAC,method='spearman'),
                         GEX_median=median(zdf$GEX),
                         GEX_mean=mean(zdf$GEX),
                         ATAC_median=median(zdf$ATAC),
                         ATAC_mean=mean(zdf$ATAC),
                         Mean=mean(c(zdf$GEX,zdf$ATAC)),
                         Median=median(c(zdf$GEX,zdf$ATAC)))
            }))
zdf=zdf[zdf$Var!='Residuals',]
zdfa=reshape2::melt(zdf[,c('Var','pcc','scc','Mean'),],id.vars=c('Var','Mean'))
pdf('concordance_GEX_ATAC.pdf',width=3,height = 3)
print(ggplot(zdf[,c('Var','pcc','scc','Mean'),],aes(x=Mean,y=scc))+
        geom_point(size=2)+scale_x_log10(limits=c(1e-4,1))+
        geom_smooth(method='lm',alpha=0.2)+
        geom_text_repel(aes(x=Mean,y=scc,label=Var))+
        #scale_y_continuous(limits=c(0,0.7))+
        theme_bw()+ggtitle(paste0('rho:',round(cor(zdf$Mean,zdf$scc,method='spearman'),2)))+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
dev.off()

