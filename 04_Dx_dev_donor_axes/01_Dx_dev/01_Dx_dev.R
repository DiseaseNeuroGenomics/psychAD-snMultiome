.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(ggplot2)
library(tidyverse)
library(scales)
library(dplyr)
library(crumblr)
library(data.table)
library(qvalue)
library(variancePartition)
####
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

setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate3/')
genes=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/metadata/cellranger_gene.csv')
genes[duplicated(genes$gene_name),'gene_name']=paste0(genes[duplicated(genes$gene_name),'gene_name'],'-1')
genes$mt = grepl('^MT-',genes$gene_name)
genes$ribo = grepl('^RPL',genes$gene_name) |grepl('^RPS',genes$gene_name) # annotate the group of mitochondrial genes as 'mt'
genes$hb=grepl('^HB',genes$gene_name)
genes$protein_coding = genes$gene_type=='protein_coding'
genes$filtered=genes$protein_coding & (!genes$ribo)& (!genes$mt) & genes$seqnames%in%paste0('chr',c(1:22,'X'))
#genes=genes[genes$filtered,]
rownames(genes)=genes$gene_name
#############################################################################################
lgb=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ML_cellstate3/lgb_anno.csv')
lgb$celltype=lgb$subtype
lgb[lgb$class%in%c('EN','IN'),]$celltype=lgb[lgb$class%in%c('EN','IN'),]$subclass


metaInfo=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/files/metaInfo_11102023.csv',row.names=1)
metaInfo$Cognitive_ResilienceUse=metaInfo$Cognitive_Resilience
metaInfo$Cognitive_ResilienceUse[metaInfo$DxPC1_clust3==1]=NA # resilience only for early/late AD
#######
anno=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno_05092025_fixed3.csv')
anno=anno[anno$inRNA,]
anno=anno[anno$subtype!='ependA',]
#########
Obs=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/RNA_all_obs.csv.gz')
covs=c('n_genes_by_counts','pct_counts_ribo','mito_genes')
####

Dx_comp_adj1=list()
Dx_comp_adj2=list()
Vp=data.frame()
for(cell in unique(lgb$cell)){
  lgb1=lgb[lgb$class==cell,]
  agg=as.data.frame(lgb1%>%group_by(order)%>%summarise(M=mean(Dx_cell),Dx_donor=mean(DxPC1)))
  rownames(agg)=agg$order;agg$order=NULL
  
  ncell=apply(table(lgb1$order,lgb1$Dx_cat_cell),2,c)
  cobj=crumblr(ncell)
  cobj$E=rbind(cobj$E,agg=agg[colnames(cobj),]$M)
  cobj$weights=rbind(cobj$weights,agg=1)
  covInfo=do.call(rbind,lapply(split(Obs[rownames(lgb1),covs],Obs[rownames(lgb1),'order']),colMeans))
  form1=~Age+Sex+scale(DxPC1)+scale(PMI..min.)
  form2=~Age+Sex+scale(DxPC1)+scale(PMI..min.)+scale(n_genes_by_counts)+scale(pct_counts_ribo)+scale(mito_genes)
  vform= ~scale(Age)+(1|Sex)+(1|SubID)+scale(DxPC1)+scale(PMI..min.)+(1|BrainRegion)+(1|APOE)+scale(prs_scaled_AD_Bellenguez)+scale(n_genes_by_counts)+scale(pct_counts_ribo)+scale(mito_genes)
  vp=fitExtractVarPartModel(cobj,vform,
        data.frame(covInfo[colnames(cobj),],metaInfo[colnames(cobj),]))
  vp=data.frame(vp,ID=rownames(vp),assay=cell)
  Vp=rbind(Vp,vp)
  Dx_comp_adj1[[cell]]=t(residuals(dream(cobj,form1,
                    data.frame(covInfo[colnames(cobj),],metaInfo[colnames(cobj),]))))
  Dx_comp_adj2[[cell]]=t(residuals(dream(cobj,form2,
                       data.frame(covInfo[colnames(cobj),],metaInfo[colnames(cobj),]))))
}
names(Dx_comp_adj1)=names(Dx_comp_adj2)=unique(lgb$cell)

f=unique(unlist(lapply(Dx_comp_adj1,rownames)))

z1=do.call(cbind,lapply(Dx_comp_adj1,function(x)data.frame(x)[f,4]))
z2=do.call(cbind,lapply(Dx_comp_adj2,function(x)data.frame(x)[f,4]))
rownames(z1)=rownames(z2)=f
write.csv(z1,file='Dx_dev.csv')
write.csv(z2,file='Dx_dev_covadj.csv') # covariates adjusted were used in following analysis
write.csv(Vp,file='Dx_cell_vp.csv')
vp=Vp[Vp$ID=='agg',]
plot_data <- vp |>
  transmute(
    assay,
    `Dx_donor` = scale.DxPC1.,
    `SubID` = SubID,
    `BrainRegion` = BrainRegion,
    `Genetics` = APOE + scale.prs_scaled_AD_Bellenguez.,
    `RNA_QC` = scale.n_genes_by_counts. + scale.mito_genes. + scale.pct_counts_ribo.,
    `Other` = scale.Age.+Sex+scale.PMI..min..,
    `Residual` = Residuals
  ) |>
  pivot_longer(-assay, names_to = "component", values_to = "variance") |>
  mutate(
    assay = factor(assay, levels = c("Astro",'Oligo','OPC','Micro_PVM','EN','IN')),
    component = factor(
      component,
      levels = c(
        "Dx_donor", "SubID",'BrainRegion', "Genetics", "RNA_QC",
        "Other", "Residual"
      )
    )
  )
##
component_colors <- c(
  "Dx_donor" = "#984EA3",
  "SubID" = "#3B5BA7",
  'BrainRegion'="#4DAF4A",
  'Genetics'="#E69F00",
  "RNA_QC" = "#7A9CC6",
  "Other" = "#A65628",
  "Residual" = "#B8B8B8"
)
pdf('Dx_cell_vp.pdf',width=5,height = 5)
print(ggplot(plot_data, aes(x = assay, y = variance, fill = component)) +
        geom_col(
          width = 0.78,
          position = position_stack(reverse = TRUE),
          colour = "white",
          linewidth = 0.2
        ) +
        scale_fill_manual(values = component_colors, drop = FALSE) +
        scale_y_continuous(
          labels = label_percent(accuracy = 1),
          breaks = seq(0, 1, 0.2),
          expand = expansion(mult = c(0, 0.01))
        ) +
        theme_bw()+coord_flip()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1)) 
print(ggplot(plot_data, aes(x = assay, y = variance, fill = component)) +
        geom_col(
          width = 0.78,
          position = position_stack(reverse = TRUE),
          colour = "white",
          linewidth = 0.2
        ) +
        geom_text(
          aes(label = if_else(variance >= 0.2, percent(variance, accuracy = 1), "")),
          position = position_stack(vjust = 0.5, reverse = TRUE),
          size = 3,
          colour = "white"
        ) +
        scale_fill_manual(values = component_colors, drop = FALSE) +
        scale_y_continuous(
          labels = label_percent(accuracy = 1),
          breaks = seq(0, 1, 0.2),
          expand = expansion(mult = c(0, 0.01))
        ) +
        theme_bw()+coord_flip()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1)) 
print(ggplot(plot_data, aes(x = assay, y = variance, fill = component)) +
        geom_col(
          width = 0.78,
          position = position_stack(reverse = TRUE),
          colour = "white",
          linewidth = 0.2
        ) +
        geom_text(
          aes(label = if_else(variance >= 0.2, percent(variance, accuracy = 1), "")),
          position = position_stack(vjust = 0.5, reverse = TRUE),
          size = 3,
          colour = "white"
        ) +
        scale_fill_manual(values = component_colors, drop = FALSE) +
        scale_y_continuous(
          labels = label_percent(accuracy = 1),
          breaks = seq(0, 1, 0.2),
          expand = expansion(mult = c(0, 0.01))
        ) +
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1)) 
dev.off()

