.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(dreamlet)
library(dplyr)
library(ggplot2)
library(ggrepel)
genes=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/metadata/cellranger_gene.csv')
genes[duplicated(genes$gene_name),'gene_name']=paste0(genes[duplicated(genes$gene_name),'gene_name'],'-1')
genes$mt = grepl('^MT-',genes$gene_name)
genes$ribo = grepl('^RPL',genes$gene_name) |grepl('^RPS',genes$gene_name) # annotate the group of mitochondrial genes as 'mt'
genes$hb=grepl('^HB',genes$gene_name)
genes$protein_coding = genes$gene_type=='protein_coding'
genes$filtered=genes$protein_coding & (!genes$ribo)& (!genes$mt)
genes=genes[genes$filtered,]
rownames(genes)=genes$gene_name
####
Gregulon=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/GRN/files/Gene_based_regulon.csv',row.names = 1)
Gregulon$raw=Gregulon$TF
z=do.call(rbind,strsplit(gsub("_extended",'',Gregulon$TF),'_'))
Gregulon$TF=z[,1]
Gregulon$Ntarget = gsub('^.+?([0-9]+).+?$','\\1',z[,4])
Gregulon$extend=grepl('_extended',Gregulon$raw)
Gregulon$TF2gene=c('m','p')[(z[,2]=='+')+1]
Gregulon$region2gene=c('m','p')[(z[,3]=='+')+1]
Gregulon$ID=paste0(Gregulon$TF,c('_ext','')[2-Gregulon$extend],
                   '_',Gregulon$TF2gene,'_',Gregulon$region2gene,'_',
                   Gregulon$Ntarget)
Gregulon=Gregulon[Gregulon$TF!=Gregulon$target,]
Gregulon$group=paste0(Gregulon$TF2gene,'_',Gregulon$region2gene,
                      c('_ext','')[2-Gregulon$extend])
###################################################################################
PROC=list()
FIT=list()
wd='/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/class'
z=local({load(file.path(wd,'proc.Rdata'));environment()})
PROC['class']=list(z$proc)
z=local({load(file.path(wd,'all/Dx_donor.Rdata'));environment()})
FIT['class']=list(z$fit)
###################################################################################
lgb=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ML_cellstate/lgb_anno.csv',row.names=1)
anno=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno_05092025_fixed3.csv',row.names=1)
Dx_anno=unique(lgb[,c('BrainRegion2','order','Dx_cat_donor')])
rownames(Dx_anno)=Dx_anno$order
########
root='/sc/arion/projects/CommonMind/roussp01a/snmulti/GRN/files/filter'
setwd(root)
cov=c("scale(n_genes_by_counts)","scale(pct_counts_ribo)",
      "scale(mito_genes)","scale(PMI..min.)")
for(cell in c('Astro','Oligo','OPC','Micro_PVM','EN','IN')){
  E=PROC[['class']][[cell]]$E
  fit1=FIT[['class']][[cell]]
  sh=intersect(names(fit1$Amean),union(Gregulon$TF,Gregulon$target))
  tf=intersect(sh,Gregulon$TF)
  target=intersect(sh,Gregulon$target)
  R=E[sh,rownames(fit1$design)]-fit1$coefficients[sh,cov]%*%t(fit1$design[,cov])
  grn1=Gregulon[,c('TF','target','ID','group')]
  grn1$ID=paste0(grn1$TF,":",grn1$target)
  cmat=cor(t(R[tf,]),t(R[target,]))
  
  sel=grn1$TF%in%tf & grn1$target%in%target
  xs=match(grn1[sel,]$TF,rownames(cmat))
  ys=match(grn1[sel,]$target,colnames(cmat))
  grn1[sel,'corr']=cmat[cbind(xs,ys)]
  grn1[sel,'t']=(grn1[sel,'corr'] / sqrt((pmax(1-grn1[sel,'corr']^2, 1e-17, na.rm = TRUE))/(ncol(R)-2)))
  grn1[sel,'p'] <- 2*pt(-abs(grn1[sel,'t']), ncol(R) - 2)
  grn1[sel,'fdr']=p.adjust(grn1[sel,'p'],method='BH')
  grn1$exp=F
  grn1[sel,]$exp=T
  mp=grepl('^p_',grn1$group)
  cc=(mp&(grn1$cor>0))|((!mp)&(grn1$cor<0))
  grn1$filter=grn1$fdr<0.05 &cc
  
  grn2=grn1[which(grn1$filter & grn1$group=='p_p'),]
  f=table(grn2$TF)
  f=f[f>=30]
  grn2=grn2[grn2$TF%in%names(f),]
  
  write.csv(grn2,file=paste0(cell,'_corr_p_filtered_grn.csv'))
}

