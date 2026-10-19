.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(dreamlet)
library(variancePartition)
library(SingleCellExperiment)
library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)
library(missMDA)
library(FactoMineR)
library(zenith)

bpparam = MulticoreParam(12)
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
split_chr <- function(df, gene_col = "gene", chr_col = "chr", seed = 1) {
  set.seed(seed)
  
  # Count unique genes per chromosome
  x <- unique(df[c(gene_col, chr_col)])
  x[[chr_col]] <- as.integer(sub("^chr", "", x[[chr_col]], ignore.case = TRUE))
  n_gene <- table(factor(x[[chr_col]], levels = 1:22))
  
  capacity <- c(5, 5, 4, 4, 4)
  round <- integer(22)
  round_nchr <- integer(5)
  round_ngene <- numeric(5)
  
  # Assign chromosomes with more genes first
  for (i in order(n_gene, decreasing = TRUE)) {
    available <- which(round_nchr < capacity)
    candidates <- available[
      round_ngene[available] == min(round_ngene[available])
    ]
    j <- candidates[sample.int(length(candidates), 1)]
    
    round[i] <- j
    round_nchr[j] <- round_nchr[j] + 1
    round_ngene[j] <- round_ngene[j] + n_gene[i]
  }
  
  result <- data.frame(
    chr = 1:22,
    n_genes = as.integer(n_gene),
    round = round
  )
  
  result[order(result$round, result$chr), ]
}

setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files')
genes=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/metadata/cellranger_gene.csv')
genes[duplicated(genes$gene_name),'gene_name']=paste0(genes[duplicated(genes$gene_name),'gene_name'],'-1')
genes$mt = grepl('^MT-',genes$gene_name)
genes$ribo = grepl('^RPL',genes$gene_name) |grepl('^RPS',genes$gene_name) # annotate the group of mitochondrial genes as 'mt'
genes$hb=grepl('^HB',genes$gene_name)
genes$protein_coding = genes$gene_type=='protein_coding'
genes$filtered=genes$protein_coding & (!genes$ribo)& (!genes$mt) & genes$seqnames%in%paste0('chr',c(1:22))
#genes=genes[genes$filtered,]
rownames(genes)=genes$gene_name

metaInfo=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/metadata/metaInfo.csv')
#
covariates=c("(1|SubID)",'BrainRegion','scale(DxPC1)',
             "scale(n_genes_by_counts)","scale(pct_counts_ribo)",
             "scale(mito_genes)","scale(mito_ribo)",
             "scale(Age)", "Sex", "scale(PMI..min.)")
form = as.formula(paste0("~  ",paste0(covariates,collapse=" + ")))
#
load('class/all/varpart.Rdata')
GEX_VP=fit
load('class/all.proc.Rdata')
GEX_Vobj=proc
load('class/all/DxPC1.Rdata')
GEX_Dobj=fit

Sel=data.frame()
for(cell in c('Astro','Micro_PVM','EN','IN','Oligo','OPC')){
  vp=GEX_VP[GEX_VP$assay==cell,]
  vp=vp[vp$gene%in%genes[genes$filtered,]$gene_name,]
  E=proc[[cell]]$E
  fit1=GEX_Dobj[[cell]]
  rowsh=names(fit1$Amean)
  colsh=rownames(fit1$design)
  cov=intersect(colnames(fit1$coefficients),
                c('(Intercept)','DxPC1','BrainRegionBM-36','BrainRegionBM-44',
                  'BrainRegionBM-9/46','scale(n_genes_by_counts)','scale(mito_genes)',
                  'scale(mito_ribo)','scale(Age)','SexMale','scale(PMI..min.)'))
  R=E[rowsh,colsh]-fit1$coefficients[rowsh,cov]%*%t(fit1$design[,cov])  
  vp=vp[order(-1*vp$SubID),]
  vp=vp[vp$gene%in%rownames(R),]
  

  #####
  geneExpr=assay(proc,cell)
  ids=colnames(geneExpr)
  Info2=dreamlet:::merge_metadata(
    colData(proc)[ids, , drop = FALSE],
    metadata(proc),
    cell,
    proc@by)
  
  Info2=Info2[colnames(geneExpr),]
  modelFit1 <- fitVarPartModel(geneExpr, form, useWeights=T,
                               weightsMatrix=geneExpr$weights,colinearityCutoff=.999999,
                               Info2 , REML = TRUE, showWarnings = FALSE,BPPARAM=bpparam)
  chunk_size <- ceiling(length(modelFit1) /12)
  chunks <- split(modelFit1, ceiling(seq_along(1:length(modelFit1)) / chunk_size))
  z=bplapply(chunks, function(x) {
    a=lapply(x,function(y){
      z=lme4::ranef(y,condVar=TRUE)[[1]]
      post_var=attributes(z)$postVar[1,1,]
      vc=as.data.frame(lme4::VarCorr(y))
      tau2 <- vc$vcov[vc$grp == "SubID" & vc$var1 == "(Intercept)"][1]
      reliability <- 1 - post_var / tau2
      reliability <- pmax(0, pmin(1, reliability))
      list(z,reliability)
    });
    a1=data.frame(do.call(cbind,lapply(a,'[[',1)));colnames(a1)=names(x)
    a2=data.frame(do.call(cbind,lapply(a,'[[',2)));dimnames(a2)=dimnames(a1)
    list(a1,a2)
  },BPPARAM = bpparam)
  names(z)=NULL
  randEf=t(do.call(cbind,lapply(z,'[[',1)))
  reliab=t(do.call(cbind,lapply(z,'[[',2)))
  reliab[is.na(reliab)]=0
  write.csv(R[as.character(vp$gene),],file=sprintf('donor/new/%s_resid.csv',cell))
  write.csv(randEf[as.character(vp$gene),],file=sprintf('donor/new/%s_RandEf.csv',cell))
  ##############top 25%
  
  sel=vp[vp$SubID>=quantile(vp$SubID,0.75) & vp$SubID>0.3,]$gene
  sel=sel[rowMeans(reliab[sel,]>0.5)>0.5] # median reliability >50%
  Sel=rbind(Sel,data.frame(cell=cell,gene=sel,threshold='top25'))

  sel=vp[vp$SubID>=quantile(vp$SubID,0.5),]$gene
  sel=sel[rowSums(reliab[sel,]>0)>0]
  Sel=rbind(Sel,data.frame(cell=cell,gene=sel,threshold='top50'))
  
  sel=vp$gene
  sel=sel[rowSums(reliab[sel,]>0)>0]
  Sel=rbind(Sel,data.frame(cell=cell,gene=sel,threshold='all'))
  ##############
}
system('pigz /sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/donor/new/*csv')
###
Sel$chr=genes[as.character(Sel$gene),]$seqnames
write.csv(Sel,file='donor/new/sel.csv',row.names=F)

zdf=do.call(rbind,
            lapply(split(Sel,Sel$threshold),function(x){
              data.frame(split_chr(genes[genes$gene_name%in%x$gene,],chr_col = 'seqnames',gene_col = 'gene_name'),
                         threshold=x$threshold[1])
            }))
write.csv(zdf,file='donor/new/chr_block.csv',row.names=F)

