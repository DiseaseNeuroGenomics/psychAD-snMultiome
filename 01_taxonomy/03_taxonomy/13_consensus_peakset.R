.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
source('/sc/arion/projects/roussp01a/pengfei/tools/myArchR.R')
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
addArchRThreads(threads = 36)
root='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/'
input='save2'
output='save3'
####################### step1 read peaks per cell type and make consensus peaks
cells=c('Adaptive','Endo','Mural','OPC','EN','IN','Micro_PVM','Astro','Oligo')
setwd(root)
peaks=GRangesList(lapply(cells,function(x){obj=loadArchRProject(x);attributes(obj)$peak_summit}))
unionPeaks <- unlist(peaks)
unionPeaks <- nonOverlappingGR(unionPeaks, by = "groupScoreQuantile", decreasing = TRUE)
unionPeaks <- lapply(split(unionPeaks, seqnames(unionPeaks)), function(x){
  mcols(x)$idx <- seq_along(x)
  x
}) %>% Reduce("c", .) %>% sortSeqlevels %>% sort
###### a new object
cellanno=read_big_csv('cell_anno_final.csv')
cellanno=cellanno[cellanno$ATAC & cellanno$class==cell,]
AD3BR=AD3BR[BiocGenerics::which(rownames(AD3BR@cellColData)%in%rownames(cellanno)),]
for(x in colnames(cellanno)){
  AD3BR@cellColData[,x]=cellanno[rownames(AD3BR@cellColData),x]
}
a=table(AD3BR@cellColData$order) # order is the sampleID
a=a[a>0]
AD3BR=AD3BR[BiocGenerics::which(AD3BR@cellColData$order%in%names(a)),] 
AD3BR@sampleColData=AD3BR@sampleColData[rownames(AD3BR@sampleColData)%in%names(a),,drop=F]
AD3BR@sampleMetadata=AD3BR@sampleMetadata[names(AD3BR@sampleMetadata)%in%names(a)]
AD3BR@peakSet <- unionPeaks
saveArchRProject(ArchRProj = AD3BR,
                 outputDirectory = file.path(root,output), load = F)
AD3BR <- loadArchRProject(path = file.path(root,output), force = FALSE, showLogo = TRUE)
AD3BR <- myaddPeakMatrix(AD3BR,matrixName='PeakMatrix')
saveRDS(AD3BR,file=file.path(root,output,'Save-ArchR-Project.rds'))
###
######## export peak matrix in h5ad
library(parallel)
library(zellkonverter)
library(anndata)
library(SingleCellExperiment)
rewrite_anndata <- function(fh){
  a <- t(ArchR:::.getMatFromArrow(ArrowFile = fh, featureDF = featureDF, 
                                  useMatrix = 'PeakMatrix', binarize = F, 
                                  useIndex = FALSE))
  
  ad=AnnData(X=as(a,'RsparseMatrix'),var=as.data.frame(featureDF[,c('score','GC','idx','seqnames','peakType')]),
             obs=data.frame(row.names=rownames(a)))
  fho=file.path('h5ad_peakMat',gsub('arrow','h5ad',basename(fh)))
  write_h5ad(ad,file=fho,compression='gzip')
}
featureDF=mcols(AD3BR@peakSet)
featureDF$seqnames=as.vector(AD3BR@peakSet@seqnames)
rownames(featureDF)=paste0(as.vector(AD3BR@peakSet@seqnames),':',
                           start(AD3BR@peakSet),'-',
                           end(AD3BR@peakSet))
if(!file.exists('h5ad_peakMat')){
  dir.exists('h5ad_peakMat')
}
z=ArchR:::.safelapply(AD3BR@sampleColData$ArrowFiles,rewrite_anndata,threads = 8)
# export gene score matrix
rewrite_anndata <- function(fh){
  a <- ArchR:::.getMatFromArrow(ArrowFile = fh, featureDF = featureDF, 
                                useMatrix = 'GeneScoreinpeakMatrix', binarize = F, 
                                useIndex = FALSE)
  
  ad=AnnData(X=as(t(a),'RsparseMatrix'),var=as.data.frame(featureDF[,c('seqnames','name')]),
             obs=data.frame(row.names=colnames(a)))
  fho=file.path('h5ad_geneMat',gsub('arrow','h5ad',basename(fh)))
  write_h5ad(ad,file=fho,compression='gzip')
}
featureDF <- ArchR:::.getFeatureDF(AD3BR@sampleColData$ArrowFiles[1], subGroup = 'GeneScoreinpeakMatrix')
rownames(featureDF)=featureDF$name
if(!file.exists('h5ad_geneMat')){
  dir.exists('h5ad_geneMat')
}
z=ArchR:::.safelapply(AD3BR@sampleColData$ArrowFiles,rewrite_anndata,threads = 8)
