.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(ArchR)
source('/sc/arion/projects/roussp01a/pengfei/tools/myArchR.R')
addArchRThreads(threads = 12)
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

args=commandArgs(trailingOnly = TRUE)
cell=args[1]
options(future.globals.maxSize = 10000 * 1024^2)
root='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/'
input='save2'
AD3BR<- loadArchRProject(path = root, force = FALSE, showLogo = TRUE)
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
saveArchRProject(ArchRProj = AD3BR,
                 utputDirectory = file.path(root,cell), load = T)
AD3BR=loadArchRProject(path = file.path(root,cell), force = FALSE, showLogo = TRUE)
AD3BR<- addGroupCoverages(ArchRProj = AD3BR, groupBy = "subtype",
                          minCells = 500, maxCells = 3000,
                          minReplicates = 2, maxReplicates =5)
AD3BR <- addReproduciblePeakSet(
  ArchRProj = AD3BR, 
  groupBy = "subtype", 
  pathToMacs2 = pathToMacs2,
  additionalParams = "--nomodel",
  plot=F
)
AD3BR <- myaddPeakMatrix(AD3BR,matrixName='PeakMatrix') #PIC peak matrix
saveRDS(AD3BR,file=file.path(root,cell,'Save-ArchR-Project.rds'))
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


