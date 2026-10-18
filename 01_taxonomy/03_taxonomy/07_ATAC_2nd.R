.libPaths(c(.libPaths(),'/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1'))
library(ArchR)
source('/sc/arion/projects/roussp01a/pengfei/tools/myArchR.R')
addArchRThreads(threads = 12)

args=commandArgs(trailingOnly = TRUE)
cell=args[1]
options(future.globals.maxSize = 10000 * 1024^2)
root='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/save2'
setwd(root)
AD3BR<- loadArchRProject(path = root, force = FALSE, showLogo = TRUE)
read_obs <- function(path){
  h5f = H5Fopen(path)
  dat=h5f$obs
  l=sapply(dat,length)
  dat1=do.call(data.frame,lapply(dat[l==2],function(x)c(NA,x$categories)[x$codes+2]))
  dat2=do.call(data.frame,dat[l!=2])
  dat=data.frame(dat1,dat2)
  if('X_index'%in%colnames(dat)){
    rownames(dat)=dat$X_index  
  }
  if('barcodekey'%in%colnames(dat)){
    rownames(dat)=dat$barcodekey
  }
  
  dat
}
dat=read_obs('qc_01.h5ad')

dat=dat[dat$cell_type==cell,]
AD3BR=AD3BR[BiocGenerics::which(rownames(AD3BR@cellColData)%in%rownames(dat)),]
AD3BR@cellColData[,'cell_type']=dat[rownames(AD3BR@cellColData),'cell_type']
peak=readRDS(paste0(root,'/PeakCalls/',peak_name,'-reproduciblePeaks.gr.rds'))# use only it's own peaks
################filter by cell
dir.create(paste0('fine_cluster/',cell),recursive = T)
setwd(paste0('fine_cluster/',cell))
a=table(AD3BR@cellColData$order[AD3BR@cellColData[,cell_anno]==cell]) # order is the sampleID
a=a[a>1]
AD3BR=AD3BR[BiocGenerics::which(AD3BR@cellColData$order%in%names(a)),] # order is the sampleID
AD3BR@sampleColData=AD3BR@sampleColData[rownames(AD3BR@sampleColData)%in%names(a),,drop=F]
AD3BR@sampleMetadata=AD3BR@sampleMetadata[names(AD3BR@sampleMetadata)%in%names(a)]
AD3BR=AD3BR[BiocGenerics::which(AD3BR@cellColData[,cell_anno]==cell),]
#############################use RNA information as guidance
AD3BR <- addIterativeLSI_with_seurat(
  ArchRProj = AD3BR,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  sampleCellsPre = 200000,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1), 
    n.start = 12
  ), 
  #    LSIMethod=4,
  binarize=F,
  excludeChr='chrX',
  varFeatures = 75000, 
  dimsToUse = 1:10,
  projectCellsPre=F,
  sampleCellsFinal=min(3e5,nrow(AD3BR@cellColData)),
  inputPeak=peak # only use the peak from its cell class
)
######################
AD3BR@reducedDims$IterativeLSI$matSVD=scale(AD3BR@reducedDims$IterativeLSI$matSVD[,-1])
colnames(AD3BR@reducedDims$IterativeLSI$matSVD)=paste0('LSI',1:ncol(AD3BR@reducedDims$IterativeLSI$matSVD))
AD3BR@reducedDims$IterativeLSI$nDimensions=AD3BR@reducedDims$IterativeLSI$nDimensions-1
AD3BR@reducedDims$IterativeLSI$svd$d=AD3BR@reducedDims$IterativeLSI$svd$d[-1]
AD3BR@reducedDims$IterativeLSI$svd$u=AD3BR@reducedDims$IterativeLSI$svd$u[,-1]
AD3BR@reducedDims$IterativeLSI$svd$v=AD3BR@reducedDims$IterativeLSI$svd$v[,-1]
AD3BR@reducedDims$IterativeLSI$corToDepth$scaled=AD3BR@reducedDims$IterativeLSI$corToDepth$scaled[-1]
names(AD3BR@reducedDims$IterativeLSI$corToDepth$scaled)=paste0('LSI',1:length(AD3BR@reducedDims$IterativeLSI$corToDepth$scaled))
print('remove PC1 as it usually explained by technical factors')


AD3BR <- addHarmony(
  ArchRProj = AD3BR,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "newID"
)
library(rhdf5)
if(file.exists("meta_PCA.h5")){
  file.remove("meta_PCA.h5")
}
h5createFile("meta_PCA.h5")
h5write(AD3BR@reducedDims$Harmony$matDR, "meta_PCA.h5","Harmony_PCA")
h5write(AD3BR@reducedDims$IterativeLSI$matSVD, "meta_PCA.h5","PCA")
write.csv(AD3BR@cellColData,file='meta.csv')


