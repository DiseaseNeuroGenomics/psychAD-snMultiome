.libPaths(c(.libPaths(),'/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1'))
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rhdf5)
options(future.globals.maxSize = 10000 * 1024^2)
read_obs <- function(path){
  h5f = H5Fopen(path)
  RNA=h5f$obs
  l=sapply(RNA,length)
  dat1=do.call(data.frame,lapply(RNA[l==2],function(x)c(NA,x$categories)[x$codes+2]))
  dat2=do.call(data.frame,RNA[l!=2])
  RNA=data.frame(dat1,dat2)
  rownames(RNA)=RNA$X_index
  RNA
}
check_setwd <- function(dir){
  if(!dir.exists(dir)){
    dir.create(dir)
  }
  setwd(dir)
}
addArchRThreads(threads = 16) 
addArchRGenome("hg38")
root='/sc/arion/projects/CommonMind/roussp01a/snmulti'
check_setwd(file.path(root,'step1/ArchR/save0'))
metaInfo=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/files/metaInfo_11102023.csv',row.names=1)
metaInfo=metaInfo[metaInfo$Estimated.number.of.cells>200 & metaInfo$Estimated.number.of.cells<20000, ]
#### convert to ArchR
sample=metaInfo$newID
ArrowFiles <- createArrowFiles(
  inputFiles = metaInfo$inputFiles,
  sampleNames = metaInfo$order, ## order is the sampleID
  filterTSS = 1, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
AD3BR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "save0",
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
for(i in colnames(metaInfo)){
  AD3BR@cellColData[, i] <- metaInfo[AD3BR$Sample,i]
}
RNA=read_obs(file.path(root,'step3/files/RNA_all_raw.h5ad'))
ind=BiocGenerics::which(AD3BR$cellNames %in% rownames(RNA))
sel_meta2=c('scDblFinder_atac','scDblFinder_gex')
AD3BR@cellColData[, 'inGEX']=(AD3BR$cellNames %in% rownames(RNA))
for(x in sel_meta2){
  AD3BR@cellColData[, x] <- NA
  AD3BR@cellColData[ind,x]=RNA[AD3BR$cellNames[ind],x]
}
saveArchRProject(ArchRProj = AD3BR, 
                 outputDirectory = file.path(root,"/step3/files/ArchR/save1"), load = T)
############ use RNA class label to call peaks
RNA=read_obs('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/RNA_all_initial_cutoff_ngene_500_10000_ntotal_100000_pct_mito1_n_neighbors100_Nhvg5000_NPC30_qc_01.h5ad')
AD3BR<- loadArchRProject(path = file.path(root,input), force = FALSE, showLogo = TRUE)

AD3BR=AD3BR[BiocGenerics::which(AD3BR$scDblFinder_atac!='doublet'),]
AD3BR=AD3BR[BiocGenerics::which(AD3BR$scDblFinder_gex!='doublet'),]
AD3BR=AD3BR[BiocGenerics::which(AD3BR$TSSEnrichment>3),]
AD3BR=AD3BR[BiocGenerics::which(AD3BR$BlacklistRatio<0.01),]
AD3BR=AD3BR[BiocGenerics::which(AD3BR$cellNames %in% rownames(RNA)),]
RNA=RNA[rownames(RNA)%in%AD3BR$cellNames,]
AD3BR@cellColData[rownames(RNA),'rna_cell_type']=RNA[,'cell_type']
a=table(AD3BR@cellColData$order)
AD3BR=AD3BR[BiocGenerics::which(AD3BR@cellColData$order%in%names(a)[a>200]),]
AD3BR<- addGroupCoverages(ArchRProj = AD3BR, groupBy = "rna_cell_type",  minCells = 500, maxCells = 3000,
                          minReplicates = 2, maxReplicates =5)
pathToMacs2 <- findMacs2()
AD3BR <- addReproduciblePeakSet(
  ArchRProj = AD3BR, 
  groupBy = "rna_cell_type", 
  pathToMacs2 = pathToMacs2
)
save(AD3BR,file=file.path(root,'inintial_peakcalling.Rdata')) ### peaks only used for taxonomy
########## initial taxonomy
root='/sc/arion/projects/CommonMind/roussp01a/snmulti'
input='save1'
output='save2'
mypeak=local({load('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ArchR/inintial_peakcalling.Rdata');environment()})$AD3BR
AD3BR<- loadArchRProject(path = file.path(root,input), force = FALSE, showLogo = TRUE)
AD3BR@peakSet=mypeak@peakSet
AD3BR=AD3BR[BiocGenerics::which(AD3BR$scDblFinder_atac!='doublet'),]
AD3BR=AD3BR[BiocGenerics::which(AD3BR$scDblFinder_gex!='doublet'),]
AD3BR=AD3BR[BiocGenerics::which(AD3BR$scDblFinder_gex!=2),]
#AD3BR <- AD3BR[BiocGenerics::which(AD3BR@cellColData$DoubletEnrichment<quantile(AD3BR@cellColData$DoubletEnrichment,0.99)),]

AD3BR=AD3BR[BiocGenerics::which(AD3BR$TSSEnrichment>3),]
AD3BR=AD3BR[BiocGenerics::which(AD3BR$BlacklistRatio<0.01),]

a=table(AD3BR@cellColData$order)
AD3BR=AD3BR[BiocGenerics::which(AD3BR@cellColData$order%in%names(a)[a>200]),]
AD3BR@sampleColData=AD3BR@sampleColData[rownames(AD3BR@sampleColData)%in%names(a)[a>200],,drop=F]
AD3BR@sampleMetadata=AD3BR@sampleMetadata[names(AD3BR@sampleMetadata)%in%names(a)[a>200]]
AD3BR <- myaddPeakMatrix(AD3BR) # add peak matrix with PIC https://doi.org/10.1038/s41592-023-02103-7
AD3BR <- addIterativeLSI_with_seurat( # we customed the code from ArchR to perform harmony before variable feature selection, this avoid prioritze sample specific features as variable features. 
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
  varFeatures = 150000, 
  dimsToUse = 1:30,
  projectCellsPre=TRUE,
  sampleCellsFinal=3e5
)
AD3BR@reducedDims$IterativeLSI$matSVD=scale(AD3BR@reducedDims$IterativeLSI$matSVD[,-1]) # 1st LSI were removed
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
AD3BR <- addUMAP(
  ArchRProj = AD3BR, 
  reducedDims = "Harmony", 
  name = "UMAPH", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

fh=AD3BR@sampleColData$ArrowFiles
Ninpeak <- ArchR:::.getColSums(ArrowFiles = fh, useMatrix = 'PeakMatrix', 
                         seqnames = paste0('chr',c(1:22,'X')),
                         threads=16)
AD3BR@cellColData$Ninpeak_initial=Ninpeak[rownames(AD3BR@cellColData)]
AD3BR@cellColData$FRiP_initial=AD3BR@cellColData$Ninpeak_initial/AD3BR@cellColData$nFrags
AD3BR@cellColData$FRiP_initial[AD3BR@cellColData$FRiP_initial>1]=1
write.csv(AD3BR@cellColData[,c('nFrags','Ninpeak_initial','FRiP_initial')],
          file='FRiP.csv')

saveArchRProject(ArchRProj = AD3BR,
                 outputDirectory = file.path(root,output), load = T)
#################################################################################
setwd(file.path(root,output))

write.csv(AD3BR@cellColData,file='initial_meta.csv')

h5createFile("initial_meta_PCA.h5")
h5write(AD3BR@reducedDims$Harmony$matDR, "initial_meta_PCA.h5","Harmony_PCA")
h5write(AD3BR@reducedDims$IterativeLSI$matSVD, "initial_meta_PCA.h5","PCA")

