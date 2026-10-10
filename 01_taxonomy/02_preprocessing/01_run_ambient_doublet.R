library(DropletUtils)
metaInfo=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/files/metaInfo_11102023.csv',row.names=1)
sample=metaInfo$newID
ambient<- function(){
  library(SoupX)
  library(Seurat)
  library(DropletUtils)
  toc = Seurat::Read10X('filtered_feature_bc_matrix')
  tod = Seurat::Read10X('raw_feature_bc_matrix')
  gobj <- CreateSeuratObject(counts = toc$`Gene Expression`)
  gobj <- NormalizeData(gobj, normalization.method = "LogNormalize", scale.factor = 10000)
  gobj <- FindVariableFeatures(gobj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(gobj)
  gobj <- ScaleData(gobj, features = all.genes)
  gobj <- RunPCA(gobj, features = VariableFeatures(object = gobj))
  gobj <- FindNeighbors(gobj, dims = 1:10)
  gobj <- FindClusters(gobj, resolution = 1)
  
  sc = SoupChannel(tod$`Gene Expression`, toc$`Gene Expression`, calcSoupProfile = FALSE)
  soupProf = data.frame(row.names = rownames(toc$`Gene Expression`), 
                        est = rowSums(toc$`Gene Expression`)/sum(toc$`Gene Expression`),
                        counts = rowSums(toc$`Gene Expression`))
  sc = setSoupProfile(sc, soupProf)
  sc = setClusters(sc, gobj$seurat_clusters)
  sc  = autoEstCont(sc, doPlot=FALSE)
  out = adjustCounts(sc, roundToInt = TRUE)
  #save(out,file='SoupX_corrected_count.Rdata')
  a=read.table('filtered_feature_bc_matrix/features.tsv.gz',sep='\t')
  idmap=a$V1
  a$V2[duplicated(a$V2)]=paste0(a$V2[duplicated(a$V2)],'-1')
  names(idmap)=a$V2
  rownames(out)[!rownames(out)%in%a$V2]=gsub('.1$','-1',rownames(out)[!rownames(out)%in%a$V2])
  write10xCounts("soupX.h5", out,type='HDF5',overwrite = T,
                 gene.id=idmap[rownames(out)],gene.symbol=rownames(out))
  #
}
doublet<- function(){
  .libPaths(c(.libPaths(),'/sc/arion/projects/roussp01a/pengfei/tools/R_4_0_3'))
  library('scDblFinder')
  library(scater)
  library(GenomicRanges)
  library(rtracklayer)
  library(Seurat)
  set.seed(1)
  toc = Seurat::Read10X('filtered_feature_bc_matrix')
  sce1 = scDblFinder(SingleCellExperiment(list(counts=toc$`Gene Expression`)))
  sce2 <- scDblFinder(SingleCellExperiment(list(counts=toc$`Peaks`)), 
                      artificialDoublets=1, aggregateFeatures=TRUE, nfeatures=25, 
                      processing="normFeatures")
  # we then launch the method
  res <- read.table('amulet/MultipletProbabilities.txt',header=T,sep='\t')
  rownames(res)=res$cell_id
  res$scDblFinder.p <- 1-colData(sce2)[rownames(res), "scDblFinder.score"]
  res$ataccombined <- apply(res[,c("scDblFinder.p", "p.value")], 1, FUN=function(x){
    x[x<0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
    suppressWarnings(aggregation::fisher(x))
  })
  res$atacFDR=p.adjust(res$ataccombined,method='BH')
  res$atacclass=c("singlet","doublet")[(res$atacFDR<0.05)+1]
  scDbl=cbind(res[,c('ataccombined','atacclass')],colData(sce1)[rownames(res),c('scDblFinder.score','scDblFinder.class')])
  save(scDbl,file='barcode_scDblFinder.Rdata')
  barcodeinfo=read.csv('per_barcode_metrics.csv')
  barcodeinfo=barcodeinfo[barcodeinfo$is_cell==1,]
  rownames(barcodeinfo)=barcodeinfo$barcode
  barcodeinfo$scDblFinder_atac=scDbl[rownames(barcodeinfo),]$atacclass
  barcodeinfo$scDblFinder_gex=scDbl[rownames(barcodeinfo),]$scDblFinder.class
  save(barcodeinfo,file='barcodeinfo.Rdata')
}
script='/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/scripts/ambient.R'
dump("ambient",script)
write('ambient()',file=script,append=T)
script='/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/scripts/doublet.R'
dump("doublet",script)
write('doublet()',file=script,append=T)

script_name=file.path('/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/','scripts',paste0('ambient_doublet_sample.sh'))
write(file=script_name,paste0("
#!/bin/bash
#BSUB -J ambient_doublet[1-",length(sample),"]%100
#BSUB -W 4:00
#BSUB -n 2
#BSUB -q express
#BSUB -P acc_roussp01a
#BSUB -R rusage[mem=8000]
#BSUB -R span[hosts=1]
#BSUB -eo /sc/arion/projects/CommonMind/roussp01a/snmulti/step2/log/ambient_doublet_%I.err
#BSUB -oo /sc/arion/projects/CommonMind/roussp01a/snmulti/step2/log/ambient_doublet_%I.out
ml R/4.0.3

sample=$(sed \"${LSB_JOBINDEX}q;d\" <<< '", paste0(sample,collapse = '\n'), "')
cd /sc/arion/projects/CommonMind/roussp01a/snmulti/step1/files
cd ${sample}/outs
amulet_path=/sc/arion/projects/roussp01a/pengfei/tools/AMULET/
mkdir amulet
python ${amulet_path}/FragmentFileOverlapCounter_multi.py --maxinsertsize 900 --expectedoverlap 2 --startbases 0 --endbases 0 atac_fragments.tsv.gz per_barcode_metrics.csv ${amulet_path}/human_autosomes.txt amulet
python ${amulet_path}/AMULET.py --expectedoverlap 2 --rfilter ${amulet_path}/hg38.blacklist.bed amulet/Overlaps.txt amulet/OverlapSummary.txt amulet
Rscript /sc/arion/projects/CommonMind/roussp01a/snmulti/step2/scripts/ambient.R
Rscript /sc/arion/projects/CommonMind/roussp01a/snmulti/step2/scripts/doublet.R"))
cat(paste0('bsub < ',script_name))
