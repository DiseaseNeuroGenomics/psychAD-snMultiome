.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(dreamlet)
library(variancePartition)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(rhdf5)

library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
library(ggplot2)

library(chromVAR)
library(Matrix)
#library(chromVARmotifs)
library(motifmatchr)
library(TFBSTools)
library(future)
library(future.apply)
# chromvar used ATAC-seq peak output, should run ATAC-seq peak (OCR first)
StringToGRanges <- function(regions, sep = c("\\:", "-"), ...) {
  library(tidyr)
  ranges.df <- data.frame(ranges = regions)
  ranges.df <- separate(
    data = ranges.df,
    col = "ranges",
    sep = paste0(sep[[1]], "|", sep[[2]]),
    into = c("chr", "start", "end")
  )
  granges <- makeGRangesFromDataFrame(df = ranges.df, ...)
  names(granges)=granges$ID=regions
  return(granges)
}
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
cdmk_folder <- function(folder){
  if(!dir.exists(folder)){
    dir.create(folder,recursive = T)  
  }
  setwd(folder)
}
#########################################################
cell='all'
level='subclass'
BR='all'
OUTDIR='/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/'
prefix_peak=paste0('ATAC_',cell)
prefix=paste0('chromvar_',cell)
metaInfo_path='/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/files/metaInfo_11102023.csv'
sampleCol='order'
NCPU=8 # less CPU to save memory
options(future.globals.maxSize= 60*1024^3)
plan(multicore, workers = NCPU)
register(MulticoreParam(NCPU))
#######################################
##########read motif ###############
motif=readRDS('/sc/arion/projects/roussp01a/pengfei/ref/cisBP_human_pfms_2021.rds')
pwms=lapply(motif,toPWM)
pwmList <- do.call(PWMatrixList, pwms)
##
genes=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/metadata/cellranger_gene.csv')
genes[duplicated(genes$gene_name),'gene_name']=paste0(genes[duplicated(genes$gene_name),'gene_name'],'-1')
genes$mt = grepl('^MT-',genes$gene_name)
genes$ribo = grepl('^RPL',genes$gene_name)|grepl('^RPS',genes$gene_name) # annotate the group of mitochondrial genes as 'mt'
genes$hb=grepl('^HB',genes$gene_name)
genes$AS=grepl('-AS$',genes$gene_name)
genes$protein_coding = genes$gene_type=='protein_coding'
rownames(genes)=genes$gene_name
pwmList=pwmList[intersect(genes$gene_name,names(pwmList))]
############################################################################
setwd(file.path(OUTDIR,level))
load(paste0('pb4dreamlet_',prefix_peak,'_',level,'.Rdata')) # use ATAC-seq peak infromation
cdmk_folder(prefix)
lst=lapply(assayNames(pb),function(CT){
  mat <- assay(pb, CT)
  if(length(sel)>1){
    sel=which(colSums(mat)>2.5e5) # remove low coverage samples
    mat=mat[,sel,drop=F]
    colnames(mat)=paste0(CT,'_',colnames(mat))
    return(mat)
  }else{
    return(NULL)
  }
})
###############chromvar and stacked assay###########
names(lst)=assayNames(pb)
lst=lst[!sapply(lst,is.null)]
mat=do.call(cbind,lst)
peak=StringToGRanges(rownames(mat))
out <- SummarizedExperiment(assays = list(counts = mat), 
                            rowRanges = peak)
out <- addGCBias(out,genome = BSgenome.Hsapiens.UCSC.hg38)
bg <- getBackgroundPeaks(object = out)
motif_ix <- matchMotifs(pwmList, out, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
dev <- computeDeviations(object =out, annotations = motif_ix,
                         background_peaks=bg)
dv=deviations(dev)
w <- array(1,dim=dim(dv))
obj <- new("EList", list(E = as.matrix(dv), weights = as.matrix(w)))
resList=list(stacked=obj)
#######
info <- lapply(names(lst),function(x) {
  info <- colData(pb)
  info$stackedAssay <- x
  rownames(info) <- paste(x, rownames(info), sep = "_")
  info
})
info <- do.call(rbind, info)
metadata = as.data.frame(dreamlet:::get_metadata_aggr_means(pb))
zdf=metadata[,metadata(pb)$agg_pars$by]
colnames(zdf)=c('assay','id')
metadata[,metadata(pb)$agg_pars$by[2]]=paste0(metadata[,metadata(pb)$agg_pars$by[1]],'_',
                                              metadata[,metadata(pb)$agg_pars$by[2]])
metadata[,metadata(pb)$agg_pars$by[1]]='stacked'
metadata=data.frame(zdf,metadata)

res.proc=proc=new("dreamletProcessedData", resList, data =droplevels(as.data.frame(info)), 
                  df_details=data.frame(assay='stacked',n_retain=dim(obj$E)[2],
                                        formula='~(1|SubID)',formDropsTerms=F,
                                        n_genes=dim(obj$E)[1],n_errors=0,error_initial=F),
                  metadata = tibble::as_tibble(metadata), by = metadata(pb)$agg_pars$by)
save(res.proc,file='stacked.proc.Rdata')
form=~(1|Dx_cat_donor)+(1|Sex)+scale(Age)+scale(PMI..min.)+(1|SubID)+(1|BrainRegion)+(1|class)+(1|class:stackedAssay)+scale(TSSEnrichment)+scale(nFrags)
varPart=fitVarPart(res.proc,form, BPPARAM = bpparam)
save(varPart,file='stackedassay_varpart.Rdata')
