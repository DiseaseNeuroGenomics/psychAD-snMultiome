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
cell=commandArgs(T)[1]
level='class'
BR='all'
OUTDIR='/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/DE3/'
prefix_peak=paste0('DxATAC_',cell)
prefix=paste0('Dxchromvar_',cell)
Dx_dev_path='/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate3/Dx_dev_covadj.csv'
Dx_cell_path='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ML_cellstate3/lgb_anno_ATAC.csv'
metaInfo_path='/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/files/metaInfo_11102023.csv'
sampleCol='order'
NCPU=8 # less CPU to save memory
options(future.globals.maxSize= 60*1024^3)
plan(multicore, workers = NCPU)
register(MulticoreParam(NCPU))
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
load(paste0(prefix_peak,'/proc.Rdata'))
#chromvar
pb=pb1

lst=future_lapply(names(proc),function(CT){
  mat <- assay(pb, CT)
  sel=colnames(proc[[CT]]$E)
  if(length(sel)>0){
    mat=mat[,sel]
  }else{
    return(NULL)
  }
  mat=mat[rowSums(mat)>0,]
  peak=StringToGRanges(rownames(mat))
  out <- SummarizedExperiment(assays = list(counts = mat), 
                              rowRanges = peak)
  out <- addGCBias(out,genome = BSgenome.Hsapiens.UCSC.hg38)
  bg <- getBackgroundPeaks(object = out)
  motif_ix <- matchMotifs(pwmList, out, 
                          genome = BSgenome.Hsapiens.UCSC.hg38)
  list(out,bg,motif_ix)
})
names(lst)=assayNames(pb)
lst=lst[sapply(lst,length)==3]
resList <- lapply(lst ,function(l) {
  out=l[[1]]
  bg=l[[2]]
  motif_ix=l[[3]]
  # computing deviations
  dev <- computeDeviations(object =out, annotations = motif_ix,
                           background_peaks=bg)
  dv=deviations(dev)
  w <- array(1,dim=dim(dv))
  obj <- new("EList", list(E = as.matrix(dv), weights = as.matrix(w)))
  obj
})

#names(resList) <- assayNames(pb)
resList <- resList[!vapply(resList, is.null, FUN.VALUE = logical(1))]
res.proc=new("dreamletProcessedData", resList, data =droplevels(as.data.frame(colData(pb))), 
             metadata = dreamlet:::get_metadata_aggr_means(pb), by = metadata(pb)$agg_pars$by)
z=reshape2::melt(do.call(rbind,pb@int_colData$n_cells))
rownames(z)=with(z,paste(Var1,Var2,sep='_'))
res.proc@metadata$ncell=z[paste0(res.proc@metadata$sampleDx,'_',unlist(res.proc@metadata[,level])),]$value
proc=res.proc
cdmk_folder(file.path(OUTDIR,prefix))
save(proc,file='proc.Rdata')
######dreamlet####################
covariates=c("(1|SubID)",'BrainRegion',
             "scale(TSSEnrichment)",
             "scale(nFrags)",'scale(logncell)',
             "scale(NucleosomeRatio)","scale(BlacklistRatio)",
             "scale(Age)", "Sex", "scale(PMI..min.)")
proc@metadata$logncell=log(proc@metadata$ncell)
proc@data$Dx_cat_donor=factor(proc@data$Dx_cat_donor,levels=c('control','early','late'),ordered=T)
proc@data$Dx_cat_cell=factor(proc@data$Dx_cat_cell,levels=c('control','early','late'),ordered=T)
proc@metadata$Dx_cell=scale(proc@metadata$Dx_cell)
proc@data$Dx_donor=scale(proc@data$Dx_donor)
#Dx_dev
Dx_dev=read.csv(Dx_dev_path,row.names = 1)
a=do.call(rbind,lapply(colnames(Dx_dev),function(x){
  data.frame(cell=x,
             Dx_dev=rowMeans(scale(Dx_dev),na.rm=T),
             Dx_dev2=rowMeans(scale(Dx_dev[,setdiff(colnames(Dx_dev),x)]),na.rm=T),
             sample=rownames(Dx_dev))}))
rownames(a)=paste0(a$cell,':',a$sample)
a=a[paste0(proc@metadata$class,":",gsub('_.*','',proc@metadata$sampleDx)),]
proc@metadata$Dx_dev=a$Dx_dev2 # Dx_dev were derived from other cell types. 

#Dx_contrast
lgb=read_big_csv(Dx_cell_path)
lgb$order=lgb[,sampleCol]
z=lgb%>%group_by(class,order,Dx_cat_cell)%>%summarise(Dx_cell=mean(Dx_cell))
M=z%>%group_by(class,order)%>%summarise(M=mean(Dx_cell),N=n())
M=as.data.frame(M[M$N>1,])
rownames(M)=paste0(M$class,':',M[,sampleCol]) # order is sample ID
z$Dx_contrast0=z$Dx_cell-M[paste0(z$class,':',z[,sampleCol]),]$M
z=z%>%group_by(class)%>%mutate(Dx_contrast=scale(Dx_contrast0)[,1])
z=data.frame(z)
rownames(z)=paste0(z[,sampleCol],"_",z$Dx_cat_cell,":",z$class)
z=z[paste0(proc@metadata$sampleDx,':',proc@metadata$class),]
proc@metadata$Dx_contrast=z$Dx_contrast
form1 = as.formula(paste0("~  ", 'scale(Dx_dev)+scale(Dx_contrast)+scale(Dx_donor)', " + ", paste0(gsub('BrainRegion','(1|BrainRegion)',gsub('Sex','(1|Sex)',covariates)), collapse=" + ")))  
fit1=fitVarPart(proc,form1,BPPARAM = bpparam)
pdf('varpart_dev_contrast.pdf',width=8,height = length(unique(fit1$assay))/0.5+5)
print(plotVarPart(fit1,label.angle =45))
dev.off()
save(fit1,file='varpart_dev_contrast.Rdata')

test <- function(){
  if(is.null(ctr$contrasts)){
    form = as.formula(paste0("~  ", paste0(ctr$variable,collapse = '+'), " + ", paste0(covariate_use, collapse=" + ")))  
  }else{
    form = as.formula(paste0("~  0+", ctr$variable, " + ", paste0(covariate_use, collapse=" + ")))  
  }
  
  
  fit = dreamlet(proc, form, 
                 contrasts = ctr$contrasts,
                 BPPARAM = bpparam)
  ctorder = assayNames(fit)
  if(is.null(ctr$contrasts)){
    
    contrasts=coefNames(fit)[grepl(ctr$variable,coefNames(fit),fixed=T)|
                               grepl(paste0(rev(unlist(strsplit(ctr$variable,":"))),collapse=':'),coefNames(fit))|
                               coefNames(fit)%in%ctr$variable]
    if(any(grepl(":",ctr$variable))){
      contrasts=coefNames(fit)[grepl(':',coefNames(fit))]
      contrasts=union(contrasts,unlist(sapply(contrasts,function(x)strsplit(x,':'))))
    }  
    ##  linear (.L), the second is quadratic (.Q), the third is cubic (.C)
    df_de=do.call(rbind,
                  lapply(contrasts,function(x)fit %>%
                           topTable(coef=x, number=Inf) %>%
                           as_tibble %>% 
                           group_by(assay) %>% 
                           summarise(
                             coef=x,
                             nGenes = length(adj.P.Val), 
                             nDE = sum(adj.P.Val < 0.05),
                             pi1 = 1 - qvalue(P.Value)$pi0) %>%
                           mutate(assay = factor(assay, ctorder))))
    DEs=do.call(rbind,lapply(contrasts,function(x)data.frame(coef=x,topTable(fit,coef=x, number=Inf))))
    save(fit,df_de,file=paste0(paste0(ctr$variable,collapse='_'),'.Rdata'))
    write.csv(DEs,file=paste0(paste0(ctr$variable,collapse='_'),'.csv'))
    write.csv(df_de,file=paste0(paste0(ctr$variable,collapse='_'),'_nde.csv'))
  }else{
    df_de = fit %>%
      topTable(coef=ctr$contrasts, number=Inf) %>%
      as_tibble %>% 
      group_by(assay) %>% 
      summarise( 
        nGenes = length(adj.P.Val), 
        nDE = sum(adj.P.Val < 0.05),
        pi1 = 1 - qvalue(P.Value)$pi0) %>%
      mutate(assay = factor(assay, ctorder))
    DEs=data.frame(coef=ctr$contrasts,topTable(fit,coef=ctr$contrasts, number=Inf))
    save(fit,df_de,file=paste0(ctr$variable,'.Rdata'))
    write.csv(DEs,file=paste0(ctr$variable,'.csv'))
    write.csv(df_de,file=paste0(ctr$variable,'_nde.csv'))
  }
  
  
}

contrasts=data.frame(variable=c('Dx_contrast','Dx_dev','Dx_donor','Dx_cell'),
                     contrasts=NA)
for(x in 1:nrow(contrasts)){
  ctr=contrasts[x,]
  if(is.na(ctr$contrasts)){
    ctr=list(variable=ctr$variable,
             contrasts=NULL)  
  }else{
    ctr=list(variable=ctr$variable,
             contrasts=ctr$contrasts)  
    
  }
  covariate_use=covariates
  if(grepl("BrainRegion",ctr$variable)){
    covariate_use=setdiff(covariate_use,'(1|BrainRegion)')
  }
  if(grepl('\\:',ctr$variable)){
    covariate_use=union(unlist(strsplit(ctr$variable,":")),covariate_use)
  }
  
  try(test())
}

