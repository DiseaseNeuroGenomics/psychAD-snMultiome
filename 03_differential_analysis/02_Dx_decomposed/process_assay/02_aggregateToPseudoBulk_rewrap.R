.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(SingleCellExperiment)
library(rhdf5)
library(Matrix)
library(tidyverse)
library(dreamlet)
library(variancePartition)
library(ggplot2)
##
format_pb <- function(pb,rm_rb_mt=T){
  cnames=unique(unlist(lapply(pb,colnames)))
  rnames=unique(unlist(lapply(pb,rownames)))
  pb1=list()
  for(x in names(pb)){
    mat=array(0,dim=c(length(rnames),length(cnames)),dimnames = list(rnames,cnames))
    mat[rownames(pb[[x]]),colnames(pb[[x]])]=pb[[x]]
    pb1[[x]]=as(mat,'sparseMatrix')
  }
  if(rm_rb_mt){
    keep=!(grepl('^RPL',rnames) | grepl('^RPS',rnames) |grepl('^MT-',rnames))
    pb1=lapply(pb1,function(x)x[keep,])
  }
  return(pb1)
} # from dreamlet
myaggregateToPseudoBulk <- function(pb,Adobs,Advar,by){
  Advar=S4Vectors::DataFrame(Advar[rownames(pb[[1]]),])
  Adobs=S4Vectors::DataFrame(Adobs)
  md=list()
  md$agg_pars <- list(assay = 'counts', by = by, fun = 'sum', scale = FALSE)
  pb <- SingleCellExperiment(pb, rowData = Advar, metadata = md)
  cd <- data.frame(Adobs[, by, drop = FALSE])
  for (i in names(cd)) {
    if (is.factor(cd[[i]])) {
      cd[[i]] <- droplevels(cd[[i]])
    }
  }
  ns <- table(cd)
  if (length(by) == 2) {
    ns <- asplit(ns, 2)
    ns <- purrr::map(ns, ~ c(unclass(.)))
  } else {
    ns <- c(unclass(ns))
  }
  int_colData(pb)$n_cells <- ns
  if (length(by) == 2) {
    cd <- Adobs
    ids <- colnames(pb)
    
    # keep columns if there is a unique value in each sample
    counts <- vapply(ids, function(u) {
      m <- as.logical(match(cd[, by[2], drop = TRUE], u, nomatch = 0))
      vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
    }, numeric(ncol(Adobs)))
    cd_keep <- apply(counts, 1, function(u) all(u == 1))
    cd_keep <- setdiff(names(which(cd_keep)), by)
    
    if (length(cd_keep) != 0) {
      m <- match(ids, cd[, by[2]], nomatch = 0)
      cd <- cd[m, cd_keep, drop = FALSE]
      rownames(cd) <- ids
      colData(pb) <- cd
    }
  }
  
  # Per sample, per cell type
  # compute mean of each continuous column in colData
  
  # Get columns that are numeric (or integers), but not factor
  cols <- vapply(Adobs, function(a) is.numeric(a) & !is.factor(a), logical(1))
  cols <- names(cols[cols])
  # exclude columns in by and already in colData(pb)
  cols <- cols[!(cols %in% by)]
  cols <- cols[!cols %in% colnames(colData(pb))]
  
  # per sample, per cell type
  # report mean of value
  metadata(pb)$aggr_means <- Adobs %>%
    as_tibble() %>%
    group_by_at(vars(metadata(pb)$agg_pars$by)) %>%
    summarise_at(cols, mean, na.rm = TRUE)
  
  return(pb)
} # from dreamlet
read_h5_mat<- function(fh,name){
  z=h5read(fh,name)
  mat=t(z$block0_values)
  colnames(mat)=z$axis0
  rownames(mat)=z$axis1
  return(mat)
}
read_pd_hdf5 <- function(fh){
  name=h5ls(fh,recursive=F)$name
  ncells=list()
  Mats=list()
  for(x in name){
    Mats[[x]]=read_h5_mat(fh,paste0(x,'/count_matrix'))
    ncells[[x]]=read_h5_mat(fh,paste0(x,'/ncell'))
  }
  list(expMats=Mats,
       ncells=ncells)
}
read_obsvar <- function(path,sub='obs'){
  dat = h5read(path,sub)
  l=sapply(dat,length)
  dat1=do.call(data.frame,lapply(dat[l==2],function(x)c(NA,x$categories)[x$codes+2]))
  dat2=do.call(data.frame,dat[l!=2])
  dat=data.frame(dat1,dat2)
  #rownames(dat)=dat$X_index
  rownames(dat)=dat[,grepl('key$',colnames(dat))]
  dat
}
NCPU=16
args=commandArgs(T)
if (length(args) !=5) {
  stop("Rscript 02_processAssays_rewrap.R path/h5ad path/Dx_cell path/metaInfo path/pseudobulk pseudobulk-prefix", call. = FALSE)
}
h5ad_path=args[1]
Dx_cell_path=args[2]
metaInfo_path=args[3]
pb_folder=args[4]
prefix=args[5]
sampleCol='order'
###

allInfo=read.csv(metaInfo_path,row.names=1)
a=b=c=allInfo
a$Dx_cat_cell='control'
b$Dx_cat_cell='early'
c$Dx_cat_cell='late'
rownames(a)=a$ID=paste0(a[,sampleCol],'_',a$Dx_cat_cell) # order is sampleID
rownames(b)=b$ID=paste0(b[,sampleCol],'_',b$Dx_cat_cell) # order is sampleID
rownames(c)=c$ID=paste0(c[,sampleCol],'_',c$Dx_cat_cell) # order is sampleID
metaInfo=rbind(a,b,c)

Advar=read_obsvar(h5ad_path,'var')
Adobs=read_obsvar(h5ad_path,'obs')
anno=read.csv(Dx_cell_path,row.names=1)

sel_index=intersect(rownames(Adobs),rownames(anno))
Adobs=data.frame(Adobs[sel_index,],anno[sel_index,c('subclass','class','subtype','Dx_cell','Dx_cat_cell')])
Adobs$sampleDx=paste0(Adobs[,sampleCol],'_',Adobs$Dx_cat_cell)
setwd(pb_folder)

for(level in c('class','subclass','subtype')){
  hdf5_path=paste0(prefix,'_',level,'_pb_persample.h5')
  allmats=read_pd_hdf5(hdf5_path)
  pb=format_pb(allmats[[1]])
  pb1=myaggregateToPseudoBulk(pb,Adobs,Advar,c(level,'sampleDx'))
  N_CELLS=pb1@int_colData$n_cells[rownames(pb1@colData)]
  pb1=pb1[,rownames(pb1@colData)%in%rownames(metaInfo)]
  pb1@int_colData$n_cells=N_CELLS[rownames(pb1@colData)]
  pb1@colData=S4Vectors::DataFrame(pb1@colData,metaInfo[rownames(pb1@colData),])
  save(pb1,metaInfo,file=paste0('pb4dreamlet_',prefix,'_',level,'.Rdata'))
}
