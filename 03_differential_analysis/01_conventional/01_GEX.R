.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(dreamlet)
library(SingleCellExperiment)
library(tidyverse)
library(zenith)
library(qvalue)
library(cowplot)
library(splines)
# python process_assay/01_pb_sampleDx.py  --input h5ad --cellanno cellanno.csv --outDir OUTDIR --outPrefix RNA
# Rscript process_assay/02_processAssay_rewrap.R h5ad cellanno.csv path_metaInfo OUTDIR RNA

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
level='class'
BR='all'
OUTDIR='/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/DE3/'
prefix='RNA'
cellanno_path='/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno_05092025_fixed3.csv'
sampleCol='order'
NCPU=32
bpparam = MulticoreParam(NCPU)
setwd(file.path(OUTDIR,level))
##########################per cell type test#############################################

load(paste0('pb4dreamlet_',prefix,'_',level,'.Rdata'))
#voom
form=~(1|BrainRegion)+(1|SubID)+scale(mito_genes)+scale(n_genes_by_counts)+scale(pct_counts_ribo)+(1|Sex)+scale(Age)+scale(Dx_donor)
proc <- processAssays(pb1, form,min.count = 5, BPPARAM = bpparam)
z=reshape2::melt(do.call(rbind,pb1@int_colData$n_cells))
rownames(z)=with(z,paste(Var1,Var2,sep='_'))
proc@metadata$ncell=z[paste0(proc@metadata[,sampleCol],'_',unlist(proc@metadata[,level])),]$value
cdmk_folder(prefix)
save(proc,file='proc.Rdata')
#dreamlet
covariates=c("(1|SubID)",'BrainRegion',
             "scale(n_genes_by_counts)",
             "scale(pct_counts_ribo)",
             "scale(mito_genes)",
             "scale(Age)", "Sex", "scale(PMI..min.)")
proc@data$Dx_cat_donor=factor(proc@data$Dx_cat_donor,levels=c('control','early','late'),ordered=T)
proc@data$Dx_donor=scale(proc@data$Dx_donor)

form1 = as.formula(paste0("~  ", 'scale(Dx_donor)', " + ", paste0(gsub('BrainRegion','(1|BrainRegion)',gsub('Sex','(1|Sex)',covariates)), collapse=" + ")))  
fit1=fitVarPart(proc,form1,BPPARAM = bpparam)
pdf('varpart.pdf',width=8,height = length(unique(fit1$assay))/0.5+5)
print(plotVarPart(fit1,label.angle =45))
dev.off()
save(fit1,file='varpart.Rdata')

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

contrasts=data.frame(variable=c('Dx_donor'),
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
##########################stacked assay #######################################
bpparam = MulticoreParam(8) # use less cpu
res.proc <- processAssays(pb.stack, form, min.count = 10,min.prop=0.01, BPPARAM = bpparam)
form=~(1|Dx_cat_donor)+(1|Sex)+scale(Age)+scale(PMI..min.)+(1|SubID)+(1|BrainRegion)+(1|class)+ (1|class:stackedAssay)+scale(mito_genes)+scale(n_genes_by_counts)+scale(mito_ribo)
varPart=fitVarPart(res.proc,form2, BPPARAM = bpparam)
save(varPart,file='stackedassay_varpart.Rdata')
