.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(Matrix)
library(BiocParallel)
library(chromVAR)
#library(chromVARmotifs)
library(motifmatchr)
library(rhdf5)
library(TFBSTools)
library(ggplot2)
library(variancePartition)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)

source('/sc/arion/projects/roussp01a/pengfei/tools/myut.R')
register(MulticoreParam(24))
##
folder=commandArgs(T)[1]

#########################################################################################
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
#pwmList=pwmList[intersect(bg,names(pwmList))]
###############################################
motif_clust=read.csv('/sc/arion/projects/roussp01a/pengfei/ref/clustermotifs_output/motif_comparison_03_clusters.csv')
rownames(motif_clust)=motif_clust$motif
clust_name=read.csv('/sc/arion/projects/roussp01a/pengfei/ref/clustermotifs_output/cluster_03_custom_name.csv',row.names=1)
motif_clust$custom_name=clust_name[as.character(motif_clust$gp),]$name
setwd(paste0('/sc/arion/projects/CommonMind/roussp01a/snmulti/motif/files/class/',folder))
mat=Seurat::Read10X_h5('RUSHMR_metacell_20.h5')
Obs=read.csv('RUSHMR_metacell_20.obs.csv.gz',row.names=1)
Var=read.csv('RUSHMR_metacell_20.var.csv.gz',row.names=1)


peak=StringToGRanges(rownames(Var))
sel=rowSums(mat)>10
#sel=(overlapsAny(peak,refpeak) | rowSums(mat)/ncol(mat)/20>0.015) & rowSums(mat)>10
Var=Var[sel,]
mat=mat[sel,]
peak=peak[sel]
#mat=read.csv('WNN_embed_pb.csv.gz',row.names=1)
#cov=read.csv('WNN_embed_cov.csv.gz',row.names=1)
#rownames(info)=paste0('X',info$bin)
out <- SummarizedExperiment(assays = list(counts = mat), 
                            rowRanges = peak)
out <- addGCBias(out,genome = BSgenome.Hsapiens.UCSC.hg38)
bg <- getBackgroundPeaks(object = out)
motif_ix <- matchMotifs(pwmList, out, 
                        genome = BSgenome.Hsapiens.UCSC.hg38)
# computing deviations
dev <- computeDeviations(object =out, annotations = motif_ix,
                         background_peaks=bg)
dv=deviations(dev)
ds=deviationScores(dev)
save(dev,ds,dv,file='RUSHMR_metacell_20_chromvar.Rdata')
write.csv(ds,file=paste0('RUSHMR_metacell_20_chromvar.csv'))

