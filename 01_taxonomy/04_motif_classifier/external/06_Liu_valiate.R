.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(data.table)
library(Matrix)
library(glmnet)
library(pROC)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(scales)
read_big_csv <- function(fh,row1asrowname=T){
  library(data.table)
  dat=fread(fh)
  dat=data.frame(dat,check.names = F)
  if(row1asrowname){
    rownames(dat)=dat[,1]
    dat[,1]=NULL
  }
  return(dat)
}
plot_split_marker_heatmap <- function(
    mat1,mat2,marker_list, ord,mark_selected = TRUE) {
  study1_name = "A"; study2_name = "B"
  marker_order <- unique(unlist(marker_list[rev(ord[ord%in%names(marker_list)])],use.names = FALSE))
  mat1 <- mat1[marker_order, ord, drop = FALSE]
  mat2 <- mat2[marker_order, ord, drop = FALSE]
  
  nr <- nrow(mat1)
  nc <- ncol(mat1)
  cell_df <- expand.grid( i = seq_len(nr), j = seq_len(nc) )
  
  df1 <- pmap_dfr(cell_df, function(i, j) {
    tibble(x = c(j - 0.5, j + 0.5, j - 0.5),
           y = c(i - 0.5, i - 0.5, i + 0.5),
           id = paste0(i, "_", j, "_study1"),
           value = mat1[i, j],
           study = study1_name)})
  
  df2 <- pmap_dfr(cell_df, function(i, j) {
    tibble(x = c(j + 0.5, j + 0.5, j - 0.5),
           y = c(i - 0.5, i + 0.5, i + 0.5),
           id = paste0(i, "_", j, "_study2"),
           value = mat2[i, j],
           study = study2_name)})
  
  plot_df <- bind_rows(df1, df2)
  
  marker_df <- imap_dfr(marker_list, function(marks, clust) {
    tibble(clust=clust,mark=marks,
           rank = seq_along(marks))  }) %>%
    filter(clust %in% ord,
           mark %in% marker_order) %>%
    mutate(x = match(clust, ord),
           y = match(mark, marker_order))
  
  lim <- max(abs(c(mat1, mat2)), na.rm = TRUE)
  
  p <- ggplot(plot_df, aes(x = y, y = x)) +geom_polygon(
    aes(fill = value, group = id),color = "grey85",linewidth = 0.15) +
    scale_fill_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0,
                         limits = c(-lim, lim),oob = squish, na.value = "grey",name = "z-score") +
    scale_y_continuous(
      breaks = seq_len(nc),
      labels = ord,
      expand = c(0, 0)) +
    scale_x_continuous(
      breaks = seq_len(nr),
      labels = marker_order,
      expand = c(0, 0)) +
    coord_fixed() +theme_bw() + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          strip.background = element_rect(fill = NA))
  
  if (mark_selected) {
    p <- p +
      geom_tile(
        data = marker_df,
        aes(x = y, y = x),
        inherit.aes = FALSE,fill = NA,color = "black",linewidth = 0.35, 
        width = 0.99, height = 0.99)
  }
  return(p)
}
#cell='EN'
cell=commandArgs(T)[1]
setwd(paste0('/sc/arion/projects/CommonMind/roussp01a/snmulti/motif/files/class/',cell))
prefix='RUSHMR'
cells=cell
anno=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno_05092025_fixed3.csv')
anno1=anno[anno$inRNA & anno$class%in%cells,]
anno1=anno1[anno1$subclass!='Epend',]

scvi1=read.csv(paste0('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/scanvi_',prefix,'_',cell,'_all.csv'),row.names=1)
scvi1=scvi1[rownames(scvi1)%in%rownames(anno1),,drop=F]
qu=anno1[rownames(scvi1),]$subclass
ref=scvi1[,1] 
z=apply(table(qu,ref),2,c)
zdf=reshape::melt(z/rowSums(z))

zdf1=as.data.frame(zdf%>%group_by(ref)%>%slice_max(order_by = value, n = 1, with_ties = FALSE))
rownames(zdf1)=zdf1$ref
map1=zdf1
###################
load('subclass.Rdata')
dev=devZ
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/dendro.Rdata')
#ord=subclass_order[subclass_order%in%meta]
Mdev=read_big_csv('RUSHMR_metacell_20_chromvar.csv')
MUMAP=read_big_csv('RUSHMR_metacell_20.obs.csv.gz')
meta=gsub('_pbulk_.*','',colnames(dev))
Mmeta=map1[gsub('_pbulk_.*','',colnames(Mdev)),1]
####AUC

get_auc <- function(cls){
  feats <- final_motifs[[cls]]
  
  subfit=score_and_auc[[cls]]$fit
  lam <- if(!is.null(subfit$lambda)) subfit$lambda[length(subfit$lambda)] else 0
  pred <- predict(subfit, newx=as.matrix(t(Mdev)[, feats, drop=FALSE]), s=lam, type="response")[,1]
  coef=coef(subfit,s=lam)
  truth <- as.numeric(Mmeta == cls)
  roc_obj <- pROC::roc(truth, pred)
  auc_value <- pROC::auc(roc_obj)
  
  
  pr <- PRROC::pr.curve(scores.class0 = pred[truth == 1],
                        scores.class1 = pred[truth == 0],
                        curve = TRUE)
  list(
    fit=subfit,
    coef=coef,
    motifs   = feats,
    AUC      = auc_value,
    roc=data.frame(data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities
    )),
    pr=data.frame(
      Recall = pr$curve[,1],
      Precision = pr$curve[,2]
    ),
    PRAUC    = pr$auc.integral)
}


auc_summary=lapply(intersect(names(final_motifs),unique(Mmeta)),get_auc)
names(auc_summary)=intersect(names(final_motifs),unique(Mmeta))
roc=do.call(rbind,lapply(names(auc_summary),function(x){a=auc_summary[[x]]$roc;a$cell=x;a}))
auc=data.frame(AUC=sapply(auc_summary,function(x)x$AUC),
               cell=names(auc_summary))
pr=do.call(rbind,lapply(names(auc_summary),function(x){a=auc_summary[[x]]$pr;a$cell=x;a}))
pauc=data.frame(PRAUC=sapply(auc_summary,function(x)x$PRAUC),
                cell=names(auc_summary))
auc=auc[order(pauc$PRAUC),]
auc$y=(1:nrow(auc))/(nrow(auc)+1)
pauc=pauc[order(pauc$PRAUC),]
pauc$y=(1:nrow(pauc))/(nrow(pauc)+1)
f=table(Mmeta)
f=round(f/sum(f),3)
pauc$base=as.numeric(f[rownames(pauc)])
write.csv(pauc,file='RUSHMR_prauc.csv')
write.csv(auc,file='RUSHMR_auc.csv')
pdf('RUSHMR_subclass_motif.pdf',width=5,height=5)
print(ggplot(roc)+geom_line(aes(x=FPR,y=TPR,col=cell))+
        theme_bw()+
        geom_text(data=auc,aes(y=y,label=paste0(gsub("EN_|IN_", "", cell),': ',round(AUC,2)),col=cell),x=0.75)+
        xlab('False postive rate')+ylab('True postive rate')+
        scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))

print(ggplot(pr)+geom_line(aes(x=Recall,y=Precision,col=cell))+
        theme_bw()+
        scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))+
        geom_text(data=pauc,aes(y=y,label=paste0(gsub("EN_|IN_", "", cell),': ',round(PRAUC,2)),col=cell),x=0.25)+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
print(ggplot(pr)+geom_line(aes(x=Recall,y=Precision,col=cell))+
        theme_bw()+
        scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))+
        geom_text(data=pauc,aes(y=y,label=paste0(gsub("EN_|IN_", "", cell),': ',round(PRAUC,2),
                                                 ' (',round(base,2),')'),col=cell),x=0.25)+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
####heatmap
ord=subclass_order[subclass_order%in%meta]
######
a=model.matrix(~0+meta)
colnames(a)=gsub('meta','',colnames(a))
a=t(t(a)/colSums(a))
zmat=as.matrix(dev[unique(unlist(final_motifs)),])%*%a
zmat=t(scale(t(zmat)))
################################################################################

a=model.matrix(~0+Mmeta)
colnames(a)=gsub('Mmeta','',colnames(a))
a=t(t(a)/colSums(a))
zmat2=as.matrix(Mdev[unique(unlist(final_motifs)),])%*%a
zmat2=t(scale(t(zmat2)))

rho=cor(reshape2::melt(zmat[rownames(zmat2),colnames(zmat2)])$value,
        reshape2::melt(zmat2)$value,method='spearman',use='complete.obs')
zmat=pmax(pmin(zmat,2),-2)
zmat2=pmax(pmin(zmat2,2),-2)
#colnames(zmat)[order.dendrogram(as.dendrogram(hc2))]
print(plot_split_marker_heatmap(zmat,zmat2,final_motifs, ord)+
        ggtitle(paste0('rho: ',round(rho,3))))
dev.off()

