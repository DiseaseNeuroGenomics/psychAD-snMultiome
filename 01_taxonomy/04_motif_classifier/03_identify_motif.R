.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(data.table)
library(Matrix)
library(glmnet)
library(FSelectorRcpp)    # mRMR
library(pROC)
library(future)
library(future.apply)
library(mRMRe)
library(dplyr)
library(ggplot2)
options(future.globals.maxSize= 60*1024^3)
plan(multicore, workers = 32)
cell=commandArgs(T)[1]

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

setwd(paste0('/sc/arion/projects/CommonMind/roussp01a/snmulti/motif/files/class/',cell))
devZ=read_big_csv('metacell_20_chromvar.csv')
meta=gsub('_pbulk_.*','',colnames(devZ))


motif_clust=read.csv('/sc/arion/projects/roussp01a/pengfei/ref/clustermotifs_output/motif_comparison_03_clusters.csv')
rownames(motif_clust)=motif_clust$motif
clust_name=read.csv('/sc/arion/projects/roussp01a/pengfei/ref/clustermotifs_output/cluster_03_custom_name.csv',row.names=1)
motif_clust$custom_name=clust_name[as.character(motif_clust$gp),]$name
motif_clust=motif_clust[rownames(devZ),]

X <- t(devZ)
y <- droplevels(factor(meta))

############ merge similar motifs################################################
group_reps <- vapply(split(seq_along(motif_clust$motif), motif_clust$gp), function(ix){
  subX <- X[, ix, drop=FALSE]
  # simple score: mean absolute ANOVA F-statistic proxy via one-vs-rest t-stat
  sc <- apply(subX, 2, function(v) {
    # pooled one-vs-rest across classes: mean(|t|)
    mean(sapply(levels(y), function(cls){
      v1 <- v[y == cls]; v0 <- v[y != cls]
      if(sd(v1)==0 || sd(v0)==0) return(0)
      abs(t.test(v1, v0)$statistic)
    }))
  })
  ix[which.max(sc)]
}, integer(1))
X_rep <- X[, group_reps, drop=FALSE]

###################split train test###############
set.seed(123)
idx <- split(seq_len(nrow(X_rep)), y)                  # stratify
test_idx <- unlist(lapply(idx, function(v) sample(v, max(1, floor(0.2*length(v))))))
train_idx <- setdiff(seq_len(nrow(X_rep)), test_idx)
Xtr <- X_rep[train_idx, , drop=FALSE]
ytr <- y[train_idx]
Xte <- X_rep[test_idx, , drop=FALSE]
yte <- y[test_idx]

###########multinomial lasso (sparse)########################
cvfit <- cv.glmnet(x = as.matrix(Xtr), y = ytr, family = "multinomial",
                   alpha = 1, type.multinomial = "ungrouped",
                   nfolds = 5, parallel = FALSE)
fit <- glmnet(x = as.matrix(Xtr), y = ytr, family = "multinomial",
              alpha = 1, type.multinomial = "ungrouped")
beta_1se <- coef(fit, s = cvfit$lambda.1se)
# Extract nonzero motifs per class (one-vs-rest coefficients)
nz_by_class <- lapply(beta_1se, function(b) {
  nz <- rownames(b)[abs(as.numeric(b)) > 1e-8]
  setdiff(nz,"(Intercept)")
})

set.seed(42)
B <- 100
stable <- lapply(levels(ytr), function(cls) character(0))  # initialize list

boot_motifs <- future_lapply(seq_len(B), function(b) {
  boot_idx <- unlist(
    lapply(split(seq_along(ytr), ytr),
           function(v) sample(v, length(v), replace = TRUE))
  )
  cvb <- cv.glmnet(
    x = as.matrix(Xtr[boot_idx, , drop = FALSE]),
    y = ytr[boot_idx],
    family = "multinomial",
    alpha = 1,
    type.multinomial = "ungrouped",
    nfolds = 5,
    parallel = FALSE
  )
  fb <- glmnet(
    x = as.matrix(Xtr[boot_idx, , drop = FALSE]),
    y = ytr[boot_idx],
    family = "multinomial",
    alpha = 1,
    type.multinomial = "ungrouped"
  )
  res <- lapply(levels(ytr), function(cls) {
    beta <- coef(fb, s = cvb$lambda.1se)[[cls]]
    nz <- setdiff(rownames(beta)[abs(as.numeric(beta)) > 1e-8],"(Intercept)")
    nz
  })
  names(res) <- levels(ytr)
  return(res)   
})

stable <- lapply(levels(ytr), function(cls) {
  sel <- unlist(lapply(boot_motifs, function(m) m[[cls]]))
  counts <- table(sel)
  keep=names(counts[counts >= 0.6 * B])  # >= 60% selection frequency
  #keep=keep[sapply(keep,function(x)mean(Xtr[ytr==cls,x])>0)]#
})
names(stable) <- levels(ytr)

############### 6) mRMR to pick top  per subclass####################
pick_top_mrmr <- function(X, y_binary, k = 3) {
  # X: matrix/data.frame of motif Z-scores (meta-cells × candidate motifs)
  # y_binary: 0/1 vector (one-vs-rest label)
  df <- data.frame(X, label = as.numeric(y_binary))   # last column = target
  dd <- mRMR.data(data = df)
  fit <- mRMR.classic(data = dd,
                      target_indices = ncol(df),
                      feature_count = min(k, ncol(df)-1))
  sel_idx <- solutions(fit)[[1]]
  setdiff(colnames(df)[sel_idx], "label")
}

final_motifs <- lapply(levels(ytr), function(cls){
  cand <- intersect(stable[[cls]], colnames(Xtr))
  if (length(cand) <= 3) return(cand)
  pick_top_mrmr(Xtr[, cand, drop = FALSE], ytr == cls, k =3)
})
names(final_motifs) <- levels(ytr)

z=setdiff(unlist(final_motifs),rownames(devZ))
if(length(z)>0){
  a=gsub('\\.','-',z)
  names(a)=z
  final_motifs=lapply(final_motifs,function(x){if(names(a)%in%x){x=c(a,setdiff(x,names(a)))}else{x}})
}

score_and_auc <- lapply(names(final_motifs), function(cls){
  feats <- final_motifs[[cls]]
  if(length(feats)==0) return(list(class=cls, auc=NA))
  # simple linear score from lasso coefficients refit on train with only feats
  if(length(feats)==1){
    df_train <- data.frame(y =ytr == cls, Xtr[, feats,drop=F])
    subfit <- glm(y ~ ., data = df_train, family = binomial)
    pred <- predict(subfit, newdata = data.frame(Xte[, feats,drop=F]), type = "response")
    coef=coef(subfit)
  }else{
    subfit <- glmnet(as.matrix(Xtr[, feats, drop=FALSE]), ytr == cls, family="binomial", alpha=1)
    lam <- if(!is.null(subfit$lambda)) subfit$lambda[length(subfit$lambda)] else 0
    pred <- predict(subfit, newx=as.matrix(Xte[, feats, drop=FALSE]), s=lam, type="response")[,1]
    coef=coef(subfit,s=lam)
  }
  truth <- as.numeric(yte == cls)
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
})
names(score_and_auc)=names(final_motifs)

dev=read_big_csv('metacell_20_chromvar.csv')
UMAP=read_big_csv('metacell_20.obs.csv.gz')
meta=gsub('_pbulk_.*','',colnames(devZ))
devZ=t(scale(t(dev)))
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/dendro.Rdata')
ord=subclass_order[subclass_order%in%meta]
if(cell%in%c('Micro_PVM','Astro','Oligo','OPC')){ord=unique(meta)}
a=model.matrix(~0+meta)
colnames(a)=gsub('meta','',colnames(a))
a=t(t(a)/colSums(a))
zmat=as.matrix(devZ[unlist(final_motifs),])%*%a
ord2=unique(unlist(final_motifs[rev(ord[ord%in%names(final_motifs)])]))#rownames(zmat)[order.dendrogram(as.dendrogram(hc1))]
#colnames(zmat)[order.dendrogram(as.dendrogram(hc2))]
pdf('subclass_motif.pdf',width=5,height=5)
zdf=reshape2::melt(t(scale(t(zmat))))
zdf$Var1=factor(zdf$Var1,levels = ord2,ordered=T)
zdf$Var2=factor(zdf$Var2,levels = ord,ordered=T)
print(ggplot(zdf)+geom_tile(aes(x=Var1,y=Var2,fill=value))+xlab('')+ylab('')+
        theme_bw() +scale_fill_gradient2(low='blue',
                                         high="red")+
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0)) +
        coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        theme(legend.position="bottom"))
print(ggplot(zdf[zdf$Var2!='EN_L5_ET' & zdf$Var1%in%unlist(final_motifs[setdiff(names(final_motifs),'EN_L5_ET')]),])+
               geom_tile(aes(x=Var1,y=Var2,fill=value))+xlab('')+ylab('')+
               theme_bw() +scale_fill_gradient2(low='blue',
                                         high="red")+
        scale_y_discrete(expand = c(0, 0)) +
        scale_x_discrete(expand = c(0, 0)) +
        coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        theme(legend.position="bottom"))
roc=do.call(rbind,lapply(names(score_and_auc),function(x){a=score_and_auc[[x]]$roc;a$cell=x;a}))
auc=data.frame(AUC=sapply(score_and_auc,function(x)x$AUC),
               cell=names(score_and_auc))
pr=do.call(rbind,lapply(names(score_and_auc),function(x){a=score_and_auc[[x]]$pr;a$cell=x;a}))
pauc=data.frame(PRAUC=sapply(score_and_auc,function(x)x$PRAUC),
                cell=names(score_and_auc))
auc=auc[order(pauc$PRAUC),]
auc$y=(1:nrow(auc))/(nrow(auc)+1)
pauc=pauc[order(pauc$PRAUC),]
pauc$y=(1:nrow(pauc))/(nrow(pauc)+1)
f=table(meta)
f=round(f/sum(f),3)
pauc$base=as.numeric(f[rownames(pauc)])
write.csv(auc,file='auc.csv')
write.csv(pauc,file='prauc.csv')
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
dev.off()
save(stable,final_motifs,devZ,meta,score_and_auc,file='subclass.Rdata')


devZ=t(devZ[,rownames(UMAP)])
dev=t(dev[,rownames(UMAP)])
pdf('subclass_UMAP.pdf',width=5,height = 5)

for(x in names(final_motifs)){
  
  tf=final_motifs[[x]]
  for(a in tf){
    umap=data.frame(UMAP,TF=as.numeric(devZ[,a]))
    umap$TF=pmin(quantile(umap$TF,0.995),umap$TF)
    umap$TF=pmax(quantile(umap$TF,1-0.995),umap$TF)
    print(ggplot(umap)+geom_point(aes(x=UMAP1,y=UMAP2,color=TF),shape='.')+
            scale_color_gradient2(low='blue',high='red',mid='grey')+
            theme_bw()+ggtitle(paste0(x,' ',a))+
            theme(axis.text.x = element_text(colour = "black"),
                  axis.text.y = element_text(colour = "black"),
                  strip.background = element_rect(fill = NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),aspect.ratio = 1))
  }
  if(length(tf)>1){
    umap=data.frame(UMAP,TF=dev[,tf]%*%score_and_auc[[x]]$coef[tf,]+score_and_auc[[x]]$coef['(Intercept)',1])
    umap$TF=1 / (1 + exp(-umap$TF))
    #umap$TF=pmin(quantile(umap$TF,0.995),umap$TF)
    #umap$TF=pmax(-max(umap$TF),umap$TF)
    print(ggplot(umap)+geom_point(aes(x=UMAP1,y=UMAP2,color=TF),shape='.')+
            #scale_color_gradient2(low='blue',high='red',mid='grey')+
            scale_color_gradientn(colors=c('grey',"#22A884FF","#414487FF"),limit=c(0,1))+
            theme_bw()+ggtitle(x)+
            theme(axis.text.x = element_text(colour = "black"),
                  axis.text.y = element_text(colour = "black"),
                  strip.background = element_rect(fill = NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),aspect.ratio = 1))
  }
}
UMAP$cell=gsub('_pbulk_.*','',rownames(UMAP))

print(ggplot(UMAP)+geom_point(aes(x=UMAP1,y=UMAP2,color=cell),shape='.')+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
dev.off()


