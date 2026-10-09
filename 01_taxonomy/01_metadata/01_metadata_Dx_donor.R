.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(dreamlet)
library(ordinal)
library(ggplot2)

setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/metadata')
metaInfo=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/files/allInfo_11102023.csv',row.names=1)
full_meta=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step2/metadata/clinical_metadata.csv')

sel_IDs=c(metaInfo$SubID,
          full_meta[which(full_meta$Brain_bank=='MSSM' & full_meta$Age>59 & rowSums(full_meta[,c('SCZ_ALL','BD_ALL','PD','DLBD','FTD','Vascular','Tumor')])==0), ]$SubID)
rownames(psychAD_meta)=psychAD_meta$SubID
psychAD_meta=full_meta[rownames(full_meta)%in%sel_IDs,]
Dx_anno=c("Age","Sex",'AD',"Ethnicity","CERAD","BRAAK_AD","BRAAK_PD","CDRScore","Plq_Mn","Plq_Mn_MFG","MCI","Dementia","Cognitive_Resilience","Cognitive_and_Tau_Resilience")
sel_Dx=c("CERAD","BRAAK_AD","Plq_Mn",
         "MidPlaquesValue","MidTanglesValue","EntorPlaquesValue",
         "EntorTanglesValue","SupPlaquesValue")
###
input=psychAD_meta[rowSums(is.na(psychAD_meta[,sel_Dx]))==0,sel_Dx]
svd1=svd(scale(input))
pdf('PCA_pathonly.pdf')
plot(svd1$d^2/sum(svd1$d^2))
plot(cumsum(svd1$d^2)/sum(svd1$d^2))
dev.off()

psychAD_meta$DxPC1=NA
psychAD_meta[rownames(input),]$DxPC1=svd1$u[,1]
NAs=intersect(psychAD_meta$SubID[is.na(psychAD_meta$DxPC1)],sel_IDs)
for(sp in NAs){
  nonNA=sel_Dx[which(!is.na(psychAD_meta[sp,sel_Dx]))]
  zdf=psychAD_meta[,c('DxPC1',nonNA)]
  lm1=lm(DxPC1~.,zdf)
  dxpc1=predict(lm1,zdf[sp,])
  psychAD_meta[sp,'DxPC1']=dxpc1
  print(psychAD_meta[sp,c('DxPC1',nonNA)])
}
##### control/early/late
psychAD_meta$seed <- NA
psychAD_meta$seed[which(psychAD_meta$CERAD>1 & psychAD_meta$BRAAK_AD >2 & 
                          psychAD_meta$CERAD+psychAD_meta$BRAAK_AD>4 &
                          psychAD_meta$CERAD+psychAD_meta$BRAAK_AD<10)]  <- "early"
psychAD_meta$seed[which(psychAD_meta$CERAD==1 & psychAD_meta$BRAAK_AD<3)]  <- "control"
psychAD_meta$seed[which(psychAD_meta$CERAD>3 & psychAD_meta$BRAAK_AD>4)]  <- "late"


psychAD_meta$seed <- factor(psychAD_meta$seed, levels = c("control", "early", "late"), ordered = TRUE)
psychAD_meta$DxPC1_z <- scale(psychAD_meta$DxPC1)[,1]
fit <- clm(seed ~ DxPC1_z,
           data = psychAD_meta[!is.na(psychAD_meta$seed), ],
           link = "probit")
pr <- predict(fit, newdata = psychAD_meta[, setdiff(names(psychAD_meta), "seed")], type = "prob")$fit
psychAD_meta$P_control <- pr[, "control"]
psychAD_meta$P_early   <- pr[, "early"]
psychAD_meta$P_late    <- pr[, "late"]

labs <- c("control", "early", "late")
psychAD_meta$Dx_cat <- labs[max.col(pr, ties.method = "first")]
cutoff=c((max(psychAD_meta$DxPC1[which(psychAD_meta$Dx_cat=='control')],na.rm=T)+min(psychAD_meta$DxPC1[which(psychAD_meta$Dx_cat=='early')],na.rm=T))/2,
         (max(psychAD_meta$DxPC1[which(psychAD_meta$Dx_cat=='early')],na.rm=T)+min(psychAD_meta$DxPC1[which(psychAD_meta$Dx_cat=='late')],na.rm=T))/2)
psychAD_meta$DxPC1_clust3=as.numeric(factor(psychAD_meta$Dx_cat))
write.csv(psychAD_meta,'psychAD_meta_PC1_pathonly.csv')
write.csv(cutoff,file='psychAD_meta_cutoff.csv')

pdf('cor_nopath.pdf')
zdf=psychAD_meta[,c(sel_Dx,'DxPC1')]
cmat=cor(zdf,use='complete.obs',method='spearman')
corrplot::corrplot.mixed(round(cmat,2),order='hclust',hclust.method='ward.D')
corrplot::corrplot(round(cmat,2),order='hclust',hclust.method='ward.D',
                   method='circle')
corrplot::corrplot(round(cmat,2),order='hclust',hclust.method='ward.D',
                   method='square')
corrplot::corrplot(round(cmat,2),order='hclust',hclust.method='ward.D',
                   method='number')
dev.off()
