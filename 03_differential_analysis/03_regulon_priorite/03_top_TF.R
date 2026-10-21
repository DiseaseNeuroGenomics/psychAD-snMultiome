.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(tidyr)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(future)
library(future.apply)
library(matrixStats)
library(dplyr)
library(dreamlet)
library(SingleCellExperiment)
library(stringr)
library(clusterProfiler)
library(rrvgo)
library(fgsea)
library(AnnotationDbi)
library(GOSemSim)
options(future.globals.maxSize= 60*1024^3)
n_cores <- 8
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

plot_split_marker_heatmap <- function(mat1,mat2) {
  library(purrr)
  library(dplyr)
  library(tibble)
  library(scales)
  study1_name = "A"; study2_name = "B"
  marker_order <- rownames(mat1)
  ord <- colnames(mat1)
  
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
  
  lim <- max(abs(c(mat1, mat2)), na.rm = TRUE)
  
  p <- ggplot(plot_df, aes(x = y, y = x)) +geom_polygon(
    aes(fill = value, group = id),color = "grey85",linewidth = 0.15) +
    #scale_fill_gradient2(low="blue",high="red",name = 'sign -log10(FDR)',
    #                     limit=c(-30,30))+
    scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),
                         limits=c(-30,30),name = 'sign -log10(FDR)')+
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
  return(p)
}
plot_split_marker_heatmap2 <- function(
    mat1,
    mat2,
    marker_list=NULL,
    sig1 = NULL,
    sig2 = NULL,
    ord = colnames(mat1),
    marker_order = rownames(mat1),
    study1_name = "study1",
    study2_name = "study2",
    
    sig_text_size = 3,
    sig_text_color = "black"
) {
  show_marker_box = !is.null(marker_list)
  stopifnot(
    all(dim(mat1) == dim(mat2))
  )
  
  nr <- nrow(mat1)
  nc <- ncol(mat1)
  
  if (is.null(sig1)) {
    sig1 <- matrix("", nrow = nr, ncol = nc)
  }
  if (is.null(sig2)) {
    sig2 <- matrix("", nrow = nr, ncol = nc)
  }
  
  stopifnot(
    all(dim(sig1) == dim(mat1)),
    all(dim(sig2) == dim(mat1))
  )
  
  if (is.null(marker_order)) {
    marker_order <- seq_len(nr)
  }
  if (is.null(ord)) {
    ord <- seq_len(nc)
  }
  
  cell_df <- expand.grid(
    i = seq_len(nr),
    j = seq_len(nc)
  )
  
  df1 <- purrr::pmap_dfr(cell_df, function(i, j) {
    tibble::tibble(
      x = c(j - 0.5, j + 0.5, j - 0.5),
      y = c(i - 0.5, i - 0.5, i + 0.5),
      id = paste0(i, "_", j, "_study1"),
      value = mat1[i, j],
      study = study1_name
    )
  })
  
  df2 <- purrr::pmap_dfr(cell_df, function(i, j) {
    tibble::tibble(
      x = c(j + 0.5, j + 0.5, j - 0.5),
      y = c(i - 0.5, i + 0.5, i + 0.5),
      id = paste0(i, "_", j, "_study2"),
      value = mat2[i, j],
      study = study2_name
    )
  })
  
  plot_df <- dplyr::bind_rows(df1, df2)
  
  sig_df1 <- purrr::pmap_dfr(cell_df, function(i, j) {
    tibble::tibble(
      x = i - 1/6,
      y = j - 1/6,
      label = sig1[i, j],
      study = study1_name
    )
  })
  
  sig_df2 <- purrr::pmap_dfr(cell_df, function(i, j) {
    tibble::tibble(
      x = i + 1/6,
      y = j + 1/6,
      label = sig2[i, j],
      study = study2_name
    )
  })
  
  sig_df <- dplyr::bind_rows(sig_df1, sig_df2) %>%
    dplyr::filter(!is.na(label), label != "")
  
  if(!is.null(marker_list)){
    marker_df <- purrr::imap_dfr(marker_list, function(marks, clust) {
      tibble::tibble(
        clust = clust,
        mark = marks,
        rank = seq_along(marks)
      )
    }) %>%
      dplyr::filter(
        clust %in% ord,
        mark %in% marker_order
      ) %>%
      dplyr::mutate(
        x = match(clust, ord),
        y = match(mark, marker_order)
      )
  }
  
  lim <- max(abs(c(mat1, mat2)), na.rm = TRUE)
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = y, y = x)) +
    ggplot2::geom_polygon(
      ggplot2::aes(fill = value, group = id),
      color = "grey85",
      linewidth = 0.15
    ) +
    #ggplot2::scale_fill_gradient2(
    #  low = "blue",
    #  mid = "white",
    #  high = "red",
    #  midpoint = 0,
    #  limits = c(-lim, lim),
    #  oob = scales::squish,
    #  na.value = "grey",
    #  name = "z-score"
    #)
    scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),
                         limits = c(-lim, lim),oob = scales::squish,na.value = "grey",
                         name = "z-score")+
    ggplot2::scale_y_continuous(
      breaks = seq_len(nc),
      labels = ord,
      expand = c(0, 0)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(nr),
      labels = marker_order,
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        #vjust = 1,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(colour = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.background = ggplot2::element_rect(fill = NA)
    )
  
  if (nrow(sig_df) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = sig_df,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = sig_text_size,
        color = sig_text_color
      )
  }
  
  if (show_marker_box) {
    p <- p +
      ggplot2::geom_tile(
        data = marker_df,
        ggplot2::aes(x = y, y = x),
        inherit.aes = FALSE,
        fill = NA,
        color = "black",
        linewidth = 0.35,
        width = 0.9,
        height = 0.9
      )
  }
  
  return(p)
}
genes=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/metadata/cellranger_gene.csv')
genes[duplicated(genes$gene_name),'gene_name']=paste0(genes[duplicated(genes$gene_name),'gene_name'],'-1')
genes$mt = grepl('^MT-',genes$gene_name)
genes$ribo = grepl('^RPL',genes$gene_name) |grepl('^RPS',genes$gene_name) # annotate the group of mitochondrial genes as 'mt'
genes$hb=grepl('^HB',genes$gene_name)
genes$protein_coding = genes$gene_type=='protein_coding'
genes$filtered=genes$protein_coding & (!genes$ribo)& (!genes$mt)
genes=genes[genes$filtered,]
rownames(genes)=genes$gene_name
cellcolor=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/color_code.csv',row.names=1)
f=cellcolor$color
names(f)=rownames(cellcolor)
cellcolor=f
lgb=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ML_cellstate3/lgb_anno.csv')
lgb=lgb[lgb$cell!='Mural',]

anno=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno_05092025_fixed3.csv')
anno=anno[anno$inRNA & anno$subclass!='Epend' & anno$subclass!='EN_L5_ET' ,]
a=unique(anno[,c('class','subtype','subclass')])
rownames(a)=a$subtype
b=unique(anno[,c('class','subclass')])
rownames(b)=b$subclass
celltype=list(subtype=a,subclass=b)
setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/DE3/')
DE1=read_big_csv(paste0('class/all/Dx_cell.csv'))
##########
DF=data.frame()
Adf=data.frame()
pdf('top_TF_select.pdf',width=4,height=4)
for(cell in unique(DE1$assay)){
  chromvar=read_big_csv(sprintf('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/depre_DE3/class_chromvar/%s_all/Dx_cell.csv',cell))
  rownames(chromvar)=chromvar$ID
  de1=DE1[DE1$assay==cell,]
  rownames(de1)=de1$ID
  sh=intersect(rownames(de1),rownames(chromvar))
  a=de1[sh,c('logFC','t','P.Value','adj.P.Val')]
  b=chromvar[sh,c('logFC','t','P.Value','adj.P.Val')]
  colnames(a)=paste0('G_',c('logFC','t','P','Q'))
  colnames(b)=paste0('C_',c('logFC','t','P','Q'))
  
  adf=data.frame(TF=sh,a,b)
  DF=rbind(DF,data.frame(adf,cell=cell))
  adf$ind=adf$G_t/abs(adf$G_t)*adf$C_t/abs(adf$C_t)*sqrt(pmin(abs(adf$C_t),25))*sqrt(pmin(abs(adf$G_t),25))
  Adf=rbind(Adf,data.frame(adf,cell=cell))
  
  a=as.data.frame(adf[ adf$C_Q<0.05 & adf$G_Q<0.05,]%>%group_by(G_logFC>0)%>%top_n(5,wt=ind))
  a=a[order(-1*a$G_t),]
  
  ord=a$ID=a$TF
  a$TF=factor(a$TF,levels=a$TF)
  print(plot_split_marker_heatmap(t(matrix(a$C_t/abs(a$C_t)*pmin((-log10(a$C_Q)),30),dimnames=list(ord,cell))),
                                  t(matrix(a$G_t/abs(a$G_t)*pmin((-log10(a$G_Q)),30),dimnames=list(ord,'G')))))
  
}
dev.off()
colnames(DF)=c('TF',paste0('GEX_',c('logFC','t','P','Padj')),
               paste0('motif_',c('logFC','t','P','Padj')),'cell type')
write.csv(DF,file=paste0('agg_TFexp_motf.csv'),row.names=F,quote = F)
sel=c()
adf=Adf[Adf$C_Q<0.05 & Adf$G_Q<0.05,]
adf=adf[order(-1*adf$ind),]
for(x in c("Astro",'Oligo','OPC','Micro_PVM','EN','IN')){
  adf1=adf[adf$cell==x,]
  sel=c(sel,c(head(adf1[adf1$C_t>0,]$TF,3),
              head(adf1[adf1$C_t<0,]$TF,3)))
  
}

a=Adf[Adf$TF%in%sel,]
a$TF=factor(a$TF,levels=unique(sel))
a$cell=factor(a$cell,levels=c("Astro",'Oligo','OPC','Micro_PVM','EN','IN'))
a$C_score=a$C_t/abs(a$C_t)*pmin((-log10(a$C_Q)),30)
a$G_score=a$G_t/abs(a$G_t)*pmin((-log10(a$G_Q)),30)
a$C_sig=c('','*')[(a$C_Q<0.05)+1]
a$G_sig=c('','*')[(a$G_Q<0.05)+1]
mat1=reshape2::acast(a,TF~cell,value.var = 'C_score')
mat2=reshape2::acast(a,TF~cell,value.var = 'G_score')
mat1b=reshape2::acast(a,TF~cell,value.var = 'C_sig')
mat2b=reshape2::acast(a,TF~cell,value.var = 'G_sig')
marker_list=lapply(split(a,a$cell),function(x)x[which(x$ind>0 & x$G_Q<0.05 & x$C_Q<0.05),]$TF)
z=t(cbind(mat1,mat2));z[is.na(z)]=0
z=pmax(pmin(z,10),-10)
ord=colnames(z)[order.dendrogram(as.dendrogram(hclust(dist(t(z)),method='ward.D2')))]
pdf('top_TF_agg.pdf')
print(plot_split_marker_heatmap2(mat1,mat2,marker_list=NULL,sig1 = mat1b,sig2 = mat2b)  )
print(plot_split_marker_heatmap2(mat1[ord,],mat2[ord,],marker_list=NULL,sig1 = mat1b[ord,],sig2 = mat2b[ord,])  )
dev.off()

map1=bitr(rownames(genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
map1=map1[!duplicated(map1$SYMBOL),]
rownames(map1)=map1$SYMBOL
ego=data.frame()
df=data.frame()
zde=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/class/all/Dx_donor.csv')
plan(multisession, workers = n_cores)
for(cell in unique(DE1$assay)){
  de1=DE1[DE1$assay==cell,]
  grn=read_big_csv(sprintf('/sc/arion/projects/CommonMind/roussp01a/snmulti/GRN/files/filter/%s_corr_p_filtered_grn.csv',cell))
  grns=split(grn$target,grn$TF)
  gl=lapply(grns,function(x)unique(map1[rownames(map1)%in%as.character(x),]$ENTREZID))
  de2=zde[zde$assay==cell,]
  bg=unique(map1[rownames(map1)%in%as.character(de2$ID),]$ENTREZID)
  
  ego1=do.call(rbind,
               future_lapply(names(grns),function(m){
                 ego <- enrichGO(gene          = gl[[m]],
                                 universe      = bg,
                                 OrgDb         = "org.Hs.eg.db",
                                 ont           = 'BP',
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 1,
                                 qvalueCutoff  = 1,
                                 readable      = TRUE)@result
                 ego$module=m;
                 ego$cell=cell;
                 ego=ego[ego$Count>=10,]
                 ego
               }))
  ego=rbind(ego1,ego)
  
  
  f=de1$t
  names(f)=de1$ID
  en1=as.data.frame(fgsea(grns,f,minSize=10))
  rownames(en1)=en1$pathway
  adf=Adf[ Adf$C_Q<0.05 & Adf$G_Q<0.05 & Adf$ind>0 &Adf$cell==cell,] 
  rownames(adf)=adf$TF
  zdf1=data.frame(adf,en1[rownames(adf),c('pval','padj','NES')])
  zdf1$direction='up'
  zdf1[zdf1$G_logFC<0,]$direction='down'
  df=rbind(df,zdf1)
} 
write.csv(df,file='Dx_cell_DEG_GRN_gsea.csv')
write.csv(ego,file='GRN_pathway_enrichment.csv')
zdf1=as.data.frame(df[which(df$G_logFC*df$NES>0 & df$padj<0.05),]%>%group_by(cell,direction)%>%top_n(3,wt=-padj))
zdf1$cell=factor(zdf1$cell,levels=c('Astro_Epend','Oligo','OPC','Micro_PVM','EN','IN'))
zdf1=zdf1[order(zdf1$cell,zdf1$NES),]
zdf1$y=paste0(zdf1$cell,':',zdf1$TF)
zdf1$y=factor(zdf1$y,levels=zdf1$y)
zdf1$sig=c("",'.','+')[1+(zdf1$pval<0.05)+(zdf1$padj<0.05)]
#
pdf('top_TF_sel_agg.pdf')
print(ggplot(zdf1)+
        geom_point(aes(x='a',y=y,color=NES,size=-log10(padj)))+
        scale_size(range=c(0.5,3))+
        scale_color_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-max(abs(zdf1$NES)),max(abs(zdf1$NES))))+
        theme_bw()+coord_equal()+facet_wrap(~cell,nrow=1)+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              legend.position='bottom'
              #panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank()
        ))
zdf1$padj[zdf1$padj<1e-30]=1e-30
print(ggplot(zdf1)+
        geom_point(aes(x='a',y=y,color=NES,size=-log10(padj)))+
        scale_size(range=c(0.5,3))+
        scale_color_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-max(abs(zdf1$NES)),max(abs(zdf1$NES))))+
        theme_bw()+coord_equal()+facet_wrap(~cell,nrow=1)+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              legend.position='bottom'
              #panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank()
        ))
dev.off()
pdf('top_TF_exp_motif.pdf',width=4,height = 4)
for(cell in unique(zdf1$cell)){
  a=zdf1[zdf1$cell==cell,]
  rownames(a)=a$TF
  ord=a$TF
  a$TF=factor(a$TF,levels=ord)
  print(plot_split_marker_heatmap(t(matrix(a$C_t/abs(a$C_t)*pmin((-log10(a$C_Q)),30),dimnames=list(ord,cell))),
                                  t(matrix(a$G_t/abs(a$G_t)*pmin((-log10(a$G_Q)),30),dimnames=list(ord,'G')))))
}
dev.off()

