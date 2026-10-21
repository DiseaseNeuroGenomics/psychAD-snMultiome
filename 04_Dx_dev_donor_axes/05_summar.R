.libPaths(c('/sc/arion/projects/roussp01a/pengfei/tools/R_4_4_1',.libPaths()))
library(ggplot2)
library(tidyverse)
library(dplyr)
library(data.table)
library(qvalue)
library(ggrepel)
library(MOFA2)
library(variancePartition)
library(crumblr)

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
clr <- function(counts, pseudocount = 0.5) {
  if (!is.data.frame(counts) & !is.matrix(counts)) {
    counts <- matrix(counts, nrow = 1)
  }
  
  log(counts + pseudocount) - rowMeans(log(counts + pseudocount))
}

plot_fgsea_go <- function(
    ego,
    sim_mat,
    go_map,
    padj_cutoff = 0.05,
    n_terms_per_factor = 5,
    sim_cut = 0.5,
    cell_order = NULL,
    factor_order = NULL,
    nes_col = "NES"
) {
  req_cols <- c("factor", "pathway", "pval", "padj", "cell")
  miss_cols <- setdiff(req_cols, colnames(ego))
  
  if (is.vector(go_map) && !is.null(names(go_map))) {
    go_df <- tibble(
      pathway = names(go_map),
      go_name = unname(go_map)
    )
  } else if (is.data.frame(go_map)) {
    if (!all(c("pathway", "go_name") %in% colnames(go_map))) {
      stop("If go_map is a data.frame, it must contain columns: pathway, go_name")
    }
    go_df <- go_map %>% distinct(pathway, go_name)
  } else {
    stop("go_map must be a named vector or a data.frame with columns pathway and go_name")
  }
  
  ego2 <- ego %>%
    mutate(
      pathway = as.character(pathway),
      factor = as.character(factor),
      cell = as.character(cell)
    ) %>%
    left_join(go_df, by = "pathway") %>%
    mutate(go_name = ifelse(is.na(go_name), pathway, go_name))
  
  if (is.null(cell_order)) {
    cell_order <- unique(ego2$cell)
  }
  if (is.null(factor_order)) {
    factor_order <- unique(ego2$factor)
  }
  
  ego2 <- ego2 %>%
    mutate(
      factor = factor(factor, levels = factor_order),
      cell = factor(cell, levels = cell_order)
    )
  
  sig <- ego2 %>%
    filter(padj < padj_cutoff)
  
  if (nrow(sig) == 0) {
    stop("No significant pathways found at padj < ", padj_cutoff)
  }
  
  ## -----------------------------
  ## 3. Assign semantic clusters within each factor
  ## -----------------------------
  cluster_one_factor <- function(df_factor, sim_mat, sim_cut = 0.5) {
    ids_all <- unique(df_factor$pathway)
    ids_in_mat <- intersect(ids_all, rownames(sim_mat))
    
    cluster_df <- tibble(pathway = ids_all)
    
    ## pathways not found in sim_mat -> each gets its own cluster
    cluster_df$cluster <- paste0("singleton_", seq_len(nrow(cluster_df)))
    
    if (length(ids_in_mat) == 1) {
      cluster_df$cluster[match(ids_in_mat, cluster_df$pathway)] <- "C1"
      return(cluster_df)
    }
    
    if (length(ids_in_mat) > 1) {
      S <- sim_mat[ids_in_mat, ids_in_mat, drop = FALSE]
      S[is.na(S)] <- 0
      S[S < 0] <- 0
      S[S > 1] <- 1
      diag(S) <- 1
      
      hc <- hclust(as.dist(1 - S), method = "average")
      cl <- cutree(hc, h = 1 - sim_cut)
      
      cluster_df$cluster[match(ids_in_mat, cluster_df$pathway)] <-
        paste0("C", cl)
    }
    
    cluster_df
  }
  
  cluster_tbl <- sig %>%
    group_split(factor) %>%
    map_dfr(function(df) {
      fac <- unique(as.character(df$factor))
      cl <- cluster_one_factor(df, sim_mat = sim_mat, sim_cut = sim_cut)
      cl$factor <- fac
      cl
    })
  
  ## -----------------------------
  ## 4. Summarize pathways across cell types
  ## -----------------------------
  term_summary <- sig %>%
    mutate(logFDR = -log10(padj)) %>%
    left_join(cluster_tbl, by = c("factor", "pathway")) %>%
    group_by(factor, pathway, go_name, cluster) %>%
    summarise(
      n_cells_sig = n_distinct(cell),
      best_padj = min(padj, na.rm = TRUE),
      best_logFDR = max(logFDR, na.rm = TRUE),
      mean_logFDR = mean(logFDR, na.rm = TRUE),
      best_absNES = if (nes_col %in% colnames(.)) {
        max(abs(.data[[nes_col]]), na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )
  
  ## -----------------------------
  ## 5. Choose 1 representative pathway per cluster
  ## -----------------------------
  cluster_rep <- term_summary %>%
    group_by(factor, cluster) %>%
    arrange(
      desc(n_cells_sig),
      desc(best_logFDR),
      desc(best_absNES),
      best_padj,
      .by_group = TRUE
    ) %>%
    slice(1) %>%
    ungroup()
  
  ## -----------------------------
  ## 6. Choose top clusters/pathways per factor
  ## -----------------------------
  selected_terms <- cluster_rep %>%
    group_by(factor) %>%
    arrange(
      desc(n_cells_sig),
      desc(best_logFDR),
      desc(best_absNES),
      best_padj,
      .by_group = TRUE
    ) %>%
    slice_head(n = n_terms_per_factor) %>%
    mutate(term_rank = row_number()) %>%
    ungroup()
  
  ## -----------------------------
  ## 7. Build plotting table
  ## -----------------------------
  plot_df <- selected_terms %>%
    dplyr::select(factor, pathway, go_name, term_rank) %>%
    tidyr::crossing(cell = factor(cell_order, levels = cell_order)) %>%
    left_join(
      ego2 %>%
        dplyr::select(factor, pathway, cell, pval, padj, any_of(nes_col)),
      by = c("factor", "pathway", "cell")
    ) %>%
    mutate(
      padj = ifelse(is.na(padj), 1, padj),
      pval = ifelse(is.na(pval), 1, pval),
      logFDR = -log10(padj),
      logFDR = pmin(logFDR, 20),
      sig = padj < padj_cutoff
    )
  
  ## Label with wrapped GO name
  plot_df <- plot_df %>%
    mutate(
      pathway_label = str_wrap(go_name, width = 80),
      pathway_label2 = tidytext::reorder_within(
        pathway_label,
        -term_rank,
        factor
      )
    )
  
  ## -----------------------------
  ## 8. Plot
  ## -----------------------------
  col=c('white','black');names(col)=c(F,T)
  print(ggplot(
    plot_df%>% filter(sig),
    aes(x = cell, y = pathway_label2)) +
      geom_point(
        aes(size = logFDR, colour  = .data[[nes_col]]),
        alpha = 0.9) +
      facet_wrap(~ factor, scales = "free_y", ncol = 1) +
      tidytext::scale_y_reordered() +
      scale_size_continuous(
        name = expression(-log[10](FDR)),
        range = c(2, 8)
      ) +#coord_equal()+
      scale_color_gradient2(
        low = "#3B4CC0",
        mid = "white",
        high = "#B40426",
        midpoint = 0,
        name = "NES") +
      labs(x = NULL,y = NULL) +
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            strip.background = element_rect(fill = "grey95"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank()))
  p <- ggplot(
    plot_df, #%>% filter(sig),
    aes(x = cell, y = pathway_label2)) +
    geom_point(
      aes(size = logFDR, fill = .data[[nes_col]],color=sig),
      alpha = 0.9,shape=21,stroke=1) +
    facet_wrap(~ factor, scales = "free_y", ncol = 1) +
    tidytext::scale_y_reordered() +
    scale_size_continuous(
      name = expression(-log[10](FDR)),
      range = c(2, 8)
    ) +#coord_equal()+
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      name = "NES") +
    scale_color_manual(values =col)+
    labs(x = NULL,y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_rect(fill = "grey95"),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank())
  
  
  list(
    plot = p,
    selected_terms = selected_terms,
    cluster_representatives = cluster_rep,
    plot_data = plot_df
  )
}

genes=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step1/metadata/cellranger_gene.csv')
genes[duplicated(genes$gene_name),'gene_name']=paste0(genes[duplicated(genes$gene_name),'gene_name'],'-1')
genes$mt = grepl('^MT-',genes$gene_name)
genes$ribo = grepl('^RPL',genes$gene_name) |grepl('^RPS',genes$gene_name) # annotate the group of mitochondrial genes as 'mt'
genes$hb=grepl('^HB',genes$gene_name)
genes$protein_coding = genes$gene_type=='protein_coding'
genes$filtered=genes$protein_coding & (!genes$ribo)& (!genes$mt) & genes$seqnames%in%paste0('chr',c(1:22,'X'))
#genes=genes[genes$filtered,]
rownames(genes)=genes$gene_name
#####
lgb=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/ML_cellstate3/lgb_anno.csv')
lgb$celltype=lgb$subtype
lgb[lgb$class%in%c('EN','IN'),]$celltype=lgb[lgb$class%in%c('EN','IN'),]$subclass


metaInfo=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/metadata/metaInfo.csv',row.names=1)
#######cell composition #########################################################################################
anno=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/cell_anno.csv')
anno=anno[anno$inRNA,]
anno=anno[anno$subtype!='ependA',]

#####neuron loss#####

Neu=anno[anno$class%in%c('EN'),]
ncells=table(Neu$order,Neu$subclass)
ncells=apply(ncells,2,c)
df=data.frame(vul=ncells[,c('EN_L2_3_IT')],
              nvul=rowSums(ncells)-ncells[,'EN_L2_3_IT'])
cobj=crumblr(df)
vul1=cobj
topTable(dream(cobj,~scale(DxPC1)+Sex+scale(Age)+BrainRegion2+(1|SubID)+scale(PMI..min.)+scale(GWR) ,metaInfo[colnames(cobj),]),'scale(DxPC1)',Inf)
topTable(dream(cobj,~scale(DxPC1)+Sex+scale(Age)+BrainRegion2+(1|SubID)+scale(PMI..min.) ,metaInfo[colnames(cobj),]),'scale(DxPC1)',Inf)
######
Neu=lgb[lgb$class%in%c('IN'),]
#ncells=table(Neu$order,Neu$subtype%in%c('IN_SST_L2_3','IN_SST_L3'))
ncells=table(Neu$order,Neu$subtype%in%c('IN_SST_L2_3'))
colnames(ncells)=c('other','IN_SST')
ncells=apply(ncells,2,c)
df=data.frame(vul=ncells[,c('IN_SST')],
              nvul=rowSums(ncells)-ncells[,'IN_SST'])
cobj=crumblr(df)
vul2=cobj
topTable(dream(cobj,~scale(DxPC1)+Sex+scale(Age)+BrainRegion2+(1|SubID)+scale(PMI..min.)+scale(GWR) ,metaInfo[colnames(cobj),]),'scale(DxPC1)',Inf)
topTable(dream(cobj,~scale(DxPC1)+Sex+scale(Age)+BrainRegion2+(1|SubID)+scale(PMI..min.) ,metaInfo[colnames(cobj),]),'scale(DxPC1)',Inf)
#####
Neu=lgb[lgb$class%in%c('EN','IN'),]
ncells=table(Neu$order,Neu$subclass)
ncells=apply(ncells,2,c)
df=data.frame(vul=rowSums(ncells[,c('EN_L2_3_IT','IN_SST')]),
              nvul=rowSums(ncells)-rowSums(ncells[,c('EN_L2_3_IT','IN_SST')]))
cobj=crumblr(df)
vul=cobj
topTable(dream(cobj,~scale(DxPC1)+Sex+scale(Age)+BrainRegion2+(1|SubID)+scale(PMI..min.)+scale(GWR) ,metaInfo[colnames(cobj),]),'scale(DxPC1)',Inf)
topTable(dream(cobj,~scale(DxPC1)+Sex+scale(Age)+BrainRegion2+(1|SubID)+scale(PMI..min.) ,metaInfo[colnames(cobj),]),'scale(DxPC1)',Inf)
sh=intersect(colnames(vul1),colnames(vul2))
vul_df=zdf=data.frame(EN=(vul1$E['vul',sh]),IN=(vul2$E['vul',sh]),
                      Neu=(vul$E['vul',sh]))
cor(zdf,use='complete.obs',method='spearman')
lapply(split(zdf,metaInfo[sh,'BrainRegion2']),cor,method='spearman')
sh=intersect(colnames(vul1),colnames(vul2))
Vul=new("EList")
Vul$E <- rbind(vul$E['vul',sh],
               vul1$E['vul',sh],
               vul2$E['vul',sh])
Vul$weights <- rbind(vul$weights['vul',sh],
                     vul1$weights['vul',sh],
                     vul2$weights['vul',sh])
rownames(Vul$E)=rownames(Vul$weights)=c('Neu','EN','IN')
####################subtypes#############
anno$celltype=anno$subtype
anno[anno$class%in%c('EN','IN'),]$celltype=anno[anno$class%in%c('EN','IN'),]$subclass
ncells=lapply(split(anno[!anno$class%in%c('Adaptive','Endo','Mural'),],anno[!anno$class%in%c('Adaptive','Endo','Mural'),]$class),function(x)table(x$order,x$celltype))
ncells=lapply(ncells,function(x)apply(x,2,c))
z1=apply(table(anno$order,anno$class),2,c)
z2=apply(table(anno$order,anno$subclass),2,c)
ncells$class=z1
ncells$subclass=z2

z=anno[anno$class%in%c('IN','EN'),]
z=apply(table(z$order,z$subclass),2,c)
ncells$neu=z
cobjs=lapply(ncells,crumblr)
###################Gex & MOFA factors####
setwd('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/donor/new')
Sel=read.csv('sel.csv')
Sels=lapply(split(Sel,Sel$threshold),function(x){rownames(x)=paste0(x$cell,':',x$gene);x})
chr_block=read.csv('chr_block.csv')
Blubs=lapply(unique(Sel$cell),function(x)read_big_csv(sprintf('%s_RandEf.csv.gz',x)))
GEXs=lapply(unique(Sel$cell),function(x)read_big_csv(sprintf('%s_resid.csv.gz',x)))
names(Blubs)=names(GEXs)=unique(Sel$cell)
sh=unique(unlist(lapply(GEXs,colnames)))
GEXs=lapply(GEXs,function(x){
  a=array(NA,dim=c(nrow(x),length(sh)),dimnames = list(rownames(x),sh))
  a[,colnames(x)]=as.matrix(x)
  a
})
#residual matrix
z=mapply(function(x,y){rownames(x)=paste0(y,":",rownames(x));x},GEXs,names(GEXs),SIMPLIFY = F)
names(z)=NULL
GEX=do.call(rbind,z)
Obs=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/RNA_all_soupX_obs.csv.gz')
CovInfo=as.data.frame(data.frame(Obs[rownames(lgb),],lgb[,c('class','order')])%>%group_by(class,order)%>%summarise(
  n_genes_by_counts=mean(n_genes_by_counts),
  pct_counts_ribo=mean(pct_counts_ribo),
  mito_genes=mean(mito_genes)))
CovInfo=do.call(rbind,lapply(split(CovInfo,CovInfo$class),function(x)data.frame(x[,1:2],scale(x[,3:5]))))
CovInfos=lapply(list('EN','IN',c("EN",'IN'),c("Astro",'Oligo','OPC','Micro_PVM'),c("Astro",'Oligo','OPC','Micro_PVM','EN','IN')),function(x){
  z=CovInfo[CovInfo$class%in%x,]
  z=data.frame(z%>%pivot_wider(names_from=class,values_from=c(n_genes_by_counts,pct_counts_ribo,mito_genes)))
  rownames(z)=z[,1]
  z[,1]=NULL
  z[rowSums(is.na(z))==0,]
  
})
names(CovInfos)=c('EN','IN','Neu','Gli','All')
for(x in c('Neu','Gli','All')){
  z=CovInfos[[x]]
  sv1=svd(z[,-1])
  npc=which(cumsum(sv1$d^2/sum(sv1$d^2))>0.85)[1]
  z=data.frame(sv1$u[,1:npc,drop=F],row.names=rownames(z))
  colnames(z)=paste0(x,'_PC',1:npc)
  CovInfos[[x]]=z
}

###psychAD
load('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/MSSM/class/DxPC1.Rdata')
NcovInfo=do.call(rbind,lapply(c('Astro','Oligo','OPC','Immune','EN','IN'),function(cell){z=fit[[cell]]$data;
data.frame(sample=rownames(z),assay=cell,scale(z[,c('n_genes','mito_genes','ribo_genes')]))}))
NcovInfos=lapply(list('EN','IN',c("EN",'IN'),c("Astro",'Oligo','OPC','Immune'),c("Astro",'Oligo','OPC','Immune',"EN",'IN')),function(x){
  z=NcovInfo[NcovInfo$assay%in%x,]
  z=data.frame(z%>%pivot_wider(names_from=assay,values_from=c(n_genes,ribo_genes,mito_genes)))
  rownames(z)=z[,1]
  z[,1]=NULL
  z[rowSums(is.na(z))==0,]
})
names(NcovInfos)=c('EN','IN','Neu','Gli','All')
for(x in c('Neu','Gli','All')){
  z=NcovInfos[[x]]
  sv1=svd(z[,-1])
  npc=which(cumsum(sv1$d^2/sum(sv1$d^2))>0.85)[1]
  z=data.frame(sv1$u[,1:npc,drop=F],row.names=rownames(z))
  colnames(z)=paste0(x,'_PC',1:npc)
  NcovInfos[[x]]=z
}

NGexs=lapply(c('Astro','Oligo','OPC','Immune','EN','IN'),function(cell){
  R=residuals(fit[[cell]])
  rownames(R)=paste0(gsub('Immune','Micro_PVM',cell),":",rownames(R))
  R
})
f=table(unlist(lapply(NGexs,colnames)))
f=names(f)[f==length(NGexs)]
NGex=do.call(rbind,lapply(NGexs,function(x)x[,f]))
##

###################################################################
#########var explained & reproducibility ###################
nfactor=25
load(sprintf('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/files/donor/new/MOFA_top_%s_10_agg.Rdata',nfactor))
####
W=agg$consensus_weights
Fs=agg$consensus_scores
F_sel=paste0('Factor',1:3)



#0. how much variance Expalined.
pdf(sprintf('MOFA_summary_%s.pdf',nfactor),width=5,height = 5)
corrplot::corrplot(cor(Fs,use='complete.obs'),
                   col = rev(corrplot::COL2())[21:180],method  = 'square',diag = T)
zdf=factor_var
zdf$lab=signif(zdf$Mean,2)
zdf[zdf$lab<2,]$lab=''
zdf$factor=factor(zdf$factor,levels=rownames(agg$factor_reproducibility))
zdf$assay=factor(zdf$assay,levels=c('Astro','Oligo','OPC','Micro_PVM','EN','IN'))
print(ggplot(zdf)+
        geom_tile(aes(y=factor,x=assay,fill=Mean),color='grey50')+
        geom_text(aes(y=factor,x=assay,label=lab))+
        scale_fill_gradient(low='white',high=RColorBrewer::brewer.pal(9,'Greens')[8])+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))
print(ggplot(zdf)+
        geom_tile(aes(y=factor,x=assay,fill=Mean),color='grey50')+
        geom_text(aes(y=factor,x=assay,label=lab))+
        scale_fill_gradient(low='white',high=RColorBrewer::brewer.pal(9,'Greens')[8],
                            limits=c(0,max(zdf$Mean)))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))
print(ggplot(zdf[zdf$factor%in%F_sel,])+
        geom_tile(aes(y=factor,x=assay,fill=Mean),color='grey50')+
        geom_text(aes(y=factor,x=assay,label=lab))+
        scale_fill_gradient(low='white',high=RColorBrewer::brewer.pal(9,'Greens')[8])+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))

#1 concordant
#a. SubID variable genes are reproducible across cohorts. 
#b. 5-fold chr-block
Mtopother=load_model(sprintf('MOFA/MOFA_top%s_10.hdf5',setdiff(c(25,50),nfactor)))   
Fother=get_factors(Mtopother)[[1]]
z=sapply(split(agg$donor_scores_long[,3:7],agg$donor_scores_long$factor),function(x){
  a=abs(cor(x,method='spearman'));mean(a[upper.tri(a)])
})
zdf=data.frame(Factor=factor(c(names(z),colnames(Fs)),levels=colnames(Fs)),
               set=rep(c('chr-block (5-mean)',sprintf('top %s%% genes',setdiff(c(25,50),nfactor))),each=ncol(Fs)),
               value=c(z,apply(abs(cor(Fs,Fother,method='spearman')),1,max)))
zdf$lab=round(zdf$value,2)
#zdf[zdf$lab<2,]$lab=''
print(ggplot(zdf)+geom_bar(aes(x=value,y=Factor,fill=set),stat='identity',position='dodge')+
        theme_bw()+ggtitle('reproducibility')+xlab('rho')+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))

print(ggplot(zdf[zdf$Factor%in%F_sel,])+
        geom_tile(aes(y=Factor,x=set,fill=value),color='grey50')+
        geom_text(aes(y=Factor,x=set,label=lab),color='white')+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
print(ggplot(zdf)+
        geom_tile(aes(y=Factor,x=set,fill=value),color='grey50')+
        geom_text(aes(y=Factor,x=set,label=lab),color='white')+
        #scale_fill_gradientn(colors=c('white', "#FCFDBFFF","#FE9F6DFF","#DE4968FF","#8C2981FF"),limits=c(0,1))+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
print(ggplot(zdf)+
        geom_tile(aes(y=Factor,x=set,fill=value),color='grey50')+
        geom_text(aes(y=Factor,x=set,label=lab),color='white')+
        #scale_fill_gradientn(colors=c('white', "#FCFDBFFF","#FE9F6DFF","#DE4968FF","#8C2981FF"),limits=c(0,1))+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')[1:5]),limits=c(0,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
print(ggplot(zdf)+
        geom_tile(aes(y=Factor,x=set,fill=value),color='grey50')+
        geom_text(aes(y=Factor,x=set,label=lab),color='white')+
        #scale_fill_gradientn(colors=c('white', "#FCFDBFFF","#FE9F6DFF","#DE4968FF","#8C2981FF"),limits=c(0,1))+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')[1:5]),limits=c(0.7,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
#c. projection across brain region
Fp=scale(t(GEX[rownames(W),])%*%W)
Fp=Fp[rownames(Fp)%in%rownames(metaInfo[metaInfo$SubID%in%rownames(Fs),]),]
proj_corr=diag(cor(Fp,Fs[metaInfo[rownames(Fp),]$SubID,],use='complete.obs',method='spearman'))
Fp=scale(t(GEX[rownames(W),])%*%W)
sh=intersect(rownames(NGex),rownames(W))
NFp=scale(t(NGex[sh,])%*%W[sh,])
sh=intersect(rownames(NFp),rownames(Fs))
#d. projections in NPSAD
Nproj_corr=diag(cor(NFp[sh,],Fs[sh,],use='complete.obs',method='spearman'))

a=data.frame(Fp,metaInfo[rownames(Fp),c('SubID','BrainRegion2')])
a=a[!duplicated(a[,c('SubID','BrainRegion2')]),]
a=do.call(rbind,lapply(colnames(Fp),function(f){
  z=cor(reshape2::acast(a,SubID~BrainRegion2,value.var=f),use='complete.obs',method='spearman')
  z[upper.tri(z)]
}))
rownames(a)=colnames(Fp)
colnames(a)=c('PFC-PHG','STG-PHG','PFC-STG')
crossbr_corr=rowMeans(a)
zdf=data.frame(Factor=factor(c(names(proj_corr),names(Nproj_corr),names(crossbr_corr)),levels=colnames(Fs)),
               rho=c(proj_corr,Nproj_corr,crossbr_corr),
               group=factor(rep(c('project-multiome','project-NPSAD','cross brain region'),each=ncol(Fs)),
                            levels=c('project-multiome','project-NPSAD','cross brain region')))
zdf$lab=round(zdf$rho,2)
print(ggplot(zdf)+
        geom_tile(aes(y=Factor,x=group,fill=rho),color='grey50')+
        geom_text(aes(y=Factor,x=group,label=lab),color='white')+
        #scale_fill_gradientn(colors=c('white', "#FCFDBFFF","#FE9F6DFF","#DE4968FF","#8C2981FF"),limits=c(0,1))+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
print(ggplot(zdf[zdf$Factor%in%F_sel,])+
        geom_tile(aes(y=Factor,x=group,fill=rho),color='grey50')+
        geom_text(aes(y=Factor,x=group,label=lab),color='white')+
        #scale_fill_gradientn(colors=c('white', "#B40426"),limits=c(0,1))+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
#######4. Glia only projection####
z=do.call(rbind,strsplit(rownames(W),'\\:'))
Wgli=W[!z[,1]%in%c("EN",'IN'),]
Fgli=scale(t(GEX[rownames(Wgli),])%*%Wgli)
Fgli=Fgli[rownames(Fgli)%in%rownames(metaInfo[metaInfo$SubID%in%rownames(Fs),]),]
Gproj_corr=diag(cor(Fgli,Fs[metaInfo[rownames(Fgli),]$SubID,],use='complete.obs',method='spearman'))
Fgli=scale(t(GEX[rownames(Wgli),])%*%Wgli)
sh=intersect(rownames(NGex),rownames(Wgli))
NFgli=scale(t(NGex[sh,])%*%W[sh,])
sh=intersect(rownames(NFgli),rownames(Fs))
NFgli_corr=diag(cor(NFgli[sh,],Fs[sh,],use='complete.obs',method='spearman'))

a=data.frame(Fgli,metaInfo[rownames(Fgli),c('SubID','BrainRegion2')])
a=a[!duplicated(a[,c('SubID','BrainRegion2')]),]
a=do.call(rbind,lapply(colnames(Fgli),function(f){
  z=cor(reshape2::acast(a,SubID~BrainRegion2,value.var=f),use='complete.obs',method='spearman')
  z[upper.tri(z)]
}))
rownames(a)=colnames(Fgli)
colnames(a)=c('PFC-PHG','STG-PHG','PFC-STG')
Gcrossbr_corr=rowMeans(a)
zdf=data.frame(Factor=factor(c(names(Gproj_corr),names(NFgli_corr),names(Gcrossbr_corr)),levels=colnames(Fs)),
               rho=c(Gproj_corr,NFgli_corr,Gcrossbr_corr),
               group=factor(rep(c('project-multiome','project-NPSAD','cross brain region'),each=ncol(Fs)),
                            levels=c('project-multiome','project-NPSAD','cross brain region')))
zdf$lab=round(zdf$rho,2)
print(ggplot(zdf)+
        geom_tile(aes(y=Factor,x=group,fill=rho),color='grey50')+
        geom_text(aes(y=Factor,x=group,label=lab),color='white')+
        #scale_fill_gradientn(colors=c('white', "#FCFDBFFF","#FE9F6DFF","#DE4968FF","#8C2981FF"),limits=c(0,1))+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+ggtitle('glia only')+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
print(ggplot(zdf[zdf$Factor%in%F_sel,])+
        geom_tile(aes(y=Factor,x=group,fill=rho),color='grey50')+
        geom_text(aes(y=Factor,x=group,label=lab),color='white')+
        #        scale_fill_gradientn(colors=c('white', "#B40426"),limits=c(0,1))+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+ggtitle('glia only')+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
dev.off()
######## Dx_dev/MOFA cell composition correlation ##################################
Dx_dev=read_big_csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/cellstate/files/ML_cellstate3/Dx_dev_covadj.csv')
pfc=rownames(Fp)[metaInfo[rownames(Fp),]$BrainRegion2=='PFC']
pfc=Fp[pfc,]
rownames(pfc)=metaInfo[rownames(pfc),]$SubID
br='PFC'
zInfo=metaInfo[metaInfo$BrainRegion2==br,]
zInfo$Dx_dev=rowMeans(scale(data.frame(Dx_dev)[rownames(zInfo),]),na.rm=T)
zInfo=zInfo[!duplicated(zInfo$SubID),]
rownames(zInfo)=zInfo$SubID
zInfo=data.frame(zInfo,data.frame(pfc)[rownames(zInfo),])
xx0=zInfo[,c('Dx_dev','BrainRegion2',colnames(Fs))]
zdf0=do.call(rbind,lapply(c(paste0('Factor',1:3)),
                          function(x)summary(lm(as.formula(sprintf('Dx_dev~DxPC1+Age+Sex+PMI..min.+scale(%s)',x)),zInfo))$coefficient[sprintf('scale(%s)',x),]))

br='STG'
zInfo=metaInfo[metaInfo$BrainRegion2==br,]
zInfo$Dx_dev=rowMeans(scale(data.frame(Dx_dev)[rownames(zInfo),]),na.rm=T)
zInfo=zInfo[!duplicated(zInfo$SubID),]
rownames(zInfo)=zInfo$SubID
zInfo=data.frame(zInfo,data.frame(pfc)[rownames(zInfo),])
xx1=zInfo[,c('Dx_dev','BrainRegion2',colnames(Fs))]
zdf1=do.call(rbind,lapply(c(paste0('Factor',1:3)),
                         function(x)summary(lm(as.formula(sprintf('Dx_dev~DxPC1+Age+Sex+PMI..min.+scale(%s)',x)),zInfo))$coefficient[sprintf('scale(%s)',x),]))
br='PHG'
zInfo=metaInfo[metaInfo$BrainRegion2==br,]
zInfo$Dx_dev=rowMeans(scale(data.frame(Dx_dev)[rownames(zInfo),]),na.rm=T)
zInfo=zInfo[!duplicated(zInfo$SubID),]
rownames(zInfo)=zInfo$SubID
zInfo=data.frame(zInfo,data.frame(pfc)[rownames(zInfo),])
xx2=zInfo[,c('Dx_dev','BrainRegion2',colnames(Fs))]

zdf2=do.call(rbind,lapply(c(paste0('Factor',1:3)),
                          function(x)summary(lm(as.formula(sprintf('Dx_dev~DxPC1+Age+Sex+PMI..min.+scale(%s)',x)),zInfo))$coefficient[sprintf('scale(%s)',x),]))
zdf0=data.frame(zdf0)
zdf1=data.frame(zdf1)
zdf2=data.frame(zdf2)
zdf0$group='PFC'
zdf1$group='STG'
zdf2$group='PHG'
zdf=rbind(zdf0,zdf1,zdf2)
colnames(zdf)=c('FC','SE','t','P','BrainRegion')
zdf$coef=paste0('Factor',1:3)
zdf=data.frame(zdf)
zdf$q=p.adjust(zdf$P,method='BH')
write.csv(zdf,file='brainregion_factor_Dx_dev_asso.csv')
write.csv(rbind(xx0,xx1,xx2),file='brainregion_F123_Dx_dev.csv')

cmat=cor(Dx_dev,use='complete.obs',method='spearman')
pdf(sprintf('Dx_dev_cor_cell_%s.pdf',nfactor),width=4,height = 4)
corrplot::corrplot(cmat[c('Astro','Oligo','OPC','Micro_PVM','EN','IN'),c('Astro','Oligo','OPC','Micro_PVM','EN','IN')],
                   col = rev(corrplot::COL2()),type = 'upper',method  = 'square',diag = FALSE)
corrplot::corrplot.mixed(cmat[c('Astro','Oligo','OPC','Micro_PVM','EN','IN'),c('Astro','Oligo','OPC','Micro_PVM','EN','IN')],
                         upper.col = rev(corrplot::COL2()),lower.col=rev(corrplot::COL2()),upper = 'square')
dev.off()
#barplot(cor(rowMeans(Dx_dev),data.frame(Fs)[metaInfo[rownames(Dx_dev),]$SubID,],use='complete.obs'),las=2)
Dx_cell_agg=as.data.frame(lgb%>%group_by(order,class)%>%summarise(M=mean(Dx_cell))%>%pivot_wider(names_from=class,values_from=M))
rownames(Dx_cell_agg)=Dx_cell_agg[,1]
covs=c('n_genes_by_counts','pct_counts_ribo','mito_genes')
Vp=data.frame()
Dif=data.frame()
for(cell in unique(lgb$class)){
  lgb1=lgb[lgb$class==cell,]
  agg1=as.data.frame(lgb1%>%group_by(order)%>%summarise(M=mean(Dx_cell),Dx_donor=mean(DxPC1)))
  rownames(agg1)=agg1$order;agg1$order=NULL
  
  ncell=apply(table(lgb1$order,lgb1$Dx_cat_cell),2,c)
  cobj=crumblr(ncell)
  cobj$E=rbind(cobj$E,agg=scale(agg1[colnames(cobj),]$M)[,1])
  cobj$weights=rbind(cobj$weights,agg=1)
  covInfo=do.call(rbind,lapply(split(Obs[rownames(lgb1),covs],Obs[rownames(lgb1),'order']),colMeans))
  
  form=~scale(Age)+Sex+scale(DxPC1)+scale(PMI..min.)+BrainRegion+(1|SubID)+scale(n_genes_by_counts)+scale(pct_counts_ribo)+scale(mito_genes)+scale(Factor1)+scale(Factor2)+scale(Factor3)
  vform= ~scale(Age)+(1|Sex)+(1|SubID)+scale(DxPC1)+scale(PMI..min.)+(1|BrainRegion)+scale(n_genes_by_counts)+scale(pct_counts_ribo)+scale(mito_genes)+scale(Factor1)+scale(Factor2)+scale(Factor3)
  zInfo=data.frame(covInfo[colnames(cobj),],metaInfo[colnames(cobj),])
  zInfo=data.frame(data.frame(Fs)[as.character(zInfo$SubID),],zInfo)
  vp=fitExtractVarPartModel(cobj,vform,
                            zInfo)
  vp=data.frame(vp,ID=rownames(vp),assay=cell)
  Vp=rbind(Vp,vp)
  zfit=eBayes(dream(cobj,form,zInfo))
  z1=topTable(zfit,'scale(Factor1)')
  topTable(diffVar(zfit),'scale(Factor1)')
  z2=topTable(zfit,'scale(Factor2)')
  z3=topTable(zfit,'scale(Factor3)')
  z4=topTable(zfit,'scale(DxPC1)')
  Dif=rbind(Dif,
            data.frame(z1,factor='Factor1',cell=cell,coef=rownames(z1)),
            data.frame(z2,factor='Factor2',cell=cell,coef=rownames(z2)),
            data.frame(z3,factor='Factor3',cell=cell,coef=rownames(z3)),
            data.frame(z4,factor='Dx_donor',cell=cell,coef=rownames(z4)))
}
zdf=reshape2::melt(Vp[Vp$ID=='agg',c('scale.Factor1.','scale.Factor2.','scale.Factor3.','scale.DxPC1.','assay')],
                   id.vars=c('assay'))
pdf(sprintf('MOFA_Dx_dev_cellcompositoincor_%s.pdf',nfactor),width=5,height = 5)
print(ggplot(zdf)+geom_bar(aes(x=assay,fill=variable,y=value),stat='identity')+
        theme_bw()+ggtitle('Dx_cell var expalined')+
        theme(axis.text.x = element_text(angle = 90,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
####
cmat=cor(Dx_dev,use='complete.obs',method='spearman')
diag(cmat)=NA
zdf1=reshape2::melt(cmat)
zdf1$Var1=factor(zdf1$Var1,levels=c('Astro','Oligo','OPC','Micro_PVM','EN','IN'))
zdf1$Var2=factor(zdf1$Var2,levels=c('Astro','Oligo','OPC','Micro_PVM','EN','IN'))
corrplot::corrplot.mixed(cmat[c('Astro','Oligo','OPC','Micro_PVM','EN','IN'),c('Astro','Oligo','OPC','Micro_PVM','EN','IN')],
                         upper.col = rev(corrplot::COL2()),upper='square')
corrplot::corrplot(cmat[c('Astro','Oligo','OPC','Micro_PVM','EN','IN'),c('Astro','Oligo','OPC','Micro_PVM','EN','IN')],
                   col = rev(corrplot::COL2()),type = 'upper',method  = 'square',diag = FALSE)
print(ggplot(zdf1)+geom_tile(aes(x=Var1,y=Var2,fill=value))+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))

zInfo=metaInfo[metaInfo$SubID%in%rownames(Fs),]
zInfo$Dx_dev=rowMeans(scale(Dx_dev),na.rm=T)[rownames(zInfo)]
rho_dx=cor(zInfo[,c('DxPC1','CERAD','BRAAK_AD','CDRScore','Dx_dev')],data.frame(Fs)[zInfo$SubID,],use='complete.obs',method='spearman')
cmat=sapply(c('DxPC1','Dx_dev','Age'),function(x){
  cor(data.frame(Fs)[zInfo$SubID,],zInfo[, x],use='complete.obs',method='spearman')
})
rownames(cmat)=colnames(Fs)
zdf=reshape::melt(cmat)
zdf$lab=round(zdf$value,2)
zdf[abs(zdf$lab)<0.2,]$lab=''
zdf$Fs=factor(zdf$X1,levels=colnames(Fs))
zdf$coef=factor(zdf$X2,levels=c("DxPC1",'Dx_dev','Age'))
print(ggplot(zdf)+
        geom_tile(aes(y=Fs,x=coef,fill=value))+
        geom_text(aes(y=Fs,x=coef,label=lab))+
        #scale_fill_gradient2(low = "#3B4CC0",mid = "white",high = "#B40426",limits=c(-1,1),name = 'rho')+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
print(ggplot(zdf[zdf$Fs%in%F_sel,])+
        geom_tile(aes(y=Fs,x=coef,fill=value))+
        geom_text(aes(y=Fs,x=coef,label=lab))+
        #        scale_fill_gradient2(low = "#3B4CC0",mid = "white",high = "#B40426",limits=c(-1,1),name = 'rho')+
        scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
        scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
        theme_bw()+coord_equal()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black")))
for(x in names(cobjs)){
  cobj=cobjs[[x]]
  
  cmat=cor(data.frame(Fs)[metaInfo$SubID,],
           data.frame(t(cobj$E))[rownames(metaInfo),],
           use='complete.obs',method='spearman')
  
  zdf=reshape::melt(cmat)
  zdf$lab=round(zdf$value,2)
  zdf[abs(zdf$lab)<0.2,]$lab=''
  zdf$Fs=factor(zdf$X1,levels=colnames(Fs))
  print(ggplot(zdf)+
          geom_tile(aes(y=Fs,x=X2,fill=value))+
          geom_text(aes(y=Fs,x=X2,label=lab))+
          #scale_fill_gradient2(low = "#3B4CC0",mid = "white",high = "#B40426",limits=c(-1,1),name = 'rho')+
          scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
          scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
          theme_bw()+coord_equal()+ggtitle(x)+
          theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
                axis.text.y = element_text(colour = "black")))
  print(ggplot(zdf[zdf$Fs%in%F_sel,])+
          geom_tile(aes(y=Fs,x=X2,fill=value))+
          geom_text(aes(y=Fs,x=X2,label=lab))+
          #scale_fill_gradient2(low = "#3B4CC0",mid = "white",high = "#B40426",limits=c(-1,1),name = 'rho')+
          scale_fill_gradientn(colors=rev(RColorBrewer::brewer.pal(9,'RdBu')),limits=c(-1,1))+
          scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+
          theme_bw()+coord_equal()+ggtitle(x)+
          theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
                axis.text.y = element_text(colour = "black")))
}
dev.off()
###

zInfo=metaInfo[metaInfo$SubID%in%rownames(Fs),]
zInfo$Dx_dev=rowMeans(scale(Dx_dev),na.rm=T)[rownames(zInfo)]
zInfo=data.frame(data.frame(Fs)[zInfo$SubID,],zInfo)

m1=calcVarPart(lme4::lmer(Dx_dev~scale(Age)+(1|Sex)+scale(DxPC1)+scale(PMI..min.)+(1|SubID)+(1|APOE)+(1|BrainRegion2),zInfo,REML=F))
m2=calcVarPart(lme4::lmer(Dx_dev~scale(Age)+(1|Sex)+scale(DxPC1)+scale(PMI..min.)+scale(Factor1)+scale(Factor2)+scale(Factor3)+(1|SubID)+(1|APOE)+(1|BrainRegion2),zInfo,REML=F))
zdf=data.frame(vp=c(m1,m2),coef=c(names(m1),names(m2)),
               mod=rep(c('base','b+Factors'),c(length(m1),length(m2))))
write.csv(zdf,file='Dx_dev_varExp_base_F123.csv')
zdf=data.frame(vp=c(m1[c('SubID')],m1['Residuals'],1-sum(m1[c('SubID','Residuals')]),
                    c(m2['SubID'],m2['Residuals'],m2['scale(Factor1)'],m2['scale(Factor2)'],
                      m2['scale(Factor3)'],1-sum(m2[c('SubID','Residuals','scale(Factor1)','scale(Factor2)',
                                                      'scale(Factor3)')]))),
               mod=rep(c('raw','Fs'),c(3,6)),
               coef=c('SubID','Resid','other','SubID','Resid','F1','F2','F3','other'))
zdf$coef=factor(zdf$coef,levels=rev(c('SubID','F1','F2','F3','other','Resid')))
col1=c('white','grey',rev(c("#FDE725FF","#21908CFF","#472D7BFF")),'grey25')
names(col1)=levels(zdf$coef)
zdf$mod=factor(zdf$mod,levels=c("raw",'Fs'))
pdf('Dx_dev_varExp_base_F123.pdf',width=4,height = 4)
print(ggplot(zdf)+geom_bar(aes(x=mod,y=vp*100,fill=coef),stat='identity',width=0.5,color='black')+
        scale_fill_manual(values=col1)+
        scale_y_continuous(expand = c(0,0))+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 2))
dev.off()
######## cell composition change all subclass################################
Fs=data.frame(Fs,Fmeta=scale(Fs[,2])-scale(Fs[,3]))
cobj=cobjs$subclass
zInfo=metaInfo[colnames(cobj),]
zInfo=data.frame(zInfo,data.frame(Fs)[as.character(zInfo$SubID),])
zInfo$Dx_dev=rowMeans(scale(Dx_dev[,c('Astro','Oligo','OPC','Micro_PVM')]),na.rm=T)[rownames(zInfo)]
zInfo$Dx_cell=rowMeans(scale(Dx_cell_agg[,c('Astro','Oligo','OPC','Micro_PVM')]),na.rm=T)[rownames(zInfo)]

form0=~scale(Age)+Sex+BrainRegion2+scale(Dx_cell)+scale(PMI..min.)+(1|SubID)
zfit0=eBayes(dream(cobj,form0,zInfo))
z0=topTable(zfit0,'scale(Dx_cell)',number=Inf)
form1=~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+(1|SubID)
zfit1=eBayes(dream(cobj,form1,zInfo))
z1=topTable(zfit1,'scale(DxPC1)',number=Inf)
form2=~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+(1|SubID)+scale(Dx_dev)
zfit2=eBayes(dream(cobj,form2,zInfo))
z2=topTable(zfit2,'scale(Dx_dev)',number=Inf)
write.csv(z0,file=sprintf('Dx_cell_subclass_crumblr_%s.csv',nfactor))
write.csv(z1,file=sprintf('Dx_donor_subclass_crumblr_%s.csv',nfactor))
write.csv(z2,file=sprintf('Dx_dev_subclass_crumblr_%s.csv',nfactor))
#zdf=data.frame(cell=rownames(z1),Dx_donor=z1$t,Dx_dev=z2[rownames(z1),])

z1$se=z1$logFC/z1$t
z2$se=z2$logFC/z2$t
zdf=data.frame(cell=rownames(z1),
               Dx_donor=z1$logFC,donor_se=z1$se,
               Dx_dev=z2[rownames(z1),]$logFC,
               dev_se=z2[rownames(z1),]$se,
               FDR=z2[rownames(z1),]$adj.P.Val)
sig=c('','FDR sig')[(zdf$FDR<0.05)+1]
col1=c('grey40','#B40426')
names(col1)=c('','FDR sig')
colorcode=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/step3/files/color_code.csv',row.names=1)
colormap_class=colorcode[unique(anno$subclass),]
names(colormap_class)=unique(anno$subclass)
log10range=c(2,5,10,15,20,25,30)
pdf(sprintf('subclass_cell_composition_change_Dx_dev_vs_Dx_donor_%s.pdf',nfactor))
print(ggplot(zdf)+
        geom_smooth(aes(x=Dx_donor,y=Dx_dev),method='lm',linewidth = 0.7)+
        geom_point(aes(x=Dx_donor,y=Dx_dev,size=-log10(FDR),color=cell))+
        geom_linerange(aes(x=Dx_donor,ymin=Dx_dev-dev_se,ymax=Dx_dev+dev_se,color=cell),linewidth = 0.35)+
        geom_linerange(aes(y=Dx_dev,xmin=Dx_donor-donor_se,xmax=Dx_donor+donor_se,color=cell),linewidth = 0.35)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        geom_vline(xintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('logFC (Dx_donor)')+ylab('logFC (Dx_dev)')+
        ggtitle(paste0('rho ',round(cor(zdf$Dx_donor,zdf$Dx_dev,use='complete.obs',method='spearman'),2)))+
        scale_color_manual(values=colormap_class)+
        scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))])+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
print(ggplot(zdf)+
        geom_smooth(aes(x=Dx_donor,y=Dx_dev),method='lm',linewidth = 0.7)+
        geom_point(aes(x=Dx_donor,y=Dx_dev,size=-log10(FDR),color=sig))+
        geom_linerange(aes(x=Dx_donor,ymin=Dx_dev-dev_se,ymax=Dx_dev+dev_se,color=sig),linewidth = 0.35)+
        geom_linerange(aes(y=Dx_dev,xmin=Dx_donor-donor_se,xmax=Dx_donor+donor_se,color=sig),linewidth = 0.35)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        geom_vline(xintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('logFC (Dx_donor)')+ylab('logFC (Dx_dev)')+
        scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))])+
        ggtitle(paste0('rho ',round(cor(zdf$Dx_donor,zdf$Dx_dev,use='complete.obs',method='spearman'),2)))+
        geom_text_repel(aes(x=Dx_donor,y=Dx_dev,label=cell),
                        color = "black",seed = 123, box.padding = 0.45,
                        point.padding = 0.25,force = 2,force_pull = 0.1,
                        max.overlaps = Inf,
                        max.time = 5, min.segment.length = 0,
                        segment.color = "grey35",linewidth = 0.3)+
        scale_color_manual(values=col1)+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
dev.off()
pdf(sprintf('subclass_cell_composition_change_MOFA_vs_Dx_donor_%s.pdf',nfactor))
for(xx in c('Factor1',"Factor2",'Factor3','Fmeta')){
  form=as.formula(sprintf('~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+(1|SubID)+scale(%s)',xx))
  zfit=eBayes(dream(cobj,form,zInfo))
  z2=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  write.csv(z2,file=sprintf('%s_subclass_crumblr_%s.csv',xx,nfactor))
  z2$se=z2$logFC/z2$t
  zdf=data.frame(cell=rownames(z1),
                 Dx_donor=z1$logFC,donor_se=z1$se,
                 Dx_dev=z2[rownames(z1),]$logFC,
                 dev_se=z2[rownames(z1),]$se,
                 FDR=z2[rownames(z1),]$adj.P.Val)
  sig=c('','FDR sig')[(zdf$FDR<0.05)+1]
  print(ggplot(zdf)+
          geom_smooth(aes(x=Dx_donor,y=Dx_dev),method='lm',linewidth = 0.7)+
          geom_point(aes(x=Dx_donor,y=Dx_dev,size=-log10(FDR),color=cell))+
          geom_linerange(aes(x=Dx_donor,ymin=Dx_dev-dev_se,ymax=Dx_dev+dev_se,color=cell),linewidth = 0.35)+
          geom_linerange(aes(y=Dx_dev,xmin=Dx_donor-donor_se,xmax=Dx_donor+donor_se,color=cell),linewidth = 0.35)+
          geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
          geom_vline(xintercept = 0,linetype = "dotted",linewidth = 0.5)+
          xlab('logFC (Dx_donor)')+ylab(sprintf('logFC (%s)',xx))+
          ggtitle(paste0('rho ',round(cor(zdf$Dx_donor,zdf$Dx_dev,use='complete.obs',method='spearman'),2)))+
          scale_color_manual(values=colormap_class)+
          scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))])+
          theme_bw()+
          theme(axis.text.x = element_text(colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 1))
  print(ggplot(zdf)+
          geom_smooth(aes(x=Dx_donor,y=Dx_dev),method='lm',linewidth = 0.7)+
          geom_point(aes(x=Dx_donor,y=Dx_dev,size=-log10(FDR),color=sig))+
          geom_linerange(aes(x=Dx_donor,ymin=Dx_dev-dev_se,ymax=Dx_dev+dev_se,color=sig),linewidth = 0.35)+
          geom_linerange(aes(y=Dx_dev,xmin=Dx_donor-donor_se,xmax=Dx_donor+donor_se,color=sig),linewidth = 0.35)+
          geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
          geom_vline(xintercept = 0,linetype = "dotted",linewidth = 0.5)+
          xlab('logFC (Dx_donor)')+ylab(sprintf('logFC (%s)',xx))+
          ggtitle(paste0('rho ',round(cor(zdf$Dx_donor,zdf$Dx_dev,use='complete.obs',method='spearman'),2)))+
          geom_text_repel(aes(x=Dx_donor,y=Dx_dev,label=cell),
                          color = "black",seed = 123, box.padding = 0.45,
                          point.padding = 0.25,force = 2,force_pull = 0.1,
                          max.overlaps = Inf,
                          max.time = 5, min.segment.length = 0,
                          segment.color = "grey35",linewidth = 0.3)+
          scale_color_manual(values=col1)+
          scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))])+
          theme_bw()+
          theme(axis.text.x = element_text(colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 1))
}
dev.off()
if(F){
  form=~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+(1|SubID)+scale(Factor1)+scale(Factor2)+scale(Factor3)
  zfit=eBayes(dream(cobj,form,zInfo))
  z2=do.call(rbind,lapply(paste0('Factor',1:3),function(xx){
    z=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
    data.frame(coef=xx, cell=rownames(z),z)
  }))
  write.csv(z2,file=sprintf('%s_subclass_crumblr_%s.csv',xx,nfactor))
  
}
###################vul neurons only #####################################

zInfo=metaInfo[colnames(Vul),]
zInfo=data.frame(zInfo,data.frame(Fs)[as.character(zInfo$SubID),])
zInfo$Dx_dev=rowMeans(scale(Dx_dev[,c('Astro','Oligo','OPC','Micro_PVM')]),na.rm=T)[rownames(zInfo)]
zInfo$Dx_cell=rowMeans(scale(Dx_cell_agg[,c('Astro','Oligo','OPC','Micro_PVM')]),na.rm=T)[rownames(zInfo)]

form0=~scale(Age)+Sex+BrainRegion2+scale(Dx_cell)+scale(PMI..min.)+(1|SubID)
zfit=eBayes(dream(Vul,form0,zInfo))
z0=topTable(zfit,'scale(Dx_cell)',number=Inf)
z0$se=z0$logFC/z0$t

form0=~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+(1|SubID)
zfit=eBayes(dream(Vul,form0,zInfo))
z1=topTable(zfit,'scale(DxPC1)',number=Inf)
z1$se=z1$logFC/z1$t
zdf1=data.frame(cell=rownames(z1),logFC=z1$logFC,se=z1$se,FDR=z1$adj.P.Val,
                group='Dx_donor')
df=zdf1
write.csv(z0,file=sprintf('Dx_cell_vulneuron_crumblr_%s.csv',nfactor))
write.csv(z1,file=sprintf('Dx_donor_vulneuron_crumblr_%s.csv',nfactor))
pdf(sprintf('Vul_neuron_crumblr_%s.pdf',nfactor),width=5,height = 5)
for(xx in c('Dx_dev','Factor1','Factor2','Factor3','Fmeta')){
  form1=as.formula(sprintf('~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+(1|SubID)+scale(%s)',xx))  
  zfit=eBayes(dream(Vul,form1,zInfo))
  z2=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  write.csv(z2,file=sprintf('%s_vulneuron_crumblr_%s.csv',xx,nfactor))
  z2$se=z2$logFC/z2$t
  zdf2=data.frame(cell=rownames(z2),logFC=z2$logFC,se=z2$se,FDR=z2$adj.P.Val,
                  group=xx)
  df=rbind(df,zdf2)
  zdf=rbind(zdf1,zdf2)
  zdf$sig=c('','+')[(zdf$FDR<0.05)+1]
  print(ggplot(zdf)+
          geom_linerange(aes(x=group,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
          geom_point(aes(x=group,y=logFC,size=-log10(FDR)))+
          geom_text(aes(x=group,y=logFC,label=sig),color='white')+
          scale_color_manual(values=col1)+
          geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
          xlab('')+ylab('logFC')+
          scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))],limits = c(1,max(-log10(zdf$FDR))))+
          theme_bw()+facet_wrap(~cell)+
          theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 2))
  
}
df$sig=c('','+')[(df$FDR<0.05)+1]
print(ggplot(df[df$group%in%c('Dx_donor','Factor1','Factor2','Factor3'),])+
        geom_linerange(aes(x=group,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
        geom_point(aes(x=group,y=logFC,size=-log10(FDR)))+
        geom_text(aes(x=group,y=logFC,label=sig),color='white')+
        scale_color_manual(values=col1)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))],limits = c(1,max(-log10(zdf$FDR))))+
        theme_bw()+facet_wrap(~cell)+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 2))
#####
dev.off()
##### DAM####
cobj=cobjs$Micro_PVM
zInfo=metaInfo[colnames(cobj),]
zInfo=data.frame(zInfo,data.frame(Fs)[as.character(zInfo$SubID),])
zInfo$Dx_dev=rowMeans(scale(Dx_dev[,c('Astro','Oligo','OPC')]),na.rm=T)[rownames(zInfo)]
zInfo$Dx_cell=rowMeans(scale(Dx_cell_agg[,c('Astro','Oligo','OPC')]),na.rm=T)[rownames(zInfo)]
form0=~scale(Age)+Sex+BrainRegion2+scale(Dx_cell)+scale(PMI..min.)+(1|SubID)
zfit=eBayes(dream(cobj,form0,zInfo))
z0=topTable(zfit,'scale(Dx_cell)',number=Inf)
z0$se=z0$logFC/z0$t


form0=~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+(1|SubID)
zfit=eBayes(dream(cobj,form0,zInfo))
z1=topTable(zfit,'scale(DxPC1)',number=Inf)
z1$se=z1$logFC/z1$t
zdf1=data.frame(cell=rownames(z1),logFC=z1$logFC,se=z1$se,FDR=z1$adj.P.Val,
                group='Dx_donor')
df=zdf1
write.csv(z0,file=sprintf('Dx_cell_Micro_PVM_crumblr_%s.csv',nfactor))
write.csv(z1,file=sprintf('Dx_donor_Micro_PVM_crumblr_%s.csv',nfactor))
pdf(sprintf('DAM_crumblr_%s.pdf',nfactor),width=5,height = 5)
for(xx in c('Dx_dev','Factor1','Factor2','Factor3','Fmeta')){
  form1=as.formula(sprintf('~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+(1|SubID)+scale(%s)',xx))  
  zfit=eBayes(dream(cobj,form1,zInfo))
  z2=topTable(zfit,sprintf('scale(%s)',xx),number=Inf)
  write.csv(z2,file=sprintf('%s_DAM_crumblr_%s.csv',xx,nfactor))
  z2$se=z2$logFC/z2$t
  zdf2=data.frame(cell=rownames(z2),logFC=z2$logFC,se=z2$se,FDR=z2$adj.P.Val,
                  group=xx)
  df=rbind(df,zdf2)
  zdf=rbind(zdf1,zdf2)
  zdf$sig=c('','+')[(zdf$FDR<0.05)+1]
  print(ggplot(zdf[zdf$cell=='DAM',])+
          geom_linerange(aes(x=group,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
          geom_point(aes(x=group,y=logFC,size=-log10(FDR)))+
          geom_text(aes(x=group,y=logFC,label=sig),color='white')+
          scale_color_manual(values=col1)+
          geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
          xlab('')+ylab('logFC')+
          #scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))],limits = c(1,max(-log10(zdf$FDR))))+
          theme_bw()+facet_wrap(~cell)+
          theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 2))
  print(ggplot(zdf)+
          geom_linerange(aes(x=group,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
          geom_point(aes(x=group,y=logFC,size=-log10(FDR)))+
          geom_text(aes(x=group,y=logFC,label=sig),color='white')+
          scale_color_manual(values=col1)+
          geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
          xlab('')+ylab('logFC')+
          #scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))],limits = c(1,max(-log10(zdf$FDR))))+
          theme_bw()+facet_wrap(~cell)+
          theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
                axis.text.y = element_text(colour = "black"),
                strip.background = element_rect(fill = NA),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),aspect.ratio = 2))
  
}
df$sig=c('','+')[(df$FDR<0.05)+1]
print(ggplot(df[df$group%in%c('Dx_donor','Factor1','Factor2','Factor3'),])+
        geom_linerange(aes(x=group,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
        geom_point(aes(x=group,y=logFC,size=-log10(FDR)))+
        geom_text(aes(x=group,y=logFC,label=sig),color='white')+
        scale_color_manual(values=col1)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))],limits = c(1,max(-log10(zdf$FDR))))+
        theme_bw()+facet_wrap(~cell)+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 2))
print(ggplot(df[df$group%in%c('Dx_donor','Factor1','Factor2','Factor3') & df$cell=='DAM',])+
        geom_linerange(aes(x=group,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
        geom_point(aes(x=group,y=logFC,size=-log10(FDR)))+
        geom_text(aes(x=group,y=logFC,label=sig),color='white')+
        scale_color_manual(values=col1)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        scale_size(breaks=log10range[log10range<max(-log10(zdf$FDR))],limits = c(1,max(-log10(zdf$FDR))))+
        theme_bw()+facet_wrap(~cell)+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 2))
#####
dev.off()
pdf(sprintf('Reilience_%s.pdf',nfactor),width=5,height = 5)
#agg1=lgb%>%group_by(order,class)%>%summarise(Dx_cell=mean(Dx_cell))
#agg1=reshape2::acast(agg1,order~class,value.var='Dx_cell')
#z=rowMeans(scale(agg1),na.rm=T)
z=rowMeans(scale(Dx_dev),na.rm=T)
zInfo=data.frame(Dx_dev=z,Dx_cell=rowMeans(scale(Dx_cell_agg[,-1]),na.rm=T)[names(z)],
                 control=rowMeans(class_mat[,grepl('control',colnames(class_mat))])[names(z)],
                 late=rowMeans(class_mat[,grepl('late',colnames(class_mat))])[names(z)],
                 metaInfo[names(z),])

zInfo=data.frame(zInfo,Fs[as.character(zInfo$SubID),])
form=~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+scale(Cognitive_ResilienceUse)+(1|APOE)
sh=zInfo[zInfo$SubID%in%rownames(Fs),]
zfit1=dream(t(scale(Fs[sh$SubID,])),form,zInfo[rownames(sh),])
z1=topTable(eBayes(zfit1),coef='scale(Cognitive_ResilienceUse)')


write.csv(z1,file=sprintf('MOFA_Cognitive_Resilience_%s.csv',nfactor))
zdf=data.frame(z1,cell=rownames(z1),se=z1$logFC/z1$t)
zdf$sig=c('','+')[(zdf$adj.P.Val<0.05)+1]
zdf$cell=factor(zdf$cell,levels=c('Fmeta',paste0('Factor',1:10)))
print(ggplot(zdf)+
        geom_linerange(aes(x=cell,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
        geom_point(aes(x=cell,y=logFC,size=-log10(adj.P.Val)))+
        geom_text(aes(x=cell,y=logFC,label=sig),color='white')+
        #scale_color_manual(values=colormap_class)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))

print(ggplot(zdf[zdf$cell%in%c('Fmeta',paste0('Factor',1:3)),])+
        geom_linerange(aes(x=cell,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
        geom_point(aes(x=cell,y=logFC,size=-log10(adj.P.Val)))+
        geom_text(aes(x=cell,y=logFC,label=sig),color='white')+
        #scale_color_manual(values=colormap_class)+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))


####


form=~scale(Age)+Sex+BrainRegion2+scale(DxPC1)+scale(PMI..min.)+(1|SubID)+scale(Cognitive_ResilienceUse)+(1|APOE)
zfit1=dream(t(scale(Dx_dev)),form,zInfo[rownames(Dx_dev),])
z1=topTable(zfit1,coef='scale(Cognitive_ResilienceUse)')
zdf=data.frame(z1,cell=rownames(z1),se=z1$logFC/z1$t)
write.csv(z1,file=sprintf('Dx_devCognitive_Resilience_%s.csv',nfactor))
colormap_class=colorcode[unique(zdf$cell),]
names(colormap_class)=unique(zdf$cell)
zdf$sig=c('','+')[(zdf$adj.P.Val<0.05)+1]
zdf$cell=factor(zdf$cell,levels=c('Astro','Oligo','OPC','Micro_PVM','EN','IN'))
print(ggplot(zdf)+
        geom_linerange(aes(x=cell,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
        geom_point(aes(x=cell,y=logFC,size=-log10(adj.P.Val)))+
        geom_text(aes(x=cell,y=logFC,label=sig),color='white')+
        #scale_color_manual(values=colormap_class)+
        scale_size(limits=c(1,max(-log10(zdf$adj.P.Val))))+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))
print(ggplot(zdf)+
        geom_linerange(aes(x=cell,ymin=logFC-se,ymax=logFC+se),linewidth = 0.35)+
        geom_point(aes(x=cell,y=logFC,size=-log10(adj.P.Val),color=cell))+
        geom_text(aes(x=cell,y=logFC,label=sig),color='white')+
        scale_color_manual(values=colormap_class)+
        scale_size(limits=c(1,max(-log10(zdf$adj.P.Val))))+
        geom_hline(yintercept = 0,linetype = "dotted",linewidth = 0.5)+
        xlab('')+ylab('logFC')+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1,colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))
#zInfo$Dx_dev=rowMeans(scale(Dx_dev),na.rm=T)[rownames(zInfo)]
zdf=zInfo[!is.na(zInfo$Cognitive_ResilienceUse),]
zdf$group1=c('least','mid','resilience')[cut(zdf$Cognitive_ResilienceUse,c(-Inf,quantile(zdf$Cognitive_ResilienceUse,c(0.2,0.8)),Inf))]
zdf$group2=c('other','resilience')[(zdf$Cognitive_ResilienceUse>quantile(zdf$Cognitive_ResilienceUse,0.8))+1]
Dx_cutoff=read.csv('/sc/arion/projects/CommonMind/roussp01a/snmulti/DE/metadata/psychAD_meta_cutoff.csv',row.names=1)
print(ggplot(zdf)+geom_smooth(aes(x=DxPC1,y=Dx_cell,color=group2))+
        geom_vline(xintercept = Dx_cutoff$x,linetype=2)+
        xlim(min(zInfo$DxPC1),max(zInfo$DxPC1))+
        theme_bw()+
        theme(axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              strip.background = element_rect(fill = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),aspect.ratio = 1))
dev.off()


################
############
if(T){
  library(AnnotationDbi)
  library(GOSemSim)
  library(fgsea)
  library(clusterProfiler) # https://yulab-smu.top/biomedical-knowledge-mining-book/021-go.html
  library(rrvgo) # https://www.bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tidytext)
  rownames(genes)=genes$gene_name
  map1=bitr(rownames(genes[genes$filtered,]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  map1=map1[!duplicated(map1$SYMBOL),]
  rownames(map1)=map1$SYMBOL
  hsGO <- godata(annoDb = 'org.Hs.eg.db', ont="BP")
  GO_DATA <- clusterProfiler:::get_GO_data('org.Hs.eg.db', 'BP', 'ENTREZID')
  gs=GO_DATA$PATHID2EXTID
  gs=gs[sapply(gs,length)>=10 & sapply(gs,length)<=1000]
  
  a=data.frame(W,do.call(rbind,strsplit(rownames(W),'\\:')))
  a=a[a$X2%in%rownames(map1),]
  a$ID=map1[as.character(a$X2),]$ENTREZID
  ego=data.frame()
  for(x in paste0('Factor',1:3)){
    statsets=lapply(split(a,a$X1),function(z){o=z[,x];names(o)=z$ID;o})
    ego=rbind(ego,
              do.call(rbind,lapply(names(statsets),function(ii){
                data.frame(fgsea(pathways = gs, 
                                 stats    = statsets[[ii]],
                                 minSize  = 10,
                                 maxSize  = 500),cell=ii,factor=x)
              })))
    
  }
  sig=unique(ego[ego$padj<0.05,]$pathway)
  ind=expand.grid(factor(sig,levels=sig),factor(sig,levels=sig))
  ind=ind[as.numeric(ind[,1])>as.numeric(ind[,2]),]
  sims=mapply(function(x,y){goSim(x, y, semData=hsGO, measure="Wang")},
              as.character(ind[,1]),as.character(ind[,2]))
  zmat=array(1,dim=rep(length(sig),2),dimnames = list(sig,sig))
  zmat[cbind(as.character(ind[,1]),as.character(ind[,2]))]=sims
  zmat[cbind(as.character(ind[,2]),as.character(ind[,1]))]=sims
  
  go_map=GO_DATA$PATHID2NAME[names(gs)]
  pdf(sprintf('pathway_F123_%s.pdf',nfactor),width=6,height = 6)
  z=plot_fgsea_go(ego=ego,sim_mat=zmat,go_map=go_map,n_terms_per_factor=4,
                  cell_order=c('Astro','Oligo','OPC','Micro_PVM','EN','IN'))
  print(z$plot)
  dev.off()
}

