# This script is for generating tables, figures to be used by the report 

# Load libraries
library(DESeq2)
library(httr)
library(tximport)
library(plotly)
library(readxl)
library(cowplot)
library(ashr)
library(repr)
library(pheatmap)
library(RColorBrewer)
library(DT)
library(patchwork)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Ss.eg.db)
library(org.Hs.eg.db)
library(ReactomePA)
library(UpSetR)
library(tidyverse)

#########################################################################################

# Get the directories to do
fp_clus <- "../../results/Make_Clusters"
fp_deg <- "../../results/Make_DEG_Lists/"

#########################################################################################

make_ortho_df <- function(){
  # Load the pig to human gene converting file from BioMart
  fpmart <- "../../data/orthologs/"
  dfmart <- read_tsv(paste0(fpmart,"pig_mart_export.txt"),show_col_types = FALSE) %>% rename(`GeneID` = 1) %>%
    filter(`Human orthology confidence [0 low, 1 high]` > 0.5) %>%
    dplyr::select(`GeneID`,`Gene name`,`Human gene stable ID`)
  # Load in file with more human gene info
  dfmart_hs <- read_tsv(paste0(fpmart,"human_mart_export.txt"),show_col_types = FALSE) %>%
    dplyr::select(-`Gene Synonym`) %>%
    distinct() %>%
    rename(`Human gene stable ID` = 1,`Human gene name` = 2, 
           `Human gene description` = 3)
  dfgene <- left_join(dfmart,dfmart_hs, by="Human gene stable ID")
  return(dfgene)
}
dfgene <- make_ortho_df()

make_deseq2_obj <- function(deseq_type){
  fp_count <- "../../data/mapped_reads"
  fp_count1b <- "../../data/mapped_reads/run1"
  fp_count2b <- "../../data/mapped_reads/run2"
  fp_meta <- "../../data/"
  dirs1 <-  list.dirs(fp_count1b,full.names = F)
  dirs1 <- dirs1[-1]
  files1 <- file.path(paste0(fp_count1b,"/",dirs1),"quant.sf")
  dirs2 <-  list.dirs(fp_count2b,full.names = F)
  dirs2 <- dirs2[-1]
  files2 <- file.path(paste0(fp_count2b,"/",dirs2),"quant.sf")
  files <- c(files1,files2)
  # read in txt2gene file
  txt2gene <- read_tsv(paste0(fp_count,"/salmon_tx2gene.tsv"),col_names = F,show_col_types = FALSE)
  # read in metadata
  coldata <- read_csv(paste0(fp_meta,"Sample info + gender with heart.csv"),show_col_types = FALSE) %>%
    filter_all(any_vars(!is.na(.)))
  coldata <- column_to_rownames(coldata,var="SampleName")
  # order the metadata the same way the counts were read in
  dirs2b <- paste0("X",dirs2)
  coldata <- coldata[ c(dirs1,dirs2b), ]
  if (deseq_type == "AllTogether"){
    # don't do anything
  } else {
    inds_tmp <- which(coldata$Tissue==deseq_type)
    files <- files[inds_tmp]
    coldata <- coldata[inds_tmp,]
  }
  txi <- tximport(files, type="salmon", tx2gene=txt2gene)
  ddsgrp <- DESeqDataSetFromTximport(txi,
                                     colData = coldata,
                                     design = ~ Group)
  ddsgrp <- DESeq(ddsgrp)
  return(list(ddsgrp,coldata))
}

make_results <- function(deseq_obj,group1,group2,atissue){
  res <- results(deseq_obj,alpha=myLFC,contrast=c("Group",
                                                  sprintf("%s%s",group1,atissue),
                                                  sprintf("%s%s",group2,atissue)))
  resLFC <- lfcShrink(deseq_obj,
                      res = res,
                      type="ashr")
  fp_count <- "../../data/mapped_reads"
  symdf <- read_tsv(paste0(fp_count,"/salmon.merged.gene_counts_length_scaled.tsv"),show_col_types = FALSE)
  resLFC$GeneName <-  symdf$gene_name
  resLFC$GeneID <-  symdf$gene_id
  resOrderLFC <- resLFC[order(resLFC$padj),]
  return(resOrderLFC)
}

make_GO_tables <- function(sig_genes,resOrderLFC){ 
  uni_genes <- left_join(resOrderLFC,dfgene, by="GeneID") %>% drop_na(`Human gene stable ID`) %>% pull(var=`Human gene stable ID`) %>% na.omit()
  ego <- enrichGO(gene          = sig_genes,
                  universe      = uni_genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  minGSSize     = 5,
                  maxGSSize     = 500,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  if (!is.null(ego)){
    ego <- filter(ego, p.adjust < 0.05)
  }
  return(ego)
}

make_reactome_tables <- function(sig_genes,resOrderLFC){
  if (length(sig_genes)>0){
    geneids <- stack(mget(sig_genes, org.Hs.egENSEMBL2EG, ifnotfound = NA)) %>% na.omit() %>% pull(values)
  } else {
    geneids <- c()
  }
  uni_genes <- left_join(resOrderLFC,dfgene, by="GeneID") %>% drop_na(`Human gene stable ID`) %>% pull(var=`Human gene stable ID`) %>% na.omit()
  uniids <- stack(mget(uni_genes, org.Hs.egENSEMBL2EG, ifnotfound = NA)) %>% na.omit() %>% pull(values)
  eRP <- enrichPathway(gene          = geneids, 
                       universe      = uniids,
                       organism     = "human",
                       minGSSize     = 5,
                       maxGSSize     = 500,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)
  if (!is.null(eRP)){
    eRP <- filter(eRP, p.adjust < 0.05)
    eRP <- setReadable(eRP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }
  return(eRP)
}

make_kegg_tables <- function(sig_genes,resOrderLFC){
  if (length(sig_genes)>0){
    geneids <- stack(mget(sig_genes, org.Hs.egENSEMBL2EG, ifnotfound = NA)) %>% na.omit() %>% pull(values)
  } else {
    geneids <- c()
  }
  uni_genes <- left_join(resOrderLFC,dfgene, by="GeneID") %>% drop_na(`Human gene stable ID`) %>% pull(var=`Human gene stable ID`) %>% na.omit()
  uniids <- stack(mget(uni_genes, org.Hs.egENSEMBL2EG, ifnotfound = NA)) %>% na.omit() %>% pull(values)
  kegg <- enrichKEGG(gene          = geneids,
                     universe      = uniids,
                     organism      = 'hsa',
                     minGSSize     = 5,
                     maxGSSize     = 500,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2)
  if (!is.null(kegg)){
    kegg <- filter(kegg, p.adjust < 0.05)
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }
  return(kegg)
}

make_heatmap <- function(deseq_obj,group1,group2,coldata,cl_genes_sc,tiss,scale_type){
  # load a count table to get the gene IDs
  fp_count <- "../../data/mapped_reads"
  symdf <- read_tsv(paste0(fp_count,"/salmon.merged.gene_counts_length_scaled.tsv"),show_col_types = FALSE)
  # create counts assay
  myassay <- counts(deseq_obj,normalized=TRUE)
  row.names(myassay) <- symdf$gene_name
  # find gene inds
  all_genes_inds <- which(symdf$gene_name %in% cl_genes_sc)
  # subset the assay to those genes
  myassay <- myassay[all_genes_inds,]
  # get the correct columns 
  grp_names <- coldata$Group
  con1_inds <- which(coldata$Group == paste0(group1,tiss))
  con2_inds <- which(coldata$Group == paste0(group2,tiss))
  all_con_inds <- c(con1_inds,con2_inds)
  # subset the assay for correct labels
  myassay <- myassay[,all_con_inds]
  # scale the data
  if (scale_type == "Relative"){
    myassay <- t(scale(t(myassay))) %>% na.omit()
  } else if (scale_type == "CPM"){
    myassay <- log2(myassay + 1) %>% na.omit()
  }
  # make label df
  df_label <- as.data.frame(colData(deseq_obj))
  df_label  <- df_label %>% filter(Group %in% c(paste0(group1,tiss),paste0(group2,tiss)))
  df_label <- df_label[,"Group", drop=FALSE]
  # all_genes_inds <- which(symdf$gene_name %in% allgenes)
  if (dim(myassay)[1]>2){
    myrowclus <- TRUE
  } else {
    myrowclus <- FALSE
  }
  hm <- pheatmap(myassay, cluster_rows=myrowclus, show_rownames=FALSE,
                 cluster_cols=FALSE, annotation_col=df_label)
  return(hm)
}

#########################################################################################

clus_dirs <-list.dirs(fp_clus)
clus_dirs <- clus_dirs[-1]
# Main Loop
for(adir in clus_dirs){
  mybasename <- strsplit(basename(adir),"__")[[1]]
  mytissue <- mybasename[3]
  mygroups <- mybasename[2]
  mygroup1 <- strsplit(mygroups,"_")[[1]][1]
  mygroup2 <- strsplit(mygroups,"_")[[1]][2]
  if(!file.exists(paste0(adir,"/cluster_membership.tsv"))){
    next
  }
  print(adir)
  # make deseq object
  fun_ret <- make_deseq2_obj(mytissue)
  ddsgrp <- fun_ret[[1]]
  col_data <- fun_ret[[2]]
  # make results
  resOrderLFC_full <- make_results(ddsgrp,mygroup1,mygroup2,mytissue)
  resOrderLFC <- as.data.frame(resOrderLFC_full)
  # read cluster membership files
  df_mem <- read_tsv(paste0(adir,"/cluster_membership.tsv"),show_col_types = FALSE)
  # load the significant gene dataframe
  ifelse(grepl("__DownGenes__", adir, fixed = TRUE),
         df_sig <- read_csv(paste0(fp_deg,mygroup1,"_",mygroup2,"_",mytissue,"_down-regulated.csv"),show_col_types = FALSE),
         ifelse(grepl("__UpGenes__", adir, fixed = TRUE),
                df_sig <- read_csv(paste0(fp_deg,mygroup1,"_",mygroup2,"_",mytissue,"_up-regulated.csv"),show_col_types = FALSE),
                ifelse(grepl("__UpDownGenes__", adir, fixed = TRUE),
                       df_sig <- rbind(read_csv(paste0(fp_deg,mygroup1,"_",mygroup2,"_",mytissue,"_up-regulated.csv"),show_col_types = FALSE),read_csv(paste0(fp_deg,mygroup1,"_",mygroup2,"_",mytissue,"_down-regulated.csv"),show_col_types = FALSE)),
                       next)))
  # loop through the clusters
  num_clus <- dim(df_mem)[2]
  for(i2 in 1:num_clus){
    cl_genes <- df_mem %>% pull(i2) %>% na.omit()
    cl_genes_sc <- df_sig %>% filter(`Human gene stable ID` %in% cl_genes) %>% pull("GeneName") %>% na.omit() %>% unique()
    clus_num <- i2-1
    # make GO results
    ego <- make_GO_tables(cl_genes,resOrderLFC)
    if (dim(as.data.frame(ego))[1]>0 & !is.null(ego)){
      write_csv(as.data.frame(ego), paste0(adir,"/",clus_num,"_go.csv"))
      ggsave(paste0(adir,"/",clus_num,"_go-dotplot.pdf"),plot=dotplot(ego, showCategory=10),device="pdf",width=6,height=6)
      ego_simp <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
      if (dim(as.data.frame(ego_simp))[1]>0 & !is.null(ego_simp)){
        write_csv(as.data.frame(ego_simp), paste0(adir,"/",clus_num,"_go-simp-0.7.csv"))
        ggsave(paste0(adir,"/",clus_num,"_go-simp-0.7-dotplot.pdf"),plot=dotplot(ego_simp, showCategory=10),device="pdf",width=6,height=6)
      }
      ego_simp <- clusterProfiler::simplify(ego, cutoff=0.5, by="p.adjust", select_fun=min)
      if (dim(as.data.frame(ego_simp))[1]>0 & !is.null(ego_simp)){
        write_csv(as.data.frame(ego_simp), paste0(adir,"/",clus_num,"_go-simp-0.5.csv"))
        ggsave(paste0(adir,"/",clus_num,"_go-simp-0.5-dotplot.pdf"),plot=dotplot(ego_simp, showCategory=10),device="pdf",width=6,height=6)
      }
    }
    erp <- make_reactome_tables(cl_genes,resOrderLFC)
    if (dim(as.data.frame(erp))[1]>0 & !is.null(erp)){
      write_csv(as.data.frame(erp), paste0(adir,"/",clus_num,"_reactome.csv"))
      ggsave(paste0(adir,"/",clus_num,"_reactome-dotplot.pdf"),plot=dotplot(erp, showCategory=10),device="pdf",width=6,height=6)
    }
    keg <- make_kegg_tables(cl_genes,resOrderLFC)
    if (dim(as.data.frame(keg))[1]>0 & !is.null(keg)){
      write_csv(as.data.frame(keg), paste0(adir,"/",clus_num,"_kegg.csv"))
      ggsave(paste0(adir,"/",clus_num,"_kegg-dotplot.pdf"),plot=dotplot(keg, showCategory=10),device="pdf",width=6,height=6)
    }
    p_hm <- make_heatmap(ddsgrp,mygroup1,mygroup2,col_data,cl_genes_sc,mytissue,"Relative")
    ggsave(paste0(adir,"/",clus_num,"_relativeheatmap.pdf"),plot=p_hm,device="pdf",width=7,height=7)
    p_hm <- make_heatmap(ddsgrp,mygroup1,mygroup2,col_data,cl_genes_sc,mytissue,"CPM")
    ggsave(paste0(adir,"/",clus_num,"_cpmheatmap.pdf"),plot=p_hm,device="pdf",width=7,height=7)
  }
}
