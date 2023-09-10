library(DESeq2)
library(dplyr)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# DESeq2 function that takes in reads matrix and runs
# gene_mat: df of read counts from featureCounts output, keeping only counts
# colData: for each col in gene_mat, specify whether it's ctrl/exp and paired/single ended
run_deseq <- function(gene_mat, colData){
  dds <- DESeqDataSetFromMatrix(countData = gene_mat,
                                colData = colData,
                                design = ~ condition)
  keep <- rowSums(counts(dds)) >= 10
  de_results <- dds[keep,] %>% DESeq() %>% results() %>% data.frame()
  de_modified = de_results %>%
    mutate(ENSG = gsub("\\..*","", row.names(de_results))) %>%
    filter(!is.na(padj))
  
  de_modified = de_modified %>%
    mutate(Symbol = mapIds(org.Hs.eg.db, keys=de_modified$ENSG, column="SYMBOL", keytype="ENSEMBL", multiVals="first")) %>%
    mutate(Color = case_when(
      padj >= 0.01 ~ 'ns',
      log2FoldChange > 1 ~ 'inc',
      log2FoldChange < -1 ~ 'dec',
      .default = 'nfc'))
    
  return(de_modified)
}

# read in gene matrix
gene_mat = read.csv('', skip=1, sep='\t', row.names=1) %>%
  dplyr::select(c(6:11))

# set colData according to each experiment
colData = data.frame(condition=c(rep('treatment', 3), rep('control', 3)),
                     type=rep('paired-end', 6),
                     row.names=colnames(gene_mat))
de_modified = run_deseq(gene_mat, colData)

# write data as file
write.csv(x=de_modified, file='')