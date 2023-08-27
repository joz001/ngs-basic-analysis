library("DESeq2")
library('dplyr')
library('tidyverse')

# change directory to one with featureCount output read matrix
gene_mat = read.csv('/data/gras1_rnaseq/FeatureCountOutput/gras1_kd_gene_counts', skip=1, sep='\t', row.names=1) %>%
  select(c(6:11))
gene_mat = read.csv('/data/nkap_rnaseq/gene_counts/nkap_try2_stranded_gene_counts', skip=1, sep='\t', row.names = 1) %>%
  select(c(6:11))


# distinguish between treatment and control 
# IMPORTANT: currently assuming files by alphabetical order are 3 treatment and
# 3 control, may not always be the case
colData = data.frame(condition=c(rep('treatment', 3), rep('control', 3)),
                     type=rep('paired-end', 6),
                     row.names=colnames(gene_mat))
colData = data.frame(condition=c('treatment', 'control', 'control', 'control', 'treatment', 'treatment'),
                     type=rep('paired-end', 6),
                     row.names=colnames(gene_mat))

# run DESeq2 on reads matrix and save
dds <- DESeqDataSetFromMatrix(countData = gene_mat,
                              colData = colData,
                              design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
de_results <- dds[keep,] %>% DESeq() %>% results() %>% data.frame()

# write data as file
write.csv(x=de_results, file='/data/RNASeq/DESeq/nkap_de_results.csv')