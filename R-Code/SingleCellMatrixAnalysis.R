library(SingleCellExperiment)
library(tictoc)
library(dplyr)
library(ggplot2)
setwd('/data/scrna_seq/count_matrices/')
files = list.files('/data/scrna_seq/count_matrices/')

make_counts <- function(x){
  f1 = read.table(x, header=TRUE, row.names=1)
  colnames(f1) = paste(substring(x, 1, 10), colnames(f1), sep='_')
  return(f1)
}

merge_by_rownames <- function(x, y){
  merged = merge(x, y, by=0, all=TRUE) %>%
    transform(row.names = Row.names, Row.names=NULL)
  merged[is.na(merged)] = 0
  return(merged)
}

{
  tic()
  counts = lapply(files, make_counts)
  toc()
  tic()
  all_counts = Reduce(merge_by_rownames, counts) 
  toc()
}

total_counts = colSums(all_counts) %>% data.frame('sample_counts' = .)
filter_6k_transcripts = total_counts %>% filter(sample_counts >= 6000)
min_6k = all_counts[colnames(all_counts) %in% row.names(filter_6k_transcripts)]

min_6k$rowsums = rowSums(min_6k)
min_6k = min_6k[min_6k$rowsums > 0, ]

colsums = colSums(min_6k)
norm = mapply('/', min_6k, colsums)

ggplot(data = total_counts, aes(x=log10(sample_counts))) +
  geom_bar(stat='bin', binwidth=.5)
   

freq = ggplot(all_counts, aes(x=rowSums)) +
  geom_histogram(binwidth = 1) +
  xlim(0, 10000)
  #ylim(0, 200)
freq

p = pca_plot(all_counts[1:3019], '', '')
p
