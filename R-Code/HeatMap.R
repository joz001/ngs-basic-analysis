library(xlsx)
library(dplyr)
library(tidyverse)
library(reshape2)

rnaseq = read.xlsx('', sheetIndex = 1, row.names=1)
all_go = read.xlsx('', sheetIndex = 1)
genes_counts = read.csv('', sep='\t', skip=1, row.names=1) %>%
  dplyr::select(c(6:11))

colnames(genes_counts) = c('exp1', 'exp2', 'exp3', 'nc1', 'nc2', 'nc3')
merged = merge(genes_counts, rnaseq, by='row.names') %>%
  column_to_rownames('Row.names')

paths = c('')

path_heat_map <- function(go_sheet, path, rnaseq_sheet, gene_counts_sheet){
  symbol_in_go_path = go_sheet %>% 
    subset(Description %in% path) %>% 
    pull(geneID) %>% 
    strsplit('/') %>%
    unlist
  
  ensg_in_go_path = rnaseq_sheet %>%
    filter(symbol %in% symbol_in_go_path) %>%
    row.names %>%
    unlist
  
  path_gene_counts = gene_counts_sheet %>% 
    filter(row.names(gene_counts_sheet) %in% ensg_in_go_path) %>%
    merge(rnaseq_sheet['symbol'], by='row.names')
  
  rownames(path_gene_counts) = path_gene_counts[,1]
  path_gene_counts[,1] = NULL
  
  path_gene_counts$nc_mean = path_gene_counts[c(4:6)] %>% rowMeans
  path_gene_counts$exp_mean = path_gene_counts[c(1:3)] %>% rowMeans
  path_gene_counts[c(1:6)] = path_gene_counts[c(1:6)] / path_gene_counts$nc_mean
  path_gene_counts  = path_gene_counts %>% arrange(exp_mean/nc_mean)
  return(path_gene_counts)
}

plot_heat_map <- function(path, single_gene_counts){
  melted = melt(single_gene_counts[-c(8,9)][c(1:min(20, nrow(single_gene_counts))),])
  melted$symbol = factor(melted$symbol, levels = single_gene_counts$symbol)
  melted$type = factor(substr(melted$variable, 1, nchar(as.character(melted$variable))-1), c("nc", "exp"))
  
  p <- ggplot(melted, aes(x=variable, y=symbol, group=type)) +
    geom_tile(aes(fill=log2(value))) +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(
      fill = "Log2 FC",
      title = paste(path, "Heat Map"),
    ) +
    theme_minimal() +
    theme(
      panel.spacing = unit(-1, "lines"),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(size = 12)
    )  +
    facet_grid(.~type, scales = "free_x", margins = FALSE, switch = "x",
               labeller = as_labeller(c("nc" = "Control", "exp" = "Knockdown")))
  return(p)
}

num = 4
path_genes = path_heat_map(go_sheet=all_go, paths[num], rnaseq_sheet=rnaseq, gene_counts_sheet=genes_counts)
p = plot_heat_map(paths[num], path_genes)
p 
