library(dplyr)
library(data.table)
library(ggplot2)

read_in <- function(file, donor_num){
  f = read.csv(file, sep='\t', row.names=1)
  f = f[colSums(f) >= 6000]
  f = f[rowSums(f) > 0,]
  f = rbind(f, donor=donor_num)
  f = setDT(f, keep.rownames = TRUE)
  colnames(f) = c('genes', paste(colnames(f)[-1], "_", strsplit(file, '[.]')[[1]][1], sep=""))
  return(f)
}

merge_files <- function(file_1, file_2){
  return(merge(file_1, file_2, all=TRUE))
}

files = c(list.files()[1:8] %>% lapply(read_in, 'D28'),
          list.files()[9:16] %>% lapply(read_in, 'D29'),
          list.files()[17:24] %>% lapply(read_in, 'D30'),
          list.files()[25:32] %>% lapply(read_in, 'D31'))

{
  start = Sys.time()
  count_mat = Reduce(merge_files, files)
  count_mat[is.na(count_mat)] = 0
  print(Sys.time() - start)
}
