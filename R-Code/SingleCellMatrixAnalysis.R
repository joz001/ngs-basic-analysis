library(dplyr)
library(data.table)
library(ggplot2)


read_in <- function(file){
  
  f = read.csv(file, sep='\t', row.names=1)
  f = f[colSums(f) >= 6000]
  f = f[rowSums(f) > 0,]
  f = setDT(f, keep.rownames = TRUE)
  colnames(f) = c('genes', paste(colnames(f)[-1], "_", strsplit(file, '[.]')[[1]][1], sep=""))
  return(f)
}

merge_files <- function(file_1, file_2){
  return(merge(file_1, file_2, all=TRUE))
}


{
  start = Sys.time()
  files = list.files() %>% unlist %>%
    lapply(read_in)
  count_mat = Reduce(merge_files, files)
  count_mat[is.na(count_mat)] = 0
  print(Sys.time() - start)
}

norm_count_mat = count_mat[,-1]
norm_count_mat = mapply("/", norm_count_mat, colSums(norm_count_mat))
norm_count_mat = norm_count_mat * 6000
#row.names(norm_count_mat) = count_mat[,1][[1]]

pre_pca = t(norm_count_mat)
pre_pca = pre_pca[, which(apply(pre_pca, 2, var) != 0)]
pca <- prcomp(pre_pca, scale=TRUE)
pca.data <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[,1], Y = pca$x[,2])
min_pca_axis = min(pca.data$X, pca.data$Y)
max_pca_axis = max(pca.data$X, pca.data$Y)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

p <- pca.data %>% 
  ggplot(aes(x=X, y=Y, label=Sample)) + 
  geom_point(size=3) +
  labs(title = 'title',
       x = paste("PC1 - ", pca.var.per[1], "%", sep = ""), 
       y = paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_bw() + 
  scale_x_continuous(limits = c(min_pca_axis * 1.25, max_pca_axis * 1.25)) +
  scale_y_continuous(limits = c(min_pca_axis * 1.25, max_pca_axis * 1.25)) +
  theme(legend.title = element_blank(),
        legend.position = 'top',
        axis.title = element_text())

p
