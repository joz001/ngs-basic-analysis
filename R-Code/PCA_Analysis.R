library(ggplot2)
library(ggforce)


gene_mat = read.csv(path, skip=1, sep='\t', row.names=1) %>%
  dplyr::select(c(6:9))

colData = data.frame(condition=c(rep('treatment', 3), rep('control', 3)),
                     type=rep('paired-end', 6),
                     row.names=colnames(gene_mat))

pca_plot <- function(gene_mat, colData, title){
  pre_pca = t(gene_mat)
  pre_pca = pre_pca[, which(apply(pre_pca, 2, var) != 0)]
  pca <- prcomp(pre_pca, scale=TRUE)
  pca.data <- data.frame(Sample = rownames(pca$x),
                         Type = colData$condition,
                         X = pca$x[,1], Y = pca$x[,2])
  min_pca_axis = min(pca.data$X, pca.data$Y)
  max_pca_axis = max(pca.data$X, pca.data$Y)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  p <- ggplot(data = pca.data, aes(x=X, y =Y, label= Sample)) + 
    geom_point(aes(color = Type), size = 3) +
    geom_mark_ellipse(aes(color = Type, label = NULL), 
                      expand = unit(4, "mm")) +
    labs(title = title,
         x = paste("PC1 - ", pca.var.per[1], "%", sep = ""), 
         y = paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
    scale_x_continuous(limits = c(min_pca_axis * 1.25, max_pca_axis * 1.25)) +
    scale_y_continuous(limits = c(min_pca_axis * 1.25, max_pca_axis * 1.25)) +
    scale_color_manual(values = c("treatment" = "red", "control" = "blue")) +
    theme_bw() +
    theme(legend.title = element_blank())

  return(p)
}
plot = pca_plot(gene_mat, colData, '')
plot
