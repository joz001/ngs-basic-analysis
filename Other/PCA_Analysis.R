library(ggplot2)
library(ggforce)


path = ""
gene_mat = read.csv(path, skip=1, sep='\t', row.names=1) %>%
  select(c(6:11))

pre_pca = t(final)
pre_pca = pre_pca[, which(apply(pre_pca, 2, var) != 0)]
pca <- prcomp(pre_pca, scale=TRUE)

rownames(pca$x)[order(pca$x[ ,1], decreasing=TRUE)[1]]


pca_plot <- function(pca, title){
  pca.data <- data.frame(Sample = rownames(pca$x),
                         Type = c(rep("exp", 3), rep("ctl", 3)),
                         X = pca$x[,1], Y = pca$x[,2])
  min_pca_axis = min(pca.data$X, pca.data$Y)
  max_pca_axis = max(pca.data$X, pca.data$Y)
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  
  p <- ggplot(data = pca.data, aes(x=X, y =Y, label= Sample)) + 
    geom_point(aes(color = Type), size = 3) +
    geom_mark_ellipse(aes(color = Type, label = NULL), 
                      expand = unit(4, "mm")) +
    labs(
      title = title,
      x = paste("PC1 - ", pca.var.per[1], "%", sep = ""), 
      y = paste("PC2 - ", pca.var.per[2], "%", sep = ""),
    ) + 
    scale_x_continuous(limits = c(min_pca_axis * 1.25, max_pca_axis * 1.25)) +
    scale_y_continuous(limits = c(min_pca_axis * 1.25, max_pca_axis * 1.25)) +
    scale_color_manual(labels = c("exp" = "Overexpression", "ctl" = "Control"), 
                       values = c("exp" = "red", "ctl" = "blue")) +
    theme_bw() +
    theme(legend.title = element_blank()
    )
  
  return(p)
}
plot = pca_plot(pca, "")
plot



## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes


pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
scree <- pca.var/sum(pca.var)
barplot((scree[1:10]*100), main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
