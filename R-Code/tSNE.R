library(Rtsne)

norm_count_mat = count_mat[-nrow(count_mat),-1] %>%
  sapply(FUN=as.numeric) %>% data.frame
norm_count_mat = mapply("/", norm_count_mat, colSums(norm_count_mat)) * 6000

{
start = Sys.time()
tsne = Rtsne(t(norm_count_mat))
print(Sys.time() - start)
}


tsne_plot <- data.frame(x = tsne$Y[,2], 
                        y = tsne$Y[,1])

plot = tsne_plot %>%
  ggplot(aes(x=x, y=y)) +
  geom_point() +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  labs(x=element_blank(),
       y=element_blank()) +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL)
plot

donor_plot = tsne_plot %>%
  cbind('donor'=count_mat[nrow(count_mat),-1] %>% unlist) %>%
  ggplot(aes(x=x, y=y, color=donor)) +
  geom_point(size=1.5) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  labs(x=element_blank(),
       y=element_blank()) +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL)
donor_plot

# alpha gene
tsne_plot$GCG = count_mat[count_mat$genes %like% 'ENSG00000115263',-1] %>% as.numeric
# alpha + beta tf
tsne_plot$MAFB = count_mat[count_mat$genes %like% 'ENSG00000204103',-1] %>% as.numeric
# beta gene
tsne_plot$INS = count_mat[count_mat$genes %like% 'ENSG00000254647',-1] %>% as.numeric
# beta tf
tsne_plot$MAFA = count_mat[count_mat$genes %like% 'ENSG00000182759',-1] %>% as.numeric
# delta gene
tsne_plot$SST = count_mat[count_mat$genes %like% 'ENSG00000157005',-1] %>% as.numeric
# pp gene
tsne_plot$PPY = count_mat[count_mat$genes %like% 'ENSG00000108849',-1] %>% as.numeric
# acinar
tsne_plot$PRSS1 = count_mat[count_mat$genes %like% 'ENSG00000204983',-1] %>% as.numeric
# duct
tsne_plot$KRT19 = count_mat[count_mat$genes %like% 'ENSG00000171345',-1] %>% as.numeric
# mes
tsne_plot$COL1A1 = count_mat[count_mat$genes %like% 'ENSG00000108821',-1] %>% as.numeric


single_gene <- function(gene){
  p = ggplot(tsne_plot, aes(x=x,y=y))  +
    geom_point(aes(color=.data[[gene]]), size=.5) +
    scale_color_gradient(low = "#eeeeee", high = "blue") +
    theme_minimal() +
    guides(color=guide_colorbar(label.position = 'left', barwidth = .5))+
    theme(
      legend.title = element_blank(),
      legend.text = element_text(angle=90, hjust=0.5),
      plot.title = element_text(hjust=.5, face = 'italic')
    ) +
    labs(title=gene,
         x=element_blank(),
         y=element_blank()) +
    scale_x_continuous(breaks = NULL) + 
    scale_y_continuous(breaks = NULL)
  return(p)
}

# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

genes = c('GCG', 'PPY', 'INS', 'PRSS1', 'SST', 'KRT19', 'MAFA', 'MAFB')
plots = lapply(genes, single_gene)

multiplot(plotlist = plots, cols=4)
