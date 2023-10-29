library('ggplot2')
library('scales')

make_ma <- function(de_modified, title){
  ma <- de_modified %>%
    ggplot(mapping=aes(x=baseMean, y=log2FoldChange, color=Color)) +
    geom_point(alpha=0.2, size=1) +
    labs(title=title,
         x='Base Mean',
         y=expression('Log'[2]*' Fold Change')) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(trans='log10', 
                       labels=trans_format('log10', math_format(10^.x))) +
    scale_color_manual(values=c('blue', 'red', '#33333333', '#33333333')) +
    theme_minimal() +
    theme(legend.position = 'none')
  return(ma)
}
ma = make_ma(de_modified, 'GRAS1 KD MA Plot')
ma