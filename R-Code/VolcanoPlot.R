library('ggplot2')

# plot genes in volcano plot
make_volcano <- function(de_modified, title){
  volcano <- ggplot(data=de_modified, mapping=aes(x=log2FoldChange, y=-log10(padj), color=Color)) +
    geom_point(alpha=0.2, size=1) +
    ggtitle(title) +
    geom_text(label=paste('Increasing:', sum(de_modified$Color == 'inc')), x=3, y=-log10(min(de_modified$padj[de_modified$padj != 0]))) +
    geom_text(label=paste('Decreasing:', sum(de_modified$Color == 'dec')), x=-3, y=-log10(min(de_modified$padj[de_modified$padj != 0]))) +
    geom_text(label=paste('NS/NFC:', sum(de_modified$Color == 'nfc' | de_modified$Color == 'ns')), x=0, y=-log10(min(de_modified$padj[de_modified$padj != 0]))) +
    geom_vline(xintercept=c(1, -1)) +
    geom_hline(yintercept=-log10(0.01)) +
    scale_color_manual(values=c('blue', 'red', 'black', 'black')) +
    theme_minimal() +
    theme(legend.position = 'none')
  
  return(volcano)
}

volcano = make_volcano(de_modified, '')
volcano
