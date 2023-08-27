library('ggplot2')
library("org.Hs.eg.db")
library("AnnotationDbi")

# convert ENSG to Symbol
de_modified = de_results %>%
  mutate(ENSG = gsub("\\..*","", row.names(de_results))) %>%
  filter(!is.na(padj))

de_modified = de_modified %>%
  mutate(Symbol = mapIds(org.Hs.eg.db, keys=de_modified$ENSG, column="SYMBOL", keytype="ENSEMBL", multiVals="first")) %>%
  mutate(Color = case_when(
    padj >= 0.01 ~ 'ns',
    log2FoldChange > 1 ~ 'inc',
    log2FoldChange < -1 ~ 'dec',
    .default = 'nfc'))

# plot genes in volcano plot
volcano <- ggplot(data=de_modified, mapping=aes(x=log2FoldChange, y=-log10(padj), color=Color, )) +
  geom_point(alpha=0.2, size=1, ) +
  ggtitle('GRAS1 Knockdown Genes Differential Expression') +
  geom_text(label=paste('Increasing:', sum(de_modified$Color == 'inc')), x=3, y=-log10(min(de_modified$padj[de_modified$padj != 0]))) +
  geom_text(label=paste('Decreasing:', sum(de_modified$Color == 'dec')), x=-3, y=-log10(min(de_modified$padj[de_modified$padj != 0]))) +
  geom_text(label=paste('NS/NFC:', sum(de_modified$Color == 'nfc' | de_modified$Color == 'ns')), x=0, y=-log10(min(de_modified$padj[de_modified$padj != 0]))) +
  geom_vline(xintercept=c(1, -1)) +
  geom_hline(yintercept=2) +
  scale_color_manual(values=c('blue', 'red', 'black', 'black')) +
  theme_minimal() +
  theme(legend.position = 'none')

# display volcano plot
volcano
