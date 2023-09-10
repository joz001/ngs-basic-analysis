library(clusterProfiler)
library(ggplot2)
library(readxl)
library(org.Hs.eg.db)
library(xlsx)

# GO Summary File
GO_summary <- function(signif_genes){
  go_df <- enrichGO(gene = signif_genes, keyType = "ENSEMBL",
                  OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH",
                  qvalueCutoff = 0.05, readable = TRUE) %>%
    data.frame
  if (length(go_df$ID) < 1)
    return(NULL)
  return(go_df[order(go_df$p.adjust),])
}

up_go = de_modified %>%
  filter(Color == 'inc') %>%
  pull(ENSG) %>% GO_summary

down_go = de_modified %>%
  filter(Color == 'dec') %>%
  pull(ENSG) %>% GO_summary 

all_go = de_modified %>%
  filter(Color == 'inc' | Color == 'dec') %>%
  pull(ENSG) %>% GO_summary

save_as = ''
write.xlsx(x = up_go, file = save_as, sheetName = "Upregulated Genes GO", append = TRUE)
write.xlsx(x = down_go, file = save_as, sheetName = "Downregulated Genes GO", append = TRUE)
write.xlsx(x = all_go, file = save_as, sheetName = "All Genes GO", append = TRUE)


# visualize top 10 up and down pathways
go_bar <- function(up_go, down_go, title){
  merged = list(up_go[1:10,], down_go[1:10,])
  top_go_df = data.frame('Pathway' = c(merged[[1]]$Description, merged[[2]]$Description),
                         'negLogPAdj_Value' = -log10(c(merged[[1]]$p.adjust, merged[[2]]$p.adjust)),
                         'count' = c(merged[[1]]$Count, merged[[2]]$Count),
                         'direction_discrete' = factor(c(rep("Up-Regulated", nrow(merged[[1]])),
                                                         rep("Down-Regulated", nrow(merged[[2]]))),
                                                       levels = c("Up-Regulated", "Down-Regulated")))
    
  top_pathways_bar = ggplot(data=top_go_df, aes(x=reorder(Pathway, negLogPAdj_Value), 
                                                y=negLogPAdj_Value, fill=direction_discrete)) +
    facet_grid(direction_discrete ~ ., scales = "free", space = "free_x") +
    geom_bar(stat="identity", show.legend = FALSE) +
    geom_text(aes(label = count, y = negLogPAdj_Value), size = 3, hjust = -0.5) +
    coord_flip() +
    labs(title = title,
         x = element_blank(), 
         y = expression(Enrichment~(-log[10](P-Value)))) +
    scale_y_continuous(limits = c(0, max(1.2 * top_go_df$negLogPAdj_Value))) +
    scale_fill_manual(values = c("red", "blue")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 12))
  
  return(top_pathways_bar)
}
bar = go_bar(up_go, down_go, '')
bar
